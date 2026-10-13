"""Illustrative full-grid BLS maxima under half-day residual block permutations.

This is not a validated generative noise model or a pilot-wide FAP estimate.
Progress is checkpointed after every permutation; rerunning resumes completed
replicates. Maximum two worker processes. No network requests or paid APIs.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

from search import CACHE, DURATIONS, ROOT, flatten, known_for, mask_known, read_curve, scan

FROZEN_SEARCH_SHA256 = "0d1216f3dd84114e0be8d44c5030324bd108a907be47573aeca76f296253614a"
_ACTIVE_SEARCH_PATH = None


def activate_search(path, expected_hash):
    """Use explicitly selected, integrity-checked implementation in every worker.

    Worker activation also supports multiprocessing spawn, where parent globals
    are not inherited. All input files and stored null measurements are unchanged.
    """
    global _ACTIVE_SEARCH_PATH, CACHE, DURATIONS, ROOT, flatten, known_for, mask_known, read_curve, scan
    path = Path(path)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected_hash:
        raise ValueError("Selected search implementation does not match required source hash.")
    if _ACTIVE_SEARCH_PATH == str(path):
        return
    spec = importlib.util.spec_from_file_location("_null_calibration_selected_search", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    CACHE, DURATIONS, ROOT = module.CACHE, module.DURATIONS, module.ROOT
    flatten, known_for, mask_known = module.flatten, module.known_for, module.mask_known
    read_curve, scan = module.read_curve, module.scan
    _ACTIVE_SEARCH_PATH = str(path)


def prepare(record, width, edge_cut=.25):
    t, f, e, _ = read_curve(record["local_path"])
    tb, fb, eb = flatten(t, f, e, width)
    edges = np.r_[0, np.flatnonzero(np.diff(t) > .2) + 1, len(t)]
    good = np.zeros(len(tb), bool)
    for lo, hi in zip(edges[:-1], edges[1:]):
        good |= (tb >= t[lo] + edge_cut) & (tb <= t[hi - 1] - edge_cut)
    good &= mask_known(tb, known_for(record["tic"]))
    tt, ff, ee = tb[good], fb[good], eb[good]
    # Identical scalar errors for original and all shuffled realizations. This
    # intentionally avoids relocating heteroscedastic uncertainty patterns.
    error_scalar = float(np.median(ee))
    uniform = np.full(len(tt), error_scalar)
    inputs = CACHE / "null_calibration"
    inputs.mkdir(parents=True, exist_ok=True)
    path = inputs / f"tic{record['tic']}_s{record['sector']}_w{width:g}_edge{edge_cut:g}.npz"
    np.savez_compressed(path, t=tt, f=ff, e=uniform, original_e=ee)
    peak, result = scan(tt, ff, uniform, n_peaks=1)
    original_peak, _ = scan(tt, ff, ee, n_peaks=1)
    return {"key": f"tic{record['tic']}_s{record['sector']}_w{width:g}",
            "tic": record["tic"], "sector": record["sector"], "width_days": width,
            "edge_cut_days": edge_cut, "input_cache": str(path),
            "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "raw_sha256": record["sha256"], "data_uri": record["data_uri"],
            "n_points": len(tt), "constant_error_ppm": error_scalar * 1e6,
            "original_error_ratio_max_min": float(ee.max() / ee.min()),
            "observed_peak": peak[0], "original_diagonal_peak": original_peak[0],
            "period_grid_count": len(result.period), "replicates": []}


def block_indices(t, block_days):
    """Split every gap-separated stretch into blocks of at most block_days."""
    edges = np.r_[0, np.flatnonzero(np.diff(t) > .2) + 1, len(t)]
    blocks = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ids = np.floor((t[lo:hi] - t[lo]) / block_days).astype(int)
        for block in np.unique(ids):
            blocks.append(np.flatnonzero(ids == block) + lo)
    return blocks


def replicate(task):
    path, tic, sector, width, rep, seed, block_days, search_path, source_hash = task
    activate_search(search_path, source_hash)
    with np.load(path) as arrays:
        t, f, e = arrays["t"], arrays["f"], arrays["e"]
    blocks = block_indices(t, block_days)
    # Width omitted: use the same permutation for paired preprocessing widths
    # when their retained sampling/block boundaries are identical.
    rng = np.random.default_rng(np.random.SeedSequence([seed, tic, sector, rep]))
    ordering = rng.permutation(len(blocks))
    source = np.concatenate([blocks[int(i)] for i in ordering])
    shuffled = np.median(f) + (f - np.median(f))[source]
    peaks, result = scan(t, shuffled, e, n_peaks=1)
    return {"tic": tic, "sector": sector, "width_days": width,
            "replicate": rep, "seed_components": [seed, tic, sector, rep],
            "n_blocks": len(blocks), "block_order": ordering.tolist(),
            "max_peak": peaks[0], "period_grid_count": len(result.period)}


def summarize(entry):
    maxima = np.array([r["max_peak"]["snr"] for r in entry["replicates"]])
    observed = entry["observed_peak"]["snr"]
    if not len(maxima):
        return None
    exceedances = int((maxima >= observed).sum())
    return {"n": len(maxima), "median_maximum_snr": float(np.median(maxima)),
            "q05_maximum_snr": float(np.quantile(maxima, .05)),
            "q95_maximum_snr": float(np.quantile(maxima, .95)),
            "largest_maximum_snr": float(maxima.max()),
            "n_null_at_least_observed": exceedances,
            "permutation_tail_fraction_plus_one": (1 + exceedances) / (1 + len(maxima)),
            "minimum_resolvable_fraction": 1 / (1 + len(maxima))}


def checkpoint(result, out):
    for entry in result["targets"]:
        entry["summary"] = summarize(entry)
    result["updated_utc"] = datetime.now(timezone.utc).isoformat()
    temp = out.with_suffix(".tmp")
    temp.write_text(json.dumps(result, indent=2, allow_nan=False))
    temp.replace(out)


def write_report(result):
    lines = ["# Illustrative empirical BLS search maxima", "",
             "These permutations compare each observed light curve with maxima over the exact same production period/duration/epoch grid. They offer a caution about ranking statistics; they do not calibrate discovery probabilities or the full pilot search.", "",
             "Residuals are obtained with the 1-day and 2-day median detrenders, 10-minute bins, known-planet masks, and 0.25-day removal around quality-filtered raw-data boundaries/gaps. Half-day contiguous blocks are randomly reordered onto the original timestamps. Original and shuffled data use identical constant errors equal to median diagonal errors, avoiding relocation of heteroscedastic errors. The original diagonal-error maximum is also retained in the JSON for comparison.", "",
             "|TIC|Sector|Width (d)|Observed maximum S/N|Null n|Median null maximum|95th percentile|Null ≥ observed|Plus-one tail fraction|",
             "|---|---|---|---|---|---|---|---|---|"]
    for entry in result["targets"]:
        summary = entry.get("summary")
        if summary:
            lines.append(f"|{entry['tic']}|{entry['sector']}|{entry['width_days']}|"
                         f"{entry['observed_peak']['snr']:.2f}|{summary['n']}|"
                         f"{summary['median_maximum_snr']:.2f}|{summary['q95_maximum_snr']:.2f}|"
                         f"{summary['n_null_at_least_observed']}|{summary['permutation_tail_fraction_plus_one']:.3f}|")
    lines += ["", "Limitations:", "",
              "- Only selected individual sectors of two targets are examined. These fractions are not global false-alarm probabilities across stars, sectors, contiguous-sector searches, detrending choices, or selected candidates. Choosing a quiet comparison star after inspecting data also limits generalization.",
              "- With 24 permutations, the smallest plus-one fraction is 1/25 = 0.04; tail probabilities cannot support high significance. The empirical 95th percentile is unstable at this sample size.",
              "- Permuting already flattened blocks fixes the detrending realization. It does not propagate uncertainties in stellar variability removal or adaptive masking. Correlations within a half-day block survive; longer correlations and their relationship to gaps are disrupted, and block joins introduce new discontinuities.",
              "- The procedure preserves observed residual excursions, including possible real events; it is not a null guaranteed to contain only stochastic noise. It assumes approximate exchangeability of blocks and cannot simulate spacecraft artifacts tied to particular timestamps.",
              "- A scalar median uncertainty slightly changes the observed ranking relative to production heteroscedastic weights. Inference uses observed and shuffled maxima under the same scalar convention.",
              "- Widths use corresponding seed/permutation when sampling agrees. Per-width summaries include the full period grid; they do not automatically correct for selecting the maximum across both widths. Paired maxima across widths are summarized below where available.", ""]
    for pair in paired_widths(result):
        lines.append(f"TIC {pair['tic']}, sector {pair['sector']}: maximum across both widths is {pair['observed_maximum_snr']:.2f}; {pair['n_null_at_least_observed']}/{pair['n']} paired null maxima are at least this large (plus-one fraction {pair['permutation_tail_fraction_plus_one']:.3f}). This includes the two tested widths only and remains subject to the limitations above.")
        lines.append("")
    if result.get("search_provenance", {}).get("frozen_search"):
        lines += ["## Pinned search implementation", "",
                  "The original null calculations used `src/search.py` from Git commit `7ecc52c`. Its exact bytes are preserved as `src/search_null_v1.py`, SHA-256 `" + FROZEN_SEARCH_SHA256 + "`. Later production changes added time-system and odd/even diagnostic guards; the saved calibration is tied to its original implementation.", "",
                  "Resume or regenerate these inputs with `python src/null_calibration.py --frozen-search`. This explicitly selects the preserved version and verifies its source hash; completed permutations are retained without repetition. The default invocation still selects production `src/search.py` and rejects this checkpoint if its hash differs. Cached raw and null-input data remain intact, and cached inputs can be reconstructed using the pinned code and retained raw products.", ""]
    path = ROOT / "reports" / "null_calibration.md"
    path.write_text("\n".join(lines))
    validation = ROOT / "reports" / "validation.md"
    if validation.exists():
        original = validation.read_text().split("\n## Empirical full-grid null appendix")[0]
        validation.write_text(original + "\n## Empirical full-grid null appendix\n\nSee [null_calibration.md](null_calibration.md) and `results/null_calibration.json` for half-day block permutation comparisons. These illustrative maxima calibrate neither the full multi-target search nor discovery probabilities.\n")


def paired_widths(result):
    pairs = []
    identities = {(entry["tic"], entry["sector"]) for entry in result["targets"]}
    for tic, sector in sorted(identities):
        entries = [entry for entry in result["targets"] if (entry["tic"], entry["sector"]) == (tic, sector)]
        if len(entries) != 2:
            continue
        a = {r["replicate"]: r for r in entries[0]["replicates"]}
        b = {r["replicate"]: r for r in entries[1]["replicates"]}
        shared = sorted(a.keys() & b.keys())
        same_blocks = all(a[i]["block_order"] == b[i]["block_order"] for i in shared)
        same_sampling = entries[0]["n_points"] == entries[1]["n_points"]
        maximum_time_difference = None
        if same_sampling:
            with np.load(entries[0]["input_cache"]) as original_a, np.load(entries[1]["input_cache"]) as original_b:
                ta, tb = original_a["t"], original_b["t"]
                maximum_time_difference = float(np.max(np.abs(ta - tb))) * 86400
                same_sampling = bool(np.allclose(ta, tb, rtol=0, atol=1e-5))
                same_sampling &= [len(block) for block in block_indices(ta, result["block_days"])] == [len(block) for block in block_indices(tb, result["block_days"])]
        if not shared or not (same_blocks and same_sampling):
            continue
        observed = max(entry["observed_peak"]["snr"] for entry in entries)
        maxima = np.array([max(a[i]["max_peak"]["snr"], b[i]["max_peak"]["snr"]) for i in shared])
        count = int((maxima >= observed).sum())
        pairs.append({"tic": tic, "sector": sector, "n": len(shared),
                      "maximum_width_sampling_difference_seconds": maximum_time_difference,
                      "observed_maximum_snr": observed, "n_null_at_least_observed": count,
                      "permutation_tail_fraction_plus_one": (1 + count) / (1 + len(shared)),
                      "median_null_maximum_snr": float(np.median(maxima)),
                      "q95_null_maximum_snr": float(np.quantile(maxima, .95))})
    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", nargs="+", default=["278727563:94", "307210830:2"])
    parser.add_argument("--replicates", type=int, default=24)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--seed", type=int, default=260113)
    parser.add_argument("--block-days", type=float, default=.5)
    parser.add_argument("--frozen-search", action="store_true",
                        help="Use exact search.py from commit 7ecc52c, preserved as search_null_v1.py.")
    args = parser.parse_args()
    search_path = ROOT / "src" / ("search_null_v1.py" if args.frozen_search else "search.py")
    selected_hash = hashlib.sha256(search_path.read_bytes()).hexdigest()
    if args.frozen_search and selected_hash != FROZEN_SEARCH_SHA256:
        raise ValueError("Pinned search source no longer matches its original commit/hash.")
    activate_search(search_path, selected_hash)
    out = ROOT / "results" / "null_calibration.json"
    result = json.loads(out.read_text()) if out.exists() else {
        "created_utc": datetime.now(timezone.utc).isoformat(), "seed": args.seed,
        "block_days": args.block_days, "desired_replicates_per_width": args.replicates,
        "period_min_days": .7, "period_max_days": "min(13, retained baseline/2)",
        "durations_days": DURATIONS.tolist(), "frequency_factor": .7, "oversample": 10,
        "search_source_sha256": selected_hash,
        "noise_assumption": "Exchangeable pre-flattened half-day residual blocks; illustrative only",
        "targets": []}
    if result["seed"] != args.seed or result["block_days"] != args.block_days:
        raise ValueError("Checkpoint seed/block size differs; use matching options.")
    if result["search_source_sha256"] != selected_hash:
        raise ValueError("Search implementation changed; archive the existing checkpoint before a new calibration.")
    result["search_provenance"] = {"frozen_search": args.frozen_search,
                                   "selected_module": str(search_path.relative_to(ROOT)),
                                   "source_sha256": selected_hash,
                                   "origin_git_commit": "7ecc52c" if args.frozen_search else None}
    records = json.loads((ROOT / "data" / "observations.json").read_text())
    for target in args.targets:
        tic, sector = map(int, target.split(":"))
        record = next(r for r in records if r["tic"] == tic and r["sector"] == sector)
        for width in (1., 2.):
            key = f"tic{tic}_s{sector}_w{width:g}"
            existing = next((entry for entry in result["targets"] if entry["key"] == key), None)
            if existing and not Path(existing["input_cache"]).exists():
                # Reconstruct disposable inputs from the retained raw FITS and
                # saved preprocessing. Raw scientific evidence is never deleted.
                reconstructed = prepare(record, width, existing["edge_cut_days"])
                if reconstructed["raw_sha256"] != existing["raw_sha256"] or reconstructed["n_points"] != existing["n_points"]:
                    raise ValueError("Reconstructed calibration input differs from checkpoint.")
                if not np.isclose(reconstructed["observed_peak"]["snr"], existing["observed_peak"]["snr"], rtol=1e-9):
                    raise ValueError("Reconstructed calibration statistic differs from checkpoint.")
                existing["input_cache"] = reconstructed["input_cache"]
                existing["input_sha256"] = reconstructed["input_sha256"]
                checkpoint(result, out)
            elif not existing:
                result["targets"].append(prepare(record, width))
                checkpoint(result, out)
                print("PREPARED", key, flush=True)
    jobs = []
    selected = {tuple(map(int, target.split(":"))) for target in args.targets}
    for entry in result["targets"]:
        if (entry["tic"], entry["sector"]) not in selected:
            continue
        completed = {r["replicate"] for r in entry["replicates"]}
        for rep in range(args.replicates):
            if rep not in completed:
                jobs.append((entry["input_cache"], entry["tic"], entry["sector"], entry["width_days"],
                             rep, args.seed, args.block_days, str(search_path), selected_hash))
    with ProcessPoolExecutor(max_workers=max(1, min(2, args.workers))) as pool:
        futures = [pool.submit(replicate, job) for job in jobs]
        for future in as_completed(futures):
            row = future.result()
            entry = next(entry for entry in result["targets"] if
                         (entry["tic"], entry["sector"], entry["width_days"]) ==
                         (row["tic"], row["sector"], row["width_days"]))
            entry["replicates"].append(row)
            entry["replicates"].sort(key=lambda r: r["replicate"])
            checkpoint(result, out)
            print("NULL", entry["key"], row["replicate"], row["max_peak"]["snr"], flush=True)
    result["paired_width_summaries"] = paired_widths(result)
    checkpoint(result, out)
    write_report(result)


if __name__ == "__main__":
    main()
