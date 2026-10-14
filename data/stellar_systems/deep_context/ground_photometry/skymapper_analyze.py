"""Reproduce the public SkyMapper DR4 epoch coverage check from retained files.

Run with /workspace/openrouter-check-venv/bin/python followed by this path.
The ephemeris is imported from the independently established TESS source summary;
no ground-data ephemeris fitting or magnitude clipping is performed.
"""
from pathlib import Path
import hashlib
import json
import datetime
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy.utils import iers

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
iers.conf.auto_download = False
source_path = REPO / "candidates/stellar_systems/tic273879360/source_summary.json"
ephem = json.loads(source_path.read_text())["fitted_phenomenological_shape"]
raw_path = ROOT / "skymapper_dr4_epochs.csv"
df = pd.read_csv(raw_path)
master = pd.read_csv(ROOT / "skymapper_dr4_master.csv", dtype={"gaia_dr3_id1": str})
target = SkyCoord(32.252462825*u.deg, -71.323170344*u.deg)
location = EarthLocation.from_geodetic(149.0661*u.deg, -31.2733*u.deg, 1165*u.m)
midtime = Time(df.mjd_start.to_numpy() + df.exp_time.to_numpy()/172800,
               format="mjd", scale="utc", location=location)
df["bjd_tdb_mid"] = (midtime.tdb + midtime.light_travel_time(target)).jd
df["phase"] = ((df.bjd_tdb_mid - (2457000 + ephem["epoch_btjd"]))
               / ephem["period_days"]) % 1
df["primary_offset_days"] = ((df.phase + .5) % 1 - .5) * ephem["period_days"]
df["secondary_offset_days"] = ((df.phase - ephem["secondary_phase"] + .5) % 1 - .5) * ephem["period_days"]
df["in_primary_window"] = abs(df.primary_offset_days) <= ephem["primary_duration_days"]/2
df["in_secondary_window"] = abs(df.secondary_offset_days) <= ephem["secondary_duration_days"]/2
df["strict_quality"] = ((df["flags"] == 0) & (df.nimaflags == 0)
                        & (df.img_qual <= 2) & np.isfinite(df.mag_psf)
                        & (df.e_mag_psf > 0))
df["master_used"] = df.use_in_clipped == 1
df["mag_psf_minus_filter_median"] = np.nan
per_filter = {}
for filt, group in df.groupby("filter"):
    good = group.loc[group.strict_quality]
    baseline = float(good.mag_psf.median()) if len(good) else float(group.mag_psf.median())
    df.loc[group.index, "mag_psf_minus_filter_median"] = group.mag_psf - baseline
    per_filter[filt] = {
        "raw_rows": len(group), "strict_quality_rows": len(good),
        "master_used_rows": int(group.master_used.sum()),
        "strict_quality_median_psf_mag": baseline,
        "strict_quality_min_psf_mag": float(good.mag_psf.min()) if len(good) else None,
        "strict_quality_max_psf_mag": float(good.mag_psf.max()) if len(good) else None,
        "strict_quality_median_error_mag": float(good.e_mag_psf.median()) if len(good) else None,
    }
phased_path = ROOT / "skymapper_dr4_phased.csv"
df.to_csv(phased_path, index=False)
columns = {}
for table in ET.fromstring((ROOT / "skymapper_tap_metadata.xml").read_text()).findall(".//table"):
    name = table.findtext("name")
    if name in ("dr4.images", "dr4.photometry", "dr4.master"):
        columns[name] = {c.findtext("name"): {"description": c.findtext("description"),
                           "unit": c.findtext("unit")} for c in table.findall("column")}
summary = {
    "target": {"tic": "273879360", "gaia_dr3": "4692910966756908288",
               "ra_deg": target.ra.deg, "dec_deg": target.dec.deg},
    "association": {"skymapper_dr4_object_id": int(master.object_id.iloc[0]),
        "smss_j": str(master.smss_j.iloc[0]), "ra_deg": float(master.raj2000.iloc[0]),
        "dec_deg": float(master.dej2000.iloc[0]),
        "target_separation_arcsec": float(target.separation(SkyCoord(
            master.raj2000.iloc[0]*u.deg, master.dej2000.iloc[0]*u.deg)).arcsec),
        "gaia_dr3_id1": master.gaia_dr3_id1.iloc[0],
        "gaia_dr3_dist1_arcsec": float(master.gaia_dr3_dist1.iloc[0]),
        "master_flags_psf": int(master.flags_psf.iloc[0]),
        "master_self_sources_within15arcsec": int(master.cnt_self_15.iloc[0]),
        "note": "Exact Gaia crossmatch; resolved target, 65.59 arcsec from foreground TIC273879358."},
    "raw_epoch_rows": len(df), "strict_quality_rows": int(df.strict_quality.sum()),
    "master_used_rows": int(df.master_used.sum()),
    "native_mjd_start_range": [float(df.mjd_start.min()), float(df.mjd_start.max())],
    "utc_midtime_range": [str(midtime.min().iso), str(midtime.max().iso)],
    "per_filter": per_filter,
    "ephemeris_source": str(source_path.relative_to(REPO)),
    "ephemeris": {k: v for k, v in ephem.items() if k in ["period_days", "period_error_days",
        "epoch_btjd", "epoch_error_days", "secondary_phase", "secondary_phase_error_conditional",
        "primary_duration_days", "secondary_duration_days"]},
    "timing": {"native_time": "dr4.images.date AS mjd_start: MJD at start of exposure",
        "scale_assumption": "UTC; scale is not explicitly declared in retained TAP column metadata",
        "conversion": "UTC mid-exposure -> TDB plus astropy barycentric light-travel correction",
        "site": "Siding Spring Observatory: lon149.0661deg, lat−31.2733deg, elevation1165m",
        "astropy_iers_download": False},
    "quality": {"strict_selection": "flags==0, nimaflags==0, img_qual<=2, finite mag_psf, e_mag_psf>0",
        "use_in_clipped": "Retained as provenance only; no variable-star clipping used for phase test",
        "notes": "Image quality3 measurements and flagged measurements are retained and separately shown."},
    "coverage": {
        "primary_in_window_all": int(df.in_primary_window.sum()),
        "secondary_in_window_all": int(df.in_secondary_window.sum()),
        "primary_in_window_strict": int((df.in_primary_window & df.strict_quality).sum()),
        "secondary_in_window_strict": int((df.in_secondary_window & df.strict_quality).sum()),
        "nearest_primary_center_hours_all": float(abs(df.primary_offset_days).min()*24),
        "nearest_secondary_center_hours_all": float(abs(df.secondary_offset_days).min()*24),
        "primary_half_duration_hours": ephem["primary_duration_days"]*12,
        "secondary_half_duration_hours": ephem["secondary_duration_days"]*12,
        "minimum_distance_past_primary_window_minutes_all": float((abs(df.primary_offset_days).min()-ephem["primary_duration_days"]/2)*1440),
        "minimum_distance_past_secondary_window_minutes_all": float((abs(df.secondary_offset_days).min()-ephem["secondary_duration_days"]/2)*1440),
    },
    "assessment": "No epoch lies inside either predicted eclipse window, including all 50 raw measurements. This is a resolved public ground-photometry retrieval and a phase-coverage check, not independent eclipse confirmation or a useful eclipse nondetection.",
    "limitations": ["Sparse clustered observations; no independent period fit or mass constraint.",
        "TESS duration windows are used as the prediction, not ground-data selected windows.",
        "Public metadata does not explicitly state UTC. A several-minute timing-scale/site ambiguity would not bridge the nearest-window gaps.",
        "Per-band magnitudes are AB PSF measurements; filter-specific medians are removed only for plotting."],
    "metadata": columns,
}
summary_path = ROOT / "skymapper_summary.json"
summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n")
colors = dict(zip("uvgriz", ["#6a3d9a", "#54278f", "#2ca25f", "#ef3b2c", "#fb6a4a", "#a6611a"]))
fig, axes = plt.subplots(1, 3, figsize=(13,4.2), gridspec_kw={"width_ratios": [1.5,1,1]})
for filt, group in df.groupby("filter"):
    for ax, x in zip(axes, [group.phase, group.primary_offset_days*24, group.secondary_offset_days*24]):
        good = group.strict_quality.to_numpy()
        ax.errorbar(x[good], group.mag_psf_minus_filter_median[good],
                    yerr=group.e_mag_psf[good], fmt="o", ms=4, color=colors[filt],
                    label=filt if ax is axes[0] else None, alpha=.8, lw=.7)
        ax.scatter(x[~good], group.mag_psf_minus_filter_median[~good], marker="x",
                   s=20, color=colors[filt], alpha=.4)
for center,duration in [(0,ephem["primary_duration_days"]), (ephem["secondary_phase"],ephem["secondary_duration_days"])]:
    half=duration/ephem["period_days"]/2
    axes[0].axvspan(center-half, center+half, color="gray", alpha=.2)
axes[0].axvspan(1-ephem["primary_duration_days"]/ephem["period_days"]/2,1,color="gray",alpha=.2)
for ax,kind in zip(axes[1:], ["primary", "secondary"]):
    ax.axvspan(-ephem[kind+"_duration_days"]*12, ephem[kind+"_duration_days"]*12,color="gray",alpha=.2)
    ax.set_xlim(-8,8); ax.set_title(kind.capitalize()+" window")
    ax.set_xlabel("Hours from predicted eclipse center")
axes[0].set_xlim(0,1); axes[0].set_xlabel("Phase at fixed TESS ephemeris")
axes[0].set_title("Resolved SkyMapper target; 50 raw epochs")
axes[0].legend(title="Filter",ncol=3,fontsize=8)
for ax in axes:
    ax.axhline(0,color="black",lw=.5); ax.invert_yaxis(); ax.grid(alpha=.2)
axes[0].set_ylabel("PSF mag − strict-quality band median")
fig.suptitle("TIC273879360: no SkyMapper exposure samples either eclipse window",fontsize=12)
fig.text(.5,.01,"Filled circles: strict quality; crosses: other retained epochs; shaded: independently predicted eclipse windows",ha="center",fontsize=9)
fig.tight_layout(rect=(0,.035,1,.94))
plot_path = ROOT / "skymapper_phase_coverage.png"
fig.savefig(plot_path,dpi=170)
inputs = [raw_path, ROOT/"skymapper_dr4_master.csv", ROOT/"skymapper_tap_metadata.xml", source_path, Path(__file__)]
for dst in [phased_path, summary_path, plot_path]:
    manifest = {"generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "file": str(dst.relative_to(REPO)), "sha256": hashlib.sha256(dst.read_bytes()).hexdigest(),
        "bytes": dst.stat().st_size,
        "inputs": [{"file": str(p.relative_to(REPO)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs],
        "assessment": summary["assessment"]}
    Path(str(dst)+".manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
print(json.dumps({k: summary[k] for k in ["association", "raw_epoch_rows", "strict_quality_rows", "master_used_rows", "utc_midtime_range", "coverage", "assessment"]},indent=2))
