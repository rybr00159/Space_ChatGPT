"""Read-only, skeptical triage of the actual-flux survey SQLite checkpoint.

This is a prioritizer, not statistical validation: scores are not false-alarm
probabilities. Known periods, host flags, and nearby contamination are kept
separate so additional periods in known systems remain investigable.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from broad_crossmatch import CatalogIndex, numeric

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = Path('/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite')
RSUN_AU = .004650467260962157
RSUN_REARTH = 109.076


def valid(value, minimum, maximum):
    number = numeric(value)
    return number if number is not None and minimum <= number <= maximum else None


def physical_checks(metadata, peak):
    radius = valid(metadata.get('RADIUS'), .05, 30.)
    mass = valid(metadata.get('MASS', metadata.get('STELLAR_MASS')), .03, 30.)
    tmag = valid(metadata.get('TESSMAG'), 0., 25.)
    teff = valid(metadata.get('TEFF'), 1500., 50000.)
    depth = max(numeric(peak.get('depth')) or 0., numeric(peak.get('local_depth')) or 0.)
    period = numeric(peak.get('period'))
    duration = numeric(peak.get('duration_days'))
    warnings = []
    if radius is None:
        warnings.append('stellar_radius_missing_or_invalid')
    if teff is None:
        warnings.append('teff_missing_or_invalid')
    if tmag is None:
        warnings.append('tmag_missing_or_invalid')
    planet_radius = RSUN_REARTH * radius * math.sqrt(depth) if radius and 0 < depth < 1 else None
    maximum_duration = None
    if mass and radius and period and period > 0:
        semimajor_au = (mass * (period / 365.25) ** 2) ** (1 / 3)
        radius_sum = RSUN_AU * radius * (1 + math.sqrt(max(depth, 0)))
        argument = radius_sum / semimajor_au
        if 0 < argument < 1:
            maximum_duration = period / math.pi * math.asin(argument)
            if duration and duration > 1.5 * maximum_duration:
                warnings.append('longer_than_circular_central_duration')
    else:
        warnings.append('duration_density_test_unavailable_no_measured_mass')
    if planet_radius and planet_radius > 22:
        warnings.append('radius_proxy_larger_than_22_rearth')
    if depth > .02:
        warnings.append('depth_above_two_percent')
    # Deep transits can be planetary around small stars; retain physically
    # plausible cases. Unmeasured radius does not supply an M-dwarf exemption.
    deep_small_star = bool(radius and radius <= .8 and planet_radius and planet_radius <= 22 and depth <= .30)
    stellar_size = bool(planet_radius and planet_radius > 22)
    excessive_depth = bool(depth > .02 and not deep_small_star)
    return {'radius_rsun': radius, 'mass_msun': mass, 'teff_kelvin': teff, 'tmag': tmag,
            'radius_proxy_rearth': planet_radius, 'circular_central_duration_hours': maximum_duration * 24 if maximum_duration else None,
            'depth_for_proxy': depth, 'deep_small_star_retained': deep_small_star,
            'physical_eb_warning': stellar_size or excessive_depth,
            'physical_warnings': ';'.join(warnings)}


def classify(record, scan, peak, index, min_snr=8.):
    local = numeric(peak.get('local_snr')) or 0.
    events = int(peak.get('n_events') or 0)
    positive = int(peak.get('n_positive_events') or 0)
    dominance = numeric(peak.get('event_snr_max_fraction'))
    odd_even = numeric(peak.get('odd_even_sigma'))
    reasons = []
    if not peak.get('screening_pass'):
        reasons.append('initial_screen_guard_failed')
    if events < 3:
        reasons.append('fewer_than_three_sampled_events')
    if positive < 3 or positive < .7 * events:
        reasons.append('insufficient_positive_event_support')
    if local < min_snr:
        reasons.append('local_diagonal_snr_below_threshold')
    if dominance is None or dominance > .65:
        reasons.append('one_event_dominates_or_unavailable')
    if (numeric(peak.get('depth')) or 0) <= 0:
        reasons.append('nonpositive_depth')
    if reasons:
        return None, reasons
    match = index.match_signal(record['tic'], peak['period'], record.get('ra'), record.get('dec'))
    physical = physical_checks(record.get('metadata', {}), peak)
    flags = record.get('catalog_flags', {})
    # Recheck host flags against the complete index: earlier worker processes
    # may have started before optional support tables were prepared.
    known = match['known_period_match'] or match['known_harmonic_match']
    eb_host = match['is_eb_host']
    odd_warning = odd_even is not None and odd_even > 3.5
    duty = peak['duration_days'] / peak['period']
    wide_duty_warning = duty > .08
    if eb_host:
        category = 'known_eb_control'
    elif known:
        category = 'known_period_control'
    elif match['possible_nearby_period_contamination']:
        category = 'nearby_known_period_warning'
    elif physical['physical_eb_warning']:
        category = 'stellar_size_or_depth_warning'
    elif odd_warning:
        category = 'odd_even_warning'
    elif wide_duty_warning:
        category = 'wide_duty_warning'
    else:
        category = 'unrecognized_lead'
    positive_fraction = positive / events
    snr_rank = min(local, 25.)
    # Saturated S/N prevents huge unvetted eclipses dominating by amplitude.
    # Repeated support and distributed event power reward usable evidence.
    score = math.log1p(snr_rank) * (1 - .5 * dominance) * positive_fraction
    score *= 1 + .08 * min(events - 3, 10)
    if odd_warning:
        score *= .35
    if physical['physical_eb_warning']:
        score *= .25
    if match['possible_nearby_period_contamination']:
        score *= .45
    if physical['radius_proxy_rearth'] and physical['radius_proxy_rearth'] > 16:
        score *= .8
    if physical['radius_rsun'] is None:
        score *= .9
    # Very broad troughs and deep short-period variations dominate otherwise.
    # Keep these in the table, but prioritize modest-duration transit shapes.
    score *= math.sqrt(min(1., .04 / duty))
    score *= math.sqrt(min(1., .003 / max(physical['depth_for_proxy'], 1.e-12)))
    names = [r['name'] for r in match['known_matches']]
    row = {'tic': record['tic'], 'sector': record.get('metadata', {}).get('SECTOR', 1),
           'category': category, 'priority_score_unvalidated': score,
           'scan_label': scan.get('label'), 'period_days': peak['period'],
           'epoch_btjd': peak.get('t0_btjd'), 'duration_hours': peak.get('duration_days', 0) * 24,
           'depth_ppm': peak['depth'] * 1e6, 'local_depth_ppm': (numeric(peak.get('local_depth')) or 0) * 1e6,
           'snr_diagonal': peak.get('snr'), 'local_snr_diagonal': local,
           'sampled_events': events, 'positive_events': positive,
           'event_snr_max_fraction': dominance, 'odd_even_sigma_diagonal': odd_even,
           'odd_even_warning': odd_warning, 'known_period_match': match['known_period_match'],
           'transit_duty_fraction': duty, 'wide_duty_warning': wide_duty_warning,
           'known_harmonic_match': match['known_harmonic_match'], 'known_rv_period_match': match['known_rv_period_match'],
           'is_confirmed_host': match['is_confirmed_host'], 'is_toi_host': match['is_toi_host'],
           'is_ctoi_host': match['is_ctoi_host'], 'is_diamante_host': match['is_diamante_host'],
           'is_published_host': match['is_published_host'], 'is_eb_host': eb_host,
           'nearby_catalog_object': match['nearby_catalog_object'],
           'possible_nearby_period_contamination': match['possible_nearby_period_contamination'],
           'matched_names': ';'.join(names), 'ra': record.get('ra'), 'dec': record.get('dec'),
           'raw_path': record.get('path'), 'raw_sha256': record.get('sha256'),
           'source_url': record.get('url'), 'method_sha256': record.get('method_sha256'),
           'novelty_status': 'catalog_screen_only_literature_and_pixel_checks_pending',
           'false_alarm_probability': None, **physical}
    return row, []


def run(db_path, output, report, max_per_category=200, min_snr=8.):
    index = CatalogIndex()
    uri = db_path.resolve().as_uri() + '?mode=ro'
    db = sqlite3.connect(uri, uri=True)
    db.execute('PRAGMA query_only=ON')
    db.execute('BEGIN')  # Snapshot the growing WAL checkpoint without mutations.
    status_counts = dict(db.execute('SELECT status,count(*) FROM stars GROUP BY status'))
    considered = 0
    categories = Counter()
    rejections = Counter()
    shortlisted = []
    for tic, status, raw in db.execute('SELECT tic,status,result_json FROM stars WHERE result_json IS NOT NULL'):
        record = json.loads(raw)
        for scan_index, scan in enumerate(record.get('scans', [])):
            for peak_index, peak in enumerate(scan.get('peaks', [])):
                considered += 1
                row, failures = classify(record, scan, peak, index, min_snr)
                if row is None:
                    rejections.update(failures)
                    continue
                row['evidence_id'] = f'stars:{tic}/scans/{scan_index}/peaks/{peak_index}'
                row['source_sqlite'] = str(db_path)
                row['source_status'] = status
                categories[row['category']] += 1
                shortlisted.append(row)
    db.rollback()
    db.close()
    # Deduplicate close aliases from the same star within 1% across masked /
    # unmasked passes, retaining the higher-ranked evidence. Harmonics remain.
    deduplicated = []
    by_tic = {}
    for row in sorted(shortlisted, key=lambda r: -r['priority_score_unvalidated']):
        prior = by_tic.setdefault(row['tic'], [])
        if any(abs(row['period_days'] / r['period_days'] - 1) < .01 for r in prior):
            continue
        prior.append(row)
        deduplicated.append(row)
    category_order = {'unrecognized_lead': 0, 'nearby_known_period_warning': 1, 'odd_even_warning': 2,
                      'wide_duty_warning': 3, 'stellar_size_or_depth_warning': 4,
                      'known_period_control': 5, 'known_eb_control': 6}
    deduplicated.sort(key=lambda r: (category_order[r['category']], -r['priority_score_unvalidated']))
    saved = []
    count_saved = Counter()
    for row in deduplicated:
        if count_saved[row['category']] >= max_per_category:
            continue
        count_saved[row['category']] += 1
        row['rank_within_category'] = count_saved[row['category']]
        saved.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    if saved:
        with output.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(saved[0]))
            writer.writeheader()
            writer.writerows(saved)
    else:
        output.write_text('tic,category,period_days\n')
    stamp = datetime.now(timezone.utc).isoformat()
    counts_lines = '\n'.join(f'- {k}: {v:,}' for k, v in categories.items())
    rows = [r for r in saved if r['category'] == 'unrecognized_lead'][:10]
    table = '| TIC | Period d | Depth ppm | Local S/N | Events | Rp proxy Earth radii | Flags |\n|---|---:|---:|---:|---:|---:|---|\n'
    for r in rows:
        rp = f"{r['radius_proxy_rearth']:.2f}" if r['radius_proxy_rearth'] is not None else 'unknown'
        hints = []
        if r['nearby_catalog_object']:
            hints.append('nearby catalog object')
        if r['is_confirmed_host'] or r['is_toi_host'] or r['is_ctoi_host']:
            hints.append('additional period in reported host')
        if r['radius_rsun'] is None:
            hints.append('unknown stellar radius')
        table += f"| {r['tic']} | {r['period_days']:.6f} | {r['depth_ppm']:.0f} | {r['local_snr_diagonal']:.1f} | {r['sampled_events']} | {rp} | {'; '.join(hints)} |\n"
    text = f"""# Broad-survey prioritization

Snapshot: {stamp}. The SQLite source was opened in read-only/query-only mode
with a consistent read transaction; the running survey was not changed.
This report and `results/broad/shortlist.csv` can be regenerated as the survey
grows. These are screening leads, not credible candidates or discoveries.

Source: `{db_path}`. Actual checkpoint statuses: `{status_counts}`.
Examined {considered:,} stored peaks. Before 1% within-star alias deduplication:

{counts_lines}

Saved at most {max_per_category} rows per category, with exact JSON evidence
references, raw product URLs/paths/hashes, and source method hashes. Category
counts above are before deduplication/capping; CSV ranks are after these steps.

{table}

Required support: initial screening guard pass; >=3 sampled events; >=3
positive-depth events and >=70% positive event fraction; local diagonal S/N
>={min_snr}; no event supplying >65% of positive event squared S/N. Event
sampling/sideband rules come from `fast_screen.event_stats`. Odd/even
differences >3.5 formal diagonal sigma are retained in a warning category.
Durations occupying >8% of a period are retained in a wide-duty warning
category, because broad variable-star troughs dominated the initial rankings.
This heuristic can demote real short-period planets and is not a veto.
This does not account for red noise, chosen-period effects, or global trials.

Crossmatching is recomputed against all prepared support catalogs. Exact or
harmonic matches are recovered controls. EBs are a separate control group;
period-matched nearby sources are contamination warnings. Known hosts at
different periods remain eligible for follow-up. Coarse 1% period matching is
a screening approximation, not ephemeris confirmation.

Rplanet proxy = Rstar * sqrt(max(BLS depth, local depth)), without a dilution
correction or transit fit. Negative sentinel radius/temperature/magnitude
values are invalid, never physical values. Proxy radii >22 Earth radii or depth
>2% are warnings, except deep cases around measured Rstar<=0.8 solar radii
with proxy radius<=22 Earth radii and depth<=30%, which are retained. These
limits prioritize follow-up and do not establish a physical planet boundary.
They can miss diluted EBs or planets around badly characterized stars.

If a measured stellar mass is available, compare duration with the circular
central-transit maximum and flag durations >1.5 times that prediction. QLP
survey metadata currently generally omit mass, so missing-mass cases are
explicitly flagged; no inferred stellar density is silently substituted.
Eccentricity, impact parameter, exposure integration, and stellar parameter
errors prevent this proxy from being a validation test.

The priority score saturates S/N at 25 and rewards recurring positive events
with distributed power, then demotes physical/odd-even/contamination warning
cases. It also reduces the rank of duty fractions >4% and depths >3,000 ppm
without discarding these cases; this addresses amplitude/variability-dominated
follow-up queues. The score and diagonal S/N have no false-alarm-probability meaning.

Next: load each leading raw QLP product; inspect raw/alternative apertures,
refine periods with an independent detrending choice and event-level residuals,
test variability and instrumental/common-mode timing, compare all sectors and
independent SPOC/TGLC data, then perform pixel/nearby-star and literature
checks. A strong folded BLS peak alone is insufficient.

Reproduce: `python src/broad_prioritize.py --db {db_path}`.
Rejection reason counts (nonexclusive): `{dict(rejections)}`.
"""
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text)
    print(json.dumps({'snapshot_utc': stamp, 'statuses': status_counts, 'considered_peaks': considered,
                      'categories': dict(categories), 'saved': dict(count_saved),
                      'top_unrecognized': [{k: r[k] for k in ['tic', 'period_days', 'depth_ppm', 'local_snr_diagonal',
                                                            'sampled_events', 'radius_proxy_rearth', 'is_confirmed_host',
                                                            'nearby_catalog_object']} for r in rows]}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, default=DEFAULT_DB)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/broad/shortlist.csv')
    parser.add_argument('--report', type=Path, default=ROOT / 'reports/broad_prioritization.md')
    parser.add_argument('--max-per-category', type=int, default=200)
    parser.add_argument('--min-snr', type=float, default=8.)
    args = parser.parse_args()
    run(args.db, args.output, args.report, args.max_per_category, args.min_snr)
