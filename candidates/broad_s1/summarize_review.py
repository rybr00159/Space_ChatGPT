"""Combine exact 389-signal frozen-v3 review checkpoints without changing them."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FROZEN_SHA = '23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333'


def main():
    queue_path = ROOT / 'results/broad/shortlist_uncapped.csv'
    rows = list(csv.DictReader(queue_path.open()))
    eligible = {r['evidence_id']: r for r in rows if r['category'] == 'unrecognized_lead'}
    main_path, extra_path = HERE / 'batch_checkpoint.json', HERE / 'extra_signals/summary.json'
    main, extra = json.loads(main_path.read_text()), json.loads(extra_path.read_text())
    assert main['errors'] == extra['errors'] == []
    assert main['method']['source_sha256'] == extra['method']['source_sha256'] == FROZEN_SHA
    results = main['results'] + extra['results']
    keys = [r['source_evidence_id'] for r in results]
    assert len(set(keys)) == len(keys) == len(eligible) == 389
    assert set(keys) == set(eligible)
    annotations = {
        238202207: '17-sector recurrence; unequal doubled-period binary fits as well, alternating depths in16/17 sectors and Gaia unresolved-source warnings favor binary; no accepted planet.',
        269994642: 'Independent sector1/92 pixels support target-neighborhood eclipses; doubled-period binary fits as well, unequal depths and GaiaRUWE2.08 favor binary; source27 prescan excluded; no accepted planet.',
        300814564: '23 usable sectors; raw half-phase secondary 2179.7+/-62.1 ppm favors stellar eclipse. Sector37 data retained but unusable.',
        273879358: 'Four-sector recurrence; raw phase0.6937 secondary 2006+/-220 ppm and Gaia RV variability favor stellar binary/source alternative.',
        270312735: 'Six-sector recurrence; raw phase0.5126 secondary 5584+/-619 ppm favors stellar binary/source alternative.',
        431481555: 'Previously recognized pulsator; later-sector recurrence weak, exploratory outside-primary mode subtraction weakens sector-1 pseudo-transit.',
        139148478: 'Weak later-sector support; shared long variability with bright rotational neighbor. Neighbor does not robustly reproduce 0.70d target dips; source unresolved.',
        260820871: 'Three-sector shallow eclipse recurrence with strong alternating depths; independent S28 pixels trace both eclipse phases to known faint neighboring EB, not hot target; catalog period is approximately harmonic, not an asserted exact match.',
        410255802: 'Three weak events; current stellar properties missing and current TIC disposition requires review.',
        79492853: 'Three weak events; current stellar properties missing and current TIC disposition requires review.',
        419589172: 'Three broad events;9.6h fitted duration and current contamination ratio8.70 require source/variability tests.',
        290028406: 'Sector28 independent signal is weaker and joint raw support only3.24; V-shaped source/aperture warning signal has tentative half-phase dip; unresolved lower-priority follow-up.',
        140803419: 'Secondary coverage is incomplete, including a quality-gap threshold boundary; actual repeated4.58d dips and strong raw modulation merit stellar/source interpretation.',
        139199396: 'Similar0.7012d signal and highly correlated processed flux with152819299 on same detector nearly2deg away favor shared optical/instrumental contamination.',
        152819299: 'Similar0.7012d signal and highly correlated processed flux with139199396 on same detector nearly2deg away favor shared optical/instrumental contamination.',
    }
    joined, verified = [], {}
    for r in results:
        q = eligible[r['source_evidence_id']]
        dossier_path = ROOT / r['evidence_file']
        d = json.loads(dossier_path.read_text())
        assert d['vetting_method']['source_sha256'] == FROZEN_SHA
        assert int(d['tic']) == int(r['tic'])
        source = d['source']
        raw_path = Path(source['path'])
        if str(raw_path) not in verified:
            assert hashlib.sha256(raw_path.read_bytes()).hexdigest() == source['sha256'] == q['raw_sha256']
            verified[str(raw_path)] = source['sha256']
        record = dict(r)
        for field in ['radius_rsun', 'mass_msun', 'teff_kelvin', 'current_tic_contamination_ratio',
                      'stellar_catalog_flags', 'stellar_properties_source', 'stellar_source_response',
                      'stellar_source_sha256', 'ra', 'dec', 'source_url', 'raw_sha256']:
            record['current_' + field if field in ['radius_rsun', 'mass_msun', 'teff_kelvin'] else field] = q[field]
        record['curated_followup_assessment'] = annotations.get(int(r['tic']),
            'Automated warnings require interpretation; this first-sector screen does not establish a planetary or stellar classification.')
        joined.append(record)
    joined.sort(key=lambda r: (-r['ranking_score_unvalidated'], r['tic'], r['source_evidence_id']))
    out_path = HERE / 'all_signal_evidence.csv'
    with out_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(joined[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(joined)
    clean = [r for r in results if not r['warnings']]
    warning_counts = Counter(v for r in results for v in r['warnings'].split(';') if v)
    paths = [queue_path, main_path, extra_path, ROOT / 'src/broad_vet.py', Path(__file__), out_path]
    summary = dict(checked_utc=datetime.now(timezone.utc).isoformat(),
                   exact_signal_rows=len(results), unique_stars=len({r['tic'] for r in results}),
                   main_batch_unique_stars=len(main['results']), alternate_signal_rows=len(extra['results']),
                   errors=0, verified_unique_raw_fits=len(verified),
                   automatic_warning_free_signals=len(clean),
                   automatic_warning_free_stars=sorted({r['tic'] for r in clean}),
                   warning_counts=dict(warning_counts),
                   accepted_new_planet_candidates=0,
                   frozen_v3_sha256=FROZEN_SHA,
                   input_output_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                   limitations=['Automated warning categories overlap and are not false-positive classifications.',
                                'All counts refer to the frozen current-property strict uncapped queue, not the entire159505-star search or older historical dossiers.',
                                'Eight warning-free first-sector signals do not mean eight planet candidates. Multi-sector, source, physical and prior-recognition checks are separate.',
                                'Diagonal signal scores and curated ranks do not provide survey-wide false-alarm probabilities.'])
    (HERE / 'all_signal_review_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ['exact_signal_rows', 'unique_stars', 'verified_unique_raw_fits', 'automatic_warning_free_signals', 'accepted_new_planet_candidates']}))


if __name__ == '__main__':
    main()
