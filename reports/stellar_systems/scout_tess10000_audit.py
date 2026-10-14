"""Rerank stellar scouts against both vetted Kostov2025 binary tables.

The ~900,000 unvetted neural detections are deliberately not treated as known
binaries. Existing source/result tables are never changed.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
from scout_catalog_audit import parse

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/stellar_systems'


def main():
    catalogs, provenance = [], []
    for number in [3, 4]:
        path = ROOT / f'data/stellar_systems/scout_novelty/tess10000_table{number}_mrt.txt'
        table = ascii.read(path, format='cds')
        catalogs.append((number, table))
        provenance.append({'table': number, 'rows': len(table), 'path': str(path.relative_to(ROOT)),
                           'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    matches = []
    sb2_path = ROOT / 'data/stellar_systems/scout_novelty/traven2020_full.tsv'
    sb2_rows = parse(sb2_path)
    sb2_coords = SkyCoord([float(row['RAJ2000']) for row in sb2_rows] * u.deg,
                         [float(row['DEJ2000']) for row in sb2_rows] * u.deg)
    provenance.append({'catalog': 'Traven2020_GALAH_SB2', 'rows': len(sb2_rows),
                       'path': str(sb2_path.relative_to(ROOT)), 'sha256': hashlib.sha256(sb2_path.read_bytes()).hexdigest()})
    for filename in ['scout_raw_measurements.csv', 'scout_selection.csv']:
        targets = pd.read_csv(OUT / filename)
        output = []
        for _, target in targets.iterrows():
            coord = SkyCoord(target.ra * u.deg, target.dec * u.deg)
            local = []
            for number, table in catalogs:
                positions = SkyCoord(np.asarray(table['RAdeg'], float) * u.deg,
                                     np.asarray(table['DEdeg'], float) * u.deg)
                separations = coord.separation(positions).arcsec
                identifiers = np.asarray(table['TIC'], int)
                keep = (separations <= 120) | (identifiers == int(target.tic))
                for index in np.flatnonzero(keep):
                    row = table[index]
                    period = float(row['Per'] if number == 3 else row['Per-TESS'])
                    factors = [.25, .333333333333, .5, 1, 2, 3, 4]
                    ratio = min(factors, key=lambda factor: abs(target.period_days / (period * factor) - 1))
                    mismatch = abs(target.period_days / (period * ratio) - 1)
                    match = {'tic': int(target.tic), 'table': number,
                             'catalog_tic': int(identifiers[index]), 'separation_arcsec': float(separations[index]),
                             'host_match': bool(identifiers[index] == int(target.tic) or separations[index] <= 2),
                             'catalog_period_days': period, 'harmonic_factor': ratio,
                             'period_fractional_mismatch': float(mismatch),
                             'period_harmonic_match': bool(mismatch < .005)}
                    local.append(match)
            matches.extend([{**row, 'source_selection': filename} for row in local])
            excluded = any(row['host_match'] or row['period_harmonic_match'] for row in local)
            sb2_sep = coord.separation(sb2_coords).arcsec
            sb2_host = bool(np.any(sb2_sep <= 2))
            # Use retained live cone snapshots when available, never infer
            # that absence of a queried snapshot means an absent catalog row.
            vsx_paths = [ROOT / f'data/stellar_systems/scout_novelty/tic{int(target.tic)}_vsx.tsv',
                         OUT / f'scout_catalogs/vsx_tic{int(target.tic)}.tsv']
            vsx_known = False
            vsx_checks = []
            for vsx_path in vsx_paths:
                if not vsx_path.exists():
                    continue
                for row in parse(vsx_path):
                    try:
                        distance = float(coord.separation(SkyCoord(float(row['RAJ2000']) * u.deg,
                                                                   float(row['DEJ2000']) * u.deg)).arcsec)
                        old = float(row['Period']) if row.get('Period') else None
                    except (ValueError, KeyError):
                        continue
                    is_eb = row.get('Type', '').startswith('E')
                    harmonic = old is not None and any(abs(target.period_days / (old * factor) - 1) < .005
                                                       for factor in [.5, 1, 2, 3, 4])
                    vsx_known |= is_eb and (distance <= 2 or harmonic)
                    vsx_checks.append({'name': row.get('Name'), 'type': row.get('Type'),
                                       'period_days': old, 'separation_arcsec': distance,
                                       'harmonic_match': harmonic})
            item = target.to_dict()
            item.update(known_tess10000_host=any(row['host_match'] for row in local),
                        known_tess10000_nearby_period_match=any(not row['host_match'] and row['period_harmonic_match'] for row in local),
                        known_tess10000_nearby_any=any(not row['host_match'] for row in local),
                        exclude_as_new_binary_tess10000=excluded,
                        known_galah_sb2_host=sb2_host, known_vsx_eb_cone=vsx_known,
                        vsx_cone_checked=any(path.exists() for path in vsx_paths),
                        vsx_matches_json=json.dumps(vsx_checks, separators=(',', ':')),
                        tess10000_matches_json=json.dumps(local, separators=(',', ':')))
            output.append(item)
        result = pd.DataFrame(output)
        name = 'scout_current_rerank.csv' if filename == 'scout_raw_measurements.csv' else 'scout_current_selection.csv'
        result.to_csv(OUT / name, index=False)
        print(filename, 'count', len(result), 'excluded_known', int(result.exclude_as_new_binary_tess10000.sum()),
              'remaining', int((~result.exclude_as_new_binary_tess10000).sum()))
        if filename == 'scout_raw_measurements.csv':
            new = result[result.robust_repeat_scout_pass & ~result.exclude_as_new_binary_tess10000
                         & ~result.known_galah_sb2_host & ~result.known_vsx_eb_cone].copy()
            # Prior whole-table host or harmonic-matched neighbor exclusions
            # remain necessary in addition to the newly recognized 2025 catalog.
            earlier = json.loads((OUT / 'scout_expanded_catalog_matches.json').read_text())
            prior_excluded = set()
            for match in earlier:
                target_period = targets.loc[targets.tic == match['tic'], 'period_days'].iloc[0]
                old = match['period_days']
                harmonic = old is not None and any(abs(target_period / (old * factor) - 1) < .005
                                                   for factor in [.5, 1, 2, 3, 4])
                if match['exact_identifier_match'] or match['separation_arcsec'] <= 2 or harmonic:
                    prior_excluded.add(match['tic'])
            new = new[~new.tic.isin(prior_excluded)].sort_values('followup_score', ascending=False)
            # This row is a known visual double; a new subcomponent would need
            # a separate hierarchical-system experiment, not this EB ranking.
            new = new[~new.tic.isin([403150111, 394136010])]
            new.to_csv(OUT / 'scout_fresh_followup.csv', index=False)
            # Operational next queue avoids all first-round targets already
            # receiving individual novelty/source audits, even when their
            # remaining evidence is not yet sufficient to classify them.
            already_individually_reviewed = [294272014, 355544861, 38825979, 141866704, 141411137,
                                             281581467, 388181335, 394136010, 197764088,
                                             260503860, 115188534, 300652458, 52285789, 115152243,
                                             403150111, 294206969, 89018908]
            new[~new.tic.isin(already_individually_reviewed)].head(64).to_csv(OUT / 'scout_next_pilot.csv', index=False)
            print('Fresh raw repeat leads after expanded tables:', len(new))
            print(new.head(40)[['tic', 'ra', 'dec', 'period_days', 'raw_depth_ppm', 'half_phase_depth_ppm',
                                'half_phase_snr', 'parity_conditional_sigma']].to_csv(index=False))
    payload = {'sources': provenance, 'matches': matches,
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'limitations': 'Finite prior tables; nonmatch is not novelty. Coordinates unpropagated, nearby period match only contamination evidence; source pixels and other catalogs/aliases/literature remain required.'}
    (OUT / 'scout_tess10000_matches.json').write_text(json.dumps(payload, indent=2) + '\n')


if __name__ == '__main__':
    main()
