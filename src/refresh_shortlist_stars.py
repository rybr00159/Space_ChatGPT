"""Batch current public TIC properties for eligible/vetted survey targets.

QLP sector-1 headers can misclassify giants as small stars. Retrieve current
TIC rows before physical triage. Original survey records and source tables
remain unchanged. Public response bodies, requests, and hashes are preserved.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from archive import invoke
from broad_crossmatch import numeric, tic_number

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = Path('/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite')
SERVICE = 'Mast.Catalogs.Filtered.Tic'
ENDPOINT = 'https://mast.stsci.edu/api/v0/invoke'


def strict_eligible(peak, min_snr=8.):
    events = int(peak.get('n_events') or 0)
    positive = int(peak.get('n_positive_events') or 0)
    dominance = numeric(peak.get('event_snr_max_fraction'))
    return bool(peak.get('screening_pass') and events >= 3 and positive >= 3
                and positive >= .7 * events and (numeric(peak.get('local_snr')) or 0) >= min_snr
                and dominance is not None and dominance <= .65
                and (numeric(peak.get('depth')) or 0) > 0)


def select_ids(inputs, db_path, include_vetted, min_snr):
    ids, origins = set(), []
    if inputs:
        for source in inputs:
            with source.open() as stream:
                rows = list(csv.DictReader(stream))
            selected = {tic_number(r.get('tic')) for r in rows} - {None}
            ids |= selected
            origins.append({'source': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                            'input_rows': len(rows), 'distinct_tics': len(selected)})
    else:
        db = sqlite3.connect(db_path.resolve().as_uri() + '?mode=ro', uri=True)
        db.execute('PRAGMA query_only=ON'); db.execute('BEGIN')
        counts = dict(db.execute('SELECT status,count(*) FROM stars GROUP BY status'))
        method_row = db.execute("SELECT value FROM survey_meta WHERE key='method_sha256'").fetchone()
        for tic, raw in db.execute('SELECT tic,result_json FROM stars WHERE result_json IS NOT NULL'):
            record = json.loads(raw)
            if any(strict_eligible(p, min_snr) for scan in record.get('scans', []) for p in scan.get('peaks', [])):
                ids.add(int(tic))
        db.rollback(); db.close()
        origins.append({'source': str(db_path), 'query_mode': 'read-only, consistent transaction',
                        'statuses_at_selection': counts, 'rule': 'any strictly eligible peak; all warning/control categories retained',
                        'minimum_local_snr_diagonal': min_snr, 'distinct_tics': len(ids),
                        'survey_method_sha256': method_row[0] if method_row else None})
    if include_vetted:
        previous = set()
        for path in (ROOT / 'candidates/broad_s1').glob('tic*/vetting.json'):
            record = json.loads(path.read_text())
            tic = tic_number(record.get('tic'))
            if tic:
                previous.add(tic)
        extra = previous - ids
        ids |= previous
        origins.append({'source': 'candidates/broad_s1/tic*/vetting.json', 'distinct_tics': len(previous),
                        'additional_historical_tics': len(extra)})
    return sorted(ids), origins


def fetch_batch(ids, directory, batch_number):
    params = {'columns': '*', 'filters': [{'paramName': 'ID', 'values': [str(tic) for tic in ids]}]}
    request_hash = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()
    path = directory / f'batch{batch_number:04}_{request_hash[:16]}.json'
    manifest_path = path.with_suffix('.json.manifest.json')
    if path.exists():
        metadata = json.loads(manifest_path.read_text())
        if metadata.get('params') != params or metadata.get('service') != SERVICE:
            raise ValueError(f'Cached request mismatch: {path}')
        if hashlib.sha256(path.read_bytes()).hexdigest() != metadata['sha256']:
            raise ValueError(f'Cached response checksum mismatch: {path}; preserve and investigate')
        response = json.loads(path.read_text())
    else:
        response = invoke(SERVICE, params)
        if response.get('status') != 'COMPLETE':
            raise ValueError('MAST catalog query is incomplete')
        body = (json.dumps(response, indent=2) + '\n').encode()
        # A missing body with an existing manifest must reproduce that snapshot.
        prior = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
        digest = hashlib.sha256(body).hexdigest()
        if prior and prior.get('sha256') != digest:
            raise ValueError(f'Remote response changed since recorded snapshot {path}; use a new snapshot directory')
        path.write_bytes(body)
        metadata = {'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'endpoint': ENDPOINT,
                    'service': SERVICE, 'params': params, 'request_sha256': request_hash,
                    'request_format': 'json', 'request_page_size': 10000, 'request_page': 1,
                    'sha256': digest, 'bytes': len(body), 'row_count': len(response.get('data', [])),
                    'requested_tics': ids, 'notes': 'Current public TIC response; row version is not retrieval date. Gaia field is a DR2 identifier.'}
        if not prior:
            manifest_path.write_text(json.dumps(metadata, indent=2) + '\n')
    rows = response.get('data', [])
    returned = [tic_number(r.get('ID')) for r in rows]
    if any(tic not in ids for tic in returned) or len(set(returned)) != len(returned):
        raise ValueError(f'Unexpected/duplicated IDs in exact-TIC query: {path}')
    return rows, path, metadata


def measured(value, minimum=None, maximum=None):
    number = numeric(value)
    if number is None or (minimum is not None and number < minimum) or (maximum is not None and number > maximum):
        return None
    return number


def normalize(tic, row, path, metadata):
    row = row or {}
    radius = measured(row.get('rad'), .01, 10000)
    mass = measured(row.get('mass'), .005, 10000)
    teff = measured(row.get('Teff'), 100, 1.e6)
    flags = []
    if not row:
        flags.append('missing_exact_tic_response')
    if radius is None:
        flags.append('radius_missing_or_invalid')
    if mass is None:
        flags.append('mass_missing_or_invalid')
    if teff is None:
        flags.append('temperature_missing_or_invalid')
    if str(row.get('lumclass', '')).upper() == 'GIANT' or (radius is not None and radius > 3):
        flags.append('evolved_or_giant_stellar_context')
    disposition = str(row.get('disposition') or '')
    if disposition and disposition.upper() not in ('', 'NONE'):
        flags.append('tic_disposition_requires_review')
    return {'tic': tic, 'query_found': bool(row), 'current_tic_record_version': row.get('version'),
            'ra': row.get('ra'), 'dec': row.get('dec'), 'stellar_radius_rsun': radius,
            'stellar_radius_error_rsun': measured(row.get('e_rad'), 0),
            'stellar_mass_msun': mass, 'stellar_mass_error_msun': measured(row.get('e_mass'), 0),
            'current_teff_kelvin': teff, 'current_teff_error_kelvin': measured(row.get('e_Teff'), 0),
            'current_logg': measured(row.get('logg'), -2, 10), 'current_tmag': measured(row.get('Tmag'), -10, 30),
            'tic_lumclass': row.get('lumclass'), 'tic_object_type': row.get('objType'),
            'tic_disposition': row.get('disposition'), 'duplicate_id': row.get('duplicate_id'),
            'tic_contamination_ratio': measured(row.get('contratio'), 0),
            'gaia_dr2': row.get('GAIA'), 'tic_gaia_quality_flag': row.get('gaiaqflag'),
            'radius_source_flag': row.get('raddflag'), 'stellar_validity_flags': ';'.join(flags),
            'source_response': str(path), 'source_sha256': metadata['sha256'],
            'catalog_retrieved_utc': metadata.get('retrieved_utc')}


def run(args):
    ids, origins = select_ids(args.input, args.db, args.include_vetted, args.min_snr)
    directory = args.snapshot_dir or ROOT / 'data/broad_catalogs' / ('current_stars_' + datetime.now(timezone.utc).strftime('%Y%m%d'))
    directory.mkdir(parents=True, exist_ok=True)
    output = args.output or directory / 'current_stellar.csv'
    # Derived refresh outputs may evolve; retain prior snapshots before changes.
    selection_hash = hashlib.sha256(json.dumps(ids).encode()).hexdigest()
    batch_rows, found = [], 0
    print(f'Current TIC selection: {len(ids)} distinct IDs; batches of {args.batch_size}', flush=True)
    for number, start in enumerate(range(0, len(ids), args.batch_size)):
        chunk = ids[start:start + args.batch_size]
        rows, path, metadata = fetch_batch(chunk, directory, number)
        mapping = {int(row['ID']): row for row in rows}
        found += len(mapping)
        batch_rows.extend(normalize(tic, mapping.get(tic), path, metadata) for tic in chunk)
        print(f'Batch {number + 1}: {len(rows)}/{len(chunk)} rows, total {start + len(chunk)}/{len(ids)}', flush=True)
    if not batch_rows:
        raise ValueError('No selected targets; no current table written')
    import io
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fieldnames=list(batch_rows[0])); writer.writeheader(); writer.writerows(batch_rows)
    body = buffer.getvalue().encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    prior_manifest = output.with_suffix(output.suffix + '.manifest.json')
    # Preserve provenance changes even if the normalized table bytes coincide.
    unchanged_prior_manifest = prior_manifest.read_bytes() if prior_manifest.exists() and output.exists() and output.read_bytes() == body else None
    if output.exists() and output.read_bytes() != body:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        output.rename(output.with_name(output.name + '.' + stamp + '.previous'))
        prior = output.with_suffix(output.suffix + '.manifest.json')
        if prior.exists():
            prior.rename(prior.with_name(prior.name + '.' + stamp + '.previous'))
    output.write_bytes(body)
    metadata = {'created_utc': datetime.now(timezone.utc).isoformat(), 'selection_origins': origins,
                'processing_method_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'invoke_helper_sha256': hashlib.sha256((ROOT / 'src/archive.py').read_bytes()).hexdigest(),
                'selection_sha256': selection_hash, 'requested_tics': ids, 'requested_count': len(ids),
                'found_count': found, 'missing_count': len(ids) - found, 'sha256': hashlib.sha256(body).hexdigest(),
                'bytes': len(body), 'columns': list(batch_rows[0]),
                'scope': 'All strict-eligible stored peaks plus optional previously vetted stars; no all-160k catalog-completeness claim.',
                'limitations': 'TIC estimates can be absent or inaccurate. Evolved stars, missing radii, duplicate/disposition flags, and crowded targets require individual Gaia/ExoFOP checks. GAIA field is DR2; DR3 IDs can differ. Original QLP and survey records are unchanged.'}
    metadata_body = (json.dumps(metadata, indent=2) + '\n').encode()
    if unchanged_prior_manifest is not None and unchanged_prior_manifest != metadata_body:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        prior_manifest.rename(prior_manifest.with_name(prior_manifest.name + '.' + stamp + '.previous'))
    prior_manifest.write_bytes(metadata_body)
    print(json.dumps({'output': str(output), 'requested': len(ids), 'found': found,
                      'evolved_context': sum('evolved_or_giant' in r['stellar_validity_flags'] for r in batch_rows),
                      'missing_radius': sum(r['stellar_radius_rsun'] is None for r in batch_rows)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, nargs='+', help='CSV(s) with tic columns; default reads all strictly eligible database targets')
    parser.add_argument('--db', type=Path, default=DEFAULT_DB)
    parser.add_argument('--include-vetted', action='store_true', help='Also include every historical detailed-vetting TIC')
    parser.add_argument('--min-snr', type=float, default=8.)
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--snapshot-dir', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 1 <= args.batch_size <= 512:
        parser.error('Batch size must be between 1 and 512')
    run(args)
