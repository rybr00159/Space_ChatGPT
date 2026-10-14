"""Review exact uncapped signal rows skipped by the one-target batch selection.

Reproduce: python candidates/broad_s1/extra_signals/review.py
Preserves the 357-star batch checkpoint; calls the frozen v3 vetting functions.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import sqlite3
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from broad_vet import DATABASE, OUT as DOSSIERS, method_identity, vet, _ranked_summary, atomic_json

OUT = Path(__file__).resolve().parent


if __name__ == '__main__':
    queue = ROOT / 'results/broad/shortlist_uncapped.csv'
    raw = queue.read_bytes()
    rows = [r for r in csv.DictReader(raw.decode().splitlines()) if r['category'] == 'unrecognized_lead']
    primary_batch = json.loads((DOSSIERS / 'batch_checkpoint.json').read_text())
    used = {r['source_evidence_id'] for r in primary_batch['results']}
    extra = [r for r in rows if r['evidence_id'] not in used]
    identity = method_identity()
    report = dict(started_utc=datetime.now(timezone.utc).isoformat(), method=identity,
                  queue_sha256=hashlib.sha256(raw).hexdigest(), main_batch_unique_targets=primary_batch['completed'],
                  main_batch_evidence_ids=list(sorted(used)), total_unrecognized_signal_rows=len(rows),
                  planned_additional_signal_rows=len(extra), results=[], errors=[],
                  limitations='Additional peaks on an already examined star are signals, not additional unique stars.')
    db = sqlite3.connect(f'file:{DATABASE}?mode=ro', uri=True)
    for row in extra:
        tic = int(row['tic'])
        try:
            record = json.loads(db.execute('SELECT result_json FROM stars WHERE tic=?', (tic,)).fetchone()[0])
            index = row['evidence_id'].split('/')
            scan = record['scans'][int(index[2])]
            peak = scan['peaks'][int(index[4])]
            if record['sha256'] != row['raw_sha256'] or abs(peak['period'] - float(row['period_days'])) > 1e-9:
                raise ValueError('Exact queue source or ephemeris mismatch')
            destination = DOSSIERS / f'tic{tic}_p{peak["period"]:.6f}' / 'vetting.json'
            cached = json.loads(destination.read_text()) if destination.exists() else None
            valid = (cached is not None and cached.get('vetting_method') == identity and
                     cached.get('source', {}).get('sha256') == record['sha256'] and
                     cached.get('survey_label') == scan['label'] and
                     abs(cached.get('input_peak', {}).get('period', -1) - peak['period']) < 1e-10 and
                     cached.get('raw_variability_guard') is not None and
                     all('additional_phase_secondary' in cached['column_tests'].get(k, {}) for k in
                         ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG')))
            if valid:
                result = cached
                print('EXTRA_CACHED', row['evidence_id'], flush=True)
            else:
                result = vet(record, peak, scan['label'])
            report['results'].append(_ranked_summary(result, row['priority_score_unvalidated'], row['evidence_id']))
        except Exception as err:
            report['errors'].append(dict(evidence_id=row['evidence_id'], error=type(err).__name__ + ': ' + str(err)))
            print('EXTRA_ERROR', row['evidence_id'], type(err).__name__, str(err), flush=True)
        report['updated_utc'] = datetime.now(timezone.utc).isoformat()
        report['completed_additional_signal_rows'] = len(report['results'])
        atomic_json(OUT / 'summary.json', report)
    db.close()
    print('EXTRA_COMPLETE', len(report['results']), 'errors', len(report['errors']), flush=True)
