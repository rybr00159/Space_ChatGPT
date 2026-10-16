"""Apply the existing fresh catalog index to empirical-error recovery seeds."""
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
import pandas as pd
from multisector_recognition import FreshCatalogIndex, rank_unknown_signals

OUT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    with gzip.open(OUT / 'all_scans.jsonl.gz', 'rt') as stream:
        records = {(r['tic'], r['sector']):r for r in map(json.loads, stream)}
    frame = pd.read_csv(OUT / 'candidate_seeds.csv')
    frame['ra'] = [float(records[(int(r.tic), int(r.sector))]['ra']) for r in frame.itertuples()]
    frame['dec'] = [float(records[(int(r.tic), int(r.sector))]['dec']) for r in frame.itertuples()]
    frame['label'] = 'unmasked_empirical_error'
    frame['recovery_method_sha256'] = sha(ROOT / 'src/multisector_missing_errors.py')
    frame['screen_method_sha256'] = sha(ROOT / 'src/multisector_screen.py')
    index = FreshCatalogIndex()
    annotated = index.annotate_dataframe(frame)
    annotated.to_csv(OUT / 'candidate_seeds_annotated.csv', index=False)
    ranked = rank_unknown_signals(annotated, limit=15)
    ranked.to_csv(OUT / 'shallow_followup_queue.csv', index=False)
    meta = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=sha(OUT / 'candidate_seeds.csv'),
        all_scans_sha256=sha(OUT / 'all_scans.jsonl.gz'),
        helper_sha256=sha(Path(__file__)),
        recognition_method_sha256=sha(ROOT / 'src/multisector_recognition.py'),
        base_recognition_method_sha256=sha(ROOT / 'src/broad_crossmatch.py'),
        output_sha256=sha(OUT / 'candidate_seeds_annotated.csv'),
        shallow_queue_sha256=sha(OUT / 'shallow_followup_queue.csv'),
        catalogs=index.metadata(), seed_rows=len(annotated),
        seed_tics=int(annotated.tic.nunique()),
        shallow_queue_tics=int(ranked.tic.nunique()),
        scope='Finite prior-report checks; catalogue absence is not novelty. Shallow queue uses the existing <5000ppm operational ranking and does not exhaust interesting deeper signals such as TIC182517223.')
    (OUT / 'recognition_provenance.json').write_text(json.dumps(meta, indent=2) + '\n')
    print('ANNOTATED_RECOVERY', len(annotated), 'seed rows;', ranked.tic.nunique(), 'shallow queue TICs')
    print(ranked[['tic','sector','period_days','local_depth_ppm','local_snr','minimum_width_local_snr']].to_string(index=False))


if __name__ == '__main__':
    main()
