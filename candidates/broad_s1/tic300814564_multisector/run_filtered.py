"""Select usable follow-up products, then run the unchanged frozen v3 routine.

Sector 37's supplied flux/error arrays yield zero valid samples. Its original
download manifest and FITS remain unchanged. This wrapper supplies a filtered
input manifest through the module's repository-input root, retaining all maths.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
import broad_vet as bv

OUT = Path(__file__).resolve().parent
EXPECTED = '23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333'


if __name__ == '__main__':
    if hashlib.sha256((ROOT / 'src/broad_vet.py').read_bytes()).hexdigest() != EXPECTED:
        raise ValueError('Frozen numerical source changed')
    source = ROOT / 'data/broad_followup_queries/tic300814564_qlp.json'
    rows = json.loads(source.read_text())
    selected, exclusions = [], []
    for row in rows:
        if row['status'] != 'retained':
            selected.append(row)
            continue
        data, _ = bv.read_source(row)
        good = data['QUALITY'] == 0
        flat = bv.preprocess_curve(data['TIME'][good], data['KSPSAP_FLUX'][good], data['KSPSAP_FLUX_ERR'][good],
                                   config=bv.FastScreenConfig(width_days=2.))
        if len(flat['t']) < 40:
            exclusions.append(dict(row, reason='Fewer than40 valid preprocessed cadences; product retained.',
                                   n_clean=flat['n_clean'], n_preprocessed=len(flat['t'])))
            selected.append(dict(row, status='retained_excluded_from_this_analysis'))
        else:
            selected.append(row)
    selection = dict(checked_utc=datetime.now(timezone.utc).isoformat(),
                     original_manifest=str(source.relative_to(ROOT)),
                     original_manifest_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                     selected_rows=selected, exclusions=exclusions,
                     frozen_numeric_source_sha256=EXPECTED,
                     wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'input_selection.json').write_text(json.dumps(selection, indent=2))
    cache = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data'))
    input_tree = cache / 'broad/followup_input_overrides/tic300814564'
    folder = input_tree / 'data/broad_followup_queries'
    folder.mkdir(parents=True, exist_ok=True)
    src_link = input_tree / 'src'
    if not src_link.exists():
        src_link.symlink_to(ROOT / 'src', target_is_directory=True)
    (folder / 'tic300814564_qlp.json').write_text(json.dumps(selected, indent=2))
    # ROOT controls input-manifest lookup only here. OUT/DATABASE stay unchanged;
    # source symlink retains identical frozen code hashes and numerical functions.
    bv.ROOT = input_tree
    result = bv.vet_multisector(300814564, selected_period=1.8611071)
    result['followup_input_selection'] = selection
    bv.atomic_json(OUT / 'multisector.json', result)
