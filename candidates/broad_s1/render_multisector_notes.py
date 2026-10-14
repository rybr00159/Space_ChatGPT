"""Correct target-specific prose without altering numeric scientific methods.

The frozen v3 renderer inherited a CD-79 151 sentence for every target. That
neighbor applies only to TIC370247818. Record text hashes; never touch arrays.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SOURCE_SHA = '23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333'
OLD = ('The bright star CD-79 151 at about 86.6 arcsec remains a possible contaminant pending pixel analysis. ')


if __name__ == '__main__':
    if hashlib.sha256((ROOT / 'src/broad_vet.py').read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError('Frozen numeric/source version unexpectedly changed')
    changes = []
    for directory in OUT.glob('tic*_multisector'):
        path = directory / 'README.md'
        if not path.exists():
            continue
        original = path.read_text()
        if OLD not in original or directory.name == 'tic370247818_multisector':
            continue
        replacement = ('Source attribution requires this target’s own pixel and nearby-star evidence. '
                       'A contaminant identified for a different target is not evidence for this source. ')
        revised = original.replace(OLD, replacement)
        numeric = directory / 'multisector.json'
        changes.append(dict(target=directory.name, path=str(path.relative_to(ROOT)),
                            previous_text_sha256=hashlib.sha256(original.encode()).hexdigest(),
                            corrected_text_sha256=hashlib.sha256(revised.encode()).hexdigest(),
                            preserved_numeric_sha256=hashlib.sha256(numeric.read_bytes()).hexdigest() if numeric.exists() else None))
        path.write_text(revised)
    report_path = OUT / 'text_corrections.json'
    previous = json.loads(report_path.read_text()) if report_path.exists() else {'changes': []}
    previous['changes'].extend(changes)
    previous.update(checked_utc=datetime.now(timezone.utc).isoformat(), preserved_numeric_source_sha256=SOURCE_SHA,
                    reason='Correct inherited target-specific contaminant prose; no numeric methods or arrays changed.')
    report_path.write_text(json.dumps(previous, indent=2))
    print('CORRECTED_READMES', len(changes))
