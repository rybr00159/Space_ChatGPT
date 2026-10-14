"""Verify/restore the five original QLP observations for the targeted audit.

Default checks existing local originals; --download retrieves missing products
into TESS_DATA_DIR/tic71316387/observations. Immutable expected hashes are never
rewritten. Changed archive bytes remain preserved as .part and cause failure.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys

import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from archive import verify_fits_checksums


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    cache = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'tic71316387/observations'
    for row in json.loads((ROOT/'candidates/tic71316387/observations.json').read_text()):
        if row['status'] != 'retained':
            continue
        original = Path(row['path'])
        path = original if original.exists() else cache / original.name
        if not path.exists():
            if not args.download:
                print(f"Sector {row['sector']}: absent; use --download")
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            response = requests.get(row['url'], timeout=(20, 120))
            response.raise_for_status()
            temporary = path.with_suffix('.part')
            temporary.write_bytes(response.content)
            if hashlib.sha256(response.content).hexdigest() != row['sha256']:
                raise ValueError(f'Archive bytes changed; preserved {temporary}; original expected SHA retained')
            temporary.replace(path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError(f'Original source SHA mismatch: {path}')
        with fits.open(path) as hd:
            hd.verify('exception')
            verify_fits_checksums(path, hd)
        print(f"Sector {row['sector']}: verified {path.stat().st_size} bytes")


if __name__ == '__main__':
    main()
