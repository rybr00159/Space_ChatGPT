"""Verify/restore immutable QLP light curves in the stellar-phase manifests.

Default verifies retained originals. --download fetches only missing sources to
their recorded workspace paths, verifies original SHA256, and retains changed
archive bytes as .part rather than rewriting the scientific expectation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import urllib.request
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from archive import verify_fits_checksums


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download',action='store_true')
    args=parser.parse_args()
    manifests=[ROOT/'data/stellar_systems/hidden_sources/analysis_manifest.json',
               ROOT/'candidates/stellar_systems/tic197764088/observations.json']
    records={}
    for manifest in manifests:
        data=json.loads(manifest.read_text())
        for row in data.get('qlp_sources',[])if isinstance(data,dict)else data:
            if row.get('status','retained')!='retained':continue
            records[row['path']]=row
    count=0;missing=0
    for row in records.values():
        p=Path(row['path'])
        if not p.exists():
            if not args.download:
                missing+=1;continue
            p.parent.mkdir(parents=True,exist_ok=True)
            temporary=p.with_suffix('.part')
            with urllib.request.urlopen(row['url'],timeout=120)as response,temporary.open('wb')as stream:
                while chunk:=response.read(1024**2):stream.write(chunk)
            if hashlib.sha256(temporary.read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError(f'Archive changed; preserved {temporary}; original expectation retained')
            temporary.replace(p)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError(f'Raw SHA256 mismatch: {p}')
        with fits.open(p)as hd:
            hd.verify('exception');verify_fits_checksums(p,hd)
        count+=1
    print(f'Verified {count} original QLP files; {missing} absent (use --download)')


if __name__=='__main__':main()
