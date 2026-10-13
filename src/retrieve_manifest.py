"""Restore exact recorded observations without repeating target selection."""
import argparse,hashlib,json,os,urllib.request
from pathlib import Path
from archive import verify_fits_checksums
from astropy.io import fits

def restore(manifest,cache):
    for rec in json.loads(Path(manifest).read_text()):
        name=rec['data_uri'].rsplit('/',1)[-1]
        path=Path(cache)/'lightcurves'/name;path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():
            tmp=path.with_suffix('.part')
            with urllib.request.urlopen(rec['retrieval_url'],timeout=120) as response,tmp.open('wb') as out:
                while chunk:=response.read(1024*1024):out.write(chunk)
            tmp.replace(path)
        sha=hashlib.sha256(path.read_bytes()).hexdigest()
        if sha!=rec['sha256']:
            raise ValueError(f'Original SHA mismatch for {name}; preserve both versions and investigate processing changes')
        with fits.open(path) as hdul:verify_fits_checksums(path,hdul)
        print('Verified',name,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',default='data/observations.json');p.add_argument('--cache',default=os.environ.get('TESS_DATA_DIR','/workspace/tess-research-data'))
    a=p.parse_args();restore(a.manifest,a.cache)
