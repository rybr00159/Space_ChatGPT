"""Retrieve exact follow-up QLP products, retaining all successful observations."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse,hashlib,json,requests
from astropy.io import fits
from archive import CACHE,verify_fits_checksums

def retrieve(tic,sector):
    s=str(tic).zfill(16);sub='/'.join(s[i:i+4] for i in range(0,16,4));name=f'hlsp_qlp_tess_ffi_s{sector:04}-{s}_tess_v01_llc.fits'
    url=f'https://archive.stsci.edu/hlsps/qlp/s{sector:04}/{sub}/{name}'
    folder=CACHE/'broad/followup';folder.mkdir(parents=True,exist_ok=True);path=folder/name
    if not path.exists():
        r=requests.get(url,timeout=45)
        if r.status_code==404:return {'tic':tic,'sector':sector,'status':'not_found','url':url}
        r.raise_for_status();tmp=path.with_suffix('.part');tmp.write_bytes(r.content);tmp.replace(path)
    with fits.open(path) as hd:
        hd.verify('exception');checks=verify_fits_checksums(path,hd)
        h=hd[1].header
        result={'tic':tic,'sector':sector,'status':'retained','url':url,'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size,'checks':checks,'n_rows':len(hd[1].data),'metadata':{k:h.get(k,hd[0].header.get(k)) for k in ['TIMESYS','BJDREFI','BJDREFF','TESSMAG','RADIUS','TEFF','CAMERA','CCD']}}
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('tic',type=int);p.add_argument('--sectors',nargs='+',type=int,required=True);a=p.parse_args()
    out=Path(__file__).resolve().parents[1]/f'data/broad_followup_queries/tic{a.tic}_qlp.json';out.parent.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(4) as pool:rows=list(pool.map(lambda sector:retrieve(a.tic,sector),a.sectors))
    out.write_text(json.dumps(rows,indent=2));print('QLP retained sectors:',[r['sector'] for r in rows if r['status']=='retained'],flush=True)
