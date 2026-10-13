"""Cache public NASA TAP catalog tables with retrieval manifests. No authentication."""
import csv
import hashlib
import io
import json
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parent
ENDPOINT='https://exoplanetarchive.ipac.caltech.edu/TAP/sync'
PS_COLS='pl_name,hostname,tic_id,gaia_dr2_id,ra,dec,sy_pnum,sy_snum,disc_year,discoverymethod,disc_facility,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_trandep,pl_rade,pl_radeerr1,pl_radeerr2,pl_ratror,pl_ratdor,pl_orbincl,pl_orbeccen,st_teff,st_rad,st_mass,sy_dist,sy_tmag,sy_vmag,sy_jmag,sy_kmag,pl_controv_flag,pl_orbper_reflink,pl_tranmid_reflink,pl_trandur_reflink'.split(',')
QUERIES={
 'confirmed_transiting.csv': 'select '+','.join(PS_COLS)+' from pscomppars where tran_flag=1',
 'toi.csv': 'select * from toi',
 'confirmed_all.csv': 'select tran_flag,rv_flag,'+','.join(PS_COLS)+' from pscomppars',
}

def retrieve(name,query):
    path=ROOT/name
    manifest=ROOT/(name+'.manifest.json')
    if path.exists() and manifest.exists():
        expected=json.loads(manifest.read_text())['sha256']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError(f'checksum mismatch in cached {name}; preserve and investigate before reuse')
        print(f'cache (checksum valid): {name}')
        return
    url=ENDPOINT+'?'+urlencode({'query':query,'format':'csv'})
    with urlopen(url,timeout=180) as response:
        body=response.read()
        status=response.status
        ctype=response.headers.get('Content-Type')
    decoded=body.decode('utf-8-sig')
    rows=list(csv.DictReader(io.StringIO(decoded)))
    if not rows or '<' in decoded[:20]:
        raise ValueError(f'invalid table response: {decoded[:200]}')
    path.write_bytes(body)
    metadata={'retrieved_utc':datetime.now(timezone.utc).isoformat(),
              'endpoint':ENDPOINT,'query':query,'url':url,'http_status':status,
              'content_type':ctype,'bytes':len(body),'row_count':len(rows),
              'sha256':hashlib.sha256(body).hexdigest(),'columns':list(rows[0]),
              'notes':'Current public archive snapshot; retrieval date is not a guarantee of completeness.'}
    manifest.write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'{name}: {len(rows)} rows, {len(body)} bytes')

if __name__=='__main__':
    for name,q in QUERIES.items():
        retrieve(name,q)

    # ExoFOP community submissions contain periods absent from the TOI catalog.
    name='ctoi.csv'
    path=ROOT/name
    manifest=ROOT/(name+'.manifest.json')
    if path.exists() and manifest.exists():
        expected=json.loads(manifest.read_text())['sha256']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('checksum mismatch in cached CTOI table')
        print('cache (checksum valid): ctoi.csv')
    else:
        url='https://exofop.ipac.caltech.edu/tess/download_ctoi.php?output=csv'
        with urlopen(url,timeout=180) as response:
            body=response.read();status=response.status
        rows=list(csv.DictReader(io.StringIO(body.decode('utf-8-sig'))))
        if not rows or 'TIC ID' not in rows[0]:
            raise ValueError('invalid ExoFOP CTOI response')
        path.write_bytes(body)
        manifest.write_text(json.dumps({'retrieved_utc':datetime.now(timezone.utc).isoformat(),'url':url,'http_status':status,'row_count':len(rows),'bytes':len(body),'sha256':hashlib.sha256(body).hexdigest(),'columns':list(rows[0])},indent=2)+'\n')
        print(f'ctoi.csv: {len(rows)} rows, {len(body)} bytes')
