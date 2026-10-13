"""MAST query/download cache with explicit identifiers and integrity manifest."""
from pathlib import Path
import argparse, csv, hashlib, json, os, time, urllib.parse, urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data'))

def invoke(service, params):
    payload = {'service': service, 'params': params, 'format': 'json', 'pagesize': 10000, 'page': 1}
    body = urllib.parse.urlencode({'request': json.dumps(payload)}).encode()
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request('https://mast.stsci.edu/api/v0/invoke', data=body), timeout=90) as response:
                result = json.load(response)
            if result.get('status') != 'COMPLETE':
                raise RuntimeError('MAST query incomplete: ' + result.get('status', 'unknown'))
            return result
        except Exception:
            if attempt == 2: raise
            time.sleep(2 * (attempt + 1))

def query(tic):
    path = CACHE / 'queries' / f'tic{tic}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists(): return json.loads(path.read_text())
    result = invoke('Mast.Caom.Filtered', {'columns': '*', 'filters': [
        {'paramName': 'obs_collection', 'values': ['TESS']},
        {'paramName': 'target_name', 'values': [str(tic)]}]})
    result['retrieval_utc'] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(result, indent=2))
    return result

def lightcurves(tic):
    rows = query(tic)['data']
    # SPOC 2-minute products, avoiding duplicate 20-second reductions and HLSPs.
    return sorted([r for r in rows if str(r.get('dataURL','')).endswith('-s_lc.fits')],
                  key=lambda r: r.get('sequence_number') or 0)

def verify_fits_checksums(path, hdul):
    """Verify raw HDU ones-complement sums, including data when DATASUM is absent.

    Astropy 8 verify_checksum assumes zero data sum without a DATASUM card;
    SPOC FITS often has CHECKSUM but no DATASUM. Verify full raw bytes instead.
    Never alter the original headers or expected checksum.
    """
    import numpy as np
    raw = path.read_bytes()
    status = []
    for index, hdu in enumerate(hdul):
        info = hdu.fileinfo()
        start = info['hdrLoc']
        stop = info['datLoc'] + info['datSpan']
        if 'CHECKSUM' not in hdu.header:
            status.append({'hdu': index, 'status': 'absent'})
            continue
        block = raw[start:stop]
        total = int(np.frombuffer(block, dtype='>u4').sum(dtype=np.uint64))
        while total >> 32:
            total = (total & 0xffffffff) + (total >> 32)
        if total != 0xffffffff:
            raise ValueError(f'Raw FITS checksum failed: {path.name} HDU {index}')
        status.append({'hdu': index, 'status': 'verified_raw_ones_complement'})
    return status

def retrieve(tic, row):
    uri = row['dataURL']
    name = uri.rsplit('/', 1)[-1]
    path = CACHE / 'lightcurves' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    url = 'https://mast.stsci.edu/api/v0.1/Download/file?' + urllib.parse.urlencode({'uri': uri})
    if not path.exists():
        tmp = path.with_suffix('.part')
        with urllib.request.urlopen(url, timeout=120) as response, tmp.open('wb') as stream:
            while block := response.read(1024 * 1024): stream.write(block)
        tmp.replace(path)
    from astropy.io import fits
    with fits.open(path) as hdul:
        hdul.verify('exception')
        integrity = verify_fits_checksums(path, hdul)
        header = hdul[0].header
        details = {k: header.get(k) for k in ('TICID','SECTOR','CAMERA','CCD','PROCVER','DATE','TESSMAG','RADIUS','TEFF')}
    return {'tic': int(tic), 'sector': row.get('sequence_number'), 'obsid': row['obsid'],
            'obs_id': row['obs_id'], 'data_uri': uri, 'retrieval_url': url,
            'local_path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'bytes': path.stat().st_size, 'metadata': details, 'fits_checksums': integrity,
            'checked_utc': datetime.now(timezone.utc).isoformat()}

def main():
    p = argparse.ArgumentParser(); p.add_argument('tic',type=int,nargs='+'); p.add_argument('--limit',type=int,default=0); p.add_argument('--recent',action='store_true')
    a = p.parse_args()
    out = ROOT/'data'/'observations.json'
    records = json.loads(out.read_text()) if out.exists() else []
    for tic in a.tic:
        rows = lightcurves(tic)
        print('TIC',tic,'SPOC sectors:',[r.get('sequence_number') for r in rows],flush=True)
        if a.limit: rows = rows[-a.limit:] if a.recent else rows[:a.limit]
        for row in rows:
            record = retrieve(tic,row)
            records = [r for r in records if r['data_uri'] != record['data_uri']] + [record]
            out.write_text(json.dumps(records,indent=2))
            print('Cached',record['metadata']['SECTOR'],record['bytes'],flush=True)
if __name__ == '__main__': main()
