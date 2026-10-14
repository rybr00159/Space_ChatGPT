"""Verify retained stellar-phase observations, result identities and controls."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from archive import verify_fits_checksums


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    catalog_paths=[]
    for folder in [ROOT/'data/stellar_systems',ROOT/'results/stellar_systems/scout_catalogs']:
        for manifest in folder.rglob('*.manifest.json'):
            data=json.loads(manifest.read_text());expected=data.get('sha256')
            if expected is None:continue
            stem=str(manifest).removesuffix('.manifest.json')
            possibilities=[Path(data.get('local_path',stem)),Path(stem),Path(stem+'.csv'),Path(stem+'.tsv')]
            p=next((x for x in possibilities if x.exists()),None)
            if p is None:raise FileNotFoundError(manifest)
            assert sha(p)==expected,p
            catalog_paths.append(str(p))
    pixel_reports=[];raw_records={};method=sha(ROOT/'src/stellar_pixels.py')
    for folder in [ROOT/'candidates/stellar_systems/tic273879358',ROOT/'candidates/stellar_systems/tic273879360']:
        for p in folder.rglob('pixel_vetting.json'):
            data=json.loads(p.read_text());assert data['script_sha256']==method,p
            pixel_reports.append(str(p.relative_to(ROOT)))
        for p in folder.rglob('raw_manifest.json'):
            row=json.loads(p.read_text());raw_records[row['path']]=row
    for p in (ROOT/'candidates/stellar_systems/hidden_sources').rglob('raw_manifest.json'):
        row=json.loads(p.read_text());raw_records[row['path']]=row
    for p,row in raw_records.items():
        assert sha(p)==row['sha256'],p
        with fits.open(p)as hd:verify_fits_checksums(Path(p),hd)
    for sector in [27,29]:
        folder=ROOT/f'candidates/stellar_systems/tic273879358/pixels/sector{sector}/large_source_validation'
        data=json.loads((folder/'pixel_vetting.json').read_text())
        gaia=pd.read_csv(folder/'gaia_positions.csv',dtype={'source_id':str})
        star=gaia[gaia.source_id=='4692910966756908288'].iloc[0]
        for kind in ['primary','secondary']:
            xy=np.load(folder/f'{kind}_differences.npz')['bootstrap_xy']
            assert len(xy)==300
            assert np.all(np.linalg.norm(xy-[star.x,star.y],axis=1)<np.linalg.norm(xy-np.asarray(data['host_pixel']),axis=1))
    summary=json.loads((ROOT/'candidates/stellar_systems/tic273879360/source_summary.json').read_text())
    assert summary['script_sha256']==sha(ROOT/'src/stellar_background273.py')
    for row in summary['inputs']:
        assert sha(ROOT/row['report_path'])==row['report_sha256']
        assert sha(ROOT/row['arrays_path'])==row['arrays_sha256']
    assert len(summary['source_core_events'])==23
    assert sha(ROOT/'src/broad_vet.py')=='23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333'
    scout=json.loads((ROOT/'results/stellar_systems/scout_manifest.json').read_text())
    assert scout['source_sha256']==sha(ROOT/'src/stellar_scout.py')
    count=0
    for p in (ROOT/'results/stellar_systems').glob('scout_tic*.json'):
        row=json.loads(p.read_text());assert row['scout_source_sha256']==scout['source_sha256'];count+=1
    assert count==192
    hidden_sha=sha(ROOT/'scripts/stellar_hidden_sources.py')
    for p in (ROOT/'candidates/stellar_systems/hidden_sources').glob('*.json'):
        assert json.loads(p.read_text())['code_sha256']==hidden_sha,p
    record=dict(checked_utc=datetime.now(timezone.utc).isoformat(),verified_catalog_snapshots=len(catalog_paths),
                verified_catalog_paths=catalog_paths,verified_pixel_reports=pixel_reports,
                verified_raw_pixel_files=len(raw_records),raw_pixel_bytes=sum(Path(p).stat().st_size for p in raw_records),
                background_native_covered_events=23,verified_scout_results=count,
                background_attribution_control='All300drawsperphaseinS27/S29favorbackgroundoverforeground',
                frozen_broad_source='unchanged',discovery_or_mass_claim=False,paid_api_calls=0,
                verification_script_sha256=sha(Path(__file__)))
    (ROOT/'results/stellar_systems/verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(f"Verified {len(catalog_paths)}catalog snapshots, {len(raw_records)}raw pixels, {len(pixel_reports)}pixel reports, all192scout identities and source controls")


if __name__=='__main__':main()
