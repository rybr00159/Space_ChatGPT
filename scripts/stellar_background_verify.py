#!/usr/bin/env python3
"""Verify retained raw bytes, method identities and reproduced fitted arrays."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import validate_background_binary as vb


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-manifest',action='store_true',help='Checkpoint current verified methods/outputs after intentional regeneration.')
    args=parser.parse_args()
    reports={name:json.loads((vb.OUT/name).read_text())for name in
             ['physical_fit.json','physical_controls.json','alias_and_period_tests.json']}
    methods={
        'src/validate_background_binary.py':sha(ROOT/'src/validate_background_binary.py'),
        'src/background_binary_controls.py':sha(ROOT/'src/background_binary_controls.py'),
        'src/background_binary_alias_tests.py':sha(ROOT/'src/background_binary_alias_tests.py'),
        'scripts/stellar_background_verify.py':sha(Path(__file__)),
    }
    assert reports['physical_fit.json']['script_sha256']==methods['src/validate_background_binary.py']
    for name,key in [('physical_controls.json','src/background_binary_controls.py'),
                     ('alias_and_period_tests.json','src/background_binary_alias_tests.py')]:
        assert reports[name]['source_sha256']==methods[key],f'Method mismatch: {name}'
        assert reports[name]['physical_method_sha256']==methods['src/validate_background_binary.py']
    events,sources,excluded=vb.extract()
    saved=reports['physical_fit.json']
    assert len(events)==saved['n_events']and excluded==saved['excluded_events']
    model=vb.BinaryModel(events,saved['supersample'])
    residual,geometry,prediction,baseline=model.residual(np.asarray(saved['best_fit']['parameters']),return_info=True)
    arrays=np.load(vb.OUT/'raw_model_arrays.npz')
    for key,value in [('time',model.t),('model',prediction),('baseline',baseline),('standardized_residual',residual),
                      ('raw_flux',np.concatenate([e['y']for e in events])),('error',np.concatenate([e['e']for e in events]))]:
        np.testing.assert_allclose(arrays[key],value,rtol=1e-11,atol=1e-11,err_msg=key)
    assert abs(float(residual@residual)-saved['best_fit']['chi2'])<1e-5
    assert saved['best_fit']['success']and all(v['success']for v in saved['alternative_starts'])
    assert max(v['chi2']for v in saved['alternative_starts'])-min(v['chi2']for v in saved['alternative_starts'])<.1
    early=reports['alias_and_period_tests.json']['early_period_diagnostic']
    assert early['unrestricted_early_fit']['success']and early['early_fit_at_full_period']['success']
    near=[v for v in early['profile']if v['delta_chi2']<4]
    assert near and all(v['success']for v in near)
    artifacts={str(p.relative_to(ROOT)):sha(p)for p in vb.OUT.iterdir()
               if p.suffix in ['.json','.png','.npz']and p.name!='verification_manifest.json'}
    manifest_path=vb.OUT/'verification_manifest.json'
    current=dict(checked_utc=datetime.now(timezone.utc).isoformat(),methods=methods,artifacts=artifacts,
                 raw_sources=[dict(path=s['path'],sha256=s['sha256'])for s in sources],
                 reproduced_events=len(events),reproduced_cadences=model.n,chi2=float(residual@residual))
    if args.write_manifest:
        manifest_path.write_text(json.dumps(current,indent=2)+'\n')
    else:
        manifest=json.loads(manifest_path.read_text())
        assert manifest['methods']==methods and manifest['artifacts']==artifacts,'Checkpoint changed; verify intentional regeneration before refreshing.'
        assert manifest['raw_sources']==current['raw_sources']
    print(f'VERIFIED {len(sources)} raw cubes, {len(methods)} methods, {len(artifacts)} outputs; reproduced {len(events)} events/{model.n} cadences.')


if __name__=='__main__':main()
