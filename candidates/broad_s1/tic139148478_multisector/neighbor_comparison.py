"""Compare the target and bright rotating neighbor at matched public cadences.

Reproduce: python candidates/broad_s1/tic139148478_multisector/neighbor_comparison.py
This measures common variability and candidate-phase depths, not pixel origin.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys

import numpy as np
from astropy.timeseries import LombScargle

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from broad_vet import load_records, read_source, local_quadratic, phase, mad
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent


if __name__ == '__main__':
    target = load_records([139148478])[0]
    target_rows = [dict(sector=1,path=target['path'],sha256=target['sha256'],url=target['url'])]
    target_rows += [r for r in json.loads((ROOT/'data/broad_followup_queries/tic139148478_qlp.json').read_text()) if r['status']=='retained']
    neighbors = {r['sector']:r for r in json.loads((ROOT/'data/broad_followup_queries/tic139148479_qlp.json').read_text()) if r['status']=='retained'}
    dossier = json.loads((OUT/'multisector.json').read_text())
    ephemerides = {r['sector']:r['ephemeris'] for r in dossier['independent_sector_searches']}
    rows, graphs = [], []
    for row in target_rows:
        sector = row['sector']
        if sector not in neighbors:
            continue
        a,_ = read_source(row)
        b,_ = read_source(neighbors[sector])
        usea = (a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['SAP_FLUX'])&np.isfinite(a['KSPSAP_FLUX_ERR'])&(a['KSPSAP_FLUX_ERR']>0)
        useb = (b['QUALITY']==0)&np.isfinite(b['TIME'])&np.isfinite(b['SAP_FLUX'])&np.isfinite(b['KSPSAP_FLUX_ERR'])&(b['KSPSAP_FLUX_ERR']>0)
        ta,fa,ea = a['TIME'][usea],a['SAP_FLUX'][usea],a['KSPSAP_FLUX_ERR'][usea]
        tb,fb,eb = b['TIME'][useb],b['SAP_FLUX'][useb],b['KSPSAP_FLUX_ERR'][useb]
        oa,ob=np.argsort(ta),np.argsort(tb)
        ta,fa,ea=ta[oa],fa[oa],ea[oa];tb,fb,eb=tb[ob],fb[ob],eb[ob]
        index=np.searchsorted(tb,ta)
        right=np.clip(index,0,len(tb)-1);left=np.clip(index-1,0,len(tb)-1)
        nearest=np.where(abs(tb[right]-ta)<abs(tb[left]-ta),right,left)
        matched=abs(tb[nearest]-ta)<.4*np.median(np.diff(ta))
        time, ft, fn=ta[matched],fa[matched]/np.median(fa[matched]),fb[nearest[matched]]/np.median(fb[nearest[matched]])
        correlation=float(np.corrcoef(ft,fn)[0,1])
        design=np.column_stack([np.ones(len(time)),fn-1.])
        beta=np.linalg.lstsq(design,ft,rcond=None)[0]
        residual=ft-design@beta
        frequency=np.linspace(.03,.3,3001)
        power=LombScargle(tb-tb.min(),fb,eb).power(frequency)
        rotation=float(1/frequency[int(np.argmax(power))])
        eph=ephemerides[sector]
        p,epoch,d=eph['period'],eph['t0_btjd'],eph['duration_days']
        target_depth,_=local_quadratic(ta,fa,ea,p,epoch,d)
        neighbor_depth,_=local_quadratic(tb,fb,eb,p,epoch,d)
        rows.append(dict(sector=sector,n_matched=int(matched.sum()),raw_flux_correlation=correlation,
                         neighbor_to_target_relative_variation_scale=float(beta[1]),
                         target_raw_mad=mad(ft),after_neighbor_linear_template_mad=mad(residual),
                         neighbor_long_timescale_ls_period_days=rotation,
                         target_independent_candidate_ephemeris=eph,
                         target_candidate_phase_depth=target_depth,neighbor_candidate_phase_depth=neighbor_depth,
                         target_source={k:row[k] for k in ('url','path','sha256')},
                         neighbor_source={k:neighbors[sector][k] for k in ('url','path','sha256')}))
        graphs.append((sector,time,ft,fn,beta,residual,p,epoch))
        print('NEIGHBOR_QA',sector,'rawcorr',correlation,'rotationP',rotation,
              'target/neighbor_candidateSNR',target_depth['snr'],neighbor_depth['snr'],flush=True)
    payload=dict(checked_utc=datetime.now(timezone.utc).isoformat(),target_tic=139148478,neighbor_tic=139148479,
                 neighbor_name='HD217027',neighbor_separation_arcsec=51.5,known_rotation_period_days=9.50244,
                 current_target_tic_contamination_ratio=.995876,sector_reports=rows,
                 limitations=['Shared rotational variability is contamination evidence, not proof of the transit source.',
                              'Candidate phase was independently optimized per target sector; its optimized score is not a false-alarm probability.',
                              'Long-timescale LS periods are exploratory and have no quoted uncertainty.',
                              'PSF weighting differs between apertures and sectors; no physical dilution correction is inferred.'],
                 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'neighbor_comparison.json').write_text(json.dumps(payload,indent=2,allow_nan=False))
    fig,axes=plt.subplots(len(graphs),2,figsize=(11,4*len(graphs)),layout='constrained');axes=np.atleast_2d(axes)
    for ax,(sector,time,ft,fn,beta,residual,p,epoch) in zip(axes,graphs):
        ax[0].plot(time,ft,'.',ms=2,label='Target raw SAP')
        ax[0].plot(time,beta[0]+beta[1]*(fn-1),lw=.8,label='Scaled bright-neighbor template')
        ax[0].set(xlabel='BTJD',ylabel=f'Sector{sector} flux');ax[0].legend(fontsize=8)
        ax[1].plot(phase(time,p,epoch),residual,'.',ms=2,color='.5')
        ax[1].set(xlabel='Candidate phase time (days)',ylabel='Residual after neighbor linear template')
    fig.savefig(OUT/'neighbor_comparison.png',dpi=130);plt.close(fig)
