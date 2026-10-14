"""Summarize source-attributed eclipses of TIC273879360 in four TESS sectors.

Fits cadence-integrated phenomenological eclipse templates to preserved native
source-core photometry. Radius/mass/binding/novelty are not inferred by this fit.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import ROOT
from broad_vet import integrated_trapezoid
from pixel_vet_49617264 import safe

OUT=ROOT/'candidates/stellar_systems/tic273879360'
P=7.586894226811602
EPOCH=1327.898092379152
PHASE=.693683559463394


def inputs():
    result=[]
    for sector in [27,29,67,95]:
        folder=ROOT/f'candidates/stellar_systems/tic273879358/pixels/sector{sector}/large_source_validation'if sector<30 else OUT/f'pixels/sector{sector}/later_source_validation'
        report_path=folder/'pixel_vetting.json';arrays_path=folder/'aperture_photometry.npz'
        report=json.loads(report_path.read_text());name='gaia_4692910966756908288'if sector<30 else'catalog_target'
        a=np.load(arrays_path);exposure=report['raw_manifest']['exposure_days']
        for kind in ['primary','secondary']:
            t,f,e=[a[name+'_'+kind+'_'+column]for column in ['time','flux','error']]
            result.append(dict(sector=sector,kind=kind,t=t,f=f,e=e,exposure=exposure,
                               stats=report['aperture_photometry'][name]['eclipses'][kind],
                               report_path=str(report_path.relative_to(ROOT)),report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
                               arrays_path=str(arrays_path.relative_to(ROOT)),arrays_sha256=hashlib.sha256(arrays_path.read_bytes()).hexdigest()))
    return result


def shape_fit(parts):
    t=np.concatenate([r['t']for r in parts]);f=np.concatenate([r['f']for r in parts]);e=np.concatenate([r['e']for r in parts]);
    groups=np.concatenate([np.full(len(r['t']),i)for i,r in enumerate(parts)])
    exposure=np.concatenate([np.full(len(r['t']),r['exposure'])for r in parts])
    kind=np.array([r['kind']=='secondary'for r in parts])[groups]
    cycle=np.rint((t-EPOCH-np.where(kind,PHASE*P,0))/P).astype(int)
    dt=t-EPOCH-cycle*P-np.where(kind,PHASE*P,0)
    indicators=np.column_stack([groups==i for i in range(len(parts))])
    def solve(parameters):
        shift,dp,secondary_shift,d0,d1,r0,r1=parameters
        x=dt-shift-cycle*dp-np.where(kind,secondary_shift,0)
        q=np.zeros(len(x))
        for select,duration,ratio in [(~kind,d0,r0),(kind,d1,r1)]:
            q[select]=integrated_trapezoid(x[select],duration,ratio*duration,exposure[select])
        design=np.column_stack([indicators,indicators*(-q[:,None])])
        a=design/e[:,None];b=f/e
        beta=np.linalg.lstsq(a,b,rcond=None)[0]
        return beta,(b-a@beta)
    initial=[0,0,0,.17,.10,.4,.35]
    fit=least_squares(lambda par:solve(par)[1],initial,
                      bounds=([-.03,-8e-5,-.08,.08,.04,.02,.02],[.03,8e-5,.08,.30,.22,.5,.5]),
                      x_scale=[.003,1e-5,.003,.17,.1,.4,.3],diff_step=1e-4,max_nfev=200)
    beta,residual=solve(fit.x)
    dof=len(f)-len(beta)-len(fit.x)
    chi2=float(residual@residual)
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*max(1.,chi2/max(1,dof))
    err=np.sqrt(np.maximum(0,np.diag(cov)))
    p=P+fit.x[1];epoch=EPOCH+fit.x[0]
    phase=(PHASE*P+fit.x[2])/p
    return dict(success=bool(fit.success),period_days=float(p),period_error_days=float(err[1]),
                epoch_btjd=float(epoch),epoch_error_days=float(err[0]),
                secondary_phase=float(phase),secondary_phase_error_conditional=float(err[2]/p),
                primary_duration_days=float(fit.x[3]),secondary_duration_days=float(fit.x[4]),
                primary_ingress_fraction=float(fit.x[5]),secondary_ingress_fraction=float(fit.x[6]),
                per_sector_kind_depths=[dict(sector=r['sector'],kind=r['kind'],depth=float(beta[len(parts)+i]))for i,r in enumerate(parts)],
                chi2=chi2,nominal_dof=dof,n_points=len(f),parameters=fit.x.tolist(),
                covariance=cov.tolist(),
                warning='Conditional joint trapezoid fit to these same data, approximate diagonal errors. Shared light-curve shape and baseline covariance/systematics are not fully covered; not a mass or physical orbital fit.')


def main():
    OUT.mkdir(parents=True,exist_ok=True);parts=inputs();shape=shape_fit(parts)
    ephemeris=json.loads((ROOT/'data/stellar_systems/tic273879360_config.json').read_text())
    records=[]
    for row in parts:
        for event in row['stats']['events']:
            records.append(dict(sector=row['sector'],kind=row['kind'],**event))
    with (OUT/'source_event_depths.csv').open('w')as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    metadata=[{k:v for k,v in r.items()if k not in ['t','f','e']}for r in parts]
    result=dict(tic=273879360,gaia_dr3_source_id='4692910966756908288',
                checked_utc=datetime.now(timezone.utc).isoformat(),
                original_foreground_aperture_tic=273879358,
                preselected_ephemeris=dict(period_days=P,epoch_btjd=EPOCH,secondary_phase=PHASE),
                observed_sectors=[27,29,67,95],source_core_events=records,
                fitted_phenomenological_shape=shape,inputs=metadata,
                excluded_sector68=dict(reason='Target lies off science silicon, column2130.37; retained raw overscan cutout is not a science observation.',
                                       manifest='candidates/stellar_systems/tic273879360/pixels/sector68/later_source_validation/raw_manifest.json'),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                assessment='Resolved background eclipsing source on/near TIC273879360, with luminous companion strongly favored; no prior classification in listed finite checks; novelty and component masses unconfirmed.',
                caveats=['Observed1pixel-core depths are diluted and differ by camera/aperture; they are not component radii.',
                         'Native TESS PRF/WCS systematics can shift empirical centroids by a fraction of a pixel.',
                         'The catalogued source already exists in Gaia/TIC and unvetted target inventories; only stellar multiplicity/classification is potentially new.'])
    (OUT/'source_summary.json').write_text(json.dumps(safe(result),indent=2)+'\n')
    fig,axes=plt.subplots(2,2,figsize=(12,8));colors={27:'tab:blue',29:'tab:orange',67:'tab:green',95:'tab:red'}
    for r in parts:
        index=0 if r['kind']=='primary' else 1
        cycle=np.rint((r['t']-EPOCH-index*PHASE*P)/P)
        x=(r['t']-EPOCH-cycle*P-index*PHASE*P)*24
        axes[0,index].scatter(x,(1-r['f'])*100,s=2,alpha=.12,c=colors[r['sector']])
        bins=np.linspace(-4,4,65);xx=[];yy=[]
        for lo,hi in zip(bins[:-1],bins[1:]):
            val=r['f'][(x>=lo)&(x<hi)]
            if len(val)>1:xx.append((lo+hi)/2);yy.append((1-val.mean())*100)
        axes[0,index].plot(xx,yy,'.',label=f"S{r['sector']}",color=colors[r['sector']])
        for ev in r['stats']['events']:
            axes[1,0].errorbar(ev['time_btjd'],ev['depth']*100,ev['depth_error']*100,
                             fmt='o'if index==0 else's',color=colors[r['sector']],alpha=.7)
    for i in [0,1]:axes[0,i].set(xlim=(-4,4),xlabel='Hours from predicted eclipse',ylabel='Source-core flux loss (%)',title='Primary'if i==0 else'Secondary');axes[0,i].legend()
    axes[1,0].set(xlabel='BTJD (TDB)',ylabel='Fitted source-core depth (%)')
    axes[1,1].axis('off');axes[1,1].text(0,.9,f"TIC273879360 / Gaia4692910966756908288\nBackground star65.59 arcsec from foreground target\nFour usable sectors; deep primary+secondary eclipses\nFitted recurrence: {shape['period_days']:.9f} d\nSecondary phase: {shape['secondary_phase']:.5f}\nS68 excluded: off science detector\nUnrecognized EB candidate; no mass/novelty confirmation",va='top',fontsize=11)
    fig.suptitle('A background eclipsing system recovered from a bright-star aperture');fig.tight_layout();fig.savefig(OUT/'source_eclipses.png',dpi=170);plt.close(fig)
    print(json.dumps(safe(shape),indent=2))


if __name__=='__main__':main()
