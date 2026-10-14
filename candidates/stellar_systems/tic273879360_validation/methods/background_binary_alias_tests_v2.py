"""Check missing eclipses and quantify early-epoch ephemeris extrapolation.

Control phases are chosen before looking at their flux: a P/2 two-body orbit
must repeat both eclipse families at these locations. Early-only likelihood
profiles are conditional on the same extraction and binary model; the profile
is not a calibrated confidence statement with all instrumental systematics.
"""
from pathlib import Path
from datetime import datetime, timezone
import copy
import hashlib
import json
import os
import numpy as np
from scipy.optimize import least_squares
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import validate_background_binary as vb
from background_binary_controls import fixed_geometry_fit
from pixel_vet_49617264 import safe


def fit_at_period(model, seed, dp=None):
    active=np.ones(12,bool)
    if dp is not None:active[7]=False
    full=seed.copy()
    if dp is not None:full[7]=dp
    lo=vb.LOW.copy();hi=vb.HIGH.copy()
    lo[6:8]=[-.5,-.003];hi[6:8]=[.5,.003]
    def expand(v):
        p=full.copy();p[active]=v;return p
    fit=least_squares(lambda v:model.residual(expand(v)),full[active],
                      bounds=(lo[active],hi[active]),x_scale=vb.SCALE[active],
                      diff_step=1e-4,max_nfev=400)
    p=expand(fit.x);r,g,_,_=model.residual(p,return_info=True)
    return dict(parameters=p.tolist(),chi2=float(r@r),geometry=g,
                success=bool(fit.success),nfev=int(fit.nfev))


def control_depth(event, seed, eclipse_kind):
    # Evaluate the already fitted primary/secondary shape translated into the
    # control window. A linear coefficient tests signed flux loss, allowing
    # brightening as well as dimming rather than imposing positive depth.
    ev=copy.deepcopy(event)
    center=seed[6]+vb.T0
    if eclipse_kind=='secondary':
        e=np.hypot(seed[3],seed[4]);w=np.arctan2(seed[4],seed[3])
        center+=vb.phase_secondary(e,w)*(vb.P0+seed[7])
    ev['t']=center+event['x']
    model=vb.BinaryModel([ev]);f,_=model.flux(seed)
    template=(1-f)/max(1e-10,float(np.max(1-f)))
    design=np.column_stack([event['baseline_design'],-template])
    a=design/event['e'][:,None];y=event['y']/event['e']
    coeff=np.linalg.lstsq(a,y,rcond=None)[0]
    r=y-a@coeff;cov=np.linalg.pinv(a.T@a)*max(1.,float(r@r)/(len(r)-design.shape[1]))
    scale=coeff[0]
    return dict(sector=event['sector'],cycle=event['cycle'],kind=event['kind'],
                depth=float(coeff[-1]/scale),error=float(np.sqrt(cov[-1,-1])/scale),
                chi2=float(r@r),n_points=len(r),central_template_expected_loss=float(np.max(1-f)))


def main():
    fit=json.loads((vb.OUT/'physical_fit.json').read_text())
    seed=np.asarray(fit['best_fit']['parameters'])
    events,_,_=vb.extract()
    early=vb.BinaryModel([e for e in events if e['sector']in [27,29]])
    later=vb.BinaryModel([e for e in events if e['sector']in [67,95]])
    unconstrained=fit_at_period(early,seed)
    fixed=fit_at_period(early,seed,seed[7])
    best_dp=unconstrained['parameters'][7]
    grid=np.unique(np.r_[np.linspace(best_dp-.00018,best_dp+.00018,25),seed[7]])
    profile=[fit_at_period(early,np.asarray(unconstrained['parameters']),dp) for dp in grid]
    # Retain the lower converged result at the full period. The first profile
    # experiment exposed sensitivity to initialization of correlated epoch/P.
    first_fixed=copy.deepcopy(fixed)
    alternatives=[fixed]+[v for dp,v in zip(grid,profile)if dp==seed[7]]
    fixed=min([v for v in alternatives if v['success']],key=lambda v:v['chi2'])
    minimum=min([unconstrained['chi2']]+[v['chi2']for v in profile])
    intervals={}
    for threshold in [1.,3.84]:
        good=np.array([v['success']for v in profile])
        delta=np.array([v['chi2']-minimum for v in profile])
        crossings=[]
        for i in range(len(grid)-1):
            if good[i]and good[i+1]and (delta[i]-threshold)*(delta[i+1]-threshold)<0:
                crossings.append(float(vb.P0+grid[i]+(grid[i+1]-grid[i])*(threshold-delta[i])/(delta[i+1]-delta[i])))
        intervals[str(threshold)]=crossings if len(crossings)==2 else None
    fixed_prediction=fixed_geometry_fit(later,np.asarray(fixed['parameters']))
    later_refit=vb.fit_model(later,seed)
    phases=[('missing_primary',.5),('missing_secondary',(vb.PHASE+.5)%1)]
    controls,sources,excluded=vb.extract(phase_centers=phases)
    rows=[control_depth(e,seed,'primary'if e['kind']=='missing_primary'else'secondary')for e in controls]
    summaries={};rng=np.random.default_rng(713273879360)
    for kind,phase in phases:
        subset=[r for r in rows if r['kind']==kind]
        depths=np.array([r['depth']for r in subset]);err=np.array([r['error']for r in subset]);weight=err**-2
        draws=[float(np.mean(depths[rng.integers(0,len(depths),len(depths))]))for _ in range(10000)]
        summaries[kind]=dict(phase=phase,n_events=len(subset),weighted_depth=float(np.sum(weight*depths)/weight.sum()),
                             conditional_white_error=float(weight.sum()**-.5),event_mean_depth=float(depths.mean()),
                             event_bootstrap_95=np.percentile(draws,[2.5,97.5]).tolist(),
                             expected_central_loss_range=[min(r['central_template_expected_loss']for r in subset),max(r['central_template_expected_loss']for r in subset)])
    initial_controls=json.loads((vb.OUT/'physical_controls.json').read_text())
    initial_early=initial_controls['epoch_holdouts']['early_to_later']['training_fit']
    result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                physical_method_sha256=hashlib.sha256(Path(vb.__file__).read_bytes()).hexdigest(),sources=sources,
                early_period_diagnostic=dict(unrestricted_early_fit=unconstrained,early_fit_at_full_period=fixed,
                     first_fixed_period_initialization=first_fixed,
                     initial_restricted_training_success=initial_early['success'],
                     initial_restricted_training_period=initial_early['geometry']['period_days'],
                     initial_restricted_prediction_penalty=initial_controls['epoch_holdouts']['early_to_later']['conditional_prediction_penalty'],
                     full_period_training_penalty=fixed['chi2']-minimum,
                     profile=[dict(period=vb.P0+dp,delta_chi2=v['chi2']-minimum,success=v['success'])for dp,v in zip(grid,profile)],
                     interpolated_conditional_profile_crossings=intervals,
                     later_prediction_at_full_period=fixed_prediction,
                     later_refit_chi2=later_refit['chi2'],
                     later_prediction_penalty_at_full_period=fixed_prediction['chi2']-later_refit['chi2'],
                     early_primary_events=sum(e['kind']=='primary'and e['sector']in [27,29]for e in events)),
                shorter_period_controls=summaries,control_events=rows,excluded_control_events=excluded,
                limitations=['Period intervals are sampled conditional likelihood diagnostics, not all-systematics credible intervals.',
                             'Period was already measured from all sectors; fixing it is a diagnosis of the failed early-only extrapolation, not an independent prediction.',
                             'Absence of repeated eclipses tests a single ordinary P/2 binary, not every contrived multiple-system hypothesis.',
                             'Control bootstrap samples events and omits shared calibration errors.'])
    (vb.OUT/'alias_and_period_tests.json').write_text(json.dumps(safe(result),indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    axes[0].plot((grid-vb.P0)*86400,[v['chi2']-minimum for v in profile],'o-')
    axes[0].axvline(seed[7]*86400,color='black',ls='--',label='All-sector period')
    axes[0].axhline(1,color='grey');axes[0].set(xlabel='Early period offset from old P (seconds)',ylabel='Profile delta chi-square',ylim=(0,25));axes[0].legend()
    for i,(kind,_)in enumerate(phases):
        rs=[r for r in rows if r['kind']==kind]
        axes[1].errorbar(np.full(len(rs),i)+np.linspace(-.12,.12,len(rs)),[100*r['depth']for r in rs],
                         yerr=[100*r['error']for r in rs],fmt='.',alpha=.7)
    axes[1].axhline(0,color='grey');axes[1].set(xticks=[0,1],xticklabels=['Extra primary at phase 0.5','Extra secondary at phase 0.192'],ylabel='Signed eclipse template depth (%)')
    fig.tight_layout();fig.savefig(vb.OUT/'alias_and_period_tests.png',dpi=170);plt.close(fig)
    print('Early period:',unconstrained['geometry']['period_days'],'full-period penalty',result['early_period_diagnostic']['full_period_training_penalty'],flush=True)
    print('Extra eclipse tests:',summaries,flush=True)


if __name__=='__main__':main()
