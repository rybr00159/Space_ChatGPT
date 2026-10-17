#!/usr/bin/env python3
"""Bounded equal-star eclipse falsification test on original QLP events.

This tests existence of stellar alternatives; it is not a Bayesian FPP, host
identification or evidence that any particular neighboring star is eclipsing.
All source ratios refer to differential flux relative to the TIC target's
reference flux. No generic aperture contamination multiplier is applied.
"""
from pathlib import Path
import hashlib,json
from datetime import datetime,timezone
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
import batman
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

OUT=Path(__file__).resolve().parent
INDEPENDENT=OUT.parent
ROOT=OUT.parents[3]
frozen=json.loads((INDEPENDENT/'ephemeris_frozen_pre_holdout.json').read_text())
epoch,period=frozen['parameters']
true_period=2*period
windows=[]
for rec in json.loads((ROOT/'data/tic214311908/observations.json').read_text()):
    with fits.open(rec['path']) as f:
        d=f[1].data; name='DET_FLUX' if 'DET_FLUX' in d.names else 'KSPSAP_FLUX'
        m=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d[name])
        t=np.array(d['TIME'])[m]; y=np.array(d[name])[m]
    exposure=float(np.median(np.diff(t)[np.diff(t)<.05]))
    for n in range(int(np.floor((t.min()-epoch)/period)),int(np.ceil((t.max()-epoch)/period))+1):
        mid=epoch+n*period; x=t-mid;use=abs(x)<.3
        if np.sum(abs(x)<.025)<1 or use.sum()<12:continue
        tx=t[use];z=y[use]/np.median(y[use])-1;xx=x[use]
        oot=abs(xx)>.1;A=np.c_[np.ones(len(xx)),xx,xx*xx]
        beta=np.linalg.lstsq(A[oot],z[oot],rcond=None)[0]
        res=z-A@beta;sigma=1.4826*np.median(abs(res[oot]-np.median(res[oot])))
        windows.append(dict(sector=rec['sector'],cycle=n,t=tx,y=z,x=xx,X=A,sigma=sigma,exposure=exposure))
y=np.concatenate([w['y'] for w in windows]);sigma=np.concatenate([np.repeat(w['sigma'],len(w['t'])) for w in windows])
X=np.zeros((len(y),len(windows)*3));at=0
for j,w in enumerate(windows):
    X[at:at+len(w['t']),j*3:j*3+3]=w['X'];at+=len(w['t'])
Xw=X/sigma[:,None];projection=np.linalg.pinv(Xw)
def marginal(signal):
    z=(y-signal)/sigma;beta=projection@z
    return z-Xw@beta,X@beta
base_res,_=marginal(np.zeros(len(y)))
previous=json.loads((INDEPENDENT/'circular_shape_assessment.json').read_text())
assert len(y)==previous['n_points'] and len(windows)==previous['n_events']
assert abs(float(base_res@base_res)-previous['baseline_chi2'])<1e-9
np.savez_compressed(OUT/'original_window_inputs.npz',flux=y,noise_sigma=sigma,baseline_design=X,
                    time_btjd=np.concatenate([w['t'] for w in windows]),
                    event_cycle=np.concatenate([np.repeat(w['cycle'],len(w['t'])) for w in windows]),
                    sector=np.concatenate([np.repeat(w['sector'],len(w['t'])) for w in windows]))

def a_from_density(rho):
    # Each star has mass M and density rho; Kepler's law uses total mass 2M.
    return (2*6.67430e-11*(rho*1000)*(true_period*86400)**2/(3*np.pi))**(1/3)

def duration_from_b_rho(b,rho):
    a=a_from_density(rho)
    return true_period/np.pi*np.arcsin(np.sqrt(max(0,4-b*b))/np.sqrt(a*a-b*b))

results=[];signals={};limb=[.35,.25]
params=batman.TransitParams();params.t0=0;params.per=true_period;params.rp=1;params.a=a_from_density(9)
params.inc=89;params.ecc=0;params.w=90;params.u=limb;params.limb_dark='quadratic'
# Equal stars have identical primary/secondary profiles. Centering every local
# window on its conjunction preserves the 2P orbital velocity in the model.
models=[batman.TransitModel(params,w['x'],supersample_factor=15,exp_time=w['exposure']) for w in windows]
def stellar_signal(b,rho,dt,f):
    a=a_from_density(rho);params.a=a;params.inc=np.degrees(np.arccos(b/a));params.t0=dt
    return np.concatenate([(model.light_curve(params)-1)*f/2 for model in models])

flux_ratios=[.03,.05,.1,.7,1.0]
density_grid=[.1,.3,1.,3.,10.,30.]
for f in flux_ratios:
    # A few broad starts prevent accepting a poor local solution.
    b0=.01 if f<.04 else .5 if f<.075 else 1.15 if f<.2 else 1.75
    for rho in density_grid:
        def fun(p):return marginal(stellar_signal(p[0],rho,p[1],f))[0]
        fits_=[least_squares(fun,[b,0],bounds=([0,-.02],[1.999,.02]),max_nfev=120,x_scale='jac')
               for b in sorted(set([b0,.2,1.7]))]
        fit=min(fits_,key=lambda z:z.fun@z.fun);b,dt=fit.x
        rec=dict(kind='fixed_density_grid',binary_flux_ratio=f,density_g_cm3=rho,impact_parameter=float(b),
                 epoch_shift_days=float(dt),chi2=float(fit.fun@fit.fun),ndof=len(y)-X.shape[1]-2,
                 duration_T14_days=float(duration_from_b_rho(b,rho)),a_over_r=float(a_from_density(rho)),
                 grid_boundary=bool(b<.001 or b>1.998 or abs(dt)>.0199),limb_darkening=limb)
        results.append(rec)
    def fun_free(p):return marginal(stellar_signal(p[0],np.exp(p[1]),p[2],f))[0]
    ordered=sorted([r for r in results if r['binary_flux_ratio']==f],key=lambda r:r['chi2'])[:3]
    fits_=[least_squares(fun_free,[r['impact_parameter'],np.log(r['density_g_cm3']),r['epoch_shift_days']],
                        bounds=([0,np.log(.1),-.02],[1.999,np.log(30),.02]),max_nfev=200,x_scale='jac') for r in ordered]
    fit=min(fits_,key=lambda z:z.fun@z.fun);b,logrho,dt=fit.x;rho=float(np.exp(logrho))
    sig=stellar_signal(b,rho,dt,f)
    rec=dict(kind='free_density',binary_flux_ratio=f,density_g_cm3=rho,impact_parameter=float(b),
             epoch_shift_days=float(dt),chi2=float(fit.fun@fit.fun),ndof=len(y)-X.shape[1]-3,
             duration_T14_days=float(duration_from_b_rho(b,rho)),a_over_r=float(a_from_density(rho)),
             grid_boundary=bool(b<.001 or b>1.998 or rho<.101 or rho>29.9 or abs(dt)>.0199),
             limb_darkening=limb,conditional_covariance_b_logrho_dt=np.linalg.pinv(fit.jac.T@fit.jac).tolist(),
             peak_model_depth_at_reference=float(-min(stellar_signal(b,rho,dt,f))))
    results.append(rec);signals[f]=sig
    print(json.dumps(rec),flush=True)

comparison_planet=next(r for r in previous['fits'] if r['limb_darkening']==limb)
for r in results:r['delta_chi2_vs_same_limb_planet']=r['chi2']-comparison_planet['chi2']
summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             original_circular_source_sha256=previous['method_sha256'],frozen_proposal_sha256=hashlib.sha256((INDEPENDENT/'ephemeris_frozen_pre_holdout.json').read_bytes()).hexdigest(),
             no_holdout_data_accessed=True,period_of_dips_days=period,binary_orbit_days=true_period,
             input_sectors=[12,39,66,93],n_events=len(windows),n_points=len(y),baseline_chi2=float(base_res@base_res),
             compared_planet=comparison_planet,flux_ratio_definition='total equal-star binary reference flux / target TIC reference flux; observed loss=f*blocked_fraction_of_one_star/2',
             assumptions=['equal radii, masses, luminosities and limb-darkening profiles','circular orbit','identical alternate eclipses at true orbit 2P',
                          'source density free over .1–30gcm−3','constant source/reference flux ratio across sectors',
                          'same original QLP windows, errors, quadratic baselines and native exposure integration as planet fit',
                          'no generic aperture contamination multiplier','diagonal noise and unchanged QLP differential normalization'],
             limitations=['No host assignment, companion mass inference, Bayesian evidence or FPP',
                          'Flux ratio is not constrained by pixel throughput/color/differential calibration here',
                          'Equal-star alternatives are a bounded subset of possible stellar models',
                          'Density independent of source stellar colors/masses; physical source correspondence requires external constraints',
                          'Delta chi-square is conditional and not a calibrated model probability'],results=results)
(OUT/'stellar_alternatives.json').write_text(json.dumps(summary,indent=2)+'\n')
pd.DataFrame([{k:v for k,v in r.items() if k not in ['conditional_covariance_b_logrho_dt','limb_darkening']} for r in results]).to_csv(OUT/'stellar_alternatives.csv',index=False)
fig,ax=plt.subplots(figsize=(10,5))
for f in flux_ratios:
    best=next(r for r in results if r['kind']=='free_density' and r['binary_flux_ratio']==f)
    fake_x=np.linspace(-.15,.15,600)
    params.a=a_from_density(best['density_g_cm3']);params.inc=np.degrees(np.arccos(best['impact_parameter']/params.a));params.t0=best['epoch_shift_days']
    model=batman.TransitModel(params,fake_x)
    ax.plot(fake_x*24,(model.light_curve(params)-1)*f/2*1000,label=f'f={f:g}, rho={best["density_g_cm3"]:.2f}, Δχ²={best["delta_chi2_vs_same_limb_planet"]:.1f}')
ax.set_xlabel('Hours from frozen ephemeris');ax.set_ylabel('Differential flux / TIC reference, ppt');ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(OUT/'stellar_alternative_shapes.png',dpi=150);plt.close(fig)
