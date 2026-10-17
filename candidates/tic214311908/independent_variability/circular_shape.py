#!/usr/bin/env python3
"""Conditional circular-transit shape test, not validation or source assignment."""
from pathlib import Path
import json, hashlib
from datetime import datetime, timezone
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
import batman
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
frozen=json.loads((OUT/'ephemeris_frozen_pre_holdout.json').read_text())
epoch,period=frozen['parameters']
exofop=json.loads((ROOT/'data/tic214311908/catalogs/exofop.json').read_text())
stellar=next(x for x in exofop['stellar_parameters'] if 'mass' in x)
mass=float(stellar['mass']); radius=float(stellar['srad'])
rho=float(stellar['dens'])
rho_solar=1.408
a_over_r=(6.67430e-11*(rho*1000)*(period*86400)**2/(3*np.pi))**(1/3)
windows=[]
for rec in json.loads((ROOT/'data/tic214311908/observations.json').read_text()):
    with fits.open(rec['path']) as f:
        d=f[1].data; name='DET_FLUX' if 'DET_FLUX' in d.names else 'KSPSAP_FLUX'
        m=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d[name])
        t=np.array(d['TIME'])[m]; y=np.array(d[name])[m]
    exposure=float(np.median(np.diff(t)[np.diff(t)<.05]))
    for n in range(int(np.floor((t.min()-epoch)/period)),int(np.ceil((t.max()-epoch)/period))+1):
        mid=epoch+n*period; x=t-mid
        use=abs(x)<.3
        if np.sum(abs(x)<.025)<1 or use.sum()<12: continue
        tx=t[use]; z=y[use]/np.median(y[use])-1; xx=x[use]
        oot=abs(xx)>.1
        A=np.c_[np.ones(len(xx)),xx,xx*xx]
        beta=np.linalg.lstsq(A[oot],z[oot],rcond=None)[0]
        res=z-A@beta
        sigma=1.4826*np.median(abs(res[oot]-np.median(res[oot])))
        windows.append(dict(sector=rec['sector'],cycle=n,t=tx,y=z,x=xx,X=A,sigma=sigma,exposure=exposure))
ys=np.concatenate([w['y'] for w in windows]); sigmas=np.concatenate([np.repeat(w['sigma'],len(w['t'])) for w in windows])
X=np.zeros((len(ys),len(windows)*3)); at=0
for j,w in enumerate(windows):
    X[at:at+len(w['t']),j*3:j*3+3]=w['X'];at+=len(w['t'])
Xw=X/sigmas[:,None]
projection=np.linalg.pinv(Xw)
def marginal_residual(signal):
    z=(ys-signal)/sigmas
    beta=projection@z
    return z-Xw@beta,beta
base_res,_=marginal_residual(np.zeros(len(ys)))
result=dict(created_utc=datetime.now(timezone.utc).isoformat(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            fixed_period_days=period,fixed_epoch_btjd=epoch,a_over_r_fixed=a_over_r,
            TIC_mass_solar=mass,TIC_radius_solar=radius,TIC_density_g_cm3=rho,
            assumptions=['signal is on target','circular dark companion','QLP relative-depth normalization unchanged by unknown dilution',
                         'fixed TIC stellar density','fixed quadratic limb darkening tested below','independent diagonal point noise',
                         'one common radius ratio and impact parameter across all sectors'],
            limitation='This exercise tests whether a circular planet-shaped transit is geometrically possible; it does not assign the source, determine companion mass, calibrate false-alarm probability, or resolve blends.',
            n_events=len(windows),n_points=len(ys),baseline_chi2=float(base_res@base_res),fits=[])
best=None
for limb in [[.35,.25],[.5,.2],[.2,.2]]:
    params=batman.TransitParams();params.t0=epoch;params.per=period;params.rp=.13;params.a=a_over_r;params.inc=89
    params.ecc=0;params.w=90;params.u=limb;params.limb_dark='quadratic'
    models=[batman.TransitModel(params,w['t'],supersample_factor=15,exp_time=w['exposure']) for w in windows]
    def signal(p):
        k,q,dt=p
        b=q*(1+k)
        params.rp=k;params.inc=np.degrees(np.arccos(b/a_over_r));params.t0=epoch+dt
        return np.concatenate([model.light_curve(params)-1 for model in models])
    def objective(p):
        return marginal_residual(signal(p))[0]
    starts=[[.13,.55,0],[.18,.80,0],[.35,.94,0]]
    opts=[least_squares(objective,s,bounds=([.01,0,-.02],[.8,.999,.02]),x_scale='jac',max_nfev=250) for s in starts]
    fit=min(opts,key=lambda v:np.sum(v.fun**2)); k,q,dt=fit.x; b=q*(1+k)
    cov=np.linalg.pinv(fit.jac.T@fit.jac)
    duration=period/np.pi*np.arcsin(np.sqrt((1+k)**2-b*b)/np.sqrt(a_over_r*a_over_r-b*b))
    rec=dict(limb_darkening=limb,radius_ratio=float(k),impact_parameter=float(b),epoch_shift_days=float(dt),
             duration_T14_days=float(duration),radius_Earth_if_on_target_and_no_extra_dilution=float(k*radius*109.076),
             chi2=float(fit.fun@fit.fun),ndof=len(ys)-len(windows)*3-3,
             conditional_covariance_k_q_dt=cov.tolist(),parameter_at_bound=bool(k>.79 or q>.998 or abs(dt)>.019))
    result['fits'].append(rec)
    if best is None or rec['chi2']<best[0]:
        _,beta=marginal_residual(signal(fit.x))
        # X@beta is a flux baseline because Xw=X/sigma.
        best=(rec['chi2'],fit.x.copy(),signal(fit.x),X@beta,limb)
(OUT/'circular_shape_assessment.json').write_text(json.dumps(result,indent=2)+'\n')
fig,axs=plt.subplots(4,1,figsize=(10,10),sharex=True)
at=0
for w in windows:
    n=len(w['t']); base=best[3][at:at+n];sig=best[2][at:at+n];z=w['y']-base
    ax=axs[[12,39,66,93].index(w['sector'])]
    ax.plot(w['x']*24,z*1000,'.',ms=1,alpha=.35)
    order=np.argsort(w['x']);ax.plot(w['x'][order]*24,sig[order]*1000,lw=.6,color='red')
    ax.set_ylabel(f"S{w['sector']} ppt");at+=n
for ax in axs:ax.set_xlim(-3,3);ax.set_ylim(-60,45)
axs[-1].set_xlabel('Hours from frozen ephemeris')
fig.tight_layout();fig.savefig(OUT/'conditional_circular_shape.png',dpi=150);plt.close(fig)
print(json.dumps(result,indent=2))
