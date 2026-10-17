#!/usr/bin/env python3
"""Independent small forward-model checks prompted by the pinned external review.

Uses our preserved ORIGINAL four-sector QLP windows, not external light curves.
No holdout data, FPP calculation, MCMC or external scripts are executed.
"""
from pathlib import Path
import hashlib,json
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import beta
import batman

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
INP=ROOT/'candidates/tic214311908/independent_variability'
freeze=json.loads((INP/'ephemeris_frozen_pre_holdout.json').read_text());epoch,P=freeze['parameters']
source=INP/'stellar_alternatives/original_window_inputs.npz'
d=np.load(source);y=d['flux'];sigma=d['noise_sigma'];X=d['baseline_design'];time=d['time_btjd'];cycles=d['event_cycle'];sectors=d['sector']
x=time-(epoch+cycles*P);Xw=X/sigma[:,None];pinv=np.linalg.pinv(Xw)
exposures={12:1800/86400,39:600/86400,66:200/86400,93:200/86400}
groups=[np.where(sectors==s)[0] for s in [12,39,66,93]]
def marginal(signal):
    z=(y-signal)/sigma
    return z-Xw@(pinv@z)
def ar(rho,period=P,massfactor=1.):return (massfactor*6.67430e-11*(rho*1000)*(period*86400)**2/(3*np.pi))**(1/3)
params=batman.TransitParams();params.t0=0;params.per=P;params.rp=.12;params.a=ar(9.445379);params.inc=89;params.ecc=0;params.w=90
params.u=[.125,.533];params.limb_dark='quadratic'
models=[batman.TransitModel(params,x[g],supersample_factor=15,exp_time=exposures[s]) for s,g in zip([12,39,66,93],groups)]
def signal(k,b,rho,dt,period=P,scale=1,massfactor=1):
    params.per=period;params.rp=k;params.a=ar(rho,period,massfactor);params.inc=np.degrees(np.arccos(b/params.a));params.t0=dt
    result=np.zeros(len(y))
    for model,g in zip(models,groups):result[g]=(model.light_curve(params)-1)*scale
    return result
fits=[]
for label,kmin,kmax,free_density in [('planet_fixed_density',.01,.175,False),('large_body_fixed_density',.197,.4,False),
                                    ('planet_free_density',.01,.175,True),('large_body_free_density',.197,.4,True)]:
    def f(p):
        k,q,dt=p[:3];rho=np.exp(p[3]) if free_density else 9.445379
        return marginal(signal(k,q*(1+k),rho,dt))
    starts=[]
    for rho,q in [(9.445379,.55),(2,.81),(16,.2)]:
        starts.append([max(kmin+.0001,.12),q,0]+([np.log(rho)] if free_density else []))
    lb=[kmin,0,-.02]+([np.log(.3)] if free_density else []);ub=[kmax,.999,.02]+([np.log(40)] if free_density else [])
    opt=min([least_squares(f,v,bounds=(lb,ub),x_scale='jac',max_nfev=250) for v in starts],key=lambda v:v.fun@v.fun)
    k,q,dt=opt.x[:3];rho=float(np.exp(opt.x[3])) if free_density else 9.445379
    record=dict(label=label,k=float(k),b=float(q*(1+k)),rho_circular_effective_g_cm3=rho,epoch_shift_days=float(dt),
                chi2=float(opt.fun@opt.fun),ndof=len(y)-X.shape[1]-len(opt.x),conditional=True)
    if 'large_body' in label and free_density:
        ratio=rho/9.445379;record['apoastron_e_min_for_TIC_density_ignoring_companion_mass']=float((1-ratio**(2/3))/(1+ratio**(2/3)))
    fits.append(record)

# Directly test the nearest-background source's basic K-star geometry using
# generous flux variation. Stellar parameters come from the pinned source's
# bg_single/best-MS models; all photometry fitted here is our original QLP.
for label,k,rho,massfactor,period,scale_factor in [
        ('background_K_equal_stars_2P',1.,1.408*.733/(.716**3),2.,2*P,.5),
        ('background_K_M_pair_P',.5998148015020005,1.408*.8102561401626047/(.7761792981138234**3),
         (1+.4563031930234208/.8102561401626047),P,1/(1+.07658290186082775))]:
    def f(p):
        q,dt,ratio=p
        return np.r_[marginal(signal(k,q*(1+k),rho,dt,period,ratio*scale_factor,massfactor)),(ratio-.031)/.0031]
    opts=[least_squares(f,[q,0,.031],bounds=([0,-.02,.0217],[.999,.02,.0403]),x_scale='jac',max_nfev=200) for q in [.02,.5,.9]]
    opt=min(opts,key=lambda v:v.fun@v.fun);q,dt,ratio=opt.x
    fits.append(dict(label=label,k=k,b=float(q*(1+k)),rho_primary_g_cm3=rho,orbital_massfactor=massfactor,
                     reference_flux_ratio=float(ratio),epoch_shift_days=float(dt),chi2_with_flux_prior=float(opt.fun@opt.fun),
                     chi2_lightcurve=float(opt.fun[:-1]@opt.fun[:-1]),ndof_lightcurve=len(y)-X.shape[1]-3,
                     reference_flux_ratio_bound=bool(ratio>.0402 or ratio<.0218),
                     limitation='tests only circular specified-source geometries; excludes no general EB population; no external photometric/secondary penalty added'))

# Photo-eccentric sign via actual numerical model contact durations.
photoecc=[]
grid=np.linspace(-.16,.16,32001)
params.t0=0;params.per=P;params.rp=.12;params.a=ar(9.445379);params.inc=90;params.u=[0,0]
for e,w in [(0,90),(.5,90),(.5,270)]:
    params.ecc=e;params.w=w
    flux=batman.TransitModel(params,grid).light_curve(params)
    t=grid[flux<1-1e-10];dur=t.max()-t.min();photoecc.append(dict(e=e,w_degrees=w,duration_days=float(dur)))
circ=photoecc[0]['duration_days']
for v in photoecc:
    g=(1+v['e']*np.sin(np.radians(v['w_degrees'])))/np.sqrt(1-v['e']**2)
    v['density_ratio_from_duration_approx']=(circ/v['duration_days'])**3;v['g_cubed']=float(g**3)

out=dict(created_utc=datetime.now(timezone.utc).isoformat(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),n_points=len(y),n_events=len(set(cycles)),
         note='Original13QLPevents only; new comparisons prompted by external findings; diagonal errors, no FPP; fixed LD .125,.533',
         fits=fits,photoeccentric_sign_checks=photoecc,
         zero_exceedance_95pct_upper_bounds={str(n):float(beta.ppf(.95,1,n)) for n in [25,30,60]})
(OUT/'shape_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
