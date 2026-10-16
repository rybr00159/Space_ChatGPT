"""Conditional photometric characterization, not validation, of TIC 229146531.

Compare a circular dark-companion transit and doubled-period luminous equal-star
grazing binaries. Catalog density constraints are explicitly conditional on a
blended stellar SED; neither hypothesis provides a measured companion mass.
"""
from pathlib import Path
import argparse,json,os,sys
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import batman
from scipy.optimize import least_squares
import multisector_followup as f
from broad_vet import integrated_trapezoid
import matplotlib.pyplot as plt

OUT=f.OUT/'tic229146531'
PRIORS=[dict(name='TIC_unresolved_SED',density=.796092,sigma=.204618),
        dict(name='NSS_distance_scaled_assumption',density=.559,sigma=.205)]

def sources():
    s=json.loads((OUT/'batch02_followup.json').read_text());j=json.loads((OUT/'batch02_discovery_joint.json').read_text())
    ephem=dict(period=j['period_days'],epoch_btjd=j['epoch_btjd'],duration_days=j['duration_days'])
    return s,[x['source']for group in ['discovery','other_batch']for x in s[group]],ephem

def datasets(records,ephem,column):
    result=[]
    for record in records:
        a,m=f.read_source(record);good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['error'])&(a['error']>0)
        stats,local=f.local_quadratic(a['TIME'][good],a[column][good],a['error'][good],ephem['period'],ephem['epoch_btjd'],max(.05,ephem['duration_days']))
        if local is None:continue
        select=abs(f.phase(local['t'],ephem['period'],ephem['epoch_btjd']))<.12
        result.append(dict(sector=record['sector'],source=record,metadata=m,t=local['t'][select],y=local['f'][select],e=local['e'][select],cycle=local['cycle'][select],exposure=float(m['metadata']['TIMEDEL']),stats=stats))
    return result

def fit_model(data,ephem,kind,prior,limb=(.45,.20),fixed_rp=None):
    doubled=kind!='dark_P';period=ephem['period']*(2 if doubled else 1);density_factor=3*np.pi/(6.67430e-8*(period*86400)**2)
    # log(a/R1), b/(1+k), epoch shift, per-sector offsets, optional log(k), log(third light + 1).
    base=3+len(data);free_k=not doubled and fixed_rp is None;diluted=kind=='diluted_equal_binary_2P'
    lower=[np.log(2.05),0,-.012]+[-.01]*len(data);upper=[np.log(30.),.99999,.012]+[.01]*len(data)
    if free_k:lower+=[np.log(.008)];upper+=[np.log(.8)]
    if diluted:lower+=[0];upper+=[np.log(21.)]
    params=batman.TransitParams();params.t0=ephem['epoch_btjd'];params.per=period;params.rp=1 if doubled else.07;params.a=5;params.inc=88;params.ecc=0.;params.w=90.;params.u=list(limb);params.limb_dark='quadratic'
    model1=[batman.TransitModel(params,d['t'],supersample_factor=9,exp_time=d['exposure'])for d in data]
    model2=[batman.TransitModel(params,d['t'],supersample_factor=9,exp_time=d['exposure'])for d in data]if doubled else None
    def evaluate(x):
        rp=1. if doubled else (fixed_rp if fixed_rp is not None else np.exp(x[base]))
        a=np.exp(x[0]);b=x[1]*(1+rp);params.a=a;params.rp=rp;params.inc=np.degrees(np.arccos(np.clip(b/a,0,.999999)));params.t0=ephem['epoch_btjd']+x[2]
        rho=density_factor*a**3/(2 if doubled else 1);third=np.exp(x[-1])-1 if diluted else 0.;models=[]
        for index,m1 in enumerate(model1):
            params.t0=ephem['epoch_btjd']+x[2];value=m1.light_curve(params)
            if doubled:
                params.t0=ephem['epoch_btjd']+x[2]+period/2;value=(value+model2[index].light_curve(params))/2
                value=(value+third)/(1+third)
            models.append(value+x[3+index])
        return models,dict(radius_ratio=rp,a_over_R1=a,impact_parameter=b,epoch_shift_days=x[2],inclination_deg=params.inc,stellar_density_g_cm3=rho,third_light_relative_combined_binary=third,orbital_period_days=period)
    def residual(x):
        model,physical=evaluate(x)
        return np.r_[np.concatenate([(value-d['y'])/d['e']for value,d in zip(model,data)]),(physical['stellar_density_g_cm3']-prior['density'])/prior['sigma']]
    starts=[]
    for fraction in [.5,.85,.97,.995]:
        x=[np.log((prior['density']*(2 if doubled else 1)/density_factor)**(1/3)),fraction,0]+[0]*len(data)
        if free_k:x+=[np.log(.07 if fraction<.97 else .20)]
        if diluted:x+=[np.log(2.)]
        starts.append(x)
    best=min([least_squares(residual,x,bounds=(lower,upper),max_nfev=400)for x in starts],key=lambda r:np.sum(r.fun**2));model,physical=evaluate(best.x)
    photometric=float(np.sum(np.concatenate([(value-d['y'])/d['e']for value,d in zip(model,data)])**2));penalty=((physical['stellar_density_g_cm3']-prior['density'])/prior['sigma'])**2
    cov=np.linalg.pinv(best.jac.T@best.jac)*max(1.,photometric/max(1,sum(len(d['t'])for d in data)-len(best.x)))
    return dict(kind=kind,prior=prior,limb_coefficients=limb,fixed_radius_ratio=fixed_rp,physical=physical,photometric_chi2=photometric,density_penalty=penalty,total_objective=photometric+penalty,n_parameters=len(best.x),n_points=sum(len(d['t'])for d in data),parameters=best.x,diagonal_parameter_covariance=cov,optimizer_success=bool(best.success),optimizer_message=best.message,models=model)

def models_run():
    summary,records,ephem=sources();output=dict(tic=229146531,ephemeris=ephem,script_sha256=f.digest(__file__),batman_version=getattr(batman,'__version__',None),assumptions=['Circular orbits and fixed quadratic limb-darkening coefficients; no measured eccentricity of the short-period source.','TIC density is inferred from an unresolved stellar SED; the known321-dayastrometricbinary invalidates strict single-star treatment.','Alternative density.559gcm−3 scales radius by the NSSdistance ratio; it is an illustrative assumption, not a new stellar measurement.','Doubled-period equal-star model has equal masses/radii/surface brightness; diluted variant allows static extra light.','Models omitspots,reflection,beaming,ellipsoidalphysicsand detailed blending; comparisons are conditional shape tests,not model odds.','Locally fitted baselines and diagonal uncertainties do not calibrate correlated-noise significance.'],fits=[],profiles=[])
    for column in ['corrected','raw']:
        data=datasets(records,ephem,column);f.save(OUT/f'physical_{column}_data.json',dict(ephemeris=ephem,sources=[dict(source=d['source'],sector=d['sector'],metadata=d['metadata'],n_points=len(d['t']),exposure_days=d['exposure'])for d in data]))
        arrays={f's{d["sector"]}_{key}':d[key] for d in data for key in ['t','y','e','cycle']};np.savez_compressed(OUT/f'physical_{column}_arrays.npz',**arrays)
        fits=[]
        for prior in PRIORS:
            for kind in ['dark_P','equal_binary_2P','diluted_equal_binary_2P']:
                fit=fit_model(data,ephem,kind,prior);fit['column']=column;models=fit.pop('models');fits.append((fit,models));output['fits'].append(fit)
                print('PHYSICAL',column,kind,prior['name'],fit['physical'],'chi2',fit['photometric_chi2'],'density_penalty',fit['density_penalty'],flush=True)
                f.save(OUT/'physical_models.json',output)
        # Fixed-k shape profile: other parameters reoptimized; not a radius posterior.
        for rp in [.035,.05,.07,.10,.15,.25,.4,.65]:
            fit=fit_model(data,ephem,'dark_P',PRIORS[0],fixed_rp=rp);fit.pop('models');fit['column']=column;output['profiles'].append(fit)
            f.save(OUT/'physical_models.json',output)
        fig,axes=plt.subplots(len(data),1,figsize=(10,3*len(data)),squeeze=False,layout='constrained')
        for i,d in enumerate(data):
            x=f.phase(d['t'],ephem['period'],ephem['epoch_btjd'])*24;ax=axes[i,0];ax.plot(x,(d['y']-1)*1000,'.',ms=1,color='.65')
            for fit,model in fits[:3]:
                order=np.argsort(x);ax.plot(x[order],(model[i][order]-1)*1000,lw=1,label=fit['kind'])
            ax.set(xlim=(-2,2),xlabel='Hours from common ephemeris',ylabel=column+' (ppt)',title='Sector'+str(d['sector']));ax.legend(fontsize=7)
        fig.suptitle('Circular shape hypotheses under unresolved-SED density assumption; not validation');fig.savefig(OUT/f'physical_models_{column}.png',dpi=150);plt.close(fig)
    f.save(OUT/'physical_models.json',output)

def timing_run():
    summary,records,ephem=sources();rows=[]
    for d in datasets(records,ephem,'raw'):
        for cycle in np.unique(d['cycle']):
            use=d['cycle']==cycle;t=d['t'][use];y=d['y'][use];e=d['e'][use];center=ephem['epoch_btjd']+cycle*ephem['period'];x=t-center
            def residual(par):return (1+par[2]-par[1]*integrated_trapezoid(x-par[0],.045,.019,d['exposure'])-y)/e
            fit=least_squares(residual,[0,.003,0],bounds=([-.025,0,-.01],[.025,.05,.01]),max_nfev=200);cov=np.linalg.pinv(fit.jac.T@fit.jac)*max(1,np.sum(fit.fun**2)/max(1,len(t)-3));err=float(np.sqrt(max(0,cov[0,0])))
            if err>0 and err<.015 and abs(fit.x[0])<.024 and fit.x[1]/np.sqrt(max(1e-20,cov[1,1]))>3:rows.append(dict(sector=d['sector'],cycle=int(cycle),time_btjd=float(center+fit.x[0]),timing_error_days=err,depth=float(fit.x[1])))
    cycle=np.array([r['cycle']for r in rows]);times=np.array([r['time_btjd']for r in rows]);errors=np.array([r['timing_error_days']for r in rows]);design=np.column_stack([np.ones(len(rows)),cycle]);w=1/errors**2
    beta=np.linalg.lstsq(design*np.sqrt(w[:,None]),times*np.sqrt(w),rcond=None)[0];resid=times-design@beta;cov=np.linalg.pinv(design.T@(design*w[:,None]))*max(1,np.sum((resid/errors)**2)/max(1,len(rows)-2))
    rng=np.random.default_rng(229146531);draws=[]
    for _ in range(2000):
        indices=np.concatenate([rng.choice(np.flatnonzero(np.array([r['sector']for r in rows])==sector),size=sum(r['sector']==sector for r in rows),replace=True)for sector in sorted(set(r['sector']for r in rows))]);a=design[indices];ww=w[indices];b=np.linalg.lstsq(a*np.sqrt(ww[:,None]),times[indices]*np.sqrt(ww),rcond=None)[0];draws.append(b)
    draws=np.array(draws);out=dict(tic=229146531,script_sha256=f.digest(__file__),common_discovery_ephemeris=ephem,events=rows,linear_epoch_btjd=float(beta[0]),linear_period_days=float(beta[1]),diagonal_covariance=cov,bootstrap_epoch_period_quantiles=np.quantile(draws,[.00135,.16,.5,.84,.99865],axis=0),n_bootstrap=len(draws),old_sector_tests=[],warning='Event timings use a fixed trapezoid template and selected discovery sectors. Bootstrap resamples events within each sector, not all temporal/systematic errors. Phase-window searches remain conditional and require trials correction; they are not per-sector ephemeris solutions.')
    for r in summary['additional']:
        record=r['source'];arrays,meta=f.read_source(record);flat=f.flatten(arrays,'raw',1.5);reference=float(np.median(flat['t']));n=round((reference-beta[0])/beta[1]);predicted=draws[:,0]+n*draws[:,1];quantiles=np.quantile(predicted,[.00135,.5,.99865]);nominal=beta[0]+n*beta[1];shiftmax=max(abs(quantiles-nominal));offsets=np.linspace(-shiftmax,shiftmax,25);stats=[]
        for shift in offsets:
            value=f.event_stats(flat['t'],flat['f'],flat['e'],float(beta[1]),float(beta[0]+shift),.045,include_details=False);stats.append(dict(shift_days=float(shift),statistics=value))
        out['old_sector_tests'].append(dict(sector=record['sector'],source=record,epoch_prediction_btjd=quantiles,maximum_absolute_phase_offset_days=float(shiftmax),window_tests=stats))
    f.save(OUT/'event_timing_and_phase_windows.json',out);np.savez_compressed(OUT/'timing_bootstrap.npz',epoch_period=draws)
    print('TIMING',len(rows),beta,'COV',cov,'bootstrap',np.quantile(draws,[.16,.5,.84],axis=0),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--timing',action='store_true');p.add_argument('--models',action='store_true');a=p.parse_args()
    if a.timing:timing_run()
    if a.models:models_run()

if __name__=='__main__':main()
