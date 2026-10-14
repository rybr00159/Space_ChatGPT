"""Independent raw-pixel extraction and physical eclipse tests of TIC273879360.

Two luminous spheres use separate quadratic limb-darkened BATMAN occultations,
sharing a Keplerian orbit. Pixel-core third light is free by sector. Local
polynomial baselines are profiled simultaneously, not inherited from the old
trapezoid detrending. This model cannot determine stellar masses or unique radii.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
import batman
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares, brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import ROOT, verify_fits_checksums
from broad_vet import mad
from pixel_vet_49617264 import safe

OUT=ROOT/'candidates/stellar_systems/tic273879360_validation'
P0=7.5868912871546215
T0=1327.8999181878462
PHASE=.6917919273084346
SECTORS=[27,29,67,95]
RA,DEC=32.25246282511974,-71.32317034350537


def phase_secondary(e,w):
    f=np.array([np.pi/2-w,3*np.pi/2-w])
    E=2*np.arctan2(np.sqrt(1-e)*np.sin(f/2),np.sqrt(1+e)*np.cos(f/2))
    M=E-e*np.sin(E)
    return float((M[1]-M[0])%(2*np.pi)/(2*np.pi))


def input_folders():
    return [ROOT/f'candidates/stellar_systems/tic273879358/pixels/sector{s}/large_source_validation'
            if s<30 else ROOT/f'candidates/stellar_systems/tic273879360/pixels/sector{s}/later_source_validation'
            for s in SECTORS]


def extract(radius=1., degree=2, phase_centers=None):
    events=[];sources=[];exclusions=[]
    if phase_centers is None:
        phase_centers=[('primary',0.),('secondary',PHASE)]
    for sector,folder in zip(SECTORS,input_folders()):
        row=json.loads((folder/'raw_manifest.json').read_text());path=Path(row['path'])
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('Native raw SHA mismatch')
        with fits.open(path)as hd:
            checks=verify_fits_checksums(path,hd)
            h,a=hd[1].header,hd[1].data
            if hd[0].header['TIMESYS']!='TDB':raise ValueError('Unexpected time system')
            wcs=WCS(hd[2].header);xy=wcs.all_world2pix([[RA,DEC]],0)[0]
            yy,xx=np.indices(a['FLUX'].shape[1:]);mask=(xx-xy[0])**2+(yy-xy[1])**2<radius**2
            good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.all(np.isfinite(a['FLUX'])&np.isfinite(a['FLUX_ERR']),axis=(1,2))
            t=np.asarray(a['TIME'][good],float);cube=np.asarray(a['FLUX'][good],float);err=np.asarray(a['FLUX_ERR'][good],float)
            # A single low-background mask, selected without eclipse times.
            reference=np.median(cube,axis=0);blank=reference<np.percentile(reference,25)
            sky=np.median(cube[:,blank],axis=1)
            flux=np.sum(cube[:,mask],axis=1)-mask.sum()*sky
            error=np.sqrt(np.sum(err[:,mask]**2,axis=1))
            exposure=float(h['TIMEDEL'])
        scale=float(np.median(flux));flux/=scale;error/=scale
        for kind,phase in phase_centers:
            epoch=T0+phase*P0
            for cycle in range(int(np.ceil((t.min()-epoch)/P0)),int(np.floor((t.max()-epoch)/P0))+1):
                center=epoch+cycle*P0;use=abs(t-center)<.48
                x=t[use]-center;y=flux[use];e=error[use]
                outside=abs(x)>.16
                reason=None
                if len(x)<12:reason='few_cadences'
                elif (outside&(x<0)).sum()<8 or (outside&(x>0)).sum()<8:reason='baseline_coverage'
                elif (abs(x)<.09).sum()<3:reason='eclipse_coverage'
                elif np.max(np.diff(x))>max(3.05*exposure,.065):reason='quality_gap'
                if reason:
                    exclusions.append(dict(sector=sector,kind=kind,cycle=cycle,reason=reason));continue
                design=np.column_stack([x**k for k in range(degree+1)])
                baseline=np.linalg.lstsq(design[outside],y[outside],rcond=None)[0]
                floor=max(mad(y[outside]-design[outside]@baseline),float(np.median(e)))
                e=np.maximum(e,floor)
                events.append(dict(sector=sector,sector_index=SECTORS.index(sector),kind=kind,cycle=cycle,
                                   center=center,t=t[use],x=x,y=y,e=e,exposure=exposure,
                                   baseline_design=design,noise_floor=floor))
        sources.append(dict(**row,analysis_aperture_radius_pixels=radius,n_pixels=int(mask.sum()),
                            raw_checksum_verification=checks,normalization_flux_e_s=scale))
    return events,sources,exclusions


class BinaryModel:
    def __init__(self,events,supersample=15,limb1=(.45,.25),limb2=(.60,.15),period_base=P0):
        self.events=events;self.limb1=limb1;self.limb2=limb2;self.supersample=supersample
        self.period_base=period_base
        self.t=np.concatenate([v['t']for v in events]);self.sindex=np.concatenate([np.full(len(v['t']),v['sector_index'],int)for v in events])
        self.exposure=np.concatenate([np.full(len(v['t']),v['exposure'])for v in events])
        self.models={}
        for exp in np.unique(self.exposure):
            use=self.exposure==exp
            a=self.parameters(.6,18,87,.35,335,T0,P0,self.limb1)
            b=self.parameters(1/.6,18/.6,87,.35,155,T0+P0*.69,P0,self.limb2)
            self.models[exp]=(use,batman.TransitModel(a,self.t[use],supersample_factor=supersample,exp_time=exp),
                                   batman.TransitModel(b,self.t[use],supersample_factor=supersample,exp_time=exp))
        self.n=len(self.t)

    @staticmethod
    def parameters(k,alpha,inc,e,w,epoch,p,limb):
        a=batman.TransitParams();a.t0=epoch;a.per=p;a.rp=k;a.a=alpha;a.inc=inc;a.ecc=e;a.w=w;a.u=list(limb);a.limb_dark='quadratic'
        return a

    def flux(self,v,dark=False,circular=False):
        k,z,logalpha,ec,es,logq,shift,dp=v[:8]
        alpha=np.exp(logalpha);e=0 if circular else np.hypot(ec,es);w=0 if circular else np.arctan2(es,ec)
        if e>=.85:return None,dict(eccentricity=float(e),invalid=True)
        b=z*(1+k);cosinc=b*(1+e*np.sin(w))/(alpha*(1-e*e))
        if cosinc>=1 or alpha*(1-e)<=1+k:return None,dict(eccentricity=float(e),invalid=True)
        inc=np.degrees(np.arccos(cosinc));p=self.period_base+dp;epoch=T0+shift
        phase=phase_secondary(e,w);q=0 if dark else np.exp(logq)
        first=self.parameters(k,alpha,inc,e,np.degrees(w)%360,epoch,p,self.limb1)
        second=self.parameters(1/k,alpha/k,inc,e,(np.degrees(w)+180)%360,epoch+phase*p,p,self.limb2)
        f=np.empty(self.n)
        for exp,(use,m1,m2)in self.models.items():
            f[use]=(m1.light_curve(first)+q*m2.light_curve(second))/(1+q)
        dilution=np.asarray(v[8:12])[self.sindex]
        f=1-dilution*(1-f)
        geometry=dict(radius_ratio=float(k),a_over_r1=float(alpha),impact_parameter_primary=float(b),inclination_deg=float(inc),
                      eccentricity=float(e),omega_deg=float(np.degrees(w)%360),secondary_phase=float(phase),
                      light_ratio_tess=float(q),surface_brightness_ratio=float(q/(k*k)),
                      period_days=float(p),epoch_btjd=float(epoch),sector_binary_light_fraction=dict(zip(map(str,SECTORS),map(float,v[8:12]))),
                      periastron_over_r1=float(alpha*(1-e)))
        return f,geometry

    def residual(self,v,dark=False,circular=False,return_info=False):
        model,geometry=self.flux(v,dark,circular)
        if model is None:
            return np.full(self.n,1e5)
        residual=[];predictions=[];baselines=[];offset=0
        for event in self.events:
            length=len(event['t']);m=model[offset:offset+length];offset+=length
            design=event['baseline_design']*m[:,None]
            a=design/event['e'][:,None];y=event['y']/event['e']
            beta=np.linalg.lstsq(a,y,rcond=None)[0];pred=design@beta
            residual.extend((event['y']-pred)/event['e']);predictions.extend(pred);baselines.extend(event['baseline_design']@beta)
        if return_info:return np.asarray(residual),geometry,np.asarray(predictions),np.asarray(baselines)
        return np.asarray(residual)


LOW=np.array([.15,.01,np.log(5),-.65,-.65,np.log(.03),-.04,-6e-5,.35,.35,.35,.35])
HIGH=np.array([1.4,.98,np.log(38),.65,.65,np.log(1.0),.04,6e-5,1.,1.,1.,1.])
SCALE=np.array([.3,.3,.5,.2,.2,.7,.003,1e-5,.2,.2,.2,.2])


def fit_model(model,seed,dark=False,circular=False,max_nfev=180):
    active=np.ones(12,bool)
    if dark:active[5]=False
    if circular:active[3:5]=False
    fixed=seed.copy()
    if circular:fixed[3:5]=0
    def expand(v):
        full=fixed.copy();full[active]=v;return full
    fit=least_squares(lambda v:model.residual(expand(v),dark,circular),fixed[active],
                      bounds=(LOW[active],HIGH[active]),x_scale=SCALE[active],diff_step=1e-4,max_nfev=max_nfev)
    full=expand(fit.x);residual,geometry,prediction,baseline=model.residual(full,dark,circular,True)
    return dict(parameters=full.tolist(),geometry=geometry,chi2=float(residual@residual),success=bool(fit.success),nfev=fit.nfev,
                model_type='dark'if dark else('circular_binary'if circular else'eccentric_luminous_binary'))


def main():
    p=argparse.ArgumentParser();p.add_argument('--starts',type=int,default=8);p.add_argument('--supersample',type=int,default=15)
    a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    events,sources,excluded=extract();model=BinaryModel(events,a.supersample)
    seeds=[]
    for k,z,alpha,wdeg in [(.65,.52,18,-25),(.5,.45,16,-15),(.8,.62,20,-30),(.95,.65,22,-10),(.55,.5,18,10),(.7,.45,16,-40),(.9,.60,24,20),(.4,.2,15,0)][:a.starts]:
        w=np.radians(wdeg);e=brentq(lambda e:phase_secondary(e,w)-PHASE,.01,.80)
        seeds.append(np.array([k,z,np.log(alpha),e*np.cos(w),e*np.sin(w),np.log(.2),0,0,.9,.75,.86,.92]))
    trials=[]
    for i,seed in enumerate(seeds):
        result=fit_model(model,seed);trials.append(result)
        (OUT/'fit_checkpoint.json').write_text(json.dumps(safe(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),trials=trials)),indent=2)+'\n')
        print('Start',i,'chi2',result['chi2'],'geometry',result['geometry'],flush=True)
    best=min(trials,key=lambda v:v['chi2']);v=np.asarray(best['parameters'])
    dark=fit_model(model,v,dark=True);circular=fit_model(model,v,circular=True)
    residual,geometry,pred,baseline=model.residual(v,return_info=True)
    emin=brentq(lambda e:phase_secondary(e,0)-PHASE,0,.8)
    result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),tic=273879360,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sources=sources,excluded_events=excluded,
                n_events=len(events),n_points=model.n,supersample=a.supersample,best_fit=best,
                alternative_starts=trials,dark_companion_fit=dark,circular_binary_fit=circular,
                conditional_phase_only_eccentricity_lower_bound=emin,
                events=[{k:v for k,v in event.items()if k in ['sector','kind','cycle','center','exposure','noise_floor']}for event in events],
                limitations=['Model assumes one two-body eclipsing system on the identified source.',
                             'Baseline polynomials, constant sector dilution and approximate fixed limb darkening omit PRF/calibration effects.',
                             'Dimensionless radius and brightness fits are degenerate; no component masses or unique absolute radii are established.',
                             'Chi-square differences are conditional model diagnostics, not calibrated false-positive or discovery probabilities.'])
    (OUT/'physical_fit.json').write_text(json.dumps(safe(result),indent=2)+'\n')
    np.savez_compressed(OUT/'raw_model_arrays.npz',time=model.t,model=pred,baseline=baseline,standardized_residual=residual,
                        raw_flux=np.concatenate([e['y']for e in events]),error=np.concatenate([e['e']for e in events]))
    fig,axes=plt.subplots(2,2,figsize=(12,8));offset=0;colors={27:'tab:blue',29:'tab:orange',67:'tab:green',95:'tab:red'}
    for event in events:
        length=len(event['t']);use=slice(offset,offset+length);offset+=length;i=0 if event['kind']=='primary'else 1
        x=event['x']*24
        axes[0,i].scatter(x,(1-event['y']/baseline[use])*100,s=2,alpha=.25,color=colors[event['sector']])
        axes[0,i].plot(x,(1-pred[use]/baseline[use])*100,color=colors[event['sector']],alpha=.7)
        axes[1,i].scatter(x,residual[use],s=2,alpha=.25,color=colors[event['sector']])
    for i in [0,1]:
        axes[0,i].set(xlim=(-5,5),ylabel='Raw aperture flux loss (%)',title='Primary'if i==0 else'Secondary')
        axes[1,i].set(xlim=(-5,5),xlabel='Hours from old predicted eclipse',ylabel='Residual / empirical error');axes[1,i].axhline(0,color='grey')
    fig.suptitle('TIC 273879360: two luminous spheres, eccentric shared orbit');fig.tight_layout();fig.savefig(OUT/'physical_fit.png',dpi=170);plt.close(fig)
    print('Best',best,'dark_delta',dark['chi2']-best['chi2'],'circular_delta',circular['chi2']-best['chi2'],flush=True)


if __name__=='__main__':main()
