"""Falsification/sensitivity tests of the TIC273879360 physical fit.

Raw extraction alternatives, limb-darkening/exposure checks and sector-heldout
predictions are diagnostic tests, not a calibrated discovery false-alarm model.
"""
from pathlib import Path
from datetime import datetime, timezone
import copy
import hashlib
import json
import numpy as np
from scipy.optimize import least_squares
import validate_background_binary as vb
from pixel_vet_49617264 import safe


def fixed_geometry_fit(model,seed):
    # Sector light fractions are the nuisance extraction quantities; the
    # geometrical/eclipse timing hypothesis is held fixed for predictions.
    def residual(dil):
        p=seed.copy();p[8:12]=dil
        return model.residual(p)
    fit=least_squares(residual,seed[8:12],bounds=(vb.LOW[8:12],vb.HIGH[8:12]),max_nfev=50)
    p=seed.copy();p[8:12]=fit.x
    r,g,_,_=model.residual(p,return_info=True)
    return dict(chi2=float(r@r),n_points=model.n,geometry=g)


def main():
    best=json.loads((vb.OUT/'physical_fit.json').read_text())
    seed=np.asarray(best['best_fit']['parameters'])
    events,_,_=vb.extract();reference=vb.BinaryModel(events,15)
    tests={}
    for name,radius,degree,limb1,limb2 in [
        ('linear_baseline',1.,1,(.45,.25),(.60,.15)),
        ('cubic_baseline',1.,3,(.45,.25),(.60,.15)),
        ('aperture_0.8pixel',.8,2,(.45,.25),(.60,.15)),
        ('aperture_1.2pixel',1.2,2,(.45,.25),(.60,.15)),
        ('weaker_limb_darkening',1.,2,(.30,.20),(.40,.20)),
        ('stronger_limb_darkening',1.,2,(.60,.20),(.75,.10))]:
        ev,_,_=vb.extract(radius,degree)
        model=vb.BinaryModel(ev,15,limb1,limb2)
        result=vb.fit_model(model,seed,max_nfev=180)
        result.update(n_events=len(ev),n_points=model.n,radius_pixels=radius,baseline_degree=degree,limb1=limb1,limb2=limb2)
        tests[name]=result
        print(name,result['chi2'],result['geometry'],flush=True)
    high=vb.BinaryModel(events,61)
    f1,_=reference.flux(seed);f2,_=high.flux(seed)
    exposure=dict(max_absolute_model_difference_ppm=float(np.max(abs(f1-f2))*1e6),
                  rms_model_difference_ppm=float(np.sqrt(np.mean((f1-f2)**2))*1e6),
                  supersamples=[15,61],note='Two physically identical orbits integrated with different midpoint samples.')
    # Train on one epoch and evaluate geometry in the other. Extracted light
    # fractions and local baselines are free; eclipse shape/timing is held fixed.
    holdouts={}
    for name,train_sectors,test_sectors in [('early_to_later',[27,29],[67,95]),('later_to_early',[67,95],[27,29])]:
        train=[e for e in events if e['sector']in train_sectors]
        test=[e for e in events if e['sector']in test_sectors]
        tr=vb.BinaryModel(train,15);fit=vb.fit_model(tr,seed,max_nfev=180)
        p=np.asarray(fit['parameters']);te=vb.BinaryModel(test,15)
        prediction=fixed_geometry_fit(te,p)
        # For comparison, a fully refitted test model quantifies failures of
        # extrapolation. A lower chi-square is not independent evidence.
        refit=vb.fit_model(te,seed,max_nfev=180)
        holdouts[name]=dict(training_fit=fit,fixed_geometry_prediction=prediction,test_refit=refit,
                            conditional_prediction_penalty=prediction['chi2']-refit['chi2'])
        print('Holdout',name,holdouts[name]['conditional_prediction_penalty'],flush=True)
    # Sector/error cluster estimates are not extreme-tail confidence intervals;
    # bootstrap residual blocks retaining simultaneous baseline fits.
    residual,_,pred,baseline=reference.residual(seed,return_info=True)
    rng=np.random.default_rng(273879360);bootstrap=[];offsets=np.cumsum([0]+[len(e['t'])for e in events])
    for draw in range(24):
        synthetic=[]
        for i,event in enumerate(events):
            ev=copy.deepcopy(event);idx=slice(offsets[i],offsets[i+1]);r=residual[idx]*event['e']
            # 30-minute synchronized blocks at native cadences, circular within
            # each local event window. True inter-event calibration is omitted.
            block=max(1,int(round(.020833333333/event['exposure'])))
            starts=rng.integers(0,len(r),(len(r)+block-1)//block)
            indices=np.concatenate([(s+np.arange(block))%len(r)for s in starts])[:len(r)]
            ev['y']=pred[idx]+r[indices];synthetic.append(ev)
        model=vb.BinaryModel(synthetic,15);fit=vb.fit_model(model,seed,max_nfev=100)
        bootstrap.append(fit['geometry'])
    ranges={key:np.percentile([r[key]for r in bootstrap],[16,50,84]).tolist()
            for key in ['eccentricity','omega_deg','radius_ratio','light_ratio_tess','a_over_r1','inclination_deg','period_days','secondary_phase']}
    report=dict(checked_utc=datetime.now(timezone.utc).isoformat(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                physical_method_sha256=hashlib.sha256(Path(vb.__file__).read_bytes()).hexdigest(),
                extraction_and_limb_tests=tests,exposure_integration_check=exposure,
                epoch_holdouts=holdouts,residual_block_bootstrap=dict(n=24,seed=273879360,interval16_50_84=ranges,draws=bootstrap,
                     limitation='Conditional30minute residual blocks inside events; not total model/stellar/PRF uncertainty or mass confidence intervals.'),
                limitations=['Free nuisance dilution reduces sensitivity to camera amplitude changes.',
                             'Small-number event blocks and shared photons limit statistical calibration.',
                             'No independent mass or definitive novelty verification is supplied by these controls.'])
    (vb.OUT/'physical_controls.json').write_text(json.dumps(safe(report),indent=2)+'\n')
    print('Bootstrap',ranges,flush=True)


if __name__=='__main__':main()
