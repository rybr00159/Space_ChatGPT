"""Physical shape/duration falsification of TIC260820871's shallow recurrence.

Frozen earlier two-lead code is imported with an exact byte hash. No planet
mass or source attribution is inferred from this conditional shape exercise.
"""
from __future__ import annotations
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from broad_vet import integrated_trapezoid
from fit_370247818 import TransitFit, sha, save, RJUP_RSUN

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT/'src/fit_broad_leads.py'
EXPECTED = 'b27aed799c4c72ec46e89b75fd3f527be300e85d91668f729dcda20ffe99268c'
if sha(FROZEN) != EXPECTED:
    raise ValueError('Previously frozen likelihood source changed; do not silently reuse')
from fit_broad_leads import verify_source, combine, EqualBinaryFit
TIC = 260820871
INPUT = ROOT/'candidates/broad_s1/tic260820871_multisector'
OUT = ROOT/'candidates/broad_s1/tic260820871_physical'


class HotPrimaryFit(TransitFit):
    bounds = (np.array([.0005, .0001, np.log(1.001), .001, .001, -.035, -35.]),
              np.array([.6, .9999, np.log(20.), .999, .999, .035, 35.]))

    def residual(self, vector, include_priors=True):
        residual = super().residual(vector, include_priors)
        if not include_priors:
            return residual
        k, bf, loga = vector[:3]
        a, b = np.exp(loga), bf*(1+k)
        # A strong one-sided continuous penalty prevents photosphere collision
        # and nonphysical inclination. Its final value is retained separately.
        safe = a*(1-self.eccentricity)-(1+k)
        cosine = b/a*(1+self.eccentricity*np.sin(np.deg2rad(self.omega)))/(1-self.eccentricity**2)
        return np.r_[residual, max(0., .001-safe)/.0003, max(0., cosine-.999999)/.0003]


def inputs():
    record = json.loads((INPUT/'multisector.json').read_text())
    timing = record['event_timing_fit']['fitted_ephemeris']
    archive = np.load(INPUT/'multisector_arrays.npz')
    chunks, sources = [], []
    for source in record['sources']:
        verified = verify_source(source); sources.append(verified)
        s = source['sector']
        t, f, e = (archive[f's{s}_local_{name}'] for name in ('t', 'f', 'e'))
        cycle = np.round((t-timing['epoch_btjd'])/timing['period']).astype(int)
        phase = t-timing['epoch_btjd']-cycle*timing['period']
        good = np.isfinite(t)&np.isfinite(f)&np.isfinite(e)&(e>0)&(abs(phase)<.16)
        chunks.append(dict(t=t[good], f=f[good], e=e[good], cycle=cycle[good],
                           sector=np.full(good.sum(), s, int), exposure=np.full(good.sum(), verified['exposure_days_exact_TIMEDEL'])))
    star = {'radius_rsun': 2.40763, 'radius_error_rsun': .2*2.40763,
            'mass_msun': 2.85, 'mass_error_msun': .2*2.85,
            'teff_K': 10942., 'teff_error_K': 1000.,
            'catalog_rho_solar': 2.85/2.40763**3,
            'catalog_rho_error_solar': .65*2.85/2.40763**3,
            'limb_darkening_u1_mean': .20, 'limb_darkening_u1_sigma': .15,
            'limb_darkening_u2_mean': .15, 'limb_darkening_u2_sigma': .20,
            'stellar_values_source': 'Current TIC source check by parent/catalog agent; cached source linked when available.',
            'uncertainty_warning': '20% R/M and1000K widths are conservative scenario assumptions, not reported catalog uncertainties.',
            'limb_darkening_source': 'Broad approximate hot-star TESS priors, not an atmosphere-grid calculation.',
            'catalog_contaminant_to_target_flux_ratio': .00246}
    provenance = {'raw_sources': sources, 'input_array_sha256': sha(INPUT/'multisector_arrays.npz'),
                  'input_record_sha256': sha(INPUT/'multisector.json'),
                  'frozen_likelihood_source_sha256': EXPECTED,
                  'sectors': sorted(set(s['sector'] for s in sources)),
                  'method': 'Unbinned bilateral OOT-quadratic local arrays, +/-0.16d around every recurrence; exact raw TIMEDEL exposures.'}
    return combine(chunks), star, timing, provenance


def fitted_primary(data, star, timing, label, eccentricity=0., omega=90., density=True, kmax=.6):
    model = HotPrimaryFit(data, star, timing, eccentricity=eccentricity, omega=omega,
                          density_prior=density, supersample=31)
    model.bounds = (model.bounds[0].copy(), model.bounds[1].copy())
    model.bounds[1][0] = kmax
    trials = []
    for k, impact in ((.025, .7), (.04, .95), (.08, 1.04), (.25, 1.22)):
        k = min(k, kmax*.95)
        a = max(np.exp(model.loga_mean), (1+k+.01)/(1-eccentricity))
        factor = (1+eccentricity*np.sin(np.deg2rad(omega)))/(1-eccentricity**2)
        b = min(impact, .97*(1+k), .95*a/factor)
        seed = np.array([k, b/(1+k), np.log(a), .35**2, .30, 0., 0.])
        trials.append(model.fit(seed, max_nfev=160))
    vector, result = min(trials, key=lambda pair: np.sum(pair[1].fun**2))
    row = model.description(vector, result)
    margin = row['a_over_rstar']*(1-eccentricity)-(1+row['radius_ratio'])
    row.update(scenario=label, periastron_above_photospheres_Rstar=margin,
               photosphere_condition_satisfied=bool(margin>0),
               prior_and_physical_penalties_chi2=float(np.sum(model.residual(vector)**2)-row['data_chi2']))
    prediction = model.physical_model(vector)
    return row, prediction


def fit_trapezoid(data, timing):
    cycle = np.round((data['t']-timing['epoch_btjd'])/timing['period']).astype(int)
    phase = data['t']-timing['epoch_btjd']-cycle*timing['period']
    sectors, si = np.unique(data['sector'], return_inverse=True)
    w = 1/data['e']**2
    sw = np.bincount(si, weights=w)
    parity = cycle%2
    def prediction(v):
        de, do, duration, fraction, shift, dp = v
        template = integrated_trapezoid(phase-shift-cycle*dp*1e-5, duration, fraction*duration, data['exposure'])
        return 1-np.where(parity==0, de, do)*template
    def residual(v):
        pred = prediction(v)
        offsets = np.bincount(si, weights=w*(data['f']-pred))/sw
        return (data['f']-pred-offsets[si])/data['e']
    result = least_squares(residual, [.0006,.0004,.088,.398,0.,0.],
                           bounds=([0,0,.02,.02,-.035,-35],[.01,.01,.2,.499,.035,35]),
                           x_scale='jac',max_nfev=200)
    v = result.x
    row = {'even_depth': float(v[0]), 'odd_depth': float(v[1]),
           'duration_T14_hours': float(v[2]*24), 'ingress_fraction': float(v[3]),
           'epoch_btjd': float(timing['epoch_btjd']+v[4]),
           'photometric_period_days': float(timing['period']+v[5]*1e-5),
           'data_chi2': float(np.sum(result.fun**2)), 'converged': bool(result.success),
           'model': 'Free alternating-depth exposure-integrated trapezoid, empirical comparator only.'}
    return row, prediction(v)


def diluted_binary(data, star, timing):
    model = EqualBinaryFit(data, star, timing)
    model.use_density_prior = False
    # Detached equal-radius Roche lobes require a/R>~2.64;3is conservative.
    lo, hi = model.bounds[0][1:].copy(), model.bounds[1][1:].copy()
    lo[1], hi[1] = np.log(3.), np.log(30.)
    bounds = (np.r_[lo, .05, 1e-5], np.r_[hi, 10., 1.])
    parity = np.round((data['t']-timing['epoch_btjd'])/timing['period']).astype(int)%2
    def prediction(v):
        curve = model.physical_model(np.r_[1.,v[:-2]])
        j, light = v[-2:]
        factor = np.where(parity==0,2/(1+j),2*j/(1+j))
        return 1-(1-curve)*factor*light
    def residual(v):
        pred = prediction(v)
        offsets = np.bincount(model.sector_idx, weights=model.weight*(data['f']-pred))/model.sector_weight
        res = (data['f']-pred-offsets[model.sector_idx])/data['e']
        # A faint unknown binary has its own unknown density/atmosphere;
        # no target density or hot-star limb priors are transferred to it.
        return res
    trials = []
    for b,a,light in ((.5,10.,.001),(1.4,7.,.003),(1.8,4.,.02)):
        seed = np.array([b/2,np.log(a),.5**2,.4,0.,0.,.6,light])
        trials.append(least_squares(residual,seed,bounds=bounds,x_scale='jac',max_nfev=200))
    result = min(trials,key=lambda r:np.sum(r.fun**2))
    v = result.x; row = model.description(np.r_[1.,v[:-2]])
    # Keep inherited target-star radius clearly separated from source fit.
    for key in ['conditional_companion_radius_rsun','conditional_companion_radius_rjup','conditional_vector_covariance','conditional_vector_errors']:
        row.pop(key,None)
    alpha = 1.5
    row.update(scenario='Unknown-source equal-radius, unequal-brightness detached binary diluted by target aperture light',
               data_chi2=float(np.sum(result.fun**2)),prior_chi2=0.,penalized_chi2=float(np.sum(result.fun**2)),
               surface_brightness_ratio=float(v[-2]), binary_fraction_of_aperture_light=float(v[-1]),
               dilution_fraction=1-float(v[-1]),converged=bool(result.success),nfev=int(result.nfev),
               vector_with_brightness_light=v.tolist(),
               approximate_diluted_ellipsoidal_amplitude_ppm=float(alpha/row['a_over_rstar']**3*np.sin(np.deg2rad(row['inclination_degrees']))**2*v[-1]*1e6),
               ellipsoidal_assumptions='Illustrative alpha1.5,q1, spherical/weak-distortion scaling, not a fitted mass/light curve.',
               source_caveat='Actual binary radius/mass/Teff/source unconstrained; target properties deliberately not used as binary density prior. No independent radius inference.')
    return row,prediction(v)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    data,star,timing,provenance=inputs()
    rows, predictions=[], {}
    empirical, predictions['alternating_trapezoid']=fit_trapezoid(data,timing)
    save(OUT/'empirical_comparator.json',empirical)
    # First fit all observed recurrences as one identical planet per P.
    row,pred=fitted_primary(data,star,timing,'identical_primary_at_photometric_P')
    rows.append(row);predictions[row['scenario']]=pred
    even=(data['cycle']%2)==0
    primary_data={k:v[even] for k,v in data.items()}
    orbit={**timing,'period':2*timing['period']}
    # Alternation changes the physical orbit hypothesis. Fit only every-other
    # primary here, so the other eclipses do not bias its physical morphology.
    for e,w in [(0.,90.),(.15,270.),(.30,270.),(.45,270.),(.30,90.)]:
        label=f'2P_even_primary_e{e:g}_omega{w:g}'
        row,pred=fitted_primary(primary_data,star,orbit,label,eccentricity=e,omega=w,kmax=.12)
        rows.append(row);predictions[label]=pred
        save(OUT/'primary_fits.json',rows)
    row,pred=fitted_primary(primary_data,star,orbit,'2P_even_primary_free_density',density=False,kmax=.12)
    rows.append(row);predictions[row['scenario']]=pred
    binary,predictions['diluted_binary']=diluted_binary(data,star,timing)
    summary={'tic':TIC,'checked_utc':datetime.now(timezone.utc).isoformat(),'provenance':provenance,
             'stellar_scenario':star,'photometric_ephemeris':timing,'empirical_comparator':empirical,
             'primary_fits':rows,'diluted_binary_alternative':binary,
             'source_sha256':sha(Path(__file__)), 'frozen_import_sha256':EXPECTED,
             'assessment':'Recurring alternating shallow signal; source/physical interpretation unresolved. No accepted new exoplanet.',
             'limitations':['Conditional diagonal errors and protected local quadratic baselines; correlated noise uncalibrated.',
                            'Hot-star rotation/pulsation/oblateness and actual stellar uncertainties may affect models.',
                            'A low catalog contamination ratio is not proof that a sub-percent-light deep eclipsing source is absent.']}
    save(OUT/'fit_results.json',summary)
    np.savez_compressed(OUT/'fit_arrays.npz',**data,**{k+'_model':v for k,v in predictions.items()})
    phase=(data['t']-timing['epoch_btjd']+timing['period']/2)%timing['period']-timing['period']/2
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for p,color in [(0,'tab:blue'),(1,'tab:orange')]:
        use=data['cycle']%2==p
        for lo,hi in zip(np.linspace(-3.84,3.84,49)[:-1],np.linspace(-3.84,3.84,49)[1:]):
            ii=use&(phase*24>=lo)&(phase*24<hi)
            if ii.sum():
                ww=1/data['e'][ii]**2
                ax[0,p].errorbar(np.average(phase[ii]*24,weights=ww),(np.average(data['f'][ii],weights=ww)-1)*1e6,
                                yerr=1/np.sqrt(ww.sum()),fmt='o',ms=3,color=color)
        order=np.argsort(phase[use]);xx=phase[use][order]*24
        for key,c in [('alternating_trapezoid','k'),('diluted_binary','tab:purple')]:
            ax[0,p].plot(xx,(predictions[key][use][order]-1)*1e6,'.',ms=1,color=c,label=key)
        ax[0,p].set(xlabel='Hours from recurrence',ylabel='Relative flux ppm',title='Even primary' if p==0 else 'Odd secondary',ylim=(-1000,400))
        ax[0,p].legend(fontsize=7)
    xx=phase[even]*24;order=np.argsort(xx)
    for row in rows[1:]:
        pred=predictions[row['scenario']]
        ax[1,0].plot(xx[order],(pred[order]-1)*1e6,'.',ms=1,label=row['scenario'])
    ax[1,0].set(xlabel='Hours from even primary',ylabel='Physical model ppm',ylim=(-1000,100))
    ax[1,0].legend(fontsize=6)
    ax[1,1].plot(phase*24,(data['f']-predictions['diluted_binary'])*1e6,'.',ms=1,alpha=.2)
    ax[1,1].set(xlabel='Hours from recurrence',ylabel='Diluted binary residual ppm')
    fig.suptitle('TIC260820871: conditional shape tests; source and mass unmeasured')
    fig.savefig(OUT/'physical_fits.png',dpi=170);plt.close(fig)
    lines=['# TIC260820871: physical morphology and source ambiguity','',
           f"Observed photometric recurrence P={timing['period']:.9f} d has alternating depths. A free integrated trapezoid on the same native photons gives even/odd depths {empirical['even_depth']*1e6:.1f}/{empirical['odd_depth']*1e6:.1f} ppm, T14={empirical['duration_T14_hours']:.3f} h, ingress fraction={empirical['ingress_fraction']:.3f}. The physical orbital interpretation of alternating eclipses is approximately twice this period. This diagnostic does not measure companion mass.",'',
           'The nominal current TIC scenario is R=2.40763 Rsun, M=2.85 Msun, Teff=10942 K. R/M20% and temperature1000K widths are explicit conservative assumptions, not catalog errors. The circular a/R at the photometric period is about1.9. Exact native30min/10min exposure durations, raw byte hashes and TDB time references are verified. The small depth and unusually prolonged ingress require testing physical geometry rather than identifying radius with sqrt(depth).','',
           '| Primary hypothesis | k | b | T14 h | Conditional R/Rjup | Data χ² | Density/physical prior χ² | Photosphere clear |',
           '|---|---:|---:|---:|---:|---:|---:|---|']
    for r in rows:
        lines.append(f"| {r['scenario']} | {r['radius_ratio']:.5f} | {r['impact_parameter']:.4f} | {r['duration_T14_hours']:.3f} | {r['conditional_companion_radius_rjup']:.2f} | {r['data_chi2']:.1f} | {r['prior_and_physical_penalties_chi2']:.2f} | {r['photosphere_condition_satisfied']} |")
    lines+=['','Every-other primary fits use the physical2P orbit. Eccentricity scenarios test apastron/periastron velocities and enforce a(1-e)>Rstar+Rcomp; an arbitrary high-e long-duration fit that intersects the photosphere is inadmissible. The free-density scenario deliberately relaxes target stellar constraints. These are scenario comparisons with diagonal errors, not a mass validation, formal model selection, or rejection of all planetary possibilities. Even-primary χ² values use fewer photons than the all-event empirical comparator and must not be directly subtracted from it.','',
            f"A separate unknown-source detached binary uses two equal radii, unequal surface brightness, and physical2P orbit, with a/R≥3 for detached Roche geometry. It fits dataχ²={binary['data_chi2']:.1f}, compared with empirical all-event χ²={empirical['data_chi2']:.1f}; brightness ratio={binary['surface_brightness_ratio']:.3f}, binary fraction of aperture light={binary['binary_fraction_of_aperture_light']:.5f}, b={binary['impact_parameter']:.3f}, a/R={binary['a_over_rstar']:.2f}. Source density/temperature are free; target properties are not incorrectly imposed on a faint blended binary. Its illustrative diluted ellipsoidal amplitude is {binary['approximate_diluted_ellipsoidal_amplitude_ppm']:.1f} ppm (alpha1.5,q1 approximate scaling). A faint binary can therefore have undetectable broad orbital humps. The catalog contamination ratio0.00246 is a prediction, not a source-localization measurement or a guaranteed exclusion of faint deep eclipses.",'',
            'A thermal/reflection budget cannot substitute for a feasible, localized planetary primary geometry. Optical eclipse brightness of a hot or intrinsically luminous object never establishes mass. Pulsation/rotation, hot-star geometry, blended stellar eclipses, uncertain catalog radius and source attribution remain alternatives. **No credible new exoplanet is accepted from these tests.** Independent pixel localization, actual stellar constraints and source-specific spectroscopy/photometry are needed.','',
            'Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_260820871.py`. Frozen earlier two-lead source hash, original raw identifiers/hashes, fitted parameters, model arrays and diagnostic plots are preserved. No observations were deleted or paid models used.','']
    (OUT/'README.md').write_text('\n'.join(lines))
    print(json.dumps({'empirical':empirical,'primaries':[{k:r[k] for k in ['scenario','radius_ratio','impact_parameter','duration_T14_hours','data_chi2','prior_and_physical_penalties_chi2','photosphere_condition_satisfied']} for r in rows], 'binary':binary},indent=2),flush=True)


if __name__=='__main__':main()
