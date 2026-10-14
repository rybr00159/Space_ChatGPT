"""Finite-exposure physical shape tests of two recurring broad-survey leads.

These conditional fits are deliberately not planet validation. In particular,
grazing stars, brown dwarfs, unresolved dilution, uncertain stellar radii, and
eccentric velocity changes must remain explicit alternatives. No paid calls.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
import batman
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from broad_vet import local_quadratic
from fit_370247818 import TransitFit, best_fit, sha, save, RJUP_RSUN

ROOT = Path(__file__).resolve().parents[1]
SUPERSAMPLE = 31


class LeadFit(TransitFit):
    bounds = (np.array([.04, .001, np.log(2.5), .001, .001, -.035, -35.]),
              np.array([.8, .9999, np.log(35.), .999, .999, .035, 35.]))


class EqualBinaryFit(LeadFit):
    """Equal brightness/discs: two identical eclipses per physical orbit.

    A dark k=1 disc occulting one star loses half the system's light. Folded
    photons repeat every photometric P; BATMAN motion and Kepler density use
    physical orbital 2P, two equal masses and individual R=combinedR/sqrt(2).
    This is a self-consistent illustrative shape alternative, not a binary
    population inference, Roche-distorted model, or stellar evolution fit.
    """
    def __init__(self, data, prior, timing):
        individual = {**prior, 'radius_rsun': prior['radius_rsun']/np.sqrt(2),
                      'radius_error_rsun': prior['radius_error_rsun']/np.sqrt(2)}
        individual['catalog_rho_solar'] = individual['mass_msun']/individual['radius_rsun']**3
        individual['catalog_rho_error_solar'] = .35*individual['catalog_rho_solar']
        super().__init__(data, individual, {**timing, 'period': 2*timing['period']},
                         dilution=.5, companion_mass_ratio=1., supersample=SUPERSAMPLE)
        cy = np.round((data['t']-timing['epoch_btjd'])/timing['period'])
        self.fold = data['t']-timing['epoch_btjd']-cy*timing['period']
        # A fitted orbital-period increment displaces each photometric event
        # by half that increment, hence fractional cycles in the parent model.
        self.cycle = cy/2

    def binary_fit(self):
        fits_ = []
        for b in (1.3, 1.6, 1.85, 1.95):
            seed = np.array([1., b/2, self.loga_mean, .6**2, .4, 0., 0.])
            fits_.append(self.fit(seed, fixed_k=1., max_nfev=180))
        return min(fits_, key=lambda pair: np.sum(pair[1].fun**2))


def combine(chunks):
    return {k: np.concatenate([x[k] for x in chunks]) for k in chunks[0]}


def verify_source(source, pixel=False):
    path = Path(source.get('path', source.get('raw_path')))
    if sha(path) != source['sha256']:
        raise ValueError(f'Raw SHA mismatch: {path}')
    with fits.open(path, memmap=False) as hd:
        h = hd[0].header if pixel else hd[1].header
        if h['TIMESYS'] != 'TDB' or abs(h['BJDREFI'] + h.get('BJDREFF', 0) - 2457000) > 1e-8:
            raise ValueError('Unexpected astronomical time system/reference')
        exposure = float(hd[1].header['TIMEDEL'])
    return {**source, 'verified_sha256': sha(path), 'exposure_days_exact_TIMEDEL': exposure}


def star269():
    path = ROOT / 'data/broad_catalogs/tic269994642_exofop.json'
    j = json.loads(path.read_text())
    row = next(s for s in j['stellar_parameters'] if s.get('srad'))
    star = {'radius_rsun': float(row['srad']), 'radius_error_rsun': float(row['srad_e']),
            'mass_msun': float(row['mass']), 'mass_error_msun': float(row['mass_e']),
            'teff_K': float(row['teff']), 'teff_error_K': float(row['teff_e']),
            'catalog_rho_solar': float(row['mass']) / float(row['srad'])**3,
            'catalog_rho_error_solar': float(row['dens_e']) / 1.408,
            'limb_darkening_u1_mean': .50, 'limb_darkening_u1_sigma': .18,
            'limb_darkening_u2_mean': .20, 'limb_darkening_u2_sigma': .20,
            'stellar_source': str(path.relative_to(ROOT)), 'stellar_source_sha256': sha(path),
            'limb_darkening_source': 'Broad approximate TESS priors for a ~4035 K star; no atmosphere interpolation.',
            'catalog_contaminant_to_target_flux_ratio': float(j['basic_info']['tic_contamination_ratio'])}
    return star


def star238(radius='tic'):
    path = ROOT / 'data/broad_catalogs/tic238202207_tic_cone_current.json'
    row = next(s for s in json.loads(path.read_text())['data'] if s['ID'] == 238202207)
    r = float(row['rad']) if radius == 'tic' else 1.1929
    star = {'radius_rsun': r, 'radius_error_rsun': .10 * r,
            'mass_msun': float(row['mass']), 'mass_error_msun': .15 * float(row['mass']),
            'teff_K': float(row['Teff']) if radius == 'tic' else 5244.236,
            'teff_error_K': float(row['e_Teff']),
            'catalog_rho_solar': float(row['mass']) / r**3,
            'catalog_rho_error_solar': .35 * float(row['mass']) / r**3,
            'limb_darkening_u1_mean': .38, 'limb_darkening_u1_sigma': .18,
            'limb_darkening_u2_mean': .22, 'limb_darkening_u2_sigma': .20,
            'stellar_source': str(path.relative_to(ROOT)), 'stellar_source_sha256': sha(path),
            'radius_hypothesis': radius, 'radius_uncertainty_assumption': '10% scenario width, not a catalog error; TIC radius/mass errors unavailable.',
            'mass_uncertainty_assumption': '15% scenario width, not a catalog error.',
            'limb_darkening_source': 'Broad approximate TESS priors; no atmosphere interpolation.',
            'stellar_warning': 'Gaia RUWE4.128, IPDmulti_peak22%; both single-star radius solutions may be biased by multiplicity.'}
    return star


def inputs238():
    folder = ROOT / 'candidates/broad_s1/tic238202207_multisector'
    record = json.loads((folder / 'multisector.json').read_text())
    timing = record['event_timing_fit']['fitted_ephemeris']
    archive = np.load(folder / 'multisector_arrays.npz')
    chunks, sources = [], []
    for source in record['sources']:
        verified = verify_source(source)
        s = source['sector']
        t, f, e = (archive[f's{s}_local_{k}'] for k in ('t', 'f', 'e'))
        cycle = archive[f's{s}_local_cycle']
        phase = t - timing['epoch_btjd'] - np.round((t-timing['epoch_btjd']) / timing['period']) * timing['period']
        good = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0) & (abs(phase) < .16)
        chunks.append(dict(t=t[good], f=f[good], e=e[good], cycle=cycle[good],
                           sector=np.full(good.sum(), s, int), exposure=np.full(good.sum(), verified['exposure_days_exact_TIMEDEL'])))
        sources.append(verified)
    data = combine(chunks)
    provenance = {'input_array_sha256': sha(folder / 'multisector_arrays.npz'),
                  'input_record_sha256': sha(folder / 'multisector.json'), 'raw_sources': sources,
                  'method': 'Unbinned local out-of-transit quadratic normalized arrays, truncated to +/-0.16d; source broad_vet.py verified by hash.',
                  'exposures_days': sorted(np.unique(data['exposure']).tolist())}
    return data, star238(), timing, provenance


def inputs269(aperture=2.):
    folder = ROOT / 'candidates/broad_s1/tic269994642_pixels'
    chunks, sources, local_stats = [], [], []
    # S1/S92 independently measured centers imply cycle1419. Adjacent cycle
    # hypotheses change P by 0.00123d, inconsistent with individual-sector BLS.
    p0 = (3804.4520957293685 - 1326.3488754694151) / 1419
    epoch = 1326.3488754694151
    for sector in (1, 92):
        record = json.loads((folder / f'sector{sector}/pixel_vetting.json').read_text())
        raw = record['raw']
        source = verify_source(raw, pixel=True)
        sources.append(source)
        with fits.open(raw['raw_path'], memmap=False) as hd:
            dd = hd[1].data
            t, q = np.asarray(dd['TIME'], float), np.asarray(dd['QUALITY'])
            ff, ee = np.asarray(dd['FLUX'], float), np.asarray(dd['FLUX_ERR'], float)
            xx, yy = np.meshgrid(np.arange(ff.shape[2]), np.arange(ff.shape[1]))
            xy = WCS(hd[2].header).all_world2pix([[311.461584, -24.649866]], 0)[0]
            keep = (q == 0) & np.isfinite(t) & np.all(np.isfinite(ff) & np.isfinite(ee), axis=(1, 2))
            t, ff, ee = t[keep], ff[keep], ee[keep]
            ref = np.median(ff, axis=0)
            blank = ref < np.percentile(ref, 30)
            background = np.median(ff[:, blank], axis=1)
            ap = (xx-xy[0])**2 + (yy-xy[1])**2 < aperture**2
            flux = (ff - background[:, None, None])[:, ap].sum(axis=1)
            # Median-sky uncertainty assumes independent approximately
            # Gaussian blank pixels; real sky correlations remain uncovered.
            sky_var = np.pi/2*np.median(ee[:, blank]**2, axis=1)/blank.sum()
            error = np.sqrt((ee[:, ap]**2).sum(axis=1) + ap.sum()**2*sky_var)
        stat, local = local_quadratic(t, flux, error, p0, epoch, .08)
        if local is None:
            raise ValueError('No adequate bilateral event windows')
        phase = local['t'] - epoch - np.round((local['t']-epoch) / p0) * p0
        use = abs(phase) < .16
        chunks.append({**{k: v[use] for k, v in local.items()},
                       'sector': np.full(use.sum(), sector, int),
                       'exposure': np.full(use.sum(), source['exposure_days_exact_TIMEDEL'])})
        local_stats.append({'sector': sector, **stat})
    provenance = {'raw_sources': sources, 'local_stats': local_stats,
                  'input_method': 'Independent raw pixel aperture: target-centered radius2pixel, median of bottom30%reference pixels sky, QUALITY0, bilateral OOT quadratic baseline excluding +/-1.15*0.08d; no fixed transit-template subtraction or pointing decorrelation.',
                  'baseline_source_sha256': sha(ROOT / 'src/broad_vet.py'), 'aperture_radius_pixels': aperture,
                  'global_cycle_hypothesis': 1419, 'alternative_cycle_spacing_days': p0 / 1419,
                  'timing_caveat': 'Two measured-era epoch initialization; no repeated photon-derived timing prior used in likelihood.'}
    return combine(chunks), star269(), {'period': p0, 'epoch_btjd': epoch}, provenance


def multi_start(model, seed=None):
    seeds = []
    if seed is not None:
        seeds.append(np.asarray(seed))
    for k in (.15, .22, .35, .6):
        seeds.append(np.array([k, min(.97, .97/(1+k)), model.loga_mean,
                               .6**2, .4, 0., 0.]))
    fits_ = [model.fit(s, max_nfev=130) for s in seeds]
    return min(fits_, key=lambda pair: np.sum(pair[1].fun**2))


def radius_profile(model, best, path, points=38):
    grid = np.unique(np.r_[np.geomspace(.10, .8, points), .13, .15, .18, .2, .25, best[0]])
    rows = json.loads(path.read_text()) if path.exists() else []
    done = {round(r['radius_ratio'], 12): r for r in rows}
    for k in sorted(grid, key=lambda value: abs(value-best[0])):
        if round(float(k), 12) in done:
            continue
        vectors = [best]
        if done:
            vectors.append(np.asarray(min(done.values(), key=lambda r: abs(r['radius_ratio']-k))['vector']))
        fits_ = [model.fit(v, fixed_k=float(k), max_nfev=100) for v in vectors]
        vec, result = min(fits_, key=lambda pair: np.sum(pair[1].fun**2))
        r = model.description(vec)
        r.update(converged=bool(result.success), nfev=int(result.nfev))
        done[round(float(k), 12)] = r
        save(path, sorted(done.values(), key=lambda r: r['radius_ratio']))
        print('PROFILE', path.parent.name, len(done), len(grid), flush=True)
    return sorted(done.values(), key=lambda r: r['radius_ratio'])


def profile_summary(rows):
    minimum = min(r['penalized_chi2'] for r in rows)
    summary = {'minimum_penalized_chi2': minimum, 'grid_min': min(r['radius_ratio'] for r in rows),
               'grid_max': max(r['radius_ratio'] for r in rows)}
    for level in (1., 3.84):
        admitted = [r for r in rows if r['penalized_chi2']-minimum <= level]
        summary[f'delta_chi2_le_{level:g}'] = {'sampled_k_range': [min(r['radius_ratio'] for r in admitted), max(r['radius_ratio'] for r in admitted)],
            'upper_bound_reached': max(r['radius_ratio'] for r in admitted) == summary['grid_max'],
            'conditional_radius_rjup_range': [min(r['conditional_companion_radius_rjup'] for r in admitted), max(r['conditional_companion_radius_rjup'] for r in admitted)]}
    summary['warning'] = 'Conditional sampled penalized likelihood profile; no calibrated coverage, covariance or mass classification. Fixed stellar-density/limb priors and diagonal errors; grazing upper bounds are artificial.'
    return summary


def fit_alternatives(tic, data, prior, timing, best, out):
    scenarios = [('free_density', prior, {'density_prior': False}),
                 ('dilution20percent', prior, {'dilution': .2}),
                 ('dilution50percent', prior, {'dilution': .5}),
                 ('massratio0.3', prior, {'companion_mass_ratio': .3}),
                 ('eccentric0.3_periastron', prior, {'eccentricity': .3, 'omega': 90.}),
                 ('eccentric0.5_periastron', prior, {'eccentricity': .5, 'omega': 90.})]
    if tic == 238202207:
        scenarios += [('gaia_radius1.1929', star238('gaia'), {}),
                      ('gaia_radius_dilution50percent', star238('gaia'), {'dilution': .5})]
    else:
        ratio = prior['catalog_contaminant_to_target_flux_ratio']
        scenarios += [('catalog_dilution', prior, {'dilution': ratio/(1+ratio)})]
    rows = []
    for label, star, arguments in scenarios:
        model = LeadFit(data, star, timing, supersample=SUPERSAMPLE, **arguments)
        vec, result = multi_start(model, best)
        row = {'scenario': label, 'stellar_prior': star, **model.description(vec, result)}
        rows.append(row)
        save(out / 'alternative_fits.json', rows)
        print('ALTERNATIVE', tic, label, row['radius_ratio'], row['penalized_chi2'], flush=True)
    return rows


def equal_binary(data, prior, timing, out):
    model = EqualBinaryFit(data, prior, timing)
    vector, result = model.binary_fit()
    row = model.description(vector, result)
    row.update(model='Identical equal-light stars, circular physical orbit twice observed photometric period',
               individual_stellar_prior=model.prior,
               physical_orbital_period_days=row['period_days'],
               observed_photometric_period_days=row['period_days']/2,
               extra_assumptions=['Each star has the same Teff and individual radius combined-light inferred radius divided by sqrt2.',
                                  'Each star mass equals nominal catalog mass with broad density prior; this is illustrative, not an isochrone binary solution.',
                                  'Two identical limb-darkened eclipses per orbit; background star contributes exactly half the light.',
                                  'No Roche deformation, spots, mutual irradiation or Doppler beaming; local detrending removes slow orbital light terms.',
                                  'Both eclipses share geometry and are therefore indistinguishable by odd/even or secondary-at-photometric-half tests.'])
    save(out/'equal_light_binary.json', row)
    np.savez_compressed(out/'equal_light_binary_model.npz', time=data['t'], flux=data['f'], error=data['e'], model=model.physical_model(vector))
    return row


def unequal_binary(data, prior, timing, equal, out):
    """Alternating equal-radius stars with a free surface brightness ratio.

    Size/density scenario is illustrative and unchanged from the equal-light
    fit; J is the secondary/primary surface brightness, not a mass ratio.
    """
    model = EqualBinaryFit(data, prior, timing)
    parity = np.round((data['t']-timing['epoch_btjd'])/timing['period']).astype(int) % 2

    def prediction(reduced):
        physical = np.r_[1., reduced[:-1]]
        equal_curve = model.physical_model(physical)
        j = reduced[-1]
        factor = np.where(parity == 0, 2/(1+j), 2*j/(1+j))
        return 1-(1-equal_curve)*factor

    def residual(reduced):
        flux = prediction(reduced)
        offsets = np.bincount(model.sector_idx, weights=model.weight*(data['f']-flux))/model.sector_weight
        datares = (data['f']-flux-offsets[model.sector_idx])/data['e']
        priors = model.residual(np.r_[1., reduced[:-1]])[len(datares):]
        return np.r_[datares, priors]

    initial = np.r_[np.asarray(equal['vector'])[1:], 1.]
    bounds = (np.r_[model.bounds[0][1:], .3], np.r_[model.bounds[1][1:], 3.])
    result = least_squares(residual, initial, bounds=bounds, x_scale='jac', max_nfev=180,
                           ftol=1e-8, xtol=1e-8, gtol=1e-7)
    r = residual(result.x)
    row = model.description(np.r_[1., result.x[:-1]])
    covariance = np.linalg.pinv(result.jac.T@result.jac)
    row.update(model='Equal radii, unequal surface brightness, alternating eclipses at physical orbit2P',
               surface_brightness_ratio_secondary_to_primary=float(result.x[-1]),
               conditional_surface_brightness_error=float(np.sqrt(covariance[-1, -1])),
               reduced_vector=result.x.tolist(), conditional_reduced_covariance=covariance.tolist(),
               data_chi2=float(np.sum(r[:len(data['t'])]**2)),
               prior_chi2=float(np.sum(r[len(data['t']):]**2)), penalized_chi2=float(np.sum(r**2)),
               delta_data_chi2_versus_equal_light=float(equal['data_chi2']-np.sum(r[:len(data['t'])]**2)),
               physical_orbital_period_days=row['period_days'],
               observed_photometric_period_days=row['period_days']/2,
               converged=bool(result.success), nfev=int(result.nfev),
               individual_stellar_prior=model.prior,
               caveats=['No mass measurement; same radii and nominal combined-light radius/√2 approximation despite unequal brightness.',
                        'J is a TESS surface-brightness/flux ratio, not a companion mass ratio or calibrated temperature.',
                        'Conditional diagonal likelihood with same unmodeled baseline/red-noise covariance as primary fits.',
                        'Improvement is a targeted diagnostic, not survey significance or validation.'])
    save(out/'unequal_brightness_binary.json', row)
    np.savez_compressed(out/'unequal_brightness_binary_model.npz', time=data['t'], flux=data['f'], error=data['e'], parity=parity, model=prediction(result.x))
    return row


def cycle_alias_tests(data, prior, timing, best, out):
    rows = []
    difference = 3804.4520957293685-1326.3488754694151
    for cycle in (1418, 1419, 1420):
        epoch = {'period': difference/cycle, 'epoch_btjd': timing['epoch_btjd']}
        model = LeadFit(data, prior, epoch, supersample=SUPERSAMPLE)
        vec, result = multi_start(model, best['vector'])
        rows.append({'integer_cycle': cycle, 'initial_period_days': epoch['period'], **model.description(vec, result)})
        save(out/'cycle_aliases.json', rows)
    return rows


def parity_robustness(tic, data, prior, timing, equal, out):
    """Independent native-photon event-amplitude and parity diagnostics.

    Freeze the physical equal-binary shape, then fit an amplitude+constant
    separately in each event window. Bootstrap event clusters within sectors;
    the two-sector269 result cannot estimate long-term stellar systematics.
    """
    model = EqualBinaryFit(data, prior, timing)
    vector = np.asarray(equal['vector'])
    template = 1-model.physical_model(vector)
    cycles = np.round((data['t']-timing['epoch_btjd'])/timing['period']).astype(int)
    measurements = []
    for sector, cycle in sorted(set(zip(data['sector'], cycles))):
        use = (data['sector'] == sector) & (cycles == cycle)
        a = np.column_stack([np.ones(use.sum()), -template[use]])
        weight = 1/data['e'][use]**2
        covariance = np.linalg.pinv(a.T@(weight[:, None]*a))
        coef = covariance@(a.T@(weight*data['f'][use]))
        residual = data['f'][use]-a@coef
        inflation = max(1., np.sum(residual**2*weight)/max(1, use.sum()-2))
        measurements.append({'sector': int(sector), 'cycle': int(cycle), 'parity': int(cycle % 2),
                             'n_samples': int(use.sum()), 'amplitude_relative_to_equal_binary': float(coef[1]),
                             'amplitude_error_residual_inflated': float(np.sqrt(covariance[1, 1]*inflation)),
                             'residual_chi2_inflation': float(inflation)})
    with (out/'event_amplitude_measurements.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(measurements[0]))
        writer.writeheader(); writer.writerows(measurements)
    groups = []
    for sector in np.unique(data['sector']):
        rows = [r for r in measurements if r['sector'] == sector]
        d = np.array([r['amplitude_relative_to_equal_binary'] for r in rows])
        e = np.array([r['amplitude_error_residual_inflated'] for r in rows])
        parity = np.array([r['parity'] for r in rows])
        if min((parity == 0).sum(), (parity == 1).sum()) < 2:
            continue
        means, errors = [], []
        for p in (0, 1):
            ii = parity == p; w = 1/e[ii]**2
            means.append(float(np.average(d[ii], weights=w)))
            errors.append(float(max(1/np.sqrt(w.sum()), np.std(d[ii], ddof=1)/np.sqrt(ii.sum()))))
        groups.append({'sector': int(sector), 'even_amplitude': means[0], 'odd_amplitude': means[1],
                       'even_error_max_diagonal_scatter': errors[0], 'odd_error_max_diagonal_scatter': errors[1],
                       'even_minus_odd': means[0]-means[1], 'difference_error': float(np.hypot(*errors)),
                       'odd_even_snr_conditional': float((means[0]-means[1])/np.hypot(*errors)),
                       'odd_to_even_amplitude_ratio': means[1]/means[0]})
    sectors = np.array([r['sector'] for r in measurements])
    parity = np.array([r['parity'] for r in measurements])
    d = np.array([r['amplitude_relative_to_equal_binary'] for r in measurements])
    e = np.array([r['amplitude_error_residual_inflated'] for r in measurements])
    unique = np.unique(sectors)
    x = np.column_stack([(sectors == s).astype(float) for s in unique]+[(.5-parity)])
    w = 1/e**2
    fit = lambda ii: np.linalg.lstsq(x[ii]*np.sqrt(w[ii, None]), d[ii]*np.sqrt(w[ii]), rcond=None)[0]
    coef = fit(np.arange(len(d)))
    covariance = np.linalg.pinv(x.T@(w[:, None]*x))
    rng = np.random.default_rng(tic)
    draws = []
    for draw in range(3000):
        indices = np.concatenate([rng.choice(np.flatnonzero(sectors == s), (sectors == s).sum(), replace=True) for s in unique])
        # Degenerate parity/bootstrap designs carry no information and are
        # explicitly rejected, not silently treated as precise differences.
        if np.linalg.matrix_rank(x[indices]) == x.shape[1]:
            draws.append(float(fit(indices)[-1]))
    pooled_sector = [g['even_minus_odd'] for g in groups]
    pooled_error = [g['difference_error'] for g in groups]
    ww = 1/np.array(pooled_error)**2
    conservative = np.average(pooled_sector, weights=ww)
    er = max(1/np.sqrt(ww.sum()), np.std(pooled_sector, ddof=1)/np.sqrt(len(groups))) if len(groups)>1 else 1/np.sqrt(ww.sum())
    summary = {'tic': tic, 'fixed_shape': 'Equal-radius binary physical shape, exact exposures; independent amplitude+constant fit per event.',
               'n_events': len(measurements), 'sector_results': groups,
               'sector_intercepts_even_minus_odd_amplitude': float(coef[-1]),
               'conditional_diagonal_error': float(np.sqrt(covariance[-1, -1])),
               'event_bootstrap_draws': len(draws), 'within_sector_event_bootstrap_error': float(np.std(draws, ddof=1)),
               'within_sector_bootstrap_16_50_84': np.percentile(draws, [16, 50, 84]).tolist(),
               'sector_depth_difference_mean': float(conservative),
               'sector_depth_difference_error_max_diagonal_sector_scatter': float(er),
               'conditional_snr_conservative_sector_comparison': float(conservative/er),
               'events_with_positive_amplitude': int((d > 0).sum()),
               'caveats': ['A targeted diagnostic, not a global significance or independent use of the same photons.',
                          'Stellar activity, fixed-shape mismatch and baseline systematics can affect alternation.',
                          'Event bootstrap keeps these observed sectors fixed; two sectors cannot characterize long-term systematics.',
                          'Same-source physical and prior QLP diagnostics must not be multiplied as independent probabilities.']}
    save(out/'parity_robustness.json', summary)
    np.savez_compressed(out/'parity_event_bootstrap.npz', even_minus_odd_amplitude=np.asarray(draws))
    return summary


def verify_binary_model(tic, data, prior, timing, equal, out):
    """Compare folded fit with two independent full-time orbital eclipses."""
    model = EqualBinaryFit(data, prior, timing)
    vector = np.asarray(equal['vector'])
    folded = model.physical_model(vector)
    params = model.params
    independent = np.empty(len(data['t']))
    for exposure in np.unique(data['exposure']):
        use = data['exposure'] == exposure
        params.t0 = timing['epoch_btjd']+vector[-2]
        primary = batman.TransitModel(params, data['t'][use], supersample_factor=SUPERSAMPLE,
                                     exp_time=float(exposure), nthreads=1).light_curve(params)
        params.t0 += params.per/2
        secondary = batman.TransitModel(params, data['t'][use], supersample_factor=SUPERSAMPLE,
                                       exp_time=float(exposure), nthreads=1).light_curve(params)
        independent[use] = 1-.5*((1-primary)+(1-secondary))
    difference = float(np.max(abs(folded-independent)))
    if difference >= 1e-7:
        raise ValueError(f'Doubled-period binary mapping failed full-orbit check: {difference}')
    row = {'tic': tic, 'max_absolute_flux_difference': difference, 'status': 'passed',
           'method': 'Independent full native-time BATMAN primary+secondary models at physical2P, equal discs/brightness, summed lost light; no photometric phase folding.',
           'native_samples': len(data['t']), 'supersample': SUPERSAMPLE,
           'physical_orbital_period_days': equal['period_days']}
    save(out/'binary_model_verification.json', row)
    return row


def save_plot(tic, data, model, best, rows, alternatives, out):
    phase = (data['t']-best['epoch_btjd']+best['period_days']/2) % best['period_days'] - best['period_days']/2
    pred = model.physical_model(np.array(best['vector']))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for sector in np.unique(data['sector']):
        ii = data['sector'] == sector
        axes[0, 0].plot(phase[ii]*24, data['f'][ii], '.', ms=1.5, alpha=.25, label=f'S{sector}')
    bins = np.linspace(-3.84, 3.84, 65)
    for lo, hi in zip(bins[:-1], bins[1:]):
        ii = (phase*24 >= lo) & (phase*24 < hi)
        if ii.sum():
            w = 1/data['e'][ii]**2
            axes[0, 0].errorbar(np.average(phase[ii]*24, weights=w), np.average(data['f'][ii], weights=w), yerr=1/np.sqrt(w.sum()), fmt='ko', ms=3)
    # Model shown for smallest exposure; actual likelihood uses every exact exposure.
    dense = np.linspace(-3.84, 3.84, 900)
    dd = dict(t=model.epoch0+dense/24, f=np.ones(900), e=np.ones(900), cycle=np.zeros(900, int), sector=np.ones(900, int), exposure=np.full(900, np.min(data['exposure'])))
    md = LeadFit(dd, model.prior, model.timing, supersample=31)
    axes[0, 0].plot(dense, md.physical_model(np.array(best['vector'])), 'k-', label='Best shortest exposure')
    binary = json.loads((out/'equal_light_binary.json').read_text())
    bd = EqualBinaryFit(dd, model.prior, model.timing)
    # Dense dummy times correspond to a single eclipse; orbital model still
    # uses2P and individual stellar density exactly as fitted on real photons.
    axes[0, 0].plot(dense, bd.physical_model(np.array(binary['vector'])), '--', color='tab:purple', lw=2, label='Equal-light 2P binary')
    for desired in (.15, .3, .6):
        row = min(rows, key=lambda r: abs(r['radius_ratio']-desired))
        axes[0, 0].plot(dense, md.physical_model(np.array(row['vector'])), label=f"Profile k={row['radius_ratio']:.2f}")
    axes[0, 0].set(xlabel='Hours from event', ylabel='Local relative flux', xlim=(-3.84, 3.84))
    axes[0, 0].set_ylim(.95 if tic == 269994642 else .97, 1.025)
    axes[0, 0].legend(fontsize=7, ncol=3)
    axes[0, 1].plot(phase*24, (data['f']-pred)*1e3, '.', ms=2, alpha=.3)
    axes[0, 1].axhline(0, color='k', lw=.5)
    axes[0, 1].set(xlabel='Hours from event', ylabel='Residual ppt')
    low = min(r['penalized_chi2'] for r in rows)
    axes[1, 0].plot([r['radius_ratio'] for r in rows], [r['penalized_chi2']-low for r in rows], 'k.-')
    axes[1, 0].axhline(3.84, color='grey', ls='--')
    axes[1, 0].set(xlabel='Radius ratio k', ylabel='Conditional Δ penalized χ²', ylim=(0, 15))
    axes[1, 1].plot([r['radius_ratio'] for r in rows], [r['impact_parameter'] for r in rows], 'k.-', label='Profile')
    axes[1, 1].plot([.1, .8], [.9, .2], '--', color='grey', label='b=1-k grazing boundary')
    axes[1, 1].set(xlabel='Radius ratio k', ylabel='Impact parameter b')
    axes[1, 1].legend()
    fig.suptitle(f'TIC {tic}: conditional transit/eclipse shape, no mass inference')
    fig.savefig(out / 'physical_fits.png', dpi=170)
    plt.close(fig)
    np.savez_compressed(out / 'fit_arrays.npz', **data, model=pred, phase_days=phase)


def report(tic, summary, out):
    best = summary['best_circular_dark_companion']
    p = summary['radius_profile_summary']
    alts = summary['alternatives']
    assumption = ('Current TIC R=0.835257±0.0872364 Rsun, M=0.631±0.0830436 Msun, Teff=4035 K. Catalog contamination ratio 0.003582 is only a prior estimate; unresolved sources are not excluded.' if tic == 269994642 else
                  'TIC/FLAME R≈1.63 Rsun versus Gaia GSP-Phot R=1.1929 Rsun. TIC radius/mass uncertainties are absent, so 10% radius/15% mass widths are explicit scenario assumptions. Gaia RUWE=4.128, IPDmulti_peak=22% undermine a clean single-star interpretation.')
    contour = p['delta_chi2_le_3.84']
    lines = [f'# TIC {tic}: physical transit/eclipse shape assessment', '',
             'This is a conditional shape analysis, **not planet validation or a mass measurement**. Finite exposure, grazing geometry, stellar assumptions and dilution materially affect the inferred companion radius.', '',
             assumption, '',
             f"{best['n_points']} unbinned local samples from sectors {summary['provenance']['sectors']} enter the diagonal likelihood. Exact retained FITS TIMEDEL exposure values are integrated with {SUPERSAMPLE} BATMAN subsamples; 61 subsamples provide a numerical convergence check. Baselines were fitted only outside protected events. Approximate error floors and baseline variance are included, but residual correlations and baseline covariance between samples are not fully modeled. Local normalization can absorb slow orbital light variations; their absence here cannot exclude a binary.", '',
             f"Circular, no-dilution, negligible-companion-mass best fit: P={best['period_days']:.9f} d, epoch BTJD {best['epoch_btjd']:.6f}, k={best['radius_ratio']:.4f}, b={best['impact_parameter']:.4f}, a/R={best['a_over_rstar']:.3f}, T14={best['duration_T14_hours']:.3f} h, T23={best['duration_T23_hours']:.3f} h; grazing={best['grazing']}. Conditional radius={best['conditional_companion_radius_rsun']:.3f} Rsun ({best['conditional_companion_radius_rjup']:.2f} Rjup). Data χ²={best['data_chi2']:.1f}/{best['nominal_dof']} nominal dof, prior penalty={best['prior_chi2']:.2f}.", '',
             f"The sampled radius profile at Δ penalized χ²≤3.84 admits k={contour['sampled_k_range'][0]:.3f}–{contour['sampled_k_range'][1]:.3f}, or {contour['conditional_radius_rjup_range'][0]:.2f}–{contour['conditional_radius_rjup_range'][1]:.2f} Rjup for the adopted radius. Upper grid boundary reached: `{contour['upper_bound_reached']}`. This upper limit is artificial when the boundary is reached. All other physical/timing parameters are reoptimized at each k. The plotted contours include broad approximate limb-darkening and log(a/R) priors. They are not calibrated confidence intervals or Bayesian posterior probabilities; conditional Gaussian covariance is unreliable along grazing degeneracies. Stellar radius uncertainty has not been marginalized into the profile radius range. Exact numerical contours are in `fit_results.json`.", '',
             '| Scenario | k | b | Conditional Rcomp/Rjup | Data chi² | Prior chi² |',
             '|---|---:|---:|---:|---:|---:|']
    for row in alts:
        lines.append(f"| {row['scenario']} | {row['radius_ratio']:.4f} | {row['impact_parameter']:.4f} | {row['conditional_companion_radius_rjup']:.2f} | {row['data_chi2']:.1f} | {row['prior_chi2']:.2f} |")
    eb = summary['equal_light_binary_alternative']
    ub = summary['unequal_brightness_binary_alternative']
    parity = summary['parity_robustness']
    lines += ['', f"An explicit **equal-light binary** has physical orbit twice the photometric period ({eb['physical_orbital_period_days']:.8f} d) and two identical grazing eclipses. Individual radii are the combined-light inferred radius divided by √2, with two equal nominal stellar masses in Kepler's law. It fits b={eb['impact_parameter']:.3f}, a/R={eb['a_over_rstar']:.2f}, T14={eb['duration_T14_hours']:.3f} h, data χ²={eb['data_chi2']:.1f}, prior χ²={eb['prior_chi2']:.2f}. Equal primary/secondary eclipses at alternating photometric epochs evade both odd/even differences and a secondary test at photometric half-phase. Their combined photocenter can be astrometrically balanced. This alternative is conditional on equal Teff/brightness and illustrative mass/radius assumptions; it does not require detectable planetary-period secondary light. Compare its absolute data residuals and priors before any classification. No population odds or Bayes factors are inferred.", '',
              f"Allowing unequal surface brightness at fixed equal radii gives secondary/primary TESS brightness ratio J={ub['surface_brightness_ratio_secondary_to_primary']:.4f}±{ub['conditional_surface_brightness_error']:.4f} (conditional diagonal covariance), data χ²={ub['data_chi2']:.1f}. The data improvement over equal brightness is Δχ²={ub['delta_data_chi2_versus_equal_light']:.2f}. This explicitly models alternating primary/secondary depths; it is not a calibrated population probability. The individual radius combined-light/√2 approximation becomes approximate for unequal brightness, and neither luminosity nor depth proves a stellar mass.", '',
              f"A separate fixed-physical-shape amplitude+constant fit to {parity['n_events']} individual events checks whether diagonal sample errors overstate this alternation. Error inflation includes per-event residual scatter, conservative even/odd event scatter, sector-to-sector scatter, and 3,000 event-cluster bootstrap draws within observed sectors. The conservative sector comparison is {parity['conditional_snr_conservative_sector_comparison']:.2f} conditional standard errors from zero; this is a **targeted diagnostic, not search-wide significance**. Sector and event measurements are in `parity_robustness.json` and `event_amplitude_measurements.csv`. With only two sectors, the TIC 269994642 bootstrap cannot estimate unobserved long-term systematics. No independent probabilities are obtained by multiplying tests that use the same photons. The alternating depths favor a binary interpretation while activity and baseline/shape mismatch remain limitations.", '',
              'Dilution is a fraction of total aperture light; 0.2 and 0.5 are illustrative unresolved-source scenarios, not measured contaminations. A faint stellar companion can reproduce a grazing primary while its optical secondary remains weak. The total-system mass affects Kepler a/R; the mass-ratio 0.3 scenario tests that density effect without imposing a stellar mass-radius relation. Eccentric scenarios change transit velocity with omega=90 degrees and fit the primary alone; they are not orbit/secondary solutions or evidence of eccentricity. Both primary luminosity and unknown extra light make a dark-transit radius conditional. Radii of planets, brown dwarfs and faint stars overlap; mass requires other information.', '',
              'The physical fits leave unresolved transit-like leads for further falsification. Follow source-specific odd/even and all-phase secondary tests, pixel localization, astrometric multiplicity, catalogs/literature and RV/high-resolution imaging before classification. **No new exoplanet is accepted on this shape analysis.**', '',
              f"Numerical integration check: max|{SUPERSAMPLE} vs 61 subsample model|={summary['integration_convergence']['max_absolute_flux_difference']:.3g}; Δ data χ²={summary['integration_convergence']['data_chi2_difference']:.4g}. Exact source/input hashes and model arrays are preserved in this directory. Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_broad_leads.py --tic {tic}`.", '']
    if tic == 269994642:
        aliases = summary.get('integer_cycle_alias_tests', [])
        lines += ['The two observed eras initially imply 1,419 photometric cycles. Explicit 1,418/1,419/1,420-cycle fits on the same native photons give data χ² '+', '.join(f"{r['integer_cycle']}: {r['data_chi2']:.1f}" for r in aliases)+'. Adjacent aliases are therefore disfavored by the selected two-era shape fit. The formal period covariance conditions on this cycle mapping and stellar/shape assumptions; no general timing-variation model is fitted.', '']
    lines += [f"The doubled-period binary implementation also passes an independent native-time full-orbit test, adding separately calculated primary and secondary occultations: max flux discrepancy {summary['binary_model_verification']['max_absolute_flux_difference']:.2g}. That verification checks the physical period/phase mapping; it does not validate a binary classification.", '']
    # Make prose readable without changing numerical identifiers.
    (out / 'README.md').write_text('\n'.join(lines))


def run(tic, points):
    out = ROOT / f'candidates/broad_s1/tic{tic}_physical'
    out.mkdir(parents=True, exist_ok=True)
    settings_path = out/'fit_config.json'
    if (out/'best_fit.json').exists():
        old_settings = json.loads(settings_path.read_text()) if settings_path.exists() else {'supersample': 15}
        if old_settings['supersample'] != SUPERSAMPLE:
            # Preserve all initial scientific evidence when refining numerical
            # integration, rather than silently reusing incompatible profiles.
            archive = out/f"preliminary_{old_settings['supersample']}_subsamples"
            archive.mkdir(exist_ok=True)
            for item in list(out.iterdir()):
                if item.is_file():
                    shutil.move(str(item), str(archive/item.name))
    save(settings_path, {'supersample': SUPERSAMPLE, 'convergence_subsample': 61})
    data, prior, timing, provenance = inputs269() if tic == 269994642 else inputs238()
    provenance['sectors'] = np.unique(data['sector']).tolist()
    provenance['n_event_windows'] = len(set(zip(data['sector'], data['cycle'])))
    digest = hashlib.sha256()
    for k in sorted(data):
        digest.update(k.encode()); digest.update(np.ascontiguousarray(data[k]).tobytes())
    provenance['scientific_input_arrays_sha256'] = digest.hexdigest()
    model = LeadFit(data, prior, timing, supersample=SUPERSAMPLE)
    previous = out / 'best_fit.json'
    if previous.exists():
        previous_record = json.loads(previous.read_text())
        if previous_record.get('input_arrays_sha256') != provenance['scientific_input_arrays_sha256']:
            raise ValueError('Inputs changed: preserve old fit and explicitly revise before reuse')
        vec = np.asarray(previous_record['best']['vector'])
        result = None
    else:
        vec, result = multi_start(model)
        save(previous, {'input_arrays_sha256': provenance['scientific_input_arrays_sha256'], 'best': model.description(vec, result)})
    best = previous_record['best'] if previous.exists() and result is None else model.description(vec, result)
    print('BEST', tic, json.dumps(best), flush=True)
    rows = radius_profile(model, vec, out / 'radius_profile.json', points)
    if (out / 'alternative_fits.json').exists() and len(json.loads((out / 'alternative_fits.json').read_text())) >= (8 if tic == 238202207 else 7):
        alternatives = json.loads((out / 'alternative_fits.json').read_text())
    else:
        alternatives = fit_alternatives(tic, data, prior, timing, vec, out)
    binary = equal_binary(data, prior, timing, out)
    verification = verify_binary_model(tic, data, prior, timing, binary, out)
    unequal = unequal_binary(data, prior, timing, binary, out)
    parity = parity_robustness(tic, data, prior, timing, binary, out)
    refined = LeadFit(data, prior, timing, supersample=61)
    old, new = model.physical_model(vec), refined.physical_model(vec)
    convergence = {'max_absolute_flux_difference': float(np.max(abs(old-new))),
                   'data_chi2_difference': float(np.sum(refined.residual(vec, False)**2)-np.sum(model.residual(vec, False)**2)),
                   'fitted_subsamples': SUPERSAMPLE, 'comparison_subsamples': 61}
    summary = {'tic': tic, 'checked_utc': datetime.now(timezone.utc).isoformat(),
               'provenance': provenance, 'stellar_prior': prior, 'initial_ephemeris': timing,
               'best_circular_dark_companion': best, 'radius_profile_summary': profile_summary(rows),
               'equal_light_binary_alternative': binary,
               'binary_model_verification': verification,
               'unequal_brightness_binary_alternative': unequal,
               'parity_robustness': parity,
               'alternatives': alternatives, 'integration_convergence': convergence,
               'software': {pkg: version(pkg) for pkg in ['numpy', 'scipy', 'astropy', 'batman-package']},
               'source_sha256': sha(Path(__file__)), 'transit_library_sha256': sha(ROOT/'src/fit_370247818.py'),
               'assessment': 'Unresolved recurring transit-like signal; no mass classification or accepted new exoplanet.'}
    save(out / 'fit_results.json', summary)
    if tic == 269994642:
        summary['integer_cycle_alias_tests'] = cycle_alias_tests(data, prior, timing, best, out)
        save(out/'fit_results.json', summary)
    save_plot(tic, data, model, best, rows, alternatives, out)
    report(tic, summary, out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tic', type=int, choices=[269994642, 238202207], nargs='+', default=[269994642, 238202207])
    parser.add_argument('--profile-points', type=int, default=38)
    args = parser.parse_args()
    for tic in args.tic:
        run(tic, args.profile_points)


if __name__ == '__main__':
    main()
