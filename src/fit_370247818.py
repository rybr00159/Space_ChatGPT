"""Conditional, finite-exposure limb-darkened fits of TIC 370247818.

Radius ratio and impact parameter remain broad enough for grazing stellar
companions. Profiles and sector bootstrap are diagnostics, not validation or
posterior probabilities of planetary mass. Never edit the active survey.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.optimize import brentq
from importlib.metadata import version
import batman
from astropy.constants import G, M_sun, R_sun
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'candidates/broad_s1/tic370247818_multisector'
OUT = ROOT / 'candidates/broad_s1/tic370247818_physical'
SEED = 370247818
RJUP_RSUN = .102763


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def inputs():
    record = json.loads((INPUT / 'multisector.json').read_text())
    timing = record['event_timing_fit']['fitted_ephemeris']
    archive = np.load(INPUT / 'multisector_arrays.npz')
    tic_rows = json.loads((ROOT / 'data/broad_catalogs/tic370247818_tic_cone.json').read_text())['data']
    star = next(row for row in tic_rows if int(row['ID']) == SEED)
    chunks, exposures, sources = [], {}, []
    for source in record['sources']:
        path = Path(source['path'])
        if sha(path) != source['sha256']:
            raise ValueError('Raw FITS SHA mismatch; preserve and investigate')
        sector = source['sector']
        with fits.open(path, memmap=False) as hd:
            h = hd[1].header
            if h.get('TIMESYS') != 'TDB' or abs(h.get('BJDREFI', 0) + h.get('BJDREFF', 0) - 2457000) > 1e-8:
                raise ValueError('Unexpected TESS time system/reference')
            exposure = float(h['TIMEDEL'])
        exposures[str(sector)] = exposure
        t, f, e = (archive[f's{sector}_local_{name}'] for name in ('t', 'f', 'e'))
        cycle = archive[f's{sector}_local_cycle']
        good = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
        chunks.append({'t': t[good], 'f': f[good], 'e': e[good], 'cycle': cycle[good],
                       'sector': np.full(good.sum(), sector, dtype=int),
                       'exposure': np.full(good.sum(), exposure)})
        sources.append({**source, 'exposure_days_exact_TIMDEL': exposure})
    data = {key: np.concatenate([c[key] for c in chunks]) for key in chunks[0]}
    prior = {'radius_rsun': float(star['rad']), 'radius_error_rsun': float(star['e_rad']),
             'mass_msun': float(star['mass']), 'mass_error_msun': float(star['e_mass']),
             'teff_K': float(star['Teff']), 'teff_error_K': float(star['e_Teff']),
             'catalog_rho_solar': float(star['rho']), 'catalog_rho_error_solar': float(star['e_rho']),
             'limb_darkening_u1_mean': .30, 'limb_darkening_u1_sigma': .15,
             'limb_darkening_u2_mean': .25, 'limb_darkening_u2_sigma': .20,
             'limb_darkening_source': 'Broad approximate TESS quadratic priors appropriate to ~6380 K; not atmosphere-grid interpolation.',
             'stellar_source': 'Cached current MAST TIC cone result, target ID 370247818'}
    provenance = {'local_arrays_sha256': sha(INPUT / 'multisector_arrays.npz'),
                  'timing_json_sha256': sha(INPUT / 'multisector.json'),
                  'tic_cone_sha256': sha(ROOT / 'data/broad_catalogs/tic370247818_tic_cone.json'),
                  'gaia_stellar_sha256': sha(ROOT / 'data/broad_catalogs/tic370247818_gaia_stellar.csv'),
                  'raw_sources': sources, 'exposures_days': exposures,
                  'n_points': len(data['t']), 'sectors': np.unique(data['sector']).tolist(),
                  'n_event_windows': len(set(zip(data['sector'], data['cycle'])))}
    digest = hashlib.sha256()
    for key in sorted(data):
        digest.update(key.encode())
        digest.update(np.ascontiguousarray(data[key]).tobytes())
    provenance['scientific_input_arrays_sha256'] = digest.hexdigest()
    return data, prior, timing, provenance


def blackbody_photon_ratio(temperature, stellar_temperature, red_edge_um=1.1):
    """Approximate top-hat 0.6--red-edge micron photon-weighted contrast.

    No claim to an actual TESS response or model stellar atmosphere. A
    separate longest-wavelength contrast supplies a more generous comparison.
    """
    wave = np.linspace(.6e-6, red_edge_um * 1e-6, 401)
    temperature = np.atleast_1d(temperature).astype(float)
    hc_k = .014387768775039337
    argument = hc_k / (wave[None, :] * np.maximum(temperature[:, None], .01))
    photons = wave[None, :]**-4 / np.expm1(np.minimum(argument, 700.))
    star = wave**-4 / np.expm1(hc_k / (wave * stellar_temperature))
    return np.trapezoid(photons, wave, axis=1) / np.trapezoid(star, wave)


def circle_overlap_fraction(k, distance):
    if distance >= 1 + k:
        return 0.
    if distance <= abs(1 - k):
        return min(1., 1 / k**2)
    x = np.clip((distance**2 + 1 - k**2) / (2 * distance), -1, 1)
    y = np.clip((distance**2 + k**2 - 1) / (2 * distance * k), -1, 1)
    radicand = (-distance + 1 + k) * (distance + 1 - k) * (distance - 1 + k) * (distance + 1 + k)
    area = np.arccos(x) + k**2 * np.arccos(y) - .5 * np.sqrt(max(0., radicand))
    return float(area / (np.pi * k**2))


def secondary_budget(profile_rows, prior, best):
    source = json.loads((INPUT / 'multisector.json').read_text())
    phase_search = source.get('additional_phase_search', {})
    pooled = source.get('pooled_additional_phase_secondaries', {})
    if not pooled:
        return {'status': 'Awaiting independently evaluated off-phase secondary',
                'reported_initial_phase': .47546039035,
                'no_planet_classification': True}
    # Preserve the full independent measurements and their exact data source.
    budget = {'status': 'Off-phase secondary supported in independent raw/aperture reductions',
              'source_json_sha256': sha(INPUT / 'multisector.json'),
              'additional_phase_search': phase_search, 'pooled_secondary_measurements': pooled,
              'assumptions': ['Approximate blackbody star/companion, top-hat photon weighting 0.6--1.1 microns, not actual TESS response.',
                  'Steady passive irradiation; thermal surface temperature is capped by the zero-redistribution substellar temperature.',
                  'Lambert scattering couples Bond albedo A_B and geometric albedo A_g=2A_B/3.',
                  'Most generous comparison sums independent A_g<=1 reflected and zero-Bond-albedo red-edge blackbody thermal limits, although those extremes need not be jointly realizable.',
                  'Full occultation assumed for the planet-light budget; grazing would reduce eclipse depth.',
                  'Measured secondary is a local box-averaged depth; ingress/exposure loss can make the true peak deeper, so inferred uniform-disc brightness temperatures are approximate.',
                  'Circular irradiation distance and target star identity assumed; eccentric/internal-luminosity alternatives remain conditional.',
                  'An occultation measures brightness, not companion mass; the optical secondary does not by itself prove a stellar mass.'],
              'planet_secondary_budget_by_radius_ratio': [], 'uniform_companion_brightness_examples': []}
    # Explicit, generous 2 Rjup illustrative planet radius cutoff; it is not a
    # universal physical or mass classification limit.
    budget['illustrative_2_rjup_radius_ratio'] = 2 * RJUP_RSUN / prior['radius_rsun']
    bond = np.linspace(0, 1, 501)
    for row in profile_rows:
        k, ar = row['radius_ratio'], row['a_over_rstar']
        t_ss = prior['teff_K'] / np.sqrt(ar)
        t_upper = t_ss * (1 - bond)**.25
        thermal = k**2 * blackbody_photon_ratio(t_upper, prior['teff_K'])
        reflected = 2 / 3 * bond * (k / ar)**2
        coupled = thermal + reflected
        occulted = circle_overlap_fraction(k, row['impact_parameter'])
        partial = occulted * thermal + min(1., 1.5 * occulted) * reflected
        best_albedo = int(np.argmax(coupled))
        hc_k = .014387768775039337
        red = 1.1e-6
        red_edge_ratio = np.expm1(hc_k / (red * prior['teff_K'])) / np.expm1(hc_k / (red * t_ss))
        upper = (k / ar)**2 + k**2 * red_edge_ratio
        budget['planet_secondary_budget_by_radius_ratio'].append({
            'radius_ratio': k, 'a_over_rstar_from_primary_profile': ar,
            'conditional_radius_rjup': row['conditional_companion_radius_rjup'],
            'max_substellar_temperature_K': float(t_ss),
            'max_Lambert_coupled_thermal_plus_reflection_ppm': float(coupled.max() * 1e6),
            'max_coupled_Bond_albedo': float(bond[best_albedo]),
            'same_impact_geometric_disc_occulted_fraction': occulted,
            'same_impact_partial_Lambert_plus_uniform_substellar_thermal_upper_ppm': float(partial.max() * 1e6),
            'very_generous_independent_extremes_rededge_envelope_ppm': float(upper * 1e6),
            'primary_penalized_chi2': row['penalized_chi2']})
    # Under a same-impact, uniform-companion illustrative eclipse model,
    # infer the brightness temperature required for the independent raw-SAP
    # secondary. Do not assert this is a fitted atmosphere or mass estimate.
    measure = pooled.get('SAP_FLUX', pooled.get('KSPSAP_FLUX'))
    if isinstance(measure, dict):
        observed = measure.get('depth')
        observed_error = measure.get('depth_err')
    else:
        observed, observed_error = None, None
    budget['adopted_raw_secondary_depth'] = observed
    budget['adopted_raw_secondary_error'] = observed_error
    if observed is not None and observed > 0:
        for desired in (.12, .20, best['radius_ratio'], .40):
            row = min(profile_rows, key=lambda r: abs(r['radius_ratio'] - desired))
            k, b = row['radius_ratio'], row['impact_parameter']
            occulted_fraction = circle_overlap_fraction(k, b)
            if occulted_fraction <= observed:
                continue
            required_ratio = observed / (k**2 * (occulted_fraction - observed))
            contrast = lambda temperature: float(blackbody_photon_ratio(temperature, prior['teff_K'])[0]) - required_ratio
            temperature = brentq(contrast, 500., 12000.)
            total_light_ratio = k**2 * required_ratio
            budget['uniform_companion_brightness_examples'].append({
                'radius_ratio': k, 'impact_parameter_same_for_secondary': b,
                'companion_radius_rsun_conditional': row['conditional_companion_radius_rsun'],
                'geometric_companion_disc_occulted_fraction': occulted_fraction,
                'companion_to_star_out_of_eclipse_light_ratio': float(total_light_ratio),
                'required_top_hat_blackbody_brightness_temperature_K': float(temperature),
                'warning': 'Same-impact approximation, uniform companion, no atmospheric spectrum; conditional EB plausibility illustration, not a joint binary fit.'})
        planet_rows = [r for r in budget['planet_secondary_budget_by_radius_ratio']
                       if r['radius_ratio'] <= budget['illustrative_2_rjup_radius_ratio']]
        budget['illustrative_2rjup_max_budget'] = {
            'coupled_Lambert_ppm': max(r['max_Lambert_coupled_thermal_plus_reflection_ppm'] for r in planet_rows),
            'independent_extremes_rededge_ppm': max(r['very_generous_independent_extremes_rededge_envelope_ppm'] for r in planet_rows),
            'observed_raw_secondary_ppm': observed * 1e6, 'observed_raw_secondary_error_ppm': observed_error * 1e6,
            'largest_sampled_profile_k': max(r['radius_ratio'] for r in planet_rows),
            'warning': 'Maximum over sampled primary profile up to illustrative 2Rjup; not a marginalized or universal planetary luminosity bound.'}
    budget['small_e_e_cos_omega_from_phase_0_47546039035'] = float((.47546039035 - .5) * np.pi / 2)
    budget['best_primary_max_e_sin_omega_for_secondary_overlap_approx'] = float(
        (1 + best['radius_ratio'] - best['impact_parameter']) /
        (1 + best['radius_ratio'] + best['impact_parameter']))
    budget['assessment'] = 'Grazing eclipsing-binary explanation favored; no credible planet candidate classification from these conditional physical fits.'
    return budget


def limb(q1, q2):
    root = np.sqrt(q1)
    return 2 * root * q2, root * (1 - 2 * q2)


class TransitFit:
    """Seven physical/timing parameters, profiled per-sector constant baselines."""
    names = ('radius_ratio', 'impact_fraction_of_1plusk', 'log_a_over_rstar',
             'limb_q1', 'limb_q2', 'epoch_shift_days', 'period_offset_1e5days')
    bounds = (np.array([.03, .001, np.log(3.), .001, .001, -.025, -20.]),
              np.array([.5, .9999, np.log(30.), .999, .999, .025, 20.]))

    def __init__(self, data, prior, timing, eccentricity=0., omega=90.,
                 dilution=0., companion_mass_ratio=0., density_prior=True, supersample=25):
        self.data, self.prior, self.timing = data, prior, timing
        self.eccentricity, self.omega, self.dilution = eccentricity, omega, dilution
        self.mass_ratio, self.use_density_prior = companion_mass_ratio, density_prior
        self.supersample = supersample
        self.period0, self.epoch0 = timing['period'], timing['epoch_btjd']
        ar = ((G.value * prior['mass_msun'] * M_sun.value * (self.period0 * 86400)**2
               / (4 * np.pi**2)) ** (1 / 3) / (prior['radius_rsun'] * R_sun.value))
        self.loga_mean = np.log(ar * (1 + companion_mass_ratio) ** (1 / 3))
        propagated = np.hypot(prior['mass_error_msun'] / prior['mass_msun'] / 3,
                              prior['radius_error_rsun'] / prior['radius_rsun'])
        self.loga_sigma = max(propagated, prior['catalog_rho_error_solar'] / prior['catalog_rho_solar'] / 3)
        self.params = batman.TransitParams()
        self.params.t0, self.params.per = 0., self.period0
        self.params.rp, self.params.a, self.params.inc = .12, np.exp(self.loga_mean), 81.
        self.params.ecc, self.params.w = eccentricity, omega
        self.params.u, self.params.limb_dark = [.3, .25], 'quadratic'
        # Fold at the reference ephemeris but keep original cycle numbers;
        # period changes alter each event's center without rebinning photons.
        self.cycle = np.round((data['t'] - self.epoch0) / self.period0).astype(int)
        self.fold = data['t'] - self.epoch0 - self.cycle * self.period0
        self.sectors, self.sector_idx = np.unique(data['sector'], return_inverse=True)
        self.weight = 1 / data['e'] ** 2
        self.sector_weight = np.bincount(self.sector_idx, weights=self.weight)
        self.models = []
        for exp in np.unique(data['exposure']):
            ii = np.flatnonzero(data['exposure'] == exp)
            # Times passed to batman are replaced in-place in model objects
            # below to incorporate fitted long-baseline period changes.
            model = batman.TransitModel(self.params, self.fold[ii], supersample_factor=supersample,
                                        exp_time=float(exp), nthreads=1)
            self.models.append((ii, float(exp), model))

    def physical_model(self, vector):
        k, bf, loga, q1, q2, shift, dp = vector
        a, b = np.exp(loga), bf * (1 + k)
        e, w = self.eccentricity, np.deg2rad(self.omega)
        cosine = b / a * (1 + e * np.sin(w)) / (1 - e**2)
        if not (0 <= cosine < 1):
            return np.ones(len(self.fold))
        self.params.t0, self.params.per = shift, self.period0 + dp * 1e-5
        self.params.rp, self.params.a, self.params.inc = k, a, float(np.rad2deg(np.arccos(cosine)))
        self.params.u = list(limb(q1, q2))
        output = np.empty(len(self.fold))
        for ii, exposure, model in self.models:
            # batman owns flattened supersampling offsets in t_supersample.
            # Reconstruct a small model at adjusted photon times; preserves
            # the actual exposure durations with no phase bin approximation.
            times = self.fold[ii] - self.cycle[ii] * dp * 1e-5
            transit = batman.TransitModel(self.params, times, supersample_factor=self.supersample,
                                         exp_time=exposure, nthreads=1)
            output[ii] = transit.light_curve(self.params)
        return 1 - (1 - output) * (1 - self.dilution)

    def residual(self, vector, include_priors=True):
        model = self.physical_model(vector)
        offsets = np.bincount(self.sector_idx, weights=self.weight * (self.data['f'] - model)) / self.sector_weight
        residual = (self.data['f'] - model - offsets[self.sector_idx]) / self.data['e']
        if not include_priors:
            return residual
        k, bf, loga, q1, q2, shift, dp = vector
        u1, u2 = limb(q1, q2)
        priors = [(u1 - self.prior['limb_darkening_u1_mean']) / self.prior['limb_darkening_u1_sigma'],
                  (u2 - self.prior['limb_darkening_u2_mean']) / self.prior['limb_darkening_u2_sigma']]
        if self.use_density_prior:
            # Recompute Kepler a/R mean for fitted P; no period prior is added.
            mean = self.loga_mean + 2 / 3 * np.log((self.period0 + dp * 1e-5) / self.period0)
            priors.append((loga - mean) / self.loga_sigma)
        return np.r_[residual, priors]

    def fit(self, initial=None, fixed_k=None, max_nfev=180):
        if initial is None:
            initial = np.array([.13, .86, self.loga_mean, .55**2, .30 / (2 * .55), 0., 0.])
        initial = np.clip(initial, self.bounds[0] + 1e-8, self.bounds[1] - 1e-8)
        if fixed_k is None:
            result = least_squares(self.residual, initial, bounds=self.bounds, x_scale='jac',
                                   max_nfev=max_nfev, ftol=1e-8, xtol=1e-8, gtol=1e-7)
            vector = result.x
        else:
            def residual(reduced):
                return self.residual(np.r_[fixed_k, reduced])
            result = least_squares(residual, initial[1:], bounds=(self.bounds[0][1:], self.bounds[1][1:]),
                                   x_scale='jac', max_nfev=max_nfev, ftol=1e-8, xtol=1e-8, gtol=1e-7)
            vector = np.r_[fixed_k, result.x]
        return vector, result

    def description(self, vector, fit=None):
        k, bf, loga, q1, q2, shift, dp = vector
        a, b, p = np.exp(loga), bf * (1 + k), self.period0 + dp * 1e-5
        e, w = self.eccentricity, np.deg2rad(self.omega)
        inc = np.arccos(b / a * (1 + e * np.sin(w)) / (1 - e**2))
        # Circular analytic contact times. Eccentric diagnostic uses the
        # common small-duration local-velocity approximation, labeled below.
        velocity = np.sqrt(1 - e**2) / (1 + e * np.sin(w))
        dur = p / np.pi * np.arcsin(np.sqrt(max(0., (1 + k)**2 - b**2)) / (a * np.sin(inc))) * velocity
        full = (p / np.pi * np.arcsin(np.sqrt((1 - k)**2 - b**2) / (a * np.sin(inc))) * velocity
                if b < 1 - k else 0.)
        u1, u2 = limb(q1, q2)
        data_res = self.residual(vector, False)
        total_res = self.residual(vector, True)
        rcomp = k * self.prior['radius_rsun']
        output = {'vector': vector.tolist(), 'radius_ratio': float(k), 'impact_parameter': float(b),
                  'grazing': bool(b > 1 - k), 'a_over_rstar': float(a), 'inclination_degrees': float(np.rad2deg(inc)),
                  'period_days': float(p), 'epoch_btjd': float(self.epoch0 + shift),
                  'duration_T14_hours': float(dur * 24), 'duration_T23_hours': float(full * 24),
                  'u1': float(u1), 'u2': float(u2), 'conditional_companion_radius_rsun': float(rcomp),
                  'conditional_companion_radius_rjup': float(rcomp / RJUP_RSUN),
                  'data_chi2': float(np.sum(data_res**2)), 'prior_chi2': float(np.sum(total_res**2) - np.sum(data_res**2)),
                  'penalized_chi2': float(np.sum(total_res**2)), 'n_points': len(data_res),
                  'nominal_dof': len(data_res) - 7 - len(self.sectors),
                  'assumptions': {'eccentricity': e, 'omega_degrees': self.omega, 'dilution_fraction': self.dilution,
                                  'companion_mass_ratio': self.mass_ratio, 'density_prior': self.use_density_prior,
                                  'contact_duration_approximation': 'circular exact geometry' if e == 0 else 'eccentric local-velocity approximation'}}
        if fit is not None:
            covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
            uncertainties = np.sqrt(np.maximum(np.diag(covariance), 0))
            output['conditional_vector_covariance'] = covariance.tolist()
            output['conditional_vector_errors'] = uncertainties.tolist()
            output['converged'] = bool(fit.success)
            output['nfev'] = int(fit.nfev)
        return output


def best_fit(model, seeds=None):
    if seeds is None:
        seeds = []
        for k in (.09, .12, .20, .35):
            seeds.append([k, min(.97, .97 / (1 + k)), model.loga_mean, .55**2, .30 / (2 * .55), 0., 0.])
    fits_ = [model.fit(np.array(seed)) for seed in seeds]
    return min(fits_, key=lambda pair: np.sum(pair[1].fun**2))


def profiles(model, best, checkpoint):
    grid = np.unique(np.r_[np.geomspace(.055, .5, 55), .08, .09, .10, .11, .12, .13, .14, .15, best[0]])
    previous = json.loads(checkpoint.read_text()) if checkpoint.exists() else []
    done = {round(row['radius_ratio'], 12): row for row in previous}
    # Warm starts from both the best solution and neighboring profile values
    # prevent a single grazing local optimum deciding the entire profile.
    for k in sorted(grid, key=lambda value: abs(value - best[0])):
        if round(float(k), 12) in done:
            continue
        trials = [model.fit(best, fixed_k=float(k))]
        if done:
            near = min(done.values(), key=lambda row: abs(row['radius_ratio'] - k))
            trials.append(model.fit(np.array(near['vector']), fixed_k=float(k)))
        vector, fit = min(trials, key=lambda pair: np.sum(pair[1].fun**2))
        done[round(float(k), 12)] = model.description(vector)
        save(checkpoint, sorted(done.values(), key=lambda row: row['radius_ratio']))
    return sorted(done.values(), key=lambda row: row['radius_ratio'])


def bootstrap(data, prior, timing, best, count, checkpoint):
    rows = json.loads(checkpoint.read_text()) if checkpoint.exists() else []
    sectors = np.unique(data['sector'])
    for draw in range(len(rows), count):
        rng = np.random.default_rng(np.random.SeedSequence([SEED, draw]))
        chosen = rng.choice(sectors, len(sectors), replace=True)
        indices = np.concatenate([np.flatnonzero(data['sector'] == sector) for sector in chosen])
        # Keep duplicate sector occurrences distinct for their nuisance offsets.
        dd = {key: values[indices].copy() for key, values in data.items()}
        dd['sector'] = np.concatenate([np.full(np.sum(data['sector'] == sector), i, int)
                                       for i, sector in enumerate(chosen)])
        model = TransitFit(dd, prior, timing)
        vector, fit = model.fit(best, max_nfev=120)
        if not fit.success:
            vector, fit = model.fit(vector, max_nfev=180)
        rows.append({'draw': draw, 'resampled_sectors': chosen.tolist(), **model.description(vector, fit)})
        save(checkpoint, rows)
        if (draw + 1) % 20 == 0:
            print('BOOTSTRAP', draw + 1, flush=True)
    return rows


def write_report(summary):
    best = summary['circular_density_prior_best']
    profile = summary['radius_ratio_profile_summary']
    boot = summary['sector_bootstrap_summary']
    secondary = summary['secondary_light_budget']
    budget = secondary.get('illustrative_2rjup_max_budget', {})
    raw_secondary = secondary.get('pooled_secondary_measurements', {}).get('SAP_FLUX', {})
    small_secondary = secondary.get('pooled_secondary_measurements', {}).get('KSPSAP_FLUX_SML', {})
    if budget:
        secondary_text = f'''An independently evaluated eclipse occurs at phase
0.47546039, about 2.07 hours earlier than exact half-phase. It was missed by
the original half-phase-only check. Raw SAP depth is
{raw_secondary['depth'] * 1e6:.1f}±{raw_secondary['depth_err'] * 1e6:.1f} ppm
(conditional S/N {raw_secondary['snr']:.2f}), and small-aperture depth is
{small_secondary['depth'] * 1e6:.1f}±{small_secondary['depth_err'] * 1e6:.1f} ppm
(S/N {small_secondary['snr']:.2f}), consistent across all eight sectors.
The initial selected-phase statistic is not a blind-search false-alarm probability.
The small-e timing approximation implies e cos(omega)≈-0.03855 if the same
binary produces both eclipses; precise inference requires a joint orbit fit.
The detected secondary also requires geometric overlap. With the main
grazing primary geometry, the usual conjunction approximation requires
e sin(omega) ≲ {secondary['best_primary_max_e_sin_omega_for_secondary_overlap_approx']:.3f}.
The illustrative e=0.3 and e=0.6 periastron planet-like primary alternatives
predict no secondary overlap and therefore cannot describe both eclipses.
They demonstrate primary-only degeneracy, not allowed joint orbit solutions.

For a generous *illustrative* companion radius ≤2 Rjup, the sampled circular
primary profile gives a coupled Lambert reflected/thermal maximum of
{budget['coupled_Lambert_ppm']:.0f} ppm. An even more generous comparison sums
independent reflected (geometric albedo ≤1) and zero-albedo, uniform-substellar
temperature blackbody thermal extremes at the reddest assumed wavelength,
1.1 micron: {budget['independent_extremes_rededge_ppm']:.0f} ppm. These extremes
need not coexist physically; this comparison grants full occultation despite
the grazing geometry. The sampled maximum is at k={budget['largest_sampled_profile_k']:.3f},
slightly below the exact illustrative 2-Rjup cutoff, and is not a likelihood
profile jointly marginalized over orbital/stellar/atmospheric uncertainties.

The observed secondary therefore strains an ordinary passive, planet-sized
interpretation under the adopted geometry. It is readily compatible with
larger, intrinsically luminous grazing companions: illustrative uniform-disc
secondary inversions require optical brightness temperature near 3650 K over
k≈0.2--0.4. At the primary best-fit k≈0.278, only ≈14% of the companion disc
is eclipsed, with companion/star total light ratio ≈0.0093. This is a plausible
faint-star geometry, although a full luminous binary fit is still needed.

The budget assumes blackbody spectra with a 0.6--1.1-micron top-hat photon
weight, rather than an actual TESS throughput/atmosphere calculation. It couples
Bond and geometric albedo using Lambert scattering, assumes steady passive
irradiation, and uses circular separation. Internal luminosity, eccentricity,
non-Lambert scattering, unusual atmospheres, target misidentification, and
blending can change it. Two Rjup is a comparison radius, not an absolute planet
limit. Eclipse brightness and radius never establish companion mass.
The measured secondary is box-averaged: unresolved ingress and exposure
averaging can make its true peak deeper. The illustrative temperature/light
ratios are approximate and should be replaced by a joint luminous eclipse fit.
Nonetheless, the repeated off-phase secondary and grazing shape favor an
**eclipsing-binary explanation**, and no credible planetary candidate is promoted.
'''
    else:
        secondary_text = 'Off-phase secondary assessment remains pending; no planetary mass inference is justified.'
    text = f'''# TIC 370247818: physical transit/eclipse assessment

This is a conditional light-curve fit, **not a mass measurement or planet
validation**. A grazing stellar companion remains an allowed explanation.
No physical radius is inferred by taking the square root of observed depth.

Eight independent TESS sectors contribute {best['n_points']} locally normalized
samples around the recurring 3.522-day signal. Exact raw FITS TIMEDEL values
are used: 30 minutes for sectors 1/5/12/13 and 10 minutes for 27/28/32/39.
BATMAN quadratic limb darkening is numerically integrated across each exposure
with 25 sub-samples. Current TIC R=1.68401±0.08391 solar radii and
M=1.28±0.20331 solar masses anchor a conservative Gaussian log(a/R) prior;
the reported TIC density uncertainty is used if larger than independent M/R
propagation. Free, physically permitted limb-darkening coefficients have broad
approximate priors u1=0.30±0.15 and u2=0.25±0.20. These are approximate TESS
priors, not a matched atmosphere-grid calculation.

Fitted parameters include k=Rcomp/Rstar over 0.03--0.5, b/(1+k) over
0--0.9999, a/R, both limb-darkening coefficients, period, epoch, and profiled
per-sector constant baselines. The prior event-timing ephemeris initializes
the fit; it is **not reused as an independent prior** because it uses the same
photons. No transit phase binning is used in the likelihood. Errors inherited
from local fits include observational scatter and approximate baseline variance;
their remaining correlations are not represented by the diagonal likelihood.

Circular, undiluted, negligible-companion-mass best fit:

| Parameter | Value |
|---|---:|
| Period d | {best['period_days']:.9f} |
| Epoch BTJD | {best['epoch_btjd']:.6f} |
| Radius ratio k | {best['radius_ratio']:.5f} |
| Impact parameter b | {best['impact_parameter']:.5f} |
| a/Rstar | {best['a_over_rstar']:.4f} |
| Inclination deg | {best['inclination_degrees']:.4f} |
| T14 h | {best['duration_T14_hours']:.4f} |
| T23 h | {best['duration_T23_hours']:.4f} |
| u1, u2 | {best['u1']:.4f}, {best['u2']:.4f} |
| Companion radius, conditional Rsun | {best['conditional_companion_radius_rsun']:.4f} |
| Companion radius, conditional Rjup | {best['conditional_companion_radius_rjup']:.3f} |
| Data chi² / nominal dof | {best['data_chi2']:.1f} / {best['nominal_dof']} |
| Prior penalty chi² | {best['prior_chi2']:.3f} |

Grazing status is `{best['grazing']}`: b>1-k removes the fully interior phase
and makes depth a degenerate function of radius, limb darkening, and geometry.
The best-fit radius above assumes the target star is the transited source,
the stated stellar radius, no dilution, a dark companion, and circular orbit.
It is not evidence of planetary mass; planets, brown dwarfs, and faint stars
overlap in radius. The giant-planet/stellar-size k=0.12 boundary is a useful
shape comparison, not a mass or classification boundary.

Radius-ratio profile diagnostics (reoptimizing all other variables):
`{profile}`. Approximate Δpenalized-chi²≤1 and ≤3.84 ranges are conditional
profile contours, not calibrated Bayesian credible intervals or global
false-alarm probabilities. Reported penalties include stellar-density and
limb-darkening priors. Compare data chi² and prior penalties separately in the
JSON. Eight-sector bootstrap ({boot.get('draws', 0)} draws) is `{boot}`;
it probes sector sensitivity but only eight independent groups are available.
Bootstrap refits condition on the same stellar priors and use the best-fit
initialization; their narrow range is not a radius posterior and does not
replace the broad, upper-bound-reaching radius profile. Gaussian covariance
is particularly unreliable along the grazing degeneracy.

Alternative illustrative fits vary eccentricity (omega=90 deg), dilution,
stellar mass ratio, density-prior removal, and Gaia stellar density. Their
likelihoods and shapes remain in `fit_results.json`; no Bayes factors or
eccentricity posterior are claimed. A circular central transit would last
roughly twice the measured duration, supporting a high-impact interpretation
under circular geometry; eccentric velocity can also shorten a transit.
An unresolved luminous companion or contaminating star changes both depth and
stellar interpretation. Stellar-companion total mass changes Kepler a/R and
is explored separately without assuming a mass-radius relation.
Eccentric alternatives shown here fit the primary alone with omega=90 degrees;
they are illustrative velocity changes and do not jointly fit the observed
secondary phase. Exact transit/occultation timing and companion luminosity
must be included before interpreting an eccentric binary orbit.

{secondary_text}

Gaia median RV 13.65±3.68 km/s is a systemic-velocity estimate. The target has
no published amplitude in the cached table, so it cannot exclude a stellar
companion. Spectroscopy and source-localized higher-cadence transit photometry
remain needed. All fit uncertainties are conditional; no discovery claimed.

Reproduce/resume: `/workspace/openrouter-check-venv/bin/python src/fit_370247818.py`.
Checkpointed profile and bootstrap results, plots, model arrays, exact input
hashes, raw retrieval identifiers, and software versions are preserved here.
No raw datasets are deleted or downloaded; no paid calls are made.
'''
    (OUT / 'README.md').write_text(text)


def plot(data, model, best, profile_rows, alternative_rows, summary):
    vector = np.array(best['vector'])
    phase_hours = (data['t'] - best['epoch_btjd'] - np.round((data['t'] - best['epoch_btjd']) / best['period_days']) * best['period_days']) * 24
    prediction = model.physical_model(vector)
    residual = data['f'] - prediction
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    colors = {1: 'tab:blue', 5: 'tab:orange', 12: 'tab:green', 13: 'tab:red',
              27: 'tab:purple', 28: 'tab:brown', 32: 'tab:pink', 39: 'tab:gray'}
    for sector in np.unique(data['sector']):
        ii = data['sector'] == sector
        axes[0, 0].plot(phase_hours[ii], data['f'][ii], '.', ms=2, alpha=.35, color=colors[int(sector)], label=f'S{sector}')
    for lower, upper in zip(np.linspace(-3, 3, 61)[:-1], np.linspace(-3, 3, 61)[1:]):
        ii = (phase_hours >= lower) & (phase_hours < upper)
        if ii.sum():
            weights = 1 / data['e'][ii]**2
            axes[0, 0].errorbar(np.average(phase_hours[ii], weights=weights), np.average(data['f'][ii], weights=weights),
                              yerr=1 / np.sqrt(weights.sum()), fmt='ko', ms=3)
    # Illustrative exposure-integrated model curves on common photon phases.
    dense = np.linspace(-3, 3, 1000)
    dummy = {'t': model.epoch0 + dense / 24, 'f': np.ones(1000), 'e': np.ones(1000),
             'cycle': np.zeros(1000, int), 'sector': np.ones(1000, int), 'exposure': np.full(1000, 1 / 144)}
    dense_model = TransitFit(dummy, model.prior, model.timing)
    axes[0, 0].plot(dense, dense_model.physical_model(vector), 'k-', lw=2, label='Best, 10-min exposure')
    for target, color in ((.10, 'tab:blue'), (.20, 'tab:red'), (.35, 'tab:purple')):
        row = min(profile_rows, key=lambda p: abs(p['radius_ratio'] - target))
        axes[0, 0].plot(dense, dense_model.physical_model(np.array(row['vector'])), color=color,
                        lw=1.2, ls='--', label=f"Profile k={row['radius_ratio']:.3f}")
    axes[0, 0].set(xlim=(-2.5, 2.5), ylim=(.982, 1.007), xlabel='Hours from fitted midtransit', ylabel='Locally normalized flux')
    axes[0, 0].legend(ncol=3, fontsize=7)
    axes[1, 0].plot(phase_hours, residual * 1e6, '.', ms=2, alpha=.3)
    axes[1, 0].axhline(0, color='k', lw=.7)
    axes[1, 0].set(xlim=(-3, 3), ylim=(-10000, 10000), xlabel='Hours from midtransit', ylabel='Model residual ppm')
    x = [r['radius_ratio'] for r in profile_rows]
    minimum = min(r['penalized_chi2'] for r in profile_rows)
    axes[0, 1].plot(x, [r['penalized_chi2'] - minimum for r in profile_rows], 'o-', ms=3, label='Data + density/LD priors')
    data_minimum = min(r['data_chi2'] for r in profile_rows)
    axes[0, 1].plot(x, [r['data_chi2'] - data_minimum for r in profile_rows], '--', label='Data only along penalized profile')
    axes[0, 1].axhline(3.84, color='.5', ls=':')
    axes[0, 1].axvline(.12, color='tab:red', ls=':')
    axes[0, 1].set(xlabel='Radius ratio Rcomp/Rstar', ylabel='Profile Δχ²', ylim=(0, 30))
    axes[0, 1].legend(fontsize=7)
    axes[1, 1].plot(x, [r['impact_parameter'] for r in profile_rows], label='Profile b')
    axes[1, 1].plot(x, 1 - np.array(x), '--', label='Grazing boundary b=1-k')
    axes[1, 1].plot(x, 1 + np.array(x), ':', label='No-overlap boundary b=1+k')
    axes[1, 1].set(xlabel='Radius ratio', ylabel='Impact parameter')
    axes[1, 1].legend(fontsize=7)
    fig.suptitle('TIC 370247818: grazing/eclipsing shape degeneracy, conditional fits')
    fig.savefig(OUT / 'physical_fit.png', dpi=150)
    plt.close(fig)
    np.savez_compressed(OUT / 'physical_fit_arrays.npz', t=data['t'], f=data['f'], e=data['e'],
                        sector=data['sector'], exposure_days=data['exposure'], prediction=prediction,
                        phase_hours=phase_hours, profile_k=x,
                        profile_penalized_chi2=[r['penalized_chi2'] for r in profile_rows])
    light = summary.get('secondary_light_budget', {})
    budgets = light.get('planet_secondary_budget_by_radius_ratio', [])
    if budgets:
        fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
        xx = [r['radius_ratio'] for r in budgets]
        ax.plot(xx, [r['max_Lambert_coupled_thermal_plus_reflection_ppm'] for r in budgets],
                label='Coupled Lambert reflection + maximal thermal, full occultation')
        ax.plot(xx, [r['very_generous_independent_extremes_rededge_envelope_ppm'] for r in budgets],
                '--', label='Independent reflected/red-edge thermal extremes')
        ax.plot(xx, [r['same_impact_partial_Lambert_plus_uniform_substellar_thermal_upper_ppm'] for r in budgets],
                ':', label='Same-impact grazing, coupled illustrative upper envelope')
        observed = light['adopted_raw_secondary_depth'] * 1e6
        err = light['adopted_raw_secondary_error'] * 1e6
        ax.axhline(observed, color='tab:red', label='Observed raw-SAP secondary')
        ax.axhspan(observed - err, observed + err, color='tab:red', alpha=.15)
        ax.axvline(light['illustrative_2_rjup_radius_ratio'], color='.4', ls=':', label='Illustrative 2 Rjup radius')
        ax.set(xlim=(.08, .5), ylim=(50, 15000), yscale='log',
               xlabel='Companion/star radius ratio along primary profile', ylabel='Secondary depth ppm',
               title='TIC 370247818: conditional passive-planet light comparisons')
        ax.legend(fontsize=8)
        fig.savefig(OUT / 'secondary_light_budget.png', dpi=150)
        plt.close(fig)


def run(args):
    OUT.mkdir(parents=True, exist_ok=True)
    data, prior, timing, provenance = inputs()
    stamp = {'source_sha256': sha(Path(__file__)), 'provenance': provenance, 'prior': prior, 'timing_initialization': timing}
    checkpoint = OUT / 'experiment.json'
    if checkpoint.exists():
        old = json.loads(checkpoint.read_text())
        # Summary/code improvements may change source hash; changing actual
        # observations requires a separately preserved experiment directory.
        if old['provenance']['local_arrays_sha256'] != provenance['local_arrays_sha256']:
            raise ValueError('Input arrays changed; preserve old evidence and use a new experiment')
        if (old['prior'] != prior or old['timing_initialization']['period'] != timing['period']
                or old['timing_initialization']['epoch_btjd'] != timing['epoch_btjd']):
            raise ValueError('Stellar prior or timing initialization changed; preserve and version the experiment')
        if old['provenance']['timing_json_sha256'] != provenance['timing_json_sha256']:
            save(OUT / ('experiment_previous_' + old['provenance']['timing_json_sha256'][:12] + '.json'), old)
    save(checkpoint, stamp)
    model = TransitFit(data, prior, timing)
    best_path = OUT / 'best_circular.json'
    if best_path.exists():
        best = json.loads(best_path.read_text())
        vector = np.array(best['vector'])
    else:
        vector, fit = best_fit(model)
        best = model.description(vector, fit)
        save(best_path, best)
    print('BEST', json.dumps({key: best[key] for key in ('radius_ratio', 'impact_parameter', 'a_over_rstar',
          'grazing', 'duration_T14_hours', 'conditional_companion_radius_rjup', 'data_chi2', 'prior_chi2')}), flush=True)
    profile_rows = profiles(model, vector, OUT / 'radius_ratio_profile.json')
    minimum = min(row['penalized_chi2'] for row in profile_rows)
    ranges = {}
    for delta in (1., 3.84, 9.):
        allowed = [r['radius_ratio'] for r in profile_rows if r['penalized_chi2'] <= minimum + delta]
        ranges[str(delta)] = [min(allowed), max(allowed)] if allowed else None
    low = min((r for r in profile_rows if r['radius_ratio'] < .12), key=lambda r: r['penalized_chi2'])
    high = min((r for r in profile_rows if r['radius_ratio'] >= .12), key=lambda r: r['penalized_chi2'])
    profile_summary = {'conditional_k_ranges_by_delta_penalized_chi2': ranges,
                       'minimum_penalized_chi2': minimum,
                       'best_k_below_0_12': low, 'best_k_atleast_0_12': high,
                       'high_minus_low_penalized_chi2': high['penalized_chi2'] - low['penalized_chi2'],
                       'note': 'Conditional penalized likelihood profile; not probability of planet versus star.'}
    print('PROFILE_RANGES', json.dumps(ranges), 'stellar_comparison_delta', profile_summary['high_minus_low_penalized_chi2'], flush=True)
    alternatives_path = OUT / 'alternative_assumptions.json'
    alternatives = json.loads(alternatives_path.read_text()) if alternatives_path.exists() else []
    alternatives_by_label = {row['label']: row for row in alternatives}
    scenarios = [('eccentric_0.3_periastron', {'eccentricity': .3}),
                 ('eccentric_0.6_periastron', {'eccentricity': .6}),
                 ('dilution_0.25', {'dilution': .25}), ('dilution_0.50', {'dilution': .50}),
                 ('companion_mass_ratio_0.2', {'companion_mass_ratio': .2}),
                 ('companion_mass_ratio_0.5', {'companion_mass_ratio': .5}),
                 ('density_prior_removed', {'density_prior': False}), ('gaia_stellar', {})]
    for label, options in scenarios:
        if label in alternatives_by_label:
            continue
        pp = dict(prior)
        if label == 'gaia_stellar':
            pp.update(radius_rsun=1.7350514, mass_msun=1.3397111)
        alternate = TransitFit(data, pp, timing, **options)
        seeds = [vector, np.r_[.10, .7, alternate.loga_mean, vector[3:]]]
        vv, ff = best_fit(alternate, seeds)
        row = {'label': label, **alternate.description(vv, ff)}
        alternatives_by_label[label] = row
        save(alternatives_path, list(alternatives_by_label.values()))
        print('ALTERNATIVE', label, 'k', round(row['radius_ratio'], 4), 'b', round(row['impact_parameter'], 4),
              'chi2', round(row['penalized_chi2'], 2), flush=True)
    for row in alternatives_by_label.values():
        ecc = row['assumptions']['eccentricity']
        omega = np.deg2rad(row['assumptions']['omega_degrees'])
        secondary_b = row['impact_parameter'] * (1 + ecc * np.sin(omega)) / (1 - ecc * np.sin(omega))
        row['secondary_impact_parameter_conjunction_approx'] = float(secondary_b)
        row['secondary_geometric_overlap_possible_approx'] = bool(secondary_b < 1 + row['radius_ratio'])
        row['joint_eclipse_warning'] = ('Primary-only illustrative fit: predicted secondary has no geometric overlap.'
            if not row['secondary_geometric_overlap_possible_approx'] else
            'Primary-only fit: secondary phase, shape and luminosity have not been jointly modeled.')
    save(alternatives_path, list(alternatives_by_label.values()))
    boots = bootstrap(data, prior, timing, vector, args.bootstrap, OUT / 'sector_bootstrap.json')
    bsummary = {'draws': len(boots), 'converged_draws': sum(row['converged'] for row in boots)}
    for key in ('radius_ratio', 'impact_parameter', 'a_over_rstar', 'duration_T14_hours',
                'period_days', 'epoch_btjd', 'conditional_companion_radius_rjup'):
        bsummary[key + '_p16_median_p84'] = np.quantile([row[key] for row in boots], [.16, .5, .84]).tolist() if boots else None
    # Test supersampling convergence on the fitted solution, without altering
    # observation evidence or allowing unconstrained interpolation of exposures.
    dense_model = TransitFit(data, prior, timing, supersample=101)
    convergence = float(np.max(np.abs(model.physical_model(vector) - dense_model.physical_model(vector))) * 1e6)
    summary = {'tic': SEED, 'checked_utc': datetime.now(timezone.utc).isoformat(),
               **stamp, 'software_versions': {'python': platform.python_version(), 'numpy': np.__version__,
                 'scipy': scipy.__version__, 'batman': batman.__version__, 'batman-package': version('batman-package'),
                 'setuptools': version('setuptools')},
               'circular_density_prior_best': best, 'radius_ratio_profile_summary': profile_summary,
               'sector_bootstrap_summary': bsummary, 'alternative_assumptions': list(alternatives_by_label.values()),
               'supersampling_25_vs_101_max_difference_ppm': convergence,
               'secondary_light_budget': secondary_budget(profile_rows, prior, best),
               'radius_to_mass_classification': None, 'discovery_claim': False, 'paid_calls': 0,
               'limitations': ['Dark-companion model; no direct companion luminosity or mass inference.',
                   'Conditional target-star identity, stellar radius, density and limb-darkening assumptions.',
                   'Circular main fit; eccentricity/dilution scenarios are illustrative, not marginalized posteriors.',
                   'Local normalized errors share baseline uncertainty; diagonal likelihood misses correlations.',
                   'Eight-sector bootstrap is a small-group sensitivity check, not guaranteed coverage.',
                   'Gaia median RV and its error do not constrain companion orbital RV amplitude.']}
    save(OUT / 'fit_results.json', summary)
    plot(data, model, best, profile_rows, alternatives, summary)
    write_report(summary)
    print('FINAL', json.dumps({'k': best['radius_ratio'], 'b': best['impact_parameter'], 'profile_ranges': ranges,
          'bootstrap': bsummary, 'integration_max_difference_ppm': convergence}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstrap', type=int, default=100)
    run(parser.parse_args())
