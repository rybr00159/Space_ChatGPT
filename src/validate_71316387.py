"""Reproducible, targeted audit of the TIC 71316387 alternating-eclipse claim.

Uses retained public QLP FITS and the catalog ephemeris, not a new period search.
Existing frozen broad-survey code is imported read-only for local baselines,
exact exposure integration, and individual-event timing measurements.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from scipy.stats import binomtest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import broad_vet as bv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'candidates/tic71316387'
PERIOD = .6518369465110034
EPOCH = 2067.327192476
WIDTH = .08  # Baseline mask, not the fitted eclipse duration.
COLUMNS = ['SAP_FLUX', 'KSPSAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']


def safe(v):
    if isinstance(v, dict):
        return {str(k): safe(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [safe(x) for x in v]
    if isinstance(v, np.ndarray):
        return safe(v.tolist())
    if isinstance(v, np.generic):
        return safe(v.item())
    if isinstance(v, float) and not np.isfinite(v):
        return None
    return v


def load_sources():
    rows = json.loads((OUT / 'observations.json').read_text())
    sources = []
    for row in rows:
        if row['status'] != 'retained':
            continue
        path = Path(row['path'])
        if not path.exists():
            path = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'tic71316387/observations' / path.name
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Raw input SHA256 mismatch')
        with fits.open(path, memmap=False) as hd:
            bv.verify_fits_checksums(path, hd)
            h, d = hd[1].header, hd[1].data
            if h.get('TIMESYS') != 'TDB' or h.get('BJDREFI', 0) + h.get('BJDREFF', 0) != 2457000:
                raise ValueError('Unexpected time coordinate')
            aliases = {k: k.replace('KSPSAP', 'DET') if k not in d.names else k
                       for k in COLUMNS + ['KSPSAP_FLUX_ERR']}
            arrays = {k: np.asarray(d[aliases.get(k, k)], float)
                      for k in ['TIME', 'QUALITY'] + COLUMNS + ['KSPSAP_FLUX_ERR']}
            cadence = float(h['TIMEDEL'])
            row['analysis_column_mapping'] = aliases
            row['analysis_exposure_days'] = cadence
        use = (arrays['QUALITY'] == 0) & np.isfinite(arrays['TIME'])
        arrays = {k: v[use] for k, v in arrays.items()}
        valid_error = np.isfinite(arrays['KSPSAP_FLUX_ERR']) & (arrays['KSPSAP_FLUX_ERR'] > 0)
        if not valid_error.any():
            # Public sector105 has an entirely NaN uncertainty column. Derive
            # an empirical scale; each event's OOT MAD and residual scatter
            # still inflate the fitted uncertainty. Do not silently drop it.
            delta = np.diff(arrays['SAP_FLUX'])
            noise = bv.mad(delta[np.isfinite(delta)]) / np.sqrt(2)
            if not np.isfinite(noise) or noise <= 0:
                raise ValueError('Cannot estimate missing uncertainty scale')
            arrays['KSPSAP_FLUX_ERR'] = np.full(len(arrays['TIME']), noise)
            row['error_fallback'] = dict(reason='No finite positive product errors',
                                         method='MAD of first differences divided by sqrt(2)',
                                         normalized_flux_scale=float(noise))
        sources.append((row, arrays, cadence))
    return sources


def integration_controls(sources, period, epoch, shape):
    """Deterministic reconstruction at real timestamps/gaps, without real noise.

    This checks parity/cadence bookkeeping and attenuation, not completeness or
    the probability that correlated noise produces alternation.
    """
    result = {}
    for name, depths in [('equal_depth', [.004, .004]),
                         ('quoted_alternating_depths', [.002487, .005293])]:
        injected = []
        for row, arrays, exposure in sources:
            t = arrays['TIME']
            cycle = np.rint((t - epoch) / period).astype(int)
            x = t - epoch - cycle * period
            q = bv.integrated_trapezoid(x, shape['duration_days'],
                                        shape['duration_days'] * shape['ingress_fraction'], exposure)
            f = 1 - np.asarray(depths)[cycle % 2] * q
            good = (np.isfinite(arrays['SAP_FLUX']) & (arrays['SAP_FLUX'] > 0)
                    & np.isfinite(arrays['KSPSAP_FLUX_ERR']) & (arrays['KSPSAP_FLUX_ERR'] > 0))
            f[~good] = np.nan
            modified = dict(arrays)
            modified['SAP_FLUX'] = f
            modified['KSPSAP_FLUX_ERR'] = np.full(len(t), .001)
            injected.append((row, modified, exposure))
        ev = event_depths(local_data(injected, 'SAP_FLUX', period, epoch), period, epoch, shape)
        stat = parity_fit(ev)
        stat['injected_even_odd_depths'] = depths
        expected = depths[0] - depths[1]
        stat['difference_bias_ppm'] = (stat['even_minus_odd_depth'] - expected) * 1e6
        if abs(stat['difference_bias_ppm']) > 1:
            raise AssertionError('Exposure/parity control failed at >1ppm')
        result[name] = stat
    return result


def local_data(sources, column, period, epoch, degree=2):
    result = []
    for row, arrays, cadence in sources:
        _, data = bv.local_quadratic(arrays['TIME'], arrays[column],
                                    arrays['KSPSAP_FLUX_ERR'], period, epoch,
                                    WIDTH, degree=degree)
        if data is None:
            continue
        data['sector'] = np.full(len(data['t']), row['sector'], int)
        data['exposure'] = np.full(len(data['t']), cadence)
        result.append(data)
    return result


def pooled_shape(data, period, epoch):
    values = {k: np.concatenate([v[k] for v in data]) for k in data[0]}
    x = values['t'] - epoch - values['cycle'] * period
    y, e = values['f'], values['e']
    sector, parity = values['sector'], values['cycle'] % 2
    groups = [(s, k) for s in sorted(set(sector)) for k in (0, 1)]
    indicators = np.column_stack([(sector == s) & (parity == k) for s, k in groups])

    def solve(par):
        shift, duration, fraction = par
        q = bv.integrated_trapezoid(x - shift, duration, fraction * duration,
                                    values['exposure'])
        design = np.column_stack([indicators, indicators * (-q[:, None])])
        a, b = design / e[:, None], y / e
        beta = np.linalg.lstsq(a, b, rcond=None)[0]
        return beta, design, (y - design @ beta) / e

    fit = least_squares(lambda par: solve(par)[2], [0., .065, .3],
                        bounds=([-.02, .025, .03], [.02, .11, .5]),
                        diff_step=1e-4, max_nfev=200)
    beta, design, residual = solve(fit.x)
    return dict(epoch_btjd=float(epoch + fit.x[0]),
                duration_days=float(fit.x[1]), ingress_fraction=float(fit.x[2]),
                depth=float(np.mean(beta[len(groups):])),
                group_depths=[dict(sector=int(s), parity=int(k), depth=float(d))
                              for (s, k), d in zip(groups, beta[len(groups):])],
                chi2=float(residual @ residual), n_points=len(y),
                success=bool(fit.success),
                note='Shared exposure-integrated trapezoid; separate depths and offsets per sector/parity. Not a physical radius fit.')


def event_depths(data, period, epoch, shape):
    events = []
    for part in data:
        for cycle in np.unique(part['cycle']):
            use = part['cycle'] == cycle
            x = part['t'][use] - epoch - cycle * period
            y, e = part['f'][use], part['e'][use]
            q = bv.integrated_trapezoid(x, shape['duration_days'],
                                        shape['duration_days'] * shape['ingress_fraction'],
                                        part['exposure'][use])
            if len(x) < 8 or np.count_nonzero(q > .25) < 2:
                continue
            design = np.column_stack([np.ones(len(x)), x, -q])
            a, b = design / e[:, None], y / e
            beta = np.linalg.lstsq(a, b, rcond=None)[0]
            residual = b - a @ beta
            chi2 = float(residual @ residual)
            covariance = np.linalg.pinv(a.T @ a) * max(1., chi2 / (len(x) - 3))
            events.append(dict(sector=int(part['sector'][use][0]), cycle=int(cycle),
                               parity=int(cycle % 2), time_btjd=float(epoch + cycle * period),
                               depth=float(beta[2]), depth_err=float(np.sqrt(covariance[2, 2])),
                               n_points=len(x), chi2=chi2, dof=len(x) - 3,
                               exposure_days=float(part['exposure'][use][0])))
    return events


def parity_fit(events):
    sectors = sorted({v['sector'] for v in events})
    sector = np.array([v['sector'] for v in events])
    parity = np.array([1 if v['parity'] == 0 else -1 for v in events])
    y = np.array([v['depth'] for v in events])
    e = np.array([v['depth_err'] for v in events])
    design = np.column_stack([sector == s for s in sectors] + [parity])
    a, b = design / e[:, None], y / e
    beta = np.linalg.lstsq(a, b, rcond=None)[0]
    residual = b - a @ beta
    inverse = np.linalg.pinv(a.T @ a)
    chi2 = float(residual @ residual)
    dof = len(y) - len(beta)
    inflated = inverse * max(1., chi2 / max(1, dof))
    cluster = None
    if len(sectors) > 1:
        meat = np.zeros_like(inverse)
        for s in sectors:
            score = a[sector == s].T @ residual[sector == s]
            meat += np.outer(score, score)
        cluster = inverse @ meat @ inverse * len(sectors) / (len(sectors) - 1) * (len(y) - 1) / max(1, dof)
    variance = inflated[-1, -1]
    if cluster is not None:
        variance = max(variance, cluster[-1, -1])
    delta, error = 2 * beta[-1], 2 * np.sqrt(variance)
    # For depth presentation only: exposure/error-weighted sector mean plus
    # common parity term. Absolute values depend on extraction and dilution.
    sector_weights = np.array([np.sum(1 / e[sector == s] ** 2) for s in sectors])
    mean = np.average(beta[:-1], weights=sector_weights)
    return dict(n_events=len(events), n_sectors=len(sectors), sectors=sectors,
                even_depth=float(mean + beta[-1]), odd_depth=float(mean - beta[-1]),
                even_minus_odd_depth=float(delta), difference_error=float(error),
                conditional_snr=float(abs(delta) / error),
                diagonal_difference_error=float(2 * np.sqrt(inverse[-1, -1])),
                scatter_inflated_difference_error=float(2 * np.sqrt(inflated[-1, -1])),
                cluster_difference_error=None if cluster is None else float(2 * np.sqrt(cluster[-1, -1])),
                chi2=chi2, dof=dof, sector_mean_depths=dict(zip(map(str, sectors), map(float, beta[:-1]))),
                note='Common alternating-depth coefficient after sector offsets. Larger scatter-inflated or sector-cluster error; few clusters limit normal-tail interpretation. Fixed hypothesis, not search false-alarm significance.')


def adjacent_pairs(events):
    lookup = {(v['sector'], v['cycle']): v for v in events}
    differences = []
    for v in events:
        if v['parity'] != 0:
            continue
        nxt = lookup.get((v['sector'], v['cycle'] + 1))
        if nxt is not None:
            differences.append(dict(sector=v['sector'], cycle=v['cycle'],
                                    odd_minus_even=nxt['depth'] - v['depth']))
    d = np.array([v['odd_minus_even'] for v in differences])
    rng = np.random.default_rng(71316387)
    samples = rng.choice(d, size=(10000, len(d)), replace=True).mean(axis=1)
    return dict(n_pairs=len(d), mean_odd_minus_even=float(d.mean()),
                empirical_standard_error=float(d.std(ddof=1) / np.sqrt(len(d))),
                bootstrap_95_interval=np.percentile(samples, [2.5, 97.5]).tolist(),
                n_positive=int((d > 0).sum()), sign_two_sided_p=float(binomtest(int((d > 0).sum()), len(d)).pvalue),
                pairs=differences,
                note='Unweighted adjacent complete even/odd pairs: removes slow baseline/dilution drift. Pair bootstrap/sign p assume independent pairs, not global search calibration.')


def main():
    sources = load_sources()
    first = local_data(sources, 'SAP_FLUX', PERIOD, EPOCH)
    initial_shape = pooled_shape(first, PERIOD, EPOCH)
    timing = bv.fit_event_timings(first, [int(v['sector'][0]) for v in first],
                                  PERIOD, initial_shape['epoch_btjd'], initial_shape)
    fitted = timing['fitted_ephemeris']
    p, epoch = fitted['period'], fitted['epoch_btjd']
    raw_data = local_data(sources, 'SAP_FLUX', p, epoch)
    shape = pooled_shape(raw_data, p, epoch)
    epoch = shape['epoch_btjd']
    tests, events, alternative = {}, {}, {}
    final_data = None
    for column in COLUMNS:
        data = local_data(sources, column, p, epoch)
        ev = event_depths(data, p, epoch, shape)
        events[column] = ev
        tests[column] = parity_fit(ev)
        if column == 'SAP_FLUX':
            final_data = data
    for degree in (0, 1):
        ev = event_depths(local_data(sources, 'SAP_FLUX', p, epoch, degree), p, epoch, shape)
        alternative[f'local_polynomial_degree_{degree}'] = parity_fit(ev)
    raw = events['SAP_FLUX']
    per_sector = {str(s): parity_fit([v for v in raw if v['sector'] == s]) for s in tests['SAP_FLUX']['sectors']}
    subsets = {name: parity_fit([v for v in raw if v['sector'] in ids]) for name, ids in
               [('raven_sectors_1_28', [1, 28]), ('independent_later_sectors', [68, 95, 105])]}
    result = dict(tic=71316387, checked_utc=datetime.now(timezone.utc).isoformat(),
                  catalog_ephemeris=dict(period_days=PERIOD, epoch_btjd=EPOCH),
                  adopted_ephemeris=dict(recurrence_days=p, epoch_btjd=epoch,
                                         proposed_binary_period_days=2 * p),
                  initial_shape=initial_shape, timing=timing, fitted_shape=shape,
                  parity_tests=tests, raw_per_sector=per_sector, subsets=subsets,
                  alternative_baselines=alternative, adjacent_pair_test=adjacent_pairs(raw),
                  noiseless_native_cadence_controls=integration_controls(sources, p, epoch, shape),
                  leave_one_sector_out={str(s): parity_fit([v for v in raw if v['sector'] != s]) for s in tests['SAP_FLUX']['sectors']},
                  events=events, source_manifest_sha256=hashlib.sha256((OUT/'observations.json').read_bytes()).hexdigest(),
                  source_metadata=[row for row, arrays, exposure in sources],
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  imported_baseline_method=bv.method_identity(),
                  claim=dict(deep_depth_ppm=5293, shallow_depth_ppm=2487, sigma=19.6,
                             implied_difference_error_ppm=(5293-2487)/19.6),
                  limitations=['QLP SAP and corrected apertures share images and are not independent experiments.',
                               'Shared shape and ephemeris are conditional; statistical errors omit full PRF and calibration uncertainty.',
                               'Alternating eclipses favor a doubled recurrence period; they alone do not identify the eclipsing star or establish its mass.',
                               'Exact claimed depths and significance require the original extraction, event mask and noise method.'])
    (OUT/'validation.json').write_text(json.dumps(safe(result), indent=2, allow_nan=False)+'\n')
    with (OUT/'event_depths.csv').open('w') as f:
        w=csv.DictWriter(f, fieldnames=['column']+list(raw[0]));w.writeheader()
        for column, rows in events.items():
            w.writerows(dict(column=column,**r) for r in rows)
    plot(final_data, raw, p, epoch, shape, per_sector)
    print(json.dumps(safe(dict(ephemeris=result['adopted_ephemeris'], shape=shape,
                              tests=tests, subsets=subsets, per_sector=per_sector,
                              pairs={k:v for k,v in result['adjacent_pair_test'].items() if k!='pairs'})),indent=2))


def plot(data, events, p, epoch, shape, per_sector):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    colors = {0: 'tab:blue', 1: 'tab:orange'}
    for k in (0, 1):
        xx, yy = [], []
        for v in data:
            use = v['cycle'] % 2 == k
            x = (v['t'][use] - epoch - v['cycle'][use] * p) * 24
            xx.extend(x); yy.extend((1-v['f'][use])*1e6)
        x,y=np.array(xx),np.array(yy)
        axes[0,0].scatter(x,y,s=2,alpha=.13,color=colors[k])
        edges=np.linspace(-2.5,2.5,51)
        bx,by,be=[],[],[]
        for a,b in zip(edges[:-1],edges[1:]):
            vals=y[(x>=a)&(x<b)]
            if len(vals)>2:
                bx.append((a+b)/2);by.append(vals.mean());be.append(vals.std(ddof=1)/np.sqrt(len(vals)))
        axes[0,0].errorbar(bx,by,be,fmt='o',ms=3,color=colors[k],label='Even (shallow)' if k==0 else 'Odd (deep)')
        ev=[v for v in events if v['parity']==k]
        axes[0,1].errorbar([v['time_btjd'] for v in ev],[v['depth']*1e6 for v in ev],
                           [v['depth_err']*1e6 for v in ev],fmt='.',alpha=.7,color=colors[k])
    axes[0,0].set(xlim=(-2.5,2.5),ylim=(-1500,11000),xlabel='Hours from each catalog recurrence',ylabel='Locally normalized flux loss (ppm)')
    axes[0,0].legend()
    axes[0,1].set(xlabel='BTJD (TDB)',ylabel='Exposure-integrated event depth (ppm)')
    for s,stats in per_sector.items():
        axes[1,0].errorbar(int(s),-stats['even_minus_odd_depth']*1e6,stats['difference_error']*1e6,fmt='o',color='black')
    axes[1,0].axhline(0,color='grey');axes[1,0].set(xlabel='TESS sector',ylabel='Deep minus shallow depth (ppm)')
    for v in data:
        x=((v['t']-epoch)/(2*p))%1
        axes[1,1].scatter(x,(1-v['f'])*1e6,s=2,alpha=.15,color='black')
    axes[1,1].set(xlabel=f'Phase at {2*p:.9f} days',ylabel='Flux loss (ppm)',ylim=(-1500,11000))
    fig.suptitle('TIC 71316387: targeted alternating-eclipse audit (raw QLP SAP)')
    fig.tight_layout();fig.savefig(OUT/'alternating_eclipses.png',dpi=170);plt.close(fig)


if __name__ == '__main__':
    main()
