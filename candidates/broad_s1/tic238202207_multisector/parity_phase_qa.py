"""Supplemental fixed-ephemeris parity/phase diagnostics; frozen v3 unchanged.

Run: python candidates/broad_s1/tic238202207_multisector/parity_phase_qa.py
The input source hashes are verified before analysis, and the numeric v3 report
and arrays remain read-only. Alternating depth is fit after sector-dependent
depths; a double-period half eclipse alone is not evidence against a planet.
"""
from pathlib import Path
import hashlib
from datetime import datetime, timezone
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
import broad_vet as bv

HERE = Path(__file__).resolve().parent


def alternating_depth(events):
    sectors = sorted({v['sector'] for v in events})
    sector = np.array([v['sector'] for v in events])
    parity = np.array([1 if v['cycle'] % 2 == 0 else -1 for v in events])
    depth = np.array([v['depth'] for v in events])
    error = np.array([v['depth_err'] for v in events])
    # Separate sector offsets absorb dilution/camera-dependent mean depths.
    design = np.column_stack([sector == s for s in sectors] + [parity])
    a, y = design / error[:, None], depth / error
    beta = np.linalg.lstsq(a, y, rcond=None)[0]
    residual = y - a @ beta
    inverse = np.linalg.pinv(a.T @ a)
    chi2 = float(np.sum(residual ** 2))
    covariance = inverse * max(1., chi2 / max(1, len(y) - len(beta)))
    # Sector-cluster sandwich covariance allows within-sector event dependence.
    meat = np.zeros_like(inverse)
    for s in sectors:
        use = sector == s
        score = a[use].T @ residual[use]
        meat += np.outer(score, score)
    correction = len(sectors) / max(1, len(sectors) - 1) * (len(y) - 1) / max(1, len(y) - len(beta))
    cluster = inverse @ meat @ inverse * correction
    err = 2 * np.sqrt(max(covariance[-1, -1], cluster[-1, -1]))
    difference = 2 * float(beta[-1])
    return dict(n_events=len(events), n_sectors=len(sectors),
                even_minus_odd_depth=difference, difference_error=err,
                difference_conditional_snr=abs(difference) / err,
                sector_mean_depths=dict(zip(map(str, sectors), map(float, beta[:-1]))),
                coefficient_order=[f'sector_{s}_mean' for s in sectors] + ['parity_signed_half_difference'],
                coefficients=beta.tolist(), covariance=covariance.tolist(),
                sector_cluster_covariance=cluster.tolist(), chi2=chi2,
                dof=len(y) - len(beta),
                note='Fixed timing ephemeris and eclipse width; per-event local-quadratic depths. Error is the larger of inflated diagonal and sector-cluster estimates. Conditional, not survey false-alarm probability.')


def main():
    report_path = HERE / 'multisector.json'
    report = json.loads(report_path.read_text())
    bv_sha = bv.method_identity()['source_sha256']
    assert bv_sha == report['vetting_method']['source_sha256']
    ephemeris = report['event_timing_fit']['fitted_ephemeris']
    p, epoch = ephemeris['period'], ephemeris['epoch_btjd']
    width = report['joint_ephemeris']['duration_days']
    columns = ['KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']
    all_events = {col: [] for col in columns}
    integrated_events = {col: [] for col in columns}
    harmonics = []
    shape = report['approximate_pooled_shape']
    for row in report['sources']:
        arrays, checks = bv.read_source(row)
        good = arrays['QUALITY'] == 0
        t, e = arrays['TIME'][good], arrays['KSPSAP_FLUX_ERR'][good]
        for col in columns:
            stats, local = bv.local_quadratic(t, arrays[col][good], e, p, epoch, width)
            all_events[col].extend(dict(sector=row['sector'], **v) for v in stats['events'])
            if local is not None:
                exposure = float(np.median(np.diff(t[np.isfinite(t)])))
                for cycle in np.unique(local['cycle']):
                    use = local['cycle'] == cycle
                    x = local['t'][use] - epoch - cycle * p
                    template = bv.integrated_trapezoid(x, shape['duration_days'],
                                                       shape['duration_days'] * shape['ingress_fraction'], exposure)
                    design = np.column_stack([np.ones(len(x)), -template])
                    a, y = design / local['e'][use, None], local['f'][use] / local['e'][use]
                    beta = np.linalg.lstsq(a, y, rcond=None)[0]
                    chi2 = float(np.sum((y - a @ beta) ** 2))
                    covariance = np.linalg.pinv(a.T @ a) * max(1., chi2 / max(1, len(x) - 2))
                    integrated_events[col].append(dict(sector=row['sector'], cycle=int(cycle),
                                                        time=float(epoch + cycle * p),
                                                        depth=float(beta[1]), depth_err=float(np.sqrt(covariance[1, 1]))))
        for factor in [1, 2]:
            test = bv.variability_diagnostics(t, arrays['SAP_FLUX'][good], e,
                                             factor * p, epoch, width)
            harmonics.append(dict(sector=row['sector'], period_factor=factor, **test))
    parity = {col: alternating_depth(events) for col, events in all_events.items()}
    integrated_parity = {col: alternating_depth(events) for col, events in integrated_events.items()}
    raw_events = integrated_events['SAP_FLUX']
    sector_parity = {str(s): alternating_depth([v for v in raw_events if v['sector'] == s])
                     for s in sorted({v['sector'] for v in raw_events})}
    splits = {name: alternating_depth([v for v in raw_events if select(v['sector'])])
              for name, select in [('early_30min', lambda s: s < 27), ('later_10min', lambda s: s >= 27)]}
    leave_one_out = {str(s): alternating_depth([v for v in raw_events if v['sector'] != s])
                     for s in sorted({v['sector'] for v in raw_events})}
    positive = sum(v['even_minus_odd_depth'] > 0 for v in sector_parity.values())
    sign_pvalue = float(binomtest(positive, len(sector_parity), .5).pvalue)
    timing_events = [dict(sector=v['sector'], cycle=v['cycle'],
                          depth=v['time_btjd'] - epoch - v['cycle'] * p,
                          depth_err=v['time_error_days'])
                     for v in report['event_timing_fit']['events']]
    timing_parity = alternating_depth(timing_events)
    result = dict(tic=report['tic'], checked_utc=datetime.now(timezone.utc).isoformat(),
                  ephemeris=dict(period_days=p, epoch_btjd=epoch),
                  parity_tests=parity, event_depths=all_events,
                  integrated_template_parity_tests=integrated_parity,
                  integrated_template_event_depths=integrated_events,
                  raw_integrated_per_sector_parity=sector_parity,
                  raw_integrated_cadence_epoch_splits=splits,
                  raw_integrated_leave_one_sector_out=leave_one_out,
                  sector_sign_check=dict(n_positive=positive, n_sectors=len(sector_parity),
                                         binomial_two_sided_pvalue=sign_pvalue,
                                         note='Exploratory fixed-hypothesis sign test assumes independent sector signs; search selection and shared calibration can violate this assumption.'),
                  conditional_timing_parity_days=timing_parity,
                  raw_sector_orbital_harmonics=harmonics,
                  source_report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
                  frozen_v3_sha256=bv_sha,
                  helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  interpretation='Absence of alternating depth does not exclude equal-depth grazing eclipsing binaries. Folding at twice the recurrence period creates two eclipses for a true-period planet as well. Harmonics are diagnostic nuisance tests, not binary classifications.',
                  limitations=['Local-width average depths depend on exposure sampling and timing; they are not fitted physical radii.',
                               'The ephemeris and width were estimated from these data; conditional errors are not independent validation.',
                               'Within-sector correlations are only partially handled by sector-cluster covariance; catalog/source ambiguity remains.'])
    bv.atomic_json(HERE / 'parity_phase_qa.json', result)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    for color, even in [('tab:blue', True), ('tab:orange', False)]:
        events = all_events['SAP_FLUX']
        use = [v for v in events if (v['cycle'] % 2 == 0) == even]
        axes[0].errorbar([v['time'] for v in use], [v['depth'] * 1e6 for v in use],
                         [v['depth_err'] * 1e6 for v in use], fmt='.', color=color,
                         alpha=.65, label='Even' if even else 'Odd')
    axes[0].legend(); axes[0].set_ylabel('Raw local depth (ppm)')
    axes[0].set_xlabel('BTJD (TDB)')
    for factor, marker in [(1, 'o'), (2, 's')]:
        vals = [r for r in harmonics if r['period_factor'] == factor and r['status'] == 'measured']
        axes[1].errorbar([r['sector'] for r in vals], [r['orbital_semiamplitude'] * 1e6 for r in vals],
                         [r['orbital_semiamplitude_err'] * 1e6 for r in vals],
                         fmt=marker, label=f'P={factor} × recurrence')
    axes[1].set_ylabel('Raw orbital semi-amplitude (ppm)')
    axes[1].set_xlabel('TESS sector'); axes[1].legend()
    fig.suptitle('TIC 238202207: conditional alternating-depth and raw phase checks')
    fig.tight_layout(); fig.savefig(HERE / 'parity_phase_qa.png', dpi=150)
    for col, stats in parity.items():
        print(col, stats['n_events'], stats['even_minus_odd_depth'], stats['difference_error'], stats['difference_conditional_snr'])
    for col, stats in integrated_parity.items():
        print('Integrated', col, stats['n_events'], stats['even_minus_odd_depth'], stats['difference_error'], stats['difference_conditional_snr'])


if __name__ == '__main__':
    main()
