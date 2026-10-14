"""Independent QLP sector-1 alternating-depth check, frozen v3 read-only."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
import broad_vet as bv


def parity_fit(events):
    d = np.array([v['depth'] for v in events])
    e = np.array([v['depth_err'] for v in events])
    even = np.array([v['cycle'] % 2 == 0 for v in events])
    a = np.column_stack([np.ones(len(d)), np.where(even, .5, -.5)]) / e[:, None]
    beta = np.linalg.lstsq(a, d / e, rcond=None)[0]
    residual = d / e - a @ beta
    chi2 = float(np.sum(residual ** 2))
    covariance = np.linalg.pinv(a.T @ a) * max(1., chi2 / max(1, len(d) - 2))
    scatter_error = np.sqrt(sum(np.var(d[even == k], ddof=1) / (even == k).sum() for k in [0, 1]))
    error = max(float(np.sqrt(covariance[1, 1])), float(scatter_error))
    return dict(n_events=len(d), n_even=int(even.sum()), n_odd=int((~even).sum()),
                mean_depth=float(beta[0]), even_minus_odd_depth=float(beta[1]),
                difference_error=error, conditional_difference_snr=float(abs(beta[1]) / error),
                chi2=chi2, dof=len(d) - 2, coefficient_order=['mean_depth', 'even_minus_odd_depth'],
                coefficients=beta.tolist(), covariance=covariance.tolist(), events=events)


def main():
    dossier_path = HERE / 'vetting.json'
    d = json.loads(dossier_path.read_text())
    assert d['vetting_method']['source_sha256'] == bv.method_identity()['source_sha256']
    raw, checks = bv.read_source(d['source'])
    good = raw['QUALITY'] == 0
    t, e = raw['TIME'][good], raw['KSPSAP_FLUX_ERR'][good]
    p = d['refined_ephemeris']['period']
    shape = d['integrated_trapezoid']
    epoch, width = shape['epoch_btjd'], d['refined_ephemeris']['duration_days']
    exposure = float(np.median(np.diff(t[np.isfinite(t)])))
    output = {}
    for col in ['KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']:
        stats, local = bv.local_quadratic(t, raw[col][good], e, p, epoch, width)
        events = []
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
            events.append(dict(cycle=int(cycle), time_btjd=float(epoch + cycle * p),
                               depth=float(beta[1]), depth_err=float(np.sqrt(covariance[1, 1]))))
        output[col] = dict(window_average=parity_fit(stats['events']), integrated_fixed_shape=parity_fit(events))
    result = dict(tic=269994642, sector=1, checked_utc=datetime.now(timezone.utc).isoformat(),
                  period=p, epoch_btjd=epoch, exposure_days=exposure, tests=output,
                  frozen_v3_sha256=bv.method_identity()['source_sha256'],
                  source_fits=d['source'], fits_checksums=checks,
                  input_dossier_sha256=hashlib.sha256(dossier_path.read_bytes()).hexdigest(),
                  helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  interpretation='Conditional alternating depth is an unequal-binary diagnostic, not a survey-calibrated exclusion of planets.',
                  limitations=['One sector only; no independent-sector covariance is available.',
                               'Depth error takes the larger of covariance with chi-square inflation and even/odd empirical event-scatter error.',
                               'Period, eclipse shape and epoch are estimated from the same data. Stellar variability, timing/shape changes and source blending remain alternatives.',
                               'A double-period half-phase eclipse by itself occurs for any true-period recurring planet.'])
    bv.atomic_json(HERE / 'parity_qa.json', result)
    for col, tests in output.items():
        for method, r in tests.items():
            print(col, method, r['even_minus_odd_depth'] * 1e6, r['difference_error'] * 1e6, r['conditional_difference_snr'])


if __name__ == '__main__':
    main()
