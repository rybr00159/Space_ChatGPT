"""Supplemental native-time QLP parity and raw phase tests for TIC260820871."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
shared = HERE.parent / 'tic238202207_multisector/parity_phase_qa.py'
spec = importlib.util.spec_from_file_location('shared_parity', shared)
qa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qa)
bv = qa.bv


def main():
    report_path = HERE / 'multisector.json'
    report = json.loads(report_path.read_text())
    assert report['vetting_method']['source_sha256'] == bv.method_identity()['source_sha256']
    timing = report['event_timing_fit']['fitted_ephemeris']
    shape = report['approximate_pooled_shape']
    period, epoch = timing['period'], timing['epoch_btjd']
    width = report['joint_ephemeris']['duration_days']
    columns = ['KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']
    events = {col: [] for col in columns}
    harmonics = []
    for row in report['sources']:
        raw, checks = bv.read_source(row)
        good = raw['QUALITY'] == 0
        t, e = raw['TIME'][good], raw['KSPSAP_FLUX_ERR'][good]
        exposure = float(np.median(np.diff(t[np.isfinite(t)])))
        for col in columns:
            stats, local = bv.local_quadratic(t, raw[col][good], e, period, epoch, width)
            for cycle in np.unique(local['cycle']):
                use = local['cycle'] == cycle
                x = local['t'][use] - epoch - cycle * period
                template = bv.integrated_trapezoid(x, shape['duration_days'],
                                                   shape['duration_days'] * shape['ingress_fraction'], exposure)
                design = np.column_stack([np.ones(len(x)), -template])
                a, y = design / local['e'][use, None], local['f'][use] / local['e'][use]
                beta = np.linalg.lstsq(a, y, rcond=None)[0]
                chi2 = float(np.sum((y - a @ beta) ** 2))
                covariance = np.linalg.pinv(a.T @ a) * max(1., chi2 / max(1, len(x) - 2))
                events[col].append(dict(sector=row['sector'], cycle=int(cycle), time=float(epoch + cycle * period),
                                        depth=float(beta[1]), depth_err=float(np.sqrt(covariance[1, 1]))))
        for factor in [1, 2]:
            diagnostics = bv.variability_diagnostics(t, raw['SAP_FLUX'][good], e, period * factor, epoch, width)
            harmonics.append(dict(sector=row['sector'], period_factor=factor, **diagnostics))
    tests = {col: qa.alternating_depth(values) for col, values in events.items()}
    per_sector = {str(s): qa.alternating_depth([v for v in events['SAP_FLUX'] if v['sector'] == s])
                  for s in sorted({v['sector'] for v in events['SAP_FLUX']})}
    result = dict(tic=260820871, checked_utc=datetime.now(timezone.utc).isoformat(),
                  ephemeris=dict(period_days=period, epoch_btjd=epoch),
                  integrated_template_parity_tests=tests, raw_integrated_per_sector_parity=per_sector,
                  integrated_template_event_depths=events, raw_sector_orbital_harmonics=harmonics,
                  input_report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
                  frozen_v3_sha256=bv.method_identity()['source_sha256'],
                  shared_helper_sha256=hashlib.sha256(shared.read_bytes()).hexdigest(),
                  helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  limitations=['Conditional fixed period and trapezoid shape; only three sectors and two observing eras.',
                               'Sector-cluster covariance is poorly estimated with three sectors; diagonal errors are inflated and the larger error retained.',
                               'Many long-gap period aliases remain possible; no survey-wide significance or binary confirmation is inferred.'])
    bv.atomic_json(HERE / 'parity_phase_qa.json', result)
    for col, r in tests.items():
        print(col, r['n_events'], r['even_minus_odd_depth'] * 1e6, r['difference_error'] * 1e6, r['difference_conditional_snr'])
    for r in harmonics:
        print('HARMONIC', r['sector'], r['period_factor'], r.get('orbital_semiamplitude'), r.get('half_period_semiamplitude'), r.get('orbital_snr'), r.get('half_period_snr'))


if __name__ == '__main__':
    main()
