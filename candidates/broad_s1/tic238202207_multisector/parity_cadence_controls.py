"""Native-time constant-planet controls for the supplemental parity statistic.

Uses limb-darkened BATMAN curves, not the fitted trapezoid, at all original
quality-zero timestamps and native exposures. This checks whether a fixed
trapezoid plus alternating exposure phases creates the observed depth offset.
It is not a correlated-noise or survey-wide false-alarm calibration.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np
import batman

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_phase_qa as qa
bv = qa.bv


def main():
    report = json.loads((HERE / 'multisector.json').read_text())
    timing = report['event_timing_fit']['fitted_ephemeris']
    shape = report['approximate_pooled_shape']
    period, epoch = timing['period'], timing['epoch_btjd']
    width = report['joint_ephemeris']['duration_days']
    rng = np.random.default_rng(238202207)
    profiles = [(0.15, 0.98), (0.20, 1.05), (0.27, 1.12)]
    controls = []
    for radius, impact in profiles:
        # Solve a/Rstar from the observed total duration under circular geometry.
        angle = np.pi * shape['duration_days'] / period
        scaled_axis = np.sqrt(impact ** 2 + ((1 + radius) ** 2 - impact ** 2) / np.sin(angle) ** 2)
        parameters = batman.TransitParams()
        parameters.t0, parameters.per, parameters.rp = epoch, period, radius
        parameters.a, parameters.inc = scaled_axis, np.degrees(np.arccos(impact / scaled_axis))
        parameters.ecc, parameters.w = 0., 90.
        parameters.u, parameters.limb_dark = [.4, .2], 'quadratic'
        for noise in [False, True]:
            events = []
            for row in report['sources']:
                raw, checks = bv.read_source(row)
                good = (raw['QUALITY'] == 0) & np.isfinite(raw['TIME']) & np.isfinite(raw['KSPSAP_FLUX_ERR']) & (raw['KSPSAP_FLUX_ERR'] > 0)
                t, error = raw['TIME'][good], raw['KSPSAP_FLUX_ERR'][good]
                exposure = float(np.median(np.diff(t)))
                flux = batman.TransitModel(parameters, t, supersample_factor=15, exp_time=exposure).light_curve(parameters)
                if noise:
                    flux += rng.normal(size=len(t)) * error
                stats, local = bv.local_quadratic(t, flux, error, period, epoch, width)
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
                    events.append(dict(sector=row['sector'], cycle=int(cycle), depth=float(beta[1]),
                                       depth_err=float(np.sqrt(covariance[1, 1]))))
            parity = qa.alternating_depth(events)
            controls.append(dict(radius_ratio=radius, impact_parameter=impact,
                                 scaled_semimajor_axis=scaled_axis, gaussian_error_noise=noise,
                                 limb_darkening=parameters.u, parity=parity))
            print(radius, impact, noise, parity['even_minus_odd_depth'] * 1e6,
                  parity['difference_error'] * 1e6, parity['difference_conditional_snr'])
    result = dict(tic=238202207, checked_utc=datetime.now(timezone.utc).isoformat(), controls=controls,
                  model='BATMAN circular quadratic limb-darkened equal recurring planet;15-subsample exposure integration',
                  frozen_v3_sha256=bv.method_identity()['source_sha256'],
                  input_report_sha256=hashlib.sha256((HERE / 'multisector.json').read_bytes()).hexdigest(),
                  helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  conclusion='These native-time synthetic controls test cadence/template bias only. They do not establish independent-source attribution, correlated-noise significance or a binary confirmation.',
                  limitations=['Only three selected profiles and one Gaussian realization per profile.',
                               'Observed quality masks, formal errors and exposures are retained; actual stellar variability/systematics are not injected.',
                               'No significance correction for selecting this star or parity diagnostic in the broad survey.'])
    bv.atomic_json(HERE / 'parity_cadence_controls.json', result)


if __name__ == '__main__':
    main()
