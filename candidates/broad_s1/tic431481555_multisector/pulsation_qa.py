"""Test broad 5.635-day dips against out-of-event fitted stellar oscillations.

Reproduce: python candidates/broad_s1/tic431481555_multisector/pulsation_qa.py
No paid calls or downloads; uses hash-verified cached observations. Frequencies
are selected from points outside the primary windows, then predict those windows.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys

import numpy as np
from astropy.timeseries import LombScargle

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from broad_vet import read_source, load_records, preprocess_curve, FastScreenConfig, local_quadratic, phase
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent


def fit_modes(t, f, e, outside, maximum_modes=8):
    baseline = float(np.ptp(t))
    frequencies = np.arange(.2, min(20., .45 / np.median(np.diff(t))), 1 / (8 * baseline))
    chosen = []
    residual = f - np.average(f[outside], weights=1 / e[outside] ** 2)
    panels = []
    for iteration in range(maximum_modes):
        ls = LombScargle(t[outside] - t.min(), residual[outside], e[outside])
        power = ls.power(frequencies, method='fast', assume_regular_frequency=True)
        allowed = np.ones(len(frequencies), bool)
        for frequency in chosen:
            allowed &= abs(frequencies - frequency) > .7 / baseline
        j = int(np.argmax(np.where(allowed, power, -np.inf)))
        chosen.append(float(frequencies[j]))
        design = [np.ones(len(t))]
        for frequency in chosen:
            angle = 2 * np.pi * frequency * (t - t.min())
            design.extend([np.sin(angle), np.cos(angle)])
        design = np.column_stack(design)
        beta = np.linalg.lstsq(design[outside] / e[outside, None], f[outside] / e[outside], rcond=None)[0]
        model = design @ beta
        residual = f - model
        panels.append(dict(iteration=iteration + 1, frequency_per_day=float(frequencies[j]),
                           period_days=float(1 / frequencies[j]), normalized_ls_power=float(power[j]),
                           amplitude=float(np.hypot(beta[-2], beta[-1]))))
    return model, panels, frequencies


if __name__ == '__main__':
    tic = 431481555
    base = load_records([tic])[0]
    sources = [dict(sector=1, path=base['path'], sha256=base['sha256'], url=base['url'])]
    sources += [r for r in json.loads((ROOT / 'data/broad_followup_queries/tic431481555_qlp.json').read_text()) if r['status'] == 'retained']
    dossier = json.loads((OUT / 'multisector.json').read_text())
    independent = {r['sector']: r['ephemeris'] for r in dossier['independent_sector_searches']}
    reports, plotdata = [], []
    for source in sources:
        data, checks = read_source(source)
        good = data['QUALITY'] == 0
        flat = preprocess_curve(data['TIME'][good], data['SAP_FLUX'][good], data['KSPSAP_FLUX_ERR'][good],
                                config=FastScreenConfig(width_days=2.))
        t, f, e = flat['t'], flat['f'], flat['e']
        eph = independent[source['sector']]
        p, epoch = eph['period'], eph['t0_btjd']
        # Exclude the full broad primary feature, not merely its box-search core.
        half_mask = max(.27, eph['duration_days'] * 1.3)
        outside = abs(phase(t, p, epoch)) > half_mask
        model, modes, grid = fit_modes(t, f, e, outside)
        corrected = f - model + np.median(model)
        before, _ = local_quadratic(t, f, e, p, epoch, eph['duration_days'])
        after, _ = local_quadratic(t, corrected, e, p, epoch, eph['duration_days'])
        inside = abs(phase(t, p, epoch)) < .5 * eph['duration_days']
        prediction = np.corrcoef(f[inside], model[inside])[0, 1] if inside.sum() > 5 else None
        reports.append(dict(sector=source['sector'], frequencies_selected_outside_primary=modes,
                            excluded_primary_half_width_days=half_mask, n_fit=int(outside.sum()),
                            before_local=before, after_local=after,
                            held_primary_flux_model_correlation=float(prediction) if prediction is not None else None,
                            first_mode_period_days=modes[0]['period_days'], frequency_grid_count=len(grid),
                            source={k: source[k] for k in ('url', 'path', 'sha256')}, checksums=checks))
        plotdata.append((source['sector'], t, f, model, corrected, p, epoch))
        print('PULSATION_QA', source['sector'], 'first_mode_period', modes[0]['period_days'],
              'candidate local SNR before/after', before['snr'], after['snr'], flush=True)
    payload = dict(tic=tic, checked_utc=datetime.now(timezone.utc).isoformat(), sector_reports=reports,
                   previously_recognized_variable_type='Gaia/VSX DSCT|GDOR|SXPHE',
                   hypothesis='Broad 5.635-day features may be pulsation minima/beat-envelope aliases.',
                   limitations=['Eight sinusoids are an exploratory model, not proof that every dip is a pulsation.',
                                'Frequencies are selected from points outside generous candidate event windows.',
                                'Periodograms and optimized scores are not calibrated false-alarm probabilities.',
                                'A pulsating star can also host planets; this does not exclude other planetary signals.'],
                   code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'pulsation_qa.json').write_text(json.dumps(payload, indent=2, allow_nan=False))
    sector, t, f, model, corrected, p, epoch = plotdata[0]
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), layout='constrained')
    axes[0].plot(t, f, '.', ms=2, color='.6')
    axes[0].plot(t, model, lw=.8, color='tab:orange')
    axes[0].set(xlabel='BTJD', ylabel='Raw SAP normalized', title='Sector 1: oscillation model fit outside candidate events')
    ph = phase(t, p, epoch)
    axes[1].plot(ph, f, '.', ms=2, color='.5')
    axes[1].plot(ph, model, '.', ms=2, color='tab:orange')
    axes[1].set(xlabel='Phase time at independently selected long period (days)', ylabel='Before')
    axes[2].plot(ph, corrected, '.', ms=2, color='.5')
    axes[2].set(xlabel='Phase time (days)', ylabel='After exploratory oscillation subtraction')
    fig.savefig(OUT / 'pulsation_qa.png', dpi=130)
    plt.close(fig)
