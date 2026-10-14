"""Brief independent-era QA for four stale-radius broad-survey leads.

Reproduce: python candidates/broad_s1/giant_survivors/recurrence_qa.py
Reads cached, hash-verified public QLP observations only. Per-sector epochs are
free, so absence at an extrapolated imprecise sector-1 epoch cannot reject a
planet. Current catalog giant classification is separate physical evidence.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from broad_vet import (load_records, read_source, refine_ephemeris, local_quadratic,
                       variability_diagnostics, preprocess_curve, FastScreenConfig, phase)
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
TARGETS = {260352037: (.88013, 10.1065), 441386747: (7.4801, 12.1823),
           406280158: (10.1259, 10.5196), 271697430: (3.05937, 10.2269)}


def analyze(tic, selected_period, current_radius):
    base = load_records([tic])[0]
    peak = min((v for scan in base['scans'] for v in scan['peaks']), key=lambda v: abs(v['period'] - selected_period))
    sources = [dict(sector=1, path=base['path'], sha256=base['sha256'], url=base['url'], metadata=base['metadata'])]
    manifest = ROOT / f'data/broad_followup_queries/tic{tic}_qlp.json'
    sources += [r for r in json.loads(manifest.read_text()) if r['status'] == 'retained' and r['sector'] != 1]
    rows, curves = [], []
    for source in sources:
        d, checks = read_source(source)
        good = d['QUALITY'] == 0
        t, f, e = d['TIME'][good], d['KSPSAP_FLUX'][good], d['KSPSAP_FLUX_ERR'][good]
        independent = refine_ephemeris(t, f, e, peak)
        p, epoch, duration = independent['period'], independent['t0_btjd'], independent['duration_days']
        primary, data = local_quadratic(t, f, e, p, epoch, duration)
        raw_primary, _ = local_quadratic(t, d['SAP_FLUX'][good], e, p, epoch, duration)
        shape = None
        # Fixed epoch/cycle counts from sector 1 are never imposed here.
        variability = variability_diagnostics(t, d['SAP_FLUX'][good], e, p, epoch, duration)
        row = dict(sector=source['sector'], independent_local_search=independent,
                   raw_local_primary=raw_primary, processed_local_primary=primary,
                   raw_phase_variability=variability, source={k: source[k] for k in ('path', 'sha256', 'url')},
                   checksum_status=checks, n_raw=len(d['TIME']), n_quality_zero=int(good.sum()))
        rows.append(row)
        flat = preprocess_curve(t, f, e, config=FastScreenConfig(width_days=2.))
        curves.append((source['sector'], flat, independent))
    p = selected_period
    # Keplerian semimajor axis for a deliberately generous 1--3 solar-mass giant.
    # Values below the photospheric radius exclude an external circular orbit;
    # giant oscillation or blended-star hypotheses remain separate possibilities.
    grav, solar_mass, solar_radius = 6.67430e-11, 1.98840987e30, 6.957e8
    axes = [(grav * mass * solar_mass * (p * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / solar_radius for mass in (1, 2, 3)]
    payload = dict(tic=tic, checked_utc=datetime.now(timezone.utc).isoformat(), selected_period_days=selected_period,
                   selected_original_peak=peak, current_tic_v82_radius_rsun=current_radius,
                   old_qlp_radius_rsun=base['metadata'].get('RADIUS'),
                   semimajor_axes_rsun_for_mass_1_2_3_solar=axes,
                   outside_photosphere_for_3_solar_mass=axes[-1] > current_radius,
                   independent_sector_qa=rows,
                   interpretation='Giant-star variability/oscillation and contamination require precedence over a small-planet interpretation based on the obsolete QLP radius.',
                   limitations=['Independent BLS epochs were optimized; their scores have no global false-alarm calibration.',
                                'This is brief giant/control QA, not a physical planet model or exhaustive catalog verification.',
                                'No rejection based solely on an uncertain extrapolated sector-1 ephemeris.',
                                'Current radii are from the separate catalog-refresh evidence, not measured from these light curves.'],
                   code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   vetting_source_sha256=hashlib.sha256((ROOT / 'src/broad_vet.py').read_bytes()).hexdigest())
    (OUT / f'tic{tic}_qa.json').write_text(json.dumps(payload, indent=2, allow_nan=False))
    representatives = [curves[0]] + curves[1:2] + curves[-1:] if len(curves) > 2 else curves
    fig, axes_plot = plt.subplots(len(representatives), 2, figsize=(11, 3 * len(representatives)), layout='constrained')
    axes_plot = np.atleast_2d(axes_plot)
    for axes_row, (sector, flat, eph) in zip(axes_plot, representatives):
        axes_row[0].plot(flat['t'], flat['f'], '.', ms=2, color='.5')
        axes_row[0].set(xlabel='BTJD', ylabel=f'Sector {sector}: normalized processed flux')
        axes_row[1].plot(phase(flat['t'], eph['period'], eph['t0_btjd']), flat['f'], '.', ms=2, color='.5')
        axes_row[1].set(xlabel='Phase time (days), independently selected epoch', ylabel='Normalized flux',
                        xlim=(-max(.3, eph['duration_days'] * 2), max(.3, eph['duration_days'] * 2)))
    fig.suptitle(f'TIC {tic}: current R={current_radius:.2f} solar; old dwarf prior invalid')
    fig.savefig(OUT / f'tic{tic}_qa.png', dpi=130)
    plt.close(fig)
    print(tic, 'sectors', len(rows), 'independent raw local SNR range',
          min(r['raw_local_primary']['snr'] for r in rows), max(r['raw_local_primary']['snr'] for r in rows),
          'a(3M)/R', axes[-1] / current_radius, flush=True)
    return payload


if __name__ == '__main__':
    summary = [analyze(tic, *values) for tic, values in TARGETS.items()]
    lines = ['# Stale-radius giant-star survivor checks', '',
             'These four initial flag-free weak leads used obsolete QLP stellar-radius metadata. '
             'Current TIC v8.2 gives radii of 10–12 solar radii. The previously inferred small planet sizes '
             'are unsupported. Cached additional observations are preserved and were checked with '
             'independently optimized local periods and epochs.', '',
             '|TIC|Trial P (days)|Old QLP R|Current TIC R|Sectors checked|a at 3 solar masses / current R|',
             '|---|---:|---:|---:|---:|---:|']
    for d in summary:
        lines.append(f"|{d['tic']}|{d['selected_period_days']:.5f}|{d['old_qlp_radius_rsun']}|"
                     f"{d['current_tic_v82_radius_rsun']:.4f}|{len(d['independent_sector_qa'])}|"
                     f"{d['semimajor_axes_rsun_for_mass_1_2_3_solar'][-1] / d['current_tic_v82_radius_rsun']:.2f}|")
    lines += ['', 'Keplerian axes are an illustrative 1–3 solar-mass range, not fitted masses. '
              'An orbit inside the photosphere cannot produce an external planetary transit on the giant. '
              'A signal originating on a blended dwarf remains possible and requires pixel source attribution. '
              'Where independently selected broad-search peaks recur, the optimization and giant correlated '
              'noise prevent a planetary interpretation from the score alone.', '',
              'Per-sector depths, event sampling, raw phase coefficients, archive source URLs, and SHA-256 '
              'identifiers are in the QA JSON files. All original FITS remain cached outside Git. '
              'Reproduce with `python candidates/broad_s1/giant_survivors/recurrence_qa.py`.']
    (OUT / 'README.md').write_text('\n'.join(lines) + '\n')
