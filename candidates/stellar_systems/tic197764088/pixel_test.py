"""Locate both eclipses of the TIC197764088 aperture in sector28 native FFIs."""
from pathlib import Path
import hashlib
import json
import os
import sys
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import verify_fits_checksums
import pixel_vet_269994642 as primitive
import pixel_vet_49617264 as fitting
from pixel_vet_260820871 import image_pointing

OUT = Path(__file__).resolve().parent
HOST = '6410698972335936512'
NEIGHBOR = '6410698633034526976'


def main():
    observation = json.loads((OUT / 'pixel_observations.json').read_text())[0]
    raw = Path(observation['path'])
    if hashlib.sha256(raw.read_bytes()).hexdigest() != observation['sha256']:
        raise ValueError('Pixel raw checksum changed')
    phot = json.loads((OUT / 'sector28_measurements.json').read_text())
    p, t0 = phot['refined']['period_days'], phot['refined']['epoch_btjd']
    shape = phot['columns']['SAP_FLUX']['primary_shape']
    primitive.DURATION, primitive.INGRESS = shape['duration_days'], shape['ingress_fraction']
    with fits.open(raw, memmap=False) as hd:
        checks = verify_fits_checksums(raw, hd)
        h, ph, w = dict(hd[0].header), dict(hd[1].header), WCS(hd[2].header)
        data = hd[1].data
        t, f, e, q = [np.asarray(data[key]) for key in ['TIME', 'FLUX', 'FLUX_ERR', 'QUALITY']]
        ffi = np.asarray(data['FFI_FILE']).astype(str)
        good = (q == 0) & np.isfinite(t) & np.all(np.isfinite(f) & np.isfinite(e), axis=(1, 2))
        t, f, e = t[good], f[good].astype(float), e[good].astype(float)
    if h['TIMESYS'] != 'TDB' or h['TIMEREF'] != 'SOLARSYSTEM' or h['BJDREFI'] != 2457000:
        raise ValueError('Unexpected pixel timing')
    stars = pd.read_csv(ROOT / 'data/stellar_systems/scout_novelty/tic197764088_gaia_context.csv',
                        dtype={'source_id': str})
    stars[['x', 'y']] = w.all_world2pix(stars[['ra', 'dec']].to_numpy(), 0)
    host = stars.loc[stars.source_id == HOST, ['x', 'y']].to_numpy()[0]
    neighbor = stars.loc[stars.source_id == NEIGHBOR, ['x', 'y']].to_numpy()[0]
    # TESSCut rectangles can include collateral columns. Require the target to
    # lie on illuminated silicon, and preserve detector indices in provenance.
    rawx = ph['1CRV5P'] + (host[0] + 1 - ph['1CRP5P']) * ph['1CDL5P']
    rawy = ph['2CRV5P'] + (host[1] + 1 - ph['2CRP5P']) * ph['2CDL5P']
    if not 44 <= rawx < 2092 or not 0 <= rawy < 2048:
        raise ValueError('Target outside science detector')
    reference = np.median(f, axis=0)
    blank = reference < np.percentile(reference, 25)
    background = np.median(f[:, blank], axis=1)
    f -= background[:, None, None]
    dx, dy = image_pointing(f, reference - np.median(reference[blank]), host, neighbor)
    stars.to_csv(OUT / 'gaia_pixel_positions.csv', index=False)
    (OUT / 'ffi_product_identifiers.txt').write_text('\n'.join(ffi) + '\n')
    yy, xx = np.indices(reference.shape)
    phases = {}
    for label, offset in [('primary', 0), ('secondary', .5)]:
        ep = t0 + p * offset
        modes = {}
        for mode in ['quadratic', 'pointing', 'box_core', 'box_broad']:
            events = primitive.event_fit(t, f, e, background, dx, dy, p, ep, ph['TIMEDEL'], mode)
            if len(events) < 3:
                raise ValueError('Insufficient observed pixel eclipse events')
            difference, error = fitting.combine(events)
            sy, sx = np.unravel_index(np.argmax(difference / error), difference.shape)
            gf, _ = fitting.gaussfit(difference, error, [sx, sy])
            gf['host_distance_pixels'] = float(np.hypot(gf['x'] - host[0], gf['y'] - host[1]))
            gf['known_EB_neighbor_distance_pixels'] = float(np.hypot(gf['x'] - neighbor[0], gf['y'] - neighbor[1]))
            fixed = {}
            for name, position in [('host', host), ('known_EB_neighbor', neighbor)]:
                result, _ = fitting.gaussfit(difference, error, position, fixed_center=True)
                result['delta_chi2_vs_free'] = result['chi2_diagonal'] - gf['chi2_diagonal']
                fixed[name] = result
            modes[mode] = {'n_events': len(events), 'source_fit': gf, 'fixed_positions': fixed}
            np.savez_compressed(OUT / f'{label}_{mode}_difference.npz', difference=difference, error=error,
                                event_differences=np.array([event['difference'] for event in events]),
                                event_errors=np.array([event['error'] for event in events]))
            if mode == 'pointing':
                selected_events = events
                boot = fitting.boot_centroids(events, [gf['x'], gf['y']], 200)
                np.savez_compressed(OUT / f'{label}_bootstrap.npz', centroids_xy=boot)
                bootstrap = {'n': len(boot), 'median_xy': np.median(boot, axis=0),
                             'host_distance_percentiles_16_50_84': np.percentile(np.linalg.norm(boot - host, axis=1), [16, 50, 84]),
                             'neighbor_distance_percentiles_16_50_84': np.percentile(np.linalg.norm(boot - neighbor, axis=1), [16, 50, 84]),
                             'fraction_closer_neighbor_than_host': np.mean(np.linalg.norm(boot - neighbor, axis=1) < np.linalg.norm(boot - host, axis=1))}
                fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
                for ax, image, title in zip(axes, [np.log10(np.maximum(reference, 1)), difference, difference / error],
                                            ['Median image log10 e-/s', 'Eclipse difference e-/s', 'Difference / diagonal error']):
                    z = ax.imshow(image, origin='lower', cmap='viridis')
                    ax.plot(*host, 'gx', ms=11, label='TIC aperture target')
                    ax.plot(*neighbor, 'o', mfc='none', mec='orange', ms=11, label='Known EB neighbor')
                    ax.plot(gf['x'], gf['y'], 'm+', ms=14, label='Eclipse source fit')
                    ax.set_title(title)
                    fig.colorbar(z, ax=ax, shrink=.7)
                axes[0].legend(fontsize=7)
                fig.suptitle(f'TIC197764088 S28 {label}: source localization')
                fig.tight_layout()
                fig.savefig(OUT / f'{label}_localization.png', dpi=160)
                plt.close(fig)
        mask = (xx - host[0]) ** 2 + (yy - host[1]) ** 2 < 2 ** 2
        photometry, _ = fitting.photometry(selected_events, {'host': host, 'known_EB_neighbor': neighbor},
                                           reference, mask)
        for measurement in photometry:
            if measurement['aperture'] == 'spoc_optimal':
                measurement['aperture'] = 'host_r2_reference'
        phases[label] = {'epoch_btjd': ep, 'modes': modes, 'bootstrap': bootstrap, 'aperture_photometry': photometry}
    provenance = {**observation, 'fits_checksums': checks, 'wcs': dict(w.to_header()),
                  'processing_version': h['PROCVER'], 'data_release': h['DATA_REL'],
                  'camera': h['CAMERA'], 'ccd': h['CCD'], 'time_system': h['TIMESYS'],
                  'time_reference': h['TIMEREF'], 'bjdref': h['BJDREFI'] + h['BJDREFF'],
                  'exposure_days': ph['TIMEDEL'], 'physical_host_column': rawx, 'physical_host_row': rawy}
    result = {'host_xy': host, 'known_EB_neighbor_xy': neighbor, 'period_days': p,
              'provenance': provenance, 'phases': phases,
              'limitations': 'Gaussian localization and approximate image-pointing basis; no calibrated TESS PRF/WCS/proper-motion systematics or posterior. Both phases and multiple templates tested independently.'}
    (OUT / 'pixel_summary.json').write_text(json.dumps(fitting.safe(result), indent=2) + '\n')
    print(json.dumps(fitting.safe({label: phase['bootstrap'] for label, phase in phases.items()}), indent=2))


if __name__ == '__main__':
    main()
