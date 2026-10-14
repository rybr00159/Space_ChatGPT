"""Native sector-28 spatial tests of the two TIC71316387 eclipse depths.

Empirical integrated-Gaussian centroids and event/block bootstraps diagnose
resolved contamination; they cannot resolve an unresolved companion or replace
calibrated TESS PRF/scene fitting. Raw cutouts and their expected hashes remain
outside Git. Re-run with the repository Python environment.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import io
import json
import os
import sys
import urllib.parse
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import verify_fits_checksums
from pixel_vet_49617264 import safe, gaussfit, combine

TIC = 71316387
RA, DEC = 330.887714849402, -42.8354976192099
PERIOD, EPOCH = .6518369465110034, 2067.327192476
OUT = ROOT / 'candidates/tic71316387/pixels'
CACHE = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'tic71316387/pixels'


def retrieve():
    url = 'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?' + urllib.parse.urlencode(
        {'ra': RA, 'dec': DEC, 'x': 9, 'y': 9, 'units': 'px', 'sector': 28})
    paths = list(CACHE.glob('tess-s0028-*.fits'))
    if not paths:
        CACHE.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=90) as response:
            blob = response.read(60000001)
        if len(blob) > 60000000:
            raise RuntimeError('Cutout exceeded 60 MB retrieval bound.')
        with zipfile.ZipFile(io.BytesIO(blob)) as zipped:
            for name in zipped.namelist():
                if name.endswith('.fits'):
                    path = CACHE / Path(name).name
                    if path.exists():
                        raise RuntimeError('Refusing to overwrite raw evidence.')
                    path.write_bytes(zipped.read(name))
        paths = list(CACHE.glob('tess-s0028-*.fits'))
    if len(paths) != 1:
        raise RuntimeError('Expected exactly one sector-28 cube.')
    return paths[0], url


def event_fits(t, flux, error, bg, shifts, dt, mode):
    events = []
    for cycle in range(int(np.ceil((t.min()-EPOCH)/PERIOD)),
                       int(np.floor((t.max()-EPOCH)/PERIOD))+1):
        center = EPOCH + cycle * PERIOD
        keep = abs(t-center) < .22
        tx = t[keep] - center
        duration, ingress = (.0617, .262) if cycle % 2 == 0 else (.0608, .327)
        offsets = np.linspace(-dt/2, dt/2, 31)
        shape = np.mean(np.clip((duration/2-abs(tx[:, None]+offsets)) /
                               (duration*ingress), 0, 1), axis=1)
        if (np.sum(shape > .05) < 3 or
                np.sum((shape == 0) & (tx < 0)) < 8 or
                np.sum((shape == 0) & (tx > 0)) < 8):
            continue
        cols = [np.ones(len(tx)), tx, tx**2]
        if mode in ('background', 'pointing'):
            nuisances = [bg[keep]]
            if mode == 'pointing':
                nuisances.extend([s[keep] for s in shifts])
            for cov in nuisances:
                cov = cov - np.mean(cov)
                if np.std(cov) > 1e-10:
                    cols.append(cov/np.std(cov))
        design = np.column_stack(cols + [-shape])
        inv = np.linalg.pinv(design)
        ff = flux[keep].reshape(len(tx), -1)
        coef = inv @ ff
        residual = ff - design @ coef
        sig = np.maximum(np.std(residual, axis=0, ddof=design.shape[1]),
                         np.sqrt(np.mean(error[keep].reshape(len(tx), -1)**2,
                                         axis=0)))
        differr = sig * np.sqrt(np.sum(inv[-1]**2))
        events.append({'cycle': cycle, 'parity': cycle % 2,
                       'center_btjd': center, 'n_in': int(np.sum(shape > .05)),
                       'n_out': int(np.sum(shape == 0)),
                       'difference': coef[-1].reshape(flux.shape[1:]),
                       'error': differr.reshape(flux.shape[1:]),
                       'residual': residual, 'projection': inv[-1],
                       'design_condition': np.linalg.cond(design)})
    return events


def boot_centroids(events, seed, roi, n, rng):
    result = []
    # Sample whole events, then synchronous five-cadence (~50-minute) residual
    # blocks across every pixel. This preserves pixel covariance and some red
    # noise but does not calibrate PRF or WCS systematics.
    for _ in range(n):
        draws = []
        for j in rng.integers(0, len(events), len(events)):
            event = events[j]
            residual = event['residual']
            count = len(residual)
            starts = rng.integers(0, count, (count+4)//5)
            ix = np.concatenate([(s+np.arange(5)) % count for s in starts])[:count]
            perturbation = (event['projection'] @ residual[ix]).reshape(
                event['difference'].shape)
            draws.append({'difference': event['difference']+perturbation,
                          'error': event['error']})
        difference, error = combine(draws)
        gf, _ = gaussfit(difference, error, seed, roi)
        result.append([gf['x'], gf['y']])
    return np.array(result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstrap', type=int, default=300)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    raw, url = retrieve()
    sha = hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest_path = OUT / 'raw_manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text())['sha256'] != sha:
        raise RuntimeError('Raw byte hash changed; retain mismatch for diagnosis.')
    with fits.open(raw) as hdus:
        checks = verify_fits_checksums(raw, hdus)
        header, table_header = dict(hdus[0].header), dict(hdus[1].header)
        wcs = WCS(hdus[2].header)
        data = hdus[1].data
        t = np.array(data['TIME'])
        quality = np.array(data['QUALITY'])
        flux = np.array(data['FLUX'], float)
        error = np.array(data['FLUX_ERR'], float)
        ffi = np.array(data['FFI_FILE']).astype(str)
        science = hashlib.sha256()
        for column, dtype in [('TIME', '<f8'), ('QUALITY', '<i4'),
                              ('FLUX', '<f4'), ('FLUX_ERR', '<f4')]:
            science.update(column.encode('ascii'))
            science.update(np.asarray(data[column], dtype=dtype).tobytes())
        science.update(np.asarray(ffi, dtype='S38').tobytes())
    if (header['TIMESYS'] != 'TDB' or header['TIMEREF'] != 'SOLARSYSTEM'
            or header['BJDREFI']+header['BJDREFF'] != 2457000):
        raise RuntimeError('Unexpected time reference.')
    dt = float(table_header['TIMEDEL'])
    host = wcs.all_world2pix([[RA, DEC]], 0)[0]
    physical = {
        'raw_column': table_header['1CRV5P'] +
                      (host[0]+1-table_header['1CRP5P'])*table_header['1CDL5P'],
        'raw_row': table_header['2CRV5P'] +
                   (host[1]+1-table_header['2CRP5P'])*table_header['2CDL5P']}
    if not (44 <= physical['raw_column'] < 2092 and
            0 <= physical['raw_row'] < 2048):
        raise RuntimeError('Target falls outside science silicon.')
    keep = ((quality == 0) & np.isfinite(t) &
            np.all(np.isfinite(flux) & np.isfinite(error), axis=(1, 2)))
    quality_counts = dict(zip(*np.unique(quality, return_counts=True)))
    total = len(t)
    t, flux, error = t[keep], flux[keep], error[keep]
    reference = np.median(flux, axis=0)
    yy, xx = np.indices(reference.shape)
    blank = reference < np.percentile(reference, 25)
    bg = np.median(flux[:, blank], axis=1)
    flux -= bg[:, None, None]
    refsky = reference - np.median(reference[blank])
    exclude = (xx-host[0])**2+(yy-host[1])**2 < 3.0**2
    gy, gx = np.gradient(refsky)
    reg = np.column_stack([refsky[~exclude], gx[~exclude], gy[~exclude],
                           np.ones(np.sum(~exclude))])
    coef = np.linalg.lstsq(reg, flux[:, ~exclude].T, rcond=None)[0]
    shifts = [-coef[1]/coef[0], -coef[2]/coef[0]]
    roi = (xx-host[0])**2+(yy-host[1])**2 < 3.0**2
    ref_fit, _ = gaussfit(refsky, np.ones_like(reference), host, roi)
    ref_center = np.array([ref_fit['x'], ref_fit['y']])
    modes, event_sets = {}, {}
    for mode in ['quadratic', 'background', 'pointing']:
        events = event_fits(t, flux, error, bg, shifts, dt, mode)
        event_sets[mode] = events
        modes[mode] = {}
        for parity, label in [(0, 'shallow'), (1, 'deep')]:
            selected = [event for event in events if event['parity'] == parity]
            difference, de = combine(selected)
            fit, _ = gaussfit(difference, de, host, roi)
            fit['offset_catalog_pixels'] = np.linalg.norm(
                np.array([fit['x'], fit['y']])-host)
            fit['offset_reference_pixels'] = np.linalg.norm(
                np.array([fit['x'], fit['y']])-ref_center)
            fit['sky_ra_dec'] = wcs.all_pix2world([[fit['x'], fit['y']]], 0)[0]
            modes[mode][label] = {
                'n_events': len(selected), 'source_fit': fit,
                'events': [{k: e[k] for k in ('cycle', 'center_btjd', 'n_in',
                           'n_out', 'design_condition')} for e in selected]}
            np.savez_compressed(OUT/f'difference_{mode}_{label}.npz',
                difference=difference, error=de,
                per_event_difference=np.array([e['difference'] for e in selected]),
                per_event_error=np.array([e['error'] for e in selected]))
    rng = np.random.default_rng(TIC)
    boots = {}
    # Quadratic is primary; pointing is a sensitivity check because the small
    # cutout does not contain an independent registration star.
    for parity, label in [(0, 'shallow'), (1, 'deep')]:
        selected = [e for e in event_sets['quadratic'] if e['parity'] == parity]
        fit = modes['quadratic'][label]['source_fit']
        boot = boot_centroids(selected, [fit['x'], fit['y']], roi,
                              args.bootstrap, rng)
        boots[label] = boot
        modes['quadratic'][label]['bootstrap'] = {
            'n': len(boot), 'median_xy': np.median(boot, axis=0),
            'percentile16_84_xy': np.percentile(boot, [16, 84], axis=0),
            'offset_catalog_pixels_16_50_84': np.percentile(
                np.linalg.norm(boot-host, axis=1), [16, 50, 84]),
            'offset_reference_pixels_16_50_84': np.percentile(
                np.linalg.norm(boot-ref_center, axis=1), [16, 50, 84]),
            'limitation': 'Whole-event plus synchronous five-cadence residual '
                          'bootstrap; no PRF/WCS systematic calibration.'}
    np.savez_compressed(OUT/'reference_and_centroid_bootstraps.npz',
        reference=reference, reference_sky_subtracted=refsky,
        shallow_centroid_xy=boots['shallow'], deep_centroid_xy=boots['deep'],
        target_xy=host, reference_gaussian_xy=ref_center,
        time_btjd=t, sky_e_s_pixel=bg, pointing_dx=shifts[0], pointing_dy=shifts[1])
    neighbors_path = ROOT/'candidates/tic71316387/catalog/gaia_neighbors.csv'
    if neighbors_path.exists():
        stars = pd.read_csv(neighbors_path, dtype={'source_id': str})
        stars[['x', 'y']] = wcs.all_world2pix(stars[['ra', 'dec']].to_numpy(), 0)
        stars['offset_tic_pixels'] = np.hypot(stars.x-host[0], stars.y-host[1])
        stars.to_csv(OUT/'gaia_pixel_positions.csv', index=False)
    else:
        stars = pd.DataFrame()
    roi_sensitivity = []
    nearest_neighbor = None
    if len(stars):
        non_host = stars[stars.source_id != '6571343840265900544']
        nearest_neighbor = non_host.loc[non_host.offset_tic_pixels.idxmin()]
    for label in ['shallow', 'deep']:
        arrays = np.load(OUT/f'difference_quadratic_{label}.npz')
        source = modes['quadratic'][label]['source_fit']
        fixed_tic, _ = gaussfit(arrays['difference'], arrays['error'], host,
                              roi, fixed_center=True)
        modes['quadratic'][label]['fixed_tic_center'] = {
            **fixed_tic, 'delta_chi2_diagonal_vs_free':
            fixed_tic['chi2_diagonal']-source['chi2_diagonal']}
        if nearest_neighbor is not None:
            neighbor = nearest_neighbor[['x', 'y']].to_numpy(dtype=float)
            fixed_neighbor, _ = gaussfit(arrays['difference'], arrays['error'],
                                       neighbor, roi, fixed_center=True)
            modes['quadratic'][label]['fixed_nearest_neighbor_center'] = {
                **fixed_neighbor, 'source_id': nearest_neighbor.source_id,
                'delta_chi2_diagonal_vs_free': fixed_neighbor['chi2_diagonal']-
                                              source['chi2_diagonal'],
                'offset_tic_pixels': float(nearest_neighbor.offset_tic_pixels),
                'caution': 'Diagonal-error Gaussian model comparison; not a '
                           'calibrated PRF source probability.'}
            modes['quadratic'][label]['bootstrap'][
                'fraction_closer_TIC_than_nearest_resolved_neighbor'] = float(
                    np.mean(np.linalg.norm(boots[label]-host, axis=1) <
                            np.linalg.norm(boots[label]-neighbor, axis=1)))
        for radius in [2.0, 2.5, 3.0, 3.5]:
            use = (xx-host[0])**2+(yy-host[1])**2 < radius**2
            rf, _ = gaussfit(refsky, np.ones_like(reference), host, use)
            df, _ = gaussfit(arrays['difference'], arrays['error'], host, use)
            roi_sensitivity.append({
                'label': label, 'radius_pixels': radius,
                'reference_x': rf['x'], 'reference_y': rf['y'],
                'difference_x': df['x'], 'difference_y': df['y'],
                'difference_minus_reference_pixels': float(np.hypot(
                    df['x']-rf['x'], df['y']-rf['y'])),
                'reference_width_bound_hit': rf['width_bound_hit'],
                'difference_width_bound_hit': df['width_bound_hit']})
    pd.DataFrame(roi_sensitivity).to_csv(OUT/'roi_sensitivity.csv', index=False)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for row, label in enumerate(['shallow', 'deep']):
        arrays = np.load(OUT/f'difference_quadratic_{label}.npz')
        fit = modes['quadratic'][label]['source_fit']
        for ax, image, title in zip(axes[row],
                [np.log10(np.maximum(refsky, 1)), arrays['difference'],
                 arrays['difference']/arrays['error']],
                ['Reference log10 e-/s', 'Eclipse difference e-/s',
                 'Difference / diagonal error']):
            limit = None if title.startswith('Reference') else np.percentile(abs(image), 99)
            im = ax.imshow(image, origin='lower', cmap='viridis' if limit is None else 'RdBu_r',
                           vmin=None if limit is None else -limit,
                           vmax=limit)
            ax.plot(*host, 'gx', ms=11, label='TIC position')
            ax.plot(*ref_center, 'c+', ms=10, label='Reference fit')
            ax.plot(fit['x'], fit['y'], 'm+', ms=13, label='Eclipse fit')
            if len(stars):
                ax.plot(stars['x'], stars['y'], 'o', mfc='none', mec='orange', ms=7)
            ax.set_xlim(-.5, 8.5); ax.set_ylim(-.5, 8.5)
            ax.set_title(f'{label}: {title}')
            fig.colorbar(im, ax=ax, shrink=.7)
    axes[0, 0].legend(loc='lower left', fontsize=8)
    fig.suptitle('TIC71316387 sector 28: both alternating eclipses')
    fig.tight_layout(); fig.savefig(OUT/'source_localization.png', dpi=170)
    plt.close(fig)
    (OUT/'ffi_product_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    manifest = {
        'tic': TIC, 'sector': 28, 'camera': header['CAMERA'], 'ccd': header['CCD'],
        'archive_url': url, 'raw_path': str(raw), 'sha256': sha,
        'science_array_sha256': science.hexdigest(), 'bytes': raw.stat().st_size,
        'fits_checksums': checks, 'processing_version': header['PROCVER'],
        'data_release': header['DATA_REL'], 'cube_date': header['DATE'],
        'date_obs': header['DATE-OBS'], 'date_end': header['DATE-END'],
        'time_system': header['TIMESYS'], 'time_reference': header['TIMEREF'],
        'bjdref': header['BJDREFI']+header['BJDREFF'], 'cadence_days': dt,
        'wcs_fits_keywords': dict(wcs.to_header()),
        'physical_target_pixels': physical,
        'ffi_identifiers': str((OUT/'ffi_product_identifiers.txt').relative_to(ROOT)),
        'checked_utc': datetime.now(timezone.utc).isoformat()}
    if not manifest_path.exists():
        manifest_path.write_text(json.dumps(safe(manifest), indent=2)+'\n')
    result = {
        'tic': TIC, 'sector': 28, 'checked_utc': datetime.now(timezone.utc).isoformat(),
        'source_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'ephemeris': {'period_days': PERIOD, 'epoch_btjd': EPOCH,
                      'parity_convention': 'even shallow, odd deep',
                      'shape_duration_days_even_odd': [.0617, .0608],
                      'shape_ingress_fractions_even_odd': [.262, .327]},
        'cadences_total': total, 'cadences_quality_zero_finite': len(t),
        'quality_counts': quality_counts, 'target_xy': host,
        'reference_gaussian_fit': ref_fit, 'modes': modes,
        'scientific_dependency_sha256': {
            name: hashlib.sha256((ROOT/'src'/name).read_bytes()).hexdigest()
            for name in ['archive.py', 'pixel_vet_49617264.py']},
        'science_array_hash_recipe': 'ASCII column names followed by little-endian '
            'C-order TIME f8, QUALITY i4, FLUX f4, FLUX_ERR f4, then S38 '
            'FFI_FILE strings.',
        'limitations': ['Integrated elliptical Gaussian is an empirical centroid '
                       'model, not a calibrated PRF or scene fit.',
                       'No unresolved-companion exclusion; centroid errors do '
                       'not include WCS and PRF systematics.',
                       'The small cutout lacks an independent registration '
                       'star; pointing regression is only a sensitivity check.',
                       'No source identity follows from aperture depth alone.']}
    (OUT/'pixel_vetting.json').write_text(json.dumps(safe(result), indent=2)+'\n')
    print(json.dumps(safe({label: modes['quadratic'][label]
                          for label in ('shallow', 'deep')}), indent=2), flush=True)


if __name__ == '__main__':
    main()
