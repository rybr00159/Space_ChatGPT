"""Separate, reproducible expanded screening; old sector-1 methods stay frozen.

Looser retention is a candidate-generation choice, not increased significance.
Two detrenders of the same data are correlated; scores must never be multiplied.
"""
from dataclasses import asdict
from pathlib import Path
import hashlib
import time

import numpy as np
from astropy.io import fits

from archive import verify_fits_checksums
from fast_screen import FastScreenConfig, screen_curve, preprocess_curve

CONFIG = dict(min_period=.5, max_period=20., durations=(.035, .05, .075, .10, .15),
              n_peaks=8, n_refine=32, compute_full_stats=True,
              grid_drift_durations=.8, scan_single_events=True)
WIDTHS = (.75, 1.5)


def screen_arrays(t, flux, error, quality=None, known=None, widths=WIDTHS):
    started = time.monotonic()
    scans = []
    for width in widths:
        config = dict(CONFIG, width_days=float(width))
        scan = screen_curve(t, flux, error, quality=quality, known=known, config=config)
        scan['width_days'] = float(width)
        scan['label'] = 'masked_residual' if known else 'unmasked'
        scans.append(scan)
    return dict(scans=scans, runtime_seconds=time.monotonic()-started,
                note='Correlated detrenders; statistics have no survey false-alarm calibration.')


def read_curve(path, tic=None, sector=None):
    path = Path(path)
    with fits.open(path, memmap=False) as hd:
        hd.verify('exception')
        checksum = verify_fits_checksums(path, hd)
        h, d = hd[1].header, hd[1].data
        ref = h.get('BJDREFI', 0) + h.get('BJDREFF', h.get('BJDREFR', 0))
        if h.get('TIMESYS') != 'TDB' or abs(ref-2457000) > 1e-8:
            raise ValueError('Unexpected time standard/reference')
        ident = h.get('TICID', hd[0].header.get('TICID'))
        sec = h.get('SECTOR', hd[0].header.get('SECTOR'))
        if tic is not None and int(ident) != int(tic):
            raise ValueError('Wrong target in product')
        if sector is not None and int(sec) != int(sector):
            raise ValueError('Wrong sector in product')
        column = next((k for k in ['DET_FLUX', 'KSPSAP_FLUX', 'PDCSAP_FLUX'] if k in d.names), None)
        if not column:
            raise ValueError('Unsupported photometric schema')
        arrays = [np.asarray(d[k], dtype=float) for k in ['TIME', column, column+'_ERR']]
        quality = np.asarray(d['QUALITY'])
        raw = np.asarray(d['SAP_FLUX'], dtype=float)
        metadata = {k: h.get(k, hd[0].header.get(k)) for k in
                    ['TICID','SECTOR','CAMERA','CCD','TESSMAG','TEFF','RADIUS','MASS','LOGG','PROCVER','CALIB','DETMET','TIMEDEL','TSTART','TSTOP','RA_OBJ','DEC_OBJ']}
    return dict(t=arrays[0], f=arrays[1], e=arrays[2], q=quality, raw=raw,
                column=column, metadata=metadata, fits_checksums=checksum,
                path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def analyze_file(record):
    curve = read_curve(record['path'], record['tic'], record['sector'])
    result = screen_arrays(curve['t'], curve['f'], curve['e'], curve['q'])
    result.update({k: v for k,v in record.items() if k != 'known'})
    result.update(metadata=curve['metadata'], flux_column=curve['column'],
                  sha256=curve['sha256'], fits_checksums=curve['fits_checksums'])
    if record.get('known'):
        residual = screen_arrays(curve['t'], curve['f'], curve['e'], curve['q'], record['known'])
        result['scans'].extend(residual['scans'])
    result['status'] = 'screened' if any(s['status']=='screened' for s in result['scans']) else 'insufficient_data'
    # Preserve compact prepared curves for exact independent follow-up.
    for width in WIDTHS:
        prepared = preprocess_curve(curve['t'], curve['f'], curve['e'], curve['q'], config=dict(CONFIG, width_days=width))
        dest = Path(record['prepared_dir']) / f'tic{record["tic"]}_s{record["sector"]:04d}_w{width:g}.npz'
        dest.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(dest, **{k: prepared[k] for k in ['t','f','e','segment']})
    return result


def seed_reasons(peak):
    reasons = list(peak.get('guard_failures', []))
    if peak.get('local_snr',0) < 5.5:
        reasons.append('local_snr_below_5p5_retention')
    if peak.get('n_positive_events',0) < max(2, .6*peak.get('n_events',0)):
        reasons.append('insufficient_positive_events')
    return reasons
