"""Actual-data controls for sector-by-sector transit screening.

Synthetic injections modify retained observed flux before detrending. No white
noise replacement is used. These small conditional tests are not completeness
measurements for TESS, a false-alarm calibration, or independent validation of
TIC 4206066. Original scientific inputs and old outputs are read-only.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import os
import sqlite3
import sys
import time

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from archive import verify_fits_checksums
from fast_screen import screen_curve, preprocess_curve, event_stats, phase

OUT = ROOT / 'results/multisector_controls'
REPORT = ROOT / 'candidates/multisector_controls'
CACHE = Path('/workspace/tess-research-data/multisector/controls')
SEED = 20261016
OLD_CONFIG = dict(n_refine=0, n_peaks=2, compute_full_stats=False)
CONTROL_PLANETS = [(307210830, 'L 98-59 b'), (79748331, 'TOI-1064 b'),
                   (306996324, 'TOI-776 b')]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read_observation(record):
    path = Path(record.get('local_path', record.get('path')))
    if sha(path) != record['sha256']:
        raise ValueError('Retained observation hash mismatch: ' + str(path))
    with fits.open(path, memmap=False) as hd:
        hd.verify('exception')
        checks = verify_fits_checksums(path, hd)
        h, d = hd[1].header, hd[1].data
        if h.get('TIMESYS') != 'TDB' or h.get('BJDREFI', 0) + h.get('BJDREFF', 0) != 2457000:
            raise ValueError('Unexpected observation time reference')
        tic = int(h.get('TICID', hd[0].header.get('TICID')))
        sector = int(h.get('SECTOR', hd[0].header.get('SECTOR')))
        if tic != int(record['tic']) or sector != int(record['sector']):
            raise ValueError('Observation target/sector mismatch')
        column = next(c for c in ('PDCSAP_FLUX', 'DET_FLUX', 'KSPSAP_FLUX') if c in d.names)
        arr = tuple(np.array(d[k], dtype=float) for k in ('TIME', column, column + '_ERR'))
        q = np.array(d['QUALITY'])
        provenance = dict(tic=tic, sector=sector, path=str(path), sha256=sha(path),
                          url=record.get('retrieval_url', record.get('url')), bytes=path.stat().st_size,
                          column=column, n_input=len(q), fits_checksums=checks,
                          procver=h.get('PROCVER', hd[0].header.get('PROCVER')))
    return (*arr, q), provenance


def all_observations():
    records = json.loads((ROOT / 'data/observations.json').read_text())
    for p in sorted((ROOT / 'candidates/tic4206066').glob('*manifest.json')):
        r = json.loads(p.read_text())
        if r.get('http_status') == 200 and 'sector' in r:
            records.append(r)
    return records


def method_hashes(method):
    names = ['src/multisector_controls.py', 'src/fast_screen.py']
    if method == 'new':
        names.append('src/multisector_screen.py')
    return {p: sha(ROOT / p) for p in names}


def scan(arrays, method, known=None):
    t, f, e, q = arrays
    if method == 'old':
        return [screen_curve(t, f, e, quality=q, known=known, config=OLD_CONFIG)]
    from multisector_screen import screen_arrays
    result = screen_arrays(t, f, e, quality=q, known=known)
    return result['scans']


def flat_peaks(scans):
    return [dict(variant=i, **p) for i, s in enumerate(scans) for p in s.get('peaks', [])]


def recovery(peaks, period, epoch, duration, times, threshold):
    """Require period AND predicted epoch agreement across the observed span."""
    finite = times[np.isfinite(times)]
    cycles = np.arange(np.ceil((finite.min() - epoch) / period),
                       np.floor((finite.max() - epoch) / period) + 1)
    centers = epoch + cycles * period
    matched = []
    for p in peaks:
        ratio_error = abs(p['period'] / period - 1)
        drift = float(np.max(abs(phase(centers, p['period'], p['t0_btjd'])))) if len(centers) else 1e99
        # The tolerance is fixed independently of fitted box width, so a broad
        # spurious box cannot win merely by encompassing the injected phase.
        exact = ratio_error < .003 and drift <= .75 * duration
        if exact:
            matched.append(dict(period=p['period'], epoch_btjd=p['t0_btjd'],
                                duration_days=p['duration_days'], local_snr=p['local_snr'],
                                epoch_max_error_days=drift, screening_pass=p['screening_pass'],
                                queued=bool(p['screening_pass'] and p['local_snr'] >= threshold),
                                variant=p['variant'], guard_failures=p['guard_failures']))
    return dict(found_period_and_epoch=bool(matched),
                queued=any(p['queued'] for p in matched), matches=matched)


def combine_fixed(rows):
    usable = [r for r in rows if r['n_events'] and r['local_depth_err'] and r['local_depth_err'] > 0]
    if not usable:
        return dict(n_sectors=0, n_events=0, depth_ppm=None, diagonal_snr=None)
    w = np.array([1 / r['local_depth_err']**2 for r in usable])
    d = np.array([r['local_depth'] for r in usable])
    depth = np.average(d, weights=w)
    return dict(n_sectors=len(usable), n_positive_sectors=int((d > 0).sum()),
                n_events=sum(r['n_events'] for r in usable), depth_ppm=float(depth * 1e6),
                depth_error_ppm=float(1e6 / np.sqrt(w.sum())),
                diagonal_snr=float(depth * np.sqrt(w.sum())),
                sector_depth_chi2=float(np.sum(w * (d - depth)**2)),
                sector_depth_chi2_dof=len(usable)-1,
                limitation='Published fixed ephemeris; diagonal uncertainty; no trials/FAP calibration.')


def controls(method):
    known = list(csv.DictReader((ROOT / 'data/known_planets.csv').open()))
    actual = list(csv.DictReader((ROOT / 'data/catalogs/confirmed_transiting.csv').open()))
    records = all_observations()
    specs = []
    for tic, name in CONTROL_PLANETS:
        p = next(r for r in known if r['name'] == name)
        a = next(r for r in actual if r['pl_name'] == name)
        specs.append(dict(tic=tic, name=name, period=float(p['period']),
                          epoch_btjd=float(p['epoch_bjd']) - 2457000,
                          duration_days=float(p['duration_hours']) / 24,
                          catalog_depth_ppm=float(a['pl_trandep']) * 1e4,
                          status='known_confirmed_planet_control'))
    specs.append(dict(tic=4206066, name='TIC 4206066 reported candidate', period=3.182785,
                      epoch_btjd=3998.736, duration_days=1.87 / 24,
                      status='reported_candidate_recovery_not_validation'))
    outputs, observations = [], []
    for spec in specs:
        rows = []
        others = [r for r in known if int(r['tic']) == spec['tic'] and r['name'] != spec['name']]
        for record in sorted((r for r in records if int(r['tic']) == spec['tic']), key=lambda r:r['sector']):
            arrays, prov = read_observation(record)
            observations.append(prov)
            scans = scan(arrays, method, known=others)
            peaks = flat_peaks(scans)
            fixed = []
            widths = (.75, 1.5) if method == 'new' else (1.,)
            for width in widths:
                prepared = preprocess_curve(*arrays[:3], quality=arrays[3], known=others,
                                            config=dict(width_days=width))
                s = event_stats(prepared['t'], prepared['f'], prepared['e'], spec['period'],
                                spec['epoch_btjd'], spec['duration_days'])
                fixed.append(dict(width_days=width, **s))
            row = dict(sector=record['sector'], scans=scans, fixed_ephemeris=fixed,
                       recovery=recovery(peaks, spec['period'], spec['epoch_btjd'],
                                         spec['duration_days'], arrays[0], 5.5 if method=='new' else 8.))
            rows.append(row)
            print('CONTROL', method, spec['name'], record['sector'], row['recovery']['found_period_and_epoch'],
                  row['recovery']['queued'], [round(x['local_snr'],2) for x in fixed], flush=True)
        outputs.append(dict(**spec, sectors=rows,
                            combined_fixed=[dict(width_days=w, **combine_fixed([r['fixed_ephemeris'][i] for r in rows]))
                                            for i,w in enumerate((.75,1.5) if method=='new' else (1.,))]))
        dump(OUT / f'actual_controls_{method}.json', dict(method=method, methods_sha256=method_hashes(method),
             generated_utc=datetime.now(timezone.utc).isoformat(), observations=observations, controls=outputs,
             fixed_ephemeris_limitations='Published ephemerides supply the answer; pooling is only a consistency/recovery control.'))


def backgrounds():
    # Fixed before any injection recovery results. The QLP artifact is the first
    # retained file returned by the initial filesystem probe, without SNR filtering.
    obs = all_observations()
    selected = [next(r for r in obs if r['tic']==150359842 and r['sector']==93)]
    db = Path('/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite')
    with sqlite3.connect(db.as_uri()+'?mode=ro', uri=True) as c:
        row = c.execute('SELECT result_json FROM stars WHERE tic=24293000').fetchone()
    record = json.loads(row[0])
    record['sector'] = 1
    selected.append(record)
    return selected


def integrated_box(t, period, epoch, duration, cadence):
    # Finite-exposure rectangular transit, 11-point exposure integration.
    offsets = ((np.arange(11)+.5)/11-.5) * cadence
    return np.mean(abs(phase(t[:,None] + offsets, period, epoch)) < duration/2, axis=1)


def injections(method, phases):
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    target = OUT / f'injections_{method}.jsonl'
    methods = method_hashes(method)
    config = dict(method=method, seed=SEED, phases=phases, periods=[3.18,7.3,11.13],
                  depths_ppm=[300,500,1000], duration_days=.08, methods_sha256=methods,
                  selection='Fixed TIC150359842 S93 and TIC24293000 S1, chosen before injection outcomes; no quiet-star selection.',
                  exposure_model='11-midpoint exposure-integrated rectangular box injected multiplicatively before detrending',
                  recovery='Period relative error<.003 AND maximum transit-center epoch discrepancy across sector <=.75 injected duration; queued also passes event guards and local SNR threshold.',
                  limitations=['Two backgrounds only; conditional sensitivity, not representative completeness.',
                               'Box shapes, fixed duration, no limb darkening/stellar radius distribution.',
                               'Raw retained noise and systematics remain; no generated white-noise substitution.',
                               'Old/new queue thresholds differ (8 versus5.5), report finding and queueing separately.',
                               'No trials-corrected significance or false-alarm probability.'])
    manifest = OUT / f'injections_{method}_config.json'
    if manifest.exists() and json.loads(manifest.read_text()) != config:
        raise ValueError('Existing injection experiment differs; preserve it and choose new output names')
    dump(manifest, config)
    done = set()
    if target.exists():
        for line in target.read_text().splitlines():
            done.add(json.loads(line)['case_id'])
    rng = np.random.default_rng(SEED)
    # Phase design shared across backgrounds and both methods.
    phase_grid = (np.arange(phases) + rng.uniform(.1,.9,size=phases))/phases
    observations = []
    for record in backgrounds():
        arrays, prov = read_observation(record)
        observations.append(prov)
        t, f, e, q = arrays
        valid = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (f>0) & (e>0) & (q==0)
        cadence = float(np.median(np.diff(t[valid])[np.diff(t[valid])>0]))
        baseline = float(np.nanmin(t))
        unchanged = scan(arrays, method)
        dump(OUT / f'unchanged_tic{record["tic"]}_s{record["sector"]}_{method}.json',
             dict(observation=prov, scans=unchanged))
        with target.open('a') as output:
            for period in config['periods']:
                for depth in config['depths_ppm']:
                    for j, fraction in enumerate(phase_grid):
                        case = f'{record["tic"]}:{record["sector"]}:{period}:{depth}:{j}'
                        if case in done:
                            continue
                        epoch = baseline + fraction * period
                        model = integrated_box(t,period,epoch,.08,cadence)
                        injected = f * (1-depth*1e-6*model)
                        scans = scan((t,injected,e,q),method)
                        peaks = flat_peaks(scans)
                        prepared = preprocess_curve(t,injected,e,quality=q)
                        fixed = event_stats(prepared['t'],prepared['f'],prepared['e'],period,epoch,.08)
                        result = dict(case_id=case,tic=record['tic'],sector=record['sector'],period_days=period,
                                      depth_ppm=depth,phase_index=j,epoch_btjd=float(epoch),
                                      n_exposure_overlapped=int(np.sum(model[valid]>0)),
                                      fixed_ephemeris_n_events=fixed['n_events'],
                                      fixed_ephemeris_snr=fixed['local_snr'],
                                      **recovery(peaks,period,epoch,.08,t,5.5 if method=='new' else 8.),
                                      peaks=[{k:p[k] for k in ['variant','period','t0_btjd','duration_days','local_snr','screening_pass','guard_failures','n_events','n_positive_events']} for p in peaks])
                        output.write(json.dumps(result,allow_nan=False)+'\n'); output.flush()
                        done.add(case)
                    print('INJECT',method,record['tic'],period,depth,'cases',len(done),flush=True)
    dump(OUT / f'injections_{method}_observations.json', observations)
    summarize()


def summarize():
    rows = []
    for p in sorted(OUT.glob('injections_*.jsonl')):
        method = p.stem.removeprefix('injections_')
        cases = [json.loads(line) for line in p.read_text().splitlines()]
        groups = sorted({(r['tic'],r['sector'],r['period_days'],r['depth_ppm']) for r in cases})
        for group in groups:
            g = [r for r in cases if (r['tic'],r['sector'],r['period_days'],r['depth_ppm'])==group]
            rows.append(dict(method=method,tic=group[0],sector=group[1],period_days=group[2],depth_ppm=group[3],
                             n=len(g),found=sum(r['found_period_and_epoch'] for r in g),
                             queued=sum(r['queued'] for r in g),
                             fewer_than_two_sampled_events=sum(r['fixed_ephemeris_n_events']<2 for r in g),
                             median_fixed_ephemeris_diagonal_snr=float(np.median([r['fixed_ephemeris_snr'] for r in g]))))
    if rows:
        with (OUT/'injection_summary.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n'); w.writeheader();w.writerows(rows)
    print('SUMMARY',len(rows),'groups',flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('stage',choices=['controls','injections','summarize'])
    p.add_argument('--method',choices=['old','new'],default='old')
    p.add_argument('--phases',type=int,default=16)
    a=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True); REPORT.mkdir(parents=True,exist_ok=True)
    if a.stage=='controls': controls(a.method)
    elif a.stage=='injections': injections(a.method,a.phases)
    else: summarize()


if __name__=='__main__':
    main()
