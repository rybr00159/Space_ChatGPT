"""Native TESSCut difference images for primary/secondary stellar eclipses.

The JSON configuration explicitly supplies ephemerides and template widths.
Empirical Gaussian fits diagnose resolved source contamination, not masses.
Raw pixels remain outside Git; archive/science identifiers and event evidence
are retained, and an existing expected raw SHA is never silently replaced.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import io
import json
import os
import urllib.parse
import urllib.request
import zipfile

os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import ROOT, verify_fits_checksums
from broad_vet import integrated_trapezoid
from pixel_vet_49617264 import gaussfit, combine, safe


def retrieve(config, sector):
    size = config.get('size', 11)
    cache = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'stellar_systems/pixels' / str(config['tic']) / f'sector{sector}'
    if config.get('analysis_label'):
        cache = cache / f'{size}x{size}'
    cache.mkdir(parents=True, exist_ok=True)
    url = 'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?' + urllib.parse.urlencode(
        dict(ra=config['ra'], dec=config['dec'], x=config.get('size', 11),
             y=config.get('size', 11), units='px', sector=sector))
    paths = list(cache.glob('tess-*.fits'))
    if not paths:
        with urllib.request.urlopen(url, timeout=120) as response:
            body = response.read(100000001)
        if len(body) > 100000000:
            raise RuntimeError('Cutout exceeds 100MB retrieval limit')
        with zipfile.ZipFile(io.BytesIO(body)) as z:
            for name in z.namelist():
                if name.endswith('.fits'):
                    p = cache / Path(name).name
                    if p.exists():
                        raise RuntimeError('Refusing raw overwrite')
                    p.write_bytes(z.read(name))
        paths = list(cache.glob('tess-*.fits'))
    if len(paths) != 1:
        raise RuntimeError(f'Expected one sector{sector} FITS, found {len(paths)}')
    return paths[0], url


def fit_events(t, f, e, exposure, config, kind, degree=2):
    eclipse = config['eclipses'][kind]
    p, epoch = config['period'], config['epoch_btjd'] + eclipse['phase'] * config['period']
    duration, ratio = eclipse['duration_days'], eclipse.get('ingress_fraction', .35)
    half = min(p * .20, max(duration * 2.8, .35))
    events = []
    for cycle in range(int(np.ceil((t.min()-epoch)/p)), int(np.floor((t.max()-epoch)/p))+1):
        center = epoch + p * cycle
        use = abs(t-center) < half
        x = t[use] - center
        if len(x) < 12:
            continue
        q = integrated_trapezoid(x, duration, ratio * duration, exposure)
        outside = q == 0
        inside = q > .05
        if (inside.sum() < 3 or (outside & (x < 0)).sum() < 8 or
                (outside & (x > 0)).sum() < 8 or
                np.max(np.diff(x)) > max(3.05 * exposure, .065)):
            continue
        design = np.column_stack([x ** k for k in range(degree+1)] + [-q])
        inverse = np.linalg.pinv(design)
        flux = f[use].reshape(len(x), -1)
        coef = inverse @ flux
        residual = flux - design @ coef
        scatter = np.maximum(np.std(residual, axis=0, ddof=design.shape[1]),
                             np.sqrt(np.mean(e[use].reshape(len(x), -1)**2, axis=0)))
        error = scatter * np.sqrt(np.sum(inverse[-1]**2))
        events.append(dict(cycle=cycle, time_btjd=center, n_in=int(inside.sum()),
                           n_out=int(outside.sum()), difference=coef[-1].reshape(f.shape[1:]),
                           error=error.reshape(f.shape[1:]), residual=residual,
                           projection=inverse[-1], condition=float(np.linalg.cond(design))))
    return events


def source_fit(events, roi, seed, stars, host):
    if not events:
        return None, None, None
    d, e = combine(events)
    fit, _ = gaussfit(d, e, seed, roi)
    xy = np.array([fit['x'], fit['y']])
    fit['offset_target_pixels'] = float(np.linalg.norm(xy-host))
    fixed = {}
    for label, position in stars.items():
        model, _ = gaussfit(d, e, position, roi, fixed_center=True)
        fixed[label] = dict(x=float(position[0]), y=float(position[1]),
                            chi2_diagonal=float(model['chi2_diagonal']),
                            delta_chi2_vs_free=float(model['chi2_diagonal']-fit['chi2_diagonal']))
    return dict(n_events=len(events), free_gaussian=fit, fixed_catalog_positions=fixed,
                events=[{k:v for k,v in ev.items() if k in ('cycle','time_btjd','n_in','n_out','condition')} for ev in events]), d, e


def bootstrap(events, roi, seed, n, rng):
    xy = []
    for _ in range(n):
        draw = []
        for i in rng.integers(0, len(events), len(events)):
            ev = events[i]
            residual = ev['residual']
            starts = rng.integers(0, len(residual), (len(residual)+4)//5)
            indices = np.concatenate([(s+np.arange(5))%len(residual) for s in starts])[:len(residual)]
            noise = (ev['projection'] @ residual[indices]).reshape(ev['difference'].shape)
            draw.append(dict(difference=ev['difference']+noise, error=ev['error']))
        d,e=combine(draw)
        fit,_=gaussfit(d,e,seed,roi)
        xy.append([fit['x'],fit['y']])
    return np.array(xy)


def aperture_depths(t, cube, errors, exposure, config, positions):
    """Measure diluted eclipse depths in target and candidate-source cores.

    Simultaneous local polynomial/template fits include baseline covariance.
    Error floors from event dispersion do not calibrate PRF or correlated noise.
    """
    yy, xx = np.indices(cube.shape[1:])
    results, arrays = {}, {}
    for name in config.get('photometry_positions', ['catalog_target']):
        if name not in positions:
            continue
        pos = positions[name]
        mask = (xx-pos[0])**2+(yy-pos[1])**2 < 1.0**2
        flux = np.sum(cube[:,mask],axis=1)
        error = np.sqrt(np.sum(errors[:,mask]**2,axis=1))
        results[name] = dict(n_pixels=int(mask.sum()),eclipses={})
        arrays[name+'_time'] = t
        arrays[name+'_raw_flux'] = flux
        for kind, eclipse in config['eclipses'].items():
            p=config['period'];epoch=config['epoch_btjd']+eclipse['phase']*p
            duration=eclipse['duration_days'];ratio=eclipse.get('ingress_fraction',.35)
            half=min(.2*p,max(2.8*duration,.35))
            ev=[];local_t=[];local_f=[];local_e=[]
            for cycle in range(int(np.ceil((t.min()-epoch)/p)),int(np.floor((t.max()-epoch)/p))+1):
                center=epoch+cycle*p;use=abs(t-center)<half
                x=t[use]-center;y=flux[use];ee=error[use]
                if len(x)<12:continue
                q=integrated_trapezoid(x,duration,ratio*duration,exposure)
                outside=q==0
                if ((q>.05).sum()<3 or (outside&(x<0)).sum()<8 or (outside&(x>0)).sum()<8
                        or np.max(np.diff(x))>max(3.05*exposure,.065)):continue
                design=np.column_stack([np.ones(len(x)),x,x*x,-q])
                floor=max(float(np.median(ee)),1e-6)
                a=design/np.maximum(ee,floor)[:,None];b=y/np.maximum(ee,floor)
                beta=np.linalg.lstsq(a,b,rcond=None)[0]
                residual=b-a@beta
                cov=np.linalg.pinv(a.T@a)*max(1.,float(residual@residual)/max(1,len(x)-4))
                baseline=design[:,:3]@beta[:3]
                if beta[0]<=0 or np.any(baseline<=0):continue
                depth=beta[-1]/beta[0]
                grad=np.array([-beta[-1]/beta[0]**2,0,0,1/beta[0]])
                err=np.sqrt(max(0.,float(grad@cov@grad)))
                ev.append(dict(cycle=cycle,time_btjd=center,depth=float(depth),depth_error=float(err),baseline_flux_e_s=float(beta[0])))
                local_t.extend(t[use]);local_f.extend(y/baseline);local_e.extend(ee/baseline)
            if ev:
                d=np.array([v['depth']for v in ev]);e=np.array([v['depth_error']for v in ev]);w=1/np.maximum(e,1e-10)**2
                mean=float(np.sum(w*d)/w.sum())
                err=max(float(1/np.sqrt(w.sum())),float(d.std(ddof=1)/np.sqrt(len(d)))if len(d)>1 else 0.)
                results[name]['eclipses'][kind]=dict(n_events=len(ev),depth=mean,error_max_diagonal_event_scatter=err,events=ev)
            else:results[name]['eclipses'][kind]=dict(n_events=0,depth=None)
            arrays[name+'_'+kind+'_time']=np.asarray(local_t)
            arrays[name+'_'+kind+'_flux']=np.asarray(local_f)
            arrays[name+'_'+kind+'_error']=np.asarray(local_e)
    return results,arrays


def analyze(config, config_path, sector, ndraws):
    out = ROOT / f"candidates/stellar_systems/tic{config['tic']}/pixels/sector{sector}"
    if config.get('analysis_label'):
        out = out / config['analysis_label']
    out.mkdir(parents=True, exist_ok=True)
    path, url = retrieve(config, sector)
    raw_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = out / 'raw_manifest.json'
    if manifest.exists() and json.loads(manifest.read_text())['sha256'] != raw_sha:
        raise RuntimeError('Raw retrieval differs: preserve and investigate')
    with fits.open(path) as hd:
        checks = verify_fits_checksums(path, hd)
        h, th = hd[0].header, hd[1].header
        wcs = WCS(hd[2].header)
        a = hd[1].data
        if h['TIMESYS'] != 'TDB' or h['BJDREFI'] + h['BJDREFF'] != 2457000:
            raise RuntimeError('Unexpected TDB coordinate')
        keep = ((a['QUALITY']==0)&np.isfinite(a['TIME'])&
                np.all(np.isfinite(a['FLUX']) & np.isfinite(a['FLUX_ERR']), axis=(1,2)))
        t=np.asarray(a['TIME'][keep],float)
        f=np.asarray(a['FLUX'][keep],float)
        e=np.asarray(a['FLUX_ERR'][keep],float)
        exposure=float(th['TIMEDEL'])
        ffi=np.asarray(a['FFI_FILE']).astype(str)
        science=hashlib.sha256()
        for name,dtype in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:
            science.update(name.encode());science.update(np.asarray(a[name],dtype=dtype).tobytes())
        host=wcs.all_world2pix([[config['ra'],config['dec']]],0)[0]
        column=th['1CRV5P']+(host[0]+1-th['1CRP5P'])*th['1CDL5P']
        row=th['2CRV5P']+(host[1]+1-th['2CRP5P'])*th['2CDL5P']
        raw_meta=dict(tic=config['tic'],sector=sector,url=url,path=str(path),
                      sha256=raw_sha,science_sha256=science.hexdigest(),bytes=path.stat().st_size,
                      checked_utc=datetime.now(timezone.utc).isoformat(),fits_checksums=checks,
                      native_cadences=len(a),quality_zero_finite=len(t),
                      camera=h.get('CAMERA'),ccd=h.get('CCD'),exposure_days=exposure,
                      detector_target_column=float(column),detector_target_row=float(row),
                      processing_version=h.get('PROCVER'),time_system='BTJD (TDB)',
                      preservation_note='Regenerated TessCut headers may differ; verify science arrays/FFI IDs before accepting a new raw expectation.')
        science_detector = bool(44<=column<2092 and 0<=row<2048)
        raw_meta['on_science_detector'] = science_detector
        if not science_detector:
            raw_meta['status'] = 'excluded_off_science_silicon'
            manifest.write_text(json.dumps(safe(raw_meta),indent=2)+'\n')
            (out/'ffi_identifiers.txt').write_text('\n'.join(ffi)+'\n')
            raise RuntimeError(f'Target outside science detector: column={column:.3f}, row={row:.3f}; raw retained and exclusion saved')
    manifest.write_text(json.dumps(safe(raw_meta),indent=2)+'\n')
    (out/'ffi_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    reference=np.median(f,axis=0)
    yy,xx=np.indices(reference.shape)
    blank=reference<np.percentile(reference,25)
    sky=np.median(f[:,blank],axis=1)
    f=f-sky[:,None,None]
    refsky=reference-np.median(reference[blank])
    gaia=pd.read_csv(ROOT/config['gaia_neighbors_file'],dtype={'source_id':str})
    gaia[['x','y']]=wcs.all_world2pix(gaia[['ra','dec']].to_numpy(),0)
    gaia.to_csv(out/'gaia_positions.csv',index=False)
    positions={'catalog_target':host}
    for _,star in gaia.iterrows():
        if (0<=star.x<reference.shape[1] and 0<=star.y<reference.shape[0]
                and star.phot_g_mean_mag<19 and np.hypot(star.x-host[0],star.y-host[1])>.35):
            positions['gaia_'+star.source_id]=np.array([star.x,star.y])
    radius = config.get('roi_radius_pixels', 3.5)
    roi=(xx-host[0])**2+(yy-host[1])**2<radius**2
    measurements, images = {}, {}
    rng=np.random.default_rng(config['tic']+sector)
    for kind in config['eclipses']:
        events=fit_events(t,f,e,exposure,config,kind)
        if not events:
            measurements[kind]=dict(status='insufficient_coverage',n_events=0);continue
        d,de=combine(events)
        seedxy=np.unravel_index(np.argmax(np.where(roi,d/np.maximum(de,.03),-np.inf)),d.shape)
        seed=np.array([seedxy[1],seedxy[0]],float)
        stats,d,de=source_fit(events,roi,seed,positions,host)
        center=np.array([stats['free_gaussian']['x'],stats['free_gaussian']['y']])
        boot=bootstrap(events,roi,center,ndraws,rng)
        stats['bootstrap']=dict(n=ndraws,median_xy=np.median(boot,axis=0),interval_16_84=np.percentile(boot,[16,84],axis=0),
                                note='Whole-event plus synchronized five-cadence residual blocks; excludes PRF/WCS systematics and unresolved companions.')
        stats['roi_sensitivity']={}
        for test_radius in [max(2.5,radius-.5),radius,radius+.5]:
            other=(xx-host[0])**2+(yy-host[1])**2<test_radius**2
            fit,_=gaussfit(d,de,center,other);stats['roi_sensitivity'][str(test_radius)]=fit
        linear=fit_events(t,f,e,exposure,config,kind,degree=1)
        stats['linear_baseline'],_,_=source_fit(linear,roi,center,positions,host)
        measurements[kind]=stats;images[kind]=(d,de)
        np.savez_compressed(out/f'{kind}_differences.npz',difference=d,error=de,
                            per_event_difference=np.array([v['difference']for v in events]),
                            per_event_error=np.array([v['error']for v in events]),bootstrap_xy=boot)
    report=dict(tic=config['tic'],sector=sector,checked_utc=datetime.now(timezone.utc).isoformat(),
                configuration=config,config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                raw_manifest=raw_meta,host_pixel=host,measurements=measurements,
                limitations=['Empirical integrated Gaussian is not a calibrated PRF.',
                             'Catalog positions are not proper-motion propagated.',
                             'Within-pixel covariance and search selection are incompletely calibrated.',
                             'Localization does not determine mass or prior recognition.'])
    phot, phot_arrays=aperture_depths(t,f,e,exposure,config,positions)
    report['aperture_photometry']=phot
    np.savez_compressed(out/'aperture_photometry.npz',**phot_arrays)
    (out/'pixel_vetting.json').write_text(json.dumps(safe(report),indent=2)+'\n')
    fig,axes=plt.subplots(1,len(images)+1,figsize=(5*(len(images)+1),5),squeeze=False)
    axes=axes[0];axes[0].imshow(np.log10(np.maximum(refsky,1)),origin='lower');axes[0].set_title('Median image (log flux)')
    for ax,(kind,(d,de))in zip(axes[1:],images.items()):
        ax.imshow(d,origin='lower',cmap='coolwarm');ax.set_title(kind+' eclipse difference')
        xy=measurements[kind]['free_gaussian'];ax.plot(xy['x'],xy['y'],'+',color='yellow',ms=12)
    for ax in axes:
        ax.plot(host[0],host[1],'x',color='lime',ms=10)
        for label,pos in positions.items():
            if label!='catalog_target':ax.plot(*pos,'o',mfc='none',mec='black',ms=8)
    fig.suptitle(f"TIC {config['tic']}: sector {sector} source tests");fig.tight_layout();fig.savefig(out/'localization.png',dpi=170);plt.close(fig)
    print(json.dumps(safe(dict(tic=config['tic'],sector=sector,host=host,
                              fits={k:dict(n_events=v['n_events'],free=v.get('free_gaussian'),bootstrap=v.get('bootstrap'))for k,v in measurements.items()})),indent=2),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('config',type=Path);p.add_argument('--sectors',type=int,nargs='+',required=True);p.add_argument('--bootstrap',type=int,default=200)
    a=p.parse_args();config=json.loads(a.config.read_text())
    status=[]
    for sector in a.sectors:
        try:
            analyze(config,a.config,sector,a.bootstrap)
            status.append(dict(sector=sector,status='completed'))
        except Exception as error:
            status.append(dict(sector=sector,status='failed_or_excluded',exception=type(error).__name__,reason=str(error)))
            print(json.dumps(status[-1]),flush=True)
    folder=ROOT/f"candidates/stellar_systems/tic{config['tic']}"
    folder.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    (folder/f'pixel_run_{stamp}.json').write_text(json.dumps(dict(config=str(a.config),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),results=status),indent=2)+'\n')
    if any(v['status']!='completed'for v in status):
        raise SystemExit(1)


if __name__=='__main__':main()
