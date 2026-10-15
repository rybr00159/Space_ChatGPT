#!/usr/bin/env python3
"""A fixed, image-level pilot for Gaia-missing stars near bright moving stars.

This is a search, not a discovery classifier. Every coadd point-source detection
is crossmatched against Gaia DR3, Legacy Survey DR10, and Pan-STARRS DR2. Model
residual detections are explicitly subjected to halo and morphology controls.
Original HTTP response bytes and failed responses are retained outside Git.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from io import BytesIO, StringIO
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = Path('/workspace/tess-research-data/new_stars/imaging')
OUT = ROOT / 'data/new_stars/imaging'
PLOTS = ROOT / 'candidates/new_stars/imaging'
TAP = 'https://gea.esac.esa.int/tap-server/tap/sync'
LEGACY_TAP = 'https://datalab.noirlab.edu/tap/sync'
PIXEL = .262
SIZE = 256
PILOT_QUERY = '''SELECT TOP 12 source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag
FROM gaiadr3.gaia_source WHERE ra BETWEEN 170 AND 240 AND dec BETWEEN -10 AND 25
AND phot_g_mean_mag BETWEEN 9 AND 12 AND parallax>20
AND SQRT(pmra*pmra+pmdec*pmdec)>200 ORDER BY phot_g_mean_mag'''
EXPANDED_QUERY = '''SELECT source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag
FROM gaiadr3.gaia_source WHERE ra BETWEEN 170 AND 240 AND dec BETWEEN -10 AND 25
AND phot_g_mean_mag BETWEEN 10 AND 15 AND parallax>10
AND SQRT(pmra*pmra+pmdec*pmdec)>150 ORDER BY source_id'''


def retrieve(path, url, params=None):
    """Cache response bytes and status; only valid successful bytes are reused."""
    path = RAW / path
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = path.with_suffix(path.suffix + '.json')
    if path.exists() and meta.exists():
        info = json.loads(meta.read_text())
        if info['status'] == 200 and sha256(path.read_bytes()).hexdigest() == info['sha256']:
            return path, info
    last = None
    for attempt in range(3):
        try:
            response = requests.get(url, params=params, timeout=120)
            payload = response.content
            info = dict(url=response.url, status=response.status_code,
                        content_type=response.headers.get('Content-Type'),
                        retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                        bytes=len(payload), sha256=sha256(payload).hexdigest())
            path.write_bytes(payload)
            meta.write_text(json.dumps(info, indent=2) + '\n')
            if response.status_code == 200:
                return path, info
            last = RuntimeError(f'HTTP {response.status_code}: {response.url}')
        except requests.RequestException as error:
            last = error
        time.sleep(attempt + 1)
    raise last


def query(path, endpoint, adql):
    path, info = retrieve(path, endpoint, dict(REQUEST='doQuery', LANG='ADQL',
                                             FORMAT='csv', QUERY=adql))
    table = pd.read_csv(path, dtype={'source_id': str, 'ref_id': str, 'objid': str,
                                    'objID': str, 'ls_id': str})
    if '<' in str(table.columns[0]):
        raise ValueError(f'Non-CSV TAP response: {path}')
    return table, info


def load_field(row):
    sid, ra, dec = str(row.source_id), row.ra, row.dec
    products = {}
    for layer in ['ls-dr10', 'ls-dr10-model', 'ls-dr10-resid']:
        path, info = retrieve(f'{sid}/{layer}.fits',
                              'https://www.legacysurvey.org/viewer/fits-cutout',
                              dict(ra=ra, dec=dec, layer=layer, pixscale=PIXEL,
                                   size=SIZE, bands='grz'))
        with fits.open(path, checksum=True) as hdus:
            cube = np.asarray(hdus[0].data, float)
            header = hdus[0].header.copy()
        if cube.shape != (3, SIZE, SIZE) or np.count_nonzero(cube) == 0:
            raise ValueError(f'No usable three-band coverage: {sid}/{layer}')
        products[layer] = (cube, header, info)
    radius = .015
    gq = f'''SELECT source_id,ra,dec,pmra,pmdec,phot_g_mean_mag
    FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),
    CIRCLE('ICRS',{ra},{dec},{radius}))'''
    gaia, gi = query(f'{sid}/gaia.csv', TAP, gq)
    lq = f'''SELECT ls_id,ra,dec,type,flux_g,flux_r,flux_z,psfsize_g,psfsize_r,psfsize_z,
    maskbits,ref_cat,ref_id FROM ls_dr10.tractor WHERE
    ra BETWEEN {ra-radius/np.cos(np.deg2rad(dec))} AND {ra+radius/np.cos(np.deg2rad(dec))}
    AND dec BETWEEN {dec-radius} AND {dec+radius}'''
    # Keep the failed first query (ADQL interpreted bare true as a column) in
    # legacy.csv; v2 intentionally includes overlapping-brick duplicates.
    legacy, li = query(f'{sid}/legacy_v2.csv', LEGACY_TAP, lq)
    cols = '[objID,raMean,decMean,nDetections,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag]'
    pp, pi = retrieve(f'{sid}/ps1.csv',
                      'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv',
                      dict(ra=ra, dec=dec, radius=radius, **{'nDetections.gte': 2}, columns=cols))
    ps1 = pd.read_csv(pp, dtype={'objID': str})
    if not {'raMean', 'decMean'}.issubset(ps1.columns):
        raise ValueError(f'Invalid PS1 response: {pp}')
    return row, products, gaia, legacy, ps1


def nearest(ra, dec, table, ra_col='ra', dec_col='dec'):
    if len(table) == 0:
        return np.inf, -1
    sep = np.hypot((table[ra_col].to_numpy()-ra)*np.cos(np.deg2rad(dec)),
                   table[dec_col].to_numpy()-dec)*3600
    idx = int(np.nanargmin(sep))
    return float(sep[idx]), idx


def robust_scale(values):
    values = values[np.isfinite(values)]
    return max(1.4826*np.median(np.abs(values-np.median(values))), 1e-12)


def search_field(bundle):
    row, products, gaia, legacy, ps1 = bundle
    sid = str(row.source_id)
    image = products['ls-dr10'][0]
    model = products['ls-dr10-model'][0]
    residual = products['ls-dr10-resid'][0]
    wcs = WCS(products['ls-dr10'][1]).celestial
    yc, xc = np.indices(image.shape[1:])
    center = np.array(wcs.world_to_pixel_values(row.ra, row.dec))
    rr = np.hypot(xc-center[0], yc-center[1])*PIXEL
    inner=1.5 if row.phot_g_mean_mag>=12 else 3.
    eligible = (rr > inner) & (rr < 28) & (xc > 12) & (xc < SIZE-12) & (yc > 12) & (yc < SIZE-12)
    annulus = (rr > 20) & (rr < 30)
    # Detection kernel FWHM 1.23 arcsec; NOT assumed to represent every PSF.
    smooth = np.array([gaussian_filter(z, 2) for z in image])
    rsmooth = np.array([gaussian_filter(z, 2) for z in residual])
    snr = np.array([(z-np.median(z[annulus])) / robust_scale(z[annulus]) for z in smooth])
    rsnr = np.array([(z-np.median(z[annulus])) / robust_scale(z[annulus]) for z in rsmooth])
    # Second-largest band SNR demands independent-filter recurrence.
    detection = np.sort(snr, axis=0)[1]
    rdetection = np.sort(rsnr, axis=0)[1]
    original_peaks=(detection>=6)&(detection==maximum_filter(detection,size=9))
    residual_peaks=(rdetection>=6)&(rdetection==maximum_filter(rdetection,size=9))
    peaks = eligible & (original_peaks|residual_peaks)
    ys, xs = np.where(peaks)
    records = []
    for y, x in zip(ys, xs):
        ra, dec = wcs.pixel_to_world_values(x, y)
        gd, gj = nearest(ra, dec, gaia)
        ld, lj = nearest(ra, dec, legacy)
        pdist, pj = nearest(ra, dec, ps1, 'raMean', 'decMean')
        # Repeat Gaia nearest-match after shifting reference positions over a
        # deliberately generous +/-8 yr window: proper motion is not novelty.
        motion_sep = gd
        if len(gaia):
            dx = (gaia.ra.to_numpy()-ra)*np.cos(np.deg2rad(dec))*3600
            dy = (gaia.dec.to_numpy()-dec)*3600
            vx = np.nan_to_num(gaia.pmra.to_numpy())*.001
            vy = np.nan_to_num(gaia.pmdec.to_numpy())*.001
            dt = np.clip(-(dx*vx+dy*vy)/(vx*vx+vy*vy+1e-30), -8, 8)
            motion_sep = float(np.min(np.hypot(dx+dt*vx, dy+dt*vy)))
        stamp = smooth[:, y-6:y+7, x-6:x+7].copy()
        yy, xx = np.indices(stamp.shape[1:])
        weights = np.maximum(np.median(stamp, axis=0)-np.median(stamp[:, [0, -1], :]), 0)
        norm = weights.sum()
        if norm > 0:
            mx, my = (weights*xx).sum()/norm, (weights*yy).sum()/norm
            cov = np.array([[(weights*(xx-mx)**2).sum(), (weights*(xx-mx)*(yy-my)).sum()],
                            [(weights*(xx-mx)*(yy-my)).sum(), (weights*(yy-my)**2).sum()]])/norm
            eig = np.linalg.eigvalsh(cov)
            axis_ratio = np.sqrt(max(eig[0], 0)/max(eig[1], 1e-30))
            width = np.sqrt(eig.sum()/2)*PIXEL
        else:
            axis_ratio, width = 0., np.inf
        # Halo-aware local annulus threshold is separate from blank-sky SNR.
        local = (np.hypot(xc-x,yc-y)>7) & (np.hypot(xc-x,yc-y)<16)
        local_snr = np.sort([(z[y,x]-np.median(z[local]))/robust_scale(z[local]) for z in rsmooth])[1]
        catalogued = ld < 1.5 or pdist < 1.5 or motion_sep < 1.5
        status = ('catalogued_gaia' if motion_sep < 1.5 else
                  'catalogued_deeper_survey' if catalogued else 'unmatched_image_peak')
        if not catalogued:
            if rdetection[y,x] < 6:
                status = 'residual_not_recurrent'
            elif axis_ratio < .65 or not (.35 < width < 1.0):
                status = 'nonstellar_morphology'
            elif local_snr < 6:
                status = 'halo_local_noise_failure'
            else:
                status = 'requires_independent_imaging'
        # Saturated DECam stars can leave almost uniform horizontal/vertical
        # trails that pass a small-stamp roundness test after smoothing. Detect
        # the extended trail, retaining rather than deleting these peaks.
        alignment = min(abs(x-center[0]), abs(y-center[1]))*PIXEL
        trail_fraction = 0.
        if alignment < .6:
            horizontal = abs(y-center[1]) < abs(x-center[0])
            strip = rsmooth[:, y-1:y+2, max(x-18,0):min(x+19,SIZE)] if horizontal else rsmooth[:, max(y-18,0):min(y+19,SIZE), x-1:x+2]
            means = np.median(strip,axis=1 if horizontal else 2)
            # Farther than 3 pixels from the peak: Gaussian PSFs fall steeply,
            # whereas saturated trails maintain nearly constant surface flux.
            mid=means.shape[1]//2
            wings=np.r_[max(0,mid-18):max(0,mid-7),min(mid+8,means.shape[1]):min(mid+19,means.shape[1])]
            ratios=[]
            for b in range(3):
                peak=rsmooth[b,y,x]-np.median(rsmooth[b][annulus])
                if peak>0 and len(wings): ratios.append(float(np.median(means[b,wings]-np.median(rsmooth[b][annulus]))/peak))
            if ratios: trail_fraction=float(np.median(ratios))
            if not catalogued and trail_fraction>.5:
                status='saturated_trail_artifact'
        phot=gaussian_photometry(image,x,y)
        nominal=np.array([np.nanmedian(legacy['psfsize_'+b]) for b in 'grz'])
        compact=[j for j in range(3) if phot[j]['snr']>4 and .65<phot[j]['fwhm']/nominal[j]<1.6
                 and abs(phot[j]['x']-11)<3.6 and abs(phot[j]['y']-11)<3.6]
        pairs=[np.hypot(phot[a]['x']-phot[b]['x'],phot[a]['y']-phot[b]['y'])*PIXEL
               for ai,a in enumerate(compact) for b in compact[ai+1:]]
        coherent_pair=bool(pairs and min(pairs)<.6)
        if status=='requires_independent_imaging' and not coherent_pair:
            status='incoherent_or_broad_band_profiles'
        records.append(dict(field_gaia_id=sid, x=int(x), y=int(y), ra=float(ra),dec=float(dec),
                            original_peak=bool(original_peaks[y,x]),residual_peak=bool(residual_peaks[y,x]),
                            separation_arcsec=float(rr[y,x]), snr_second_band=float(detection[y,x]),
                            residual_snr_second_band=float(rdetection[y,x]), residual_local_snr=float(local_snr),
                            axis_ratio=float(axis_ratio), width_arcsec=float(width),
                            gaia_sep_arcsec=gd, gaia_motion_sep_arcsec=motion_sep,
                            legacy_sep_arcsec=ld, legacy_type=str(legacy.iloc[lj]['type']) if lj >= 0 else '',
                            legacy_maskbits=int(legacy.iloc[lj]['maskbits']) if lj >= 0 else -1,
                            legacy_id=str(legacy.iloc[lj]['ls_id']) if lj >= 0 else '',
                            ps1_sep_arcsec=pdist,ps1_id=str(ps1.iloc[pj]['objID']) if pj >= 0 else '',
                            axis_alignment_arcsec=float(alignment),trail_wing_fraction=trail_fraction,
                            gaussian_flux_g=phot[0]['flux'],gaussian_flux_r=phot[1]['flux'],
                            gaussian_flux_z=phot[2]['flux'],
                            gaussian_fwhm_g_arcsec=phot[0]['fwhm'],gaussian_fwhm_r_arcsec=phot[1]['fwhm'],
                            gaussian_fwhm_z_arcsec=phot[2]['fwhm'],compact_band_count=len(compact),
                            compact_coherent_pair=coherent_pair,
                            gaussian_position_offset_gr_arcsec=float(np.hypot(phot[0]['x']-phot[1]['x'],phot[0]['y']-phot[1]['y'])*PIXEL),
                            status=status))
    # Bright-star halo makes a global noise threshold optimistic. Its rejected
    # peaks remain in the output; they are not quietly removed.
    summary = dict(gaia_id=sid, ra=float(row.ra), dec=float(row.dec),
                   g_mag=float(row.phot_g_mean_mag), parallax_mas=float(row.parallax),
                   pm_masyr=float(np.hypot(row.pmra,row.pmdec)),gaia_catalog_rows=len(gaia),
                   legacy_catalog_rows=len(legacy),ps1_catalog_rows=len(ps1),
                   detections=len(records), image_model_residual_max_abs=float(np.max(np.abs(image-model-residual))),
                   files=[v[2] for v in products.values()])
    # Positive control is an actually detected catalogued point source, rather
    # than an assumption that the image endpoint returned a useful sky field.
    summary['matched_psf_positive_controls'] = sum(r['legacy_sep_arcsec']<1.5 and r['legacy_type']=='PSF' for r in records)
    make_plot(sid, image, residual, records, center)
    injections=injection_controls(row,smooth,annulus,center,legacy,wcs)
    return records, summary, injections


def gaussian_photometry(image,x,y):
    """Descriptive Gaussian+plane fit; correlated coadds lack error maps here.

    The Gaussian flux is in the Legacy viewer's nanomaggy convention. Width is
    a morphology diagnostic, not a PSF classification or calibrated photometry.
    Saturated or extended structures may hit the fit bounds and remain flagged.
    """
    result=[]
    for b in range(3):
        stamp=image[b,y-11:y+12,x-11:x+12]
        yy,xx=np.indices(stamp.shape)
        noise=robust_scale(stamp)
        def fun(p):
            flux,x0,y0,sig,a,bx,by=p
            prediction=flux/(2*np.pi*sig*sig)*np.exp(-((xx-x0)**2+(yy-y0)**2)/(2*sig*sig))+a+bx*(xx-11)+by*(yy-11)
            return ((prediction-stamp)/noise).ravel()
        initial=[max(float(stamp.sum()-np.median(stamp)*stamp.size),.01),11,11,2,float(np.median(stamp)),0,0]
        fit=least_squares(fun,initial,bounds=([0,7,7,.5,-np.inf,-np.inf,-np.inf],[np.inf,15,15,5,np.inf,np.inf,np.inf]),max_nfev=120)
        p=fit.x
        try:
            covariance=np.linalg.pinv(fit.jac.T@fit.jac)
            snr=float(p[0]/np.sqrt(covariance[0,0])) if covariance[0,0]>0 else 0.
        except np.linalg.LinAlgError: snr=0.
        result.append(dict(flux=float(p[0]),x=float(p[1]),y=float(p[2]),fwhm=float(2.35482*p[3]*PIXEL),snr=snr))
    return result


def injection_controls(row,smooth,annulus,center,legacy,wcs):
    """Static-coadd, flat-AB Gaussian injections; not survey completeness.

    This measures only the matched-filter detection step, not PSF subtraction,
    image-level masking, or the artifact classifier. Use known injection flux
    changes so existing peaks cannot masquerade as recovered injected sources.
    """
    result=[]
    scales=np.array([robust_scale(z[annulus]) for z in smooth])
    medians=np.array([np.median(z[annulus]) for z in smooth])
    fwhm=np.array([np.nanmedian(legacy['psfsize_'+b]) for b in 'grz'])
    sigmas=np.sqrt((fwhm/(2.35482*PIXEL))**2+2**2)
    for sep in [4,8,16,24]:
        for angle in [22.5,112.5,202.5,292.5]:
            x=center[0]+sep/PIXEL*np.cos(np.deg2rad(angle))
            y=center[1]+sep/PIXEL*np.sin(np.deg2rad(angle))
            ra,dec=wcs.pixel_to_world_values(x,y)
            distance,_=nearest(ra,dec,legacy)
            if distance<2: continue
            ix,iy=int(round(x)),int(round(y))
            yy,xx=np.indices((13,13));dx=xx+ix-6-x;dy=yy+iy-6-y
            original=smooth[:,iy-6:iy+7,ix-6:ix+7]
            for mag in [17,19,21,23]:
                flux=10**((22.5-mag)/2.5)
                signal=np.array([flux/(2*np.pi*s*s)*np.exp(-(dx*dx+dy*dy)/(2*s*s)) for s in sigmas])
                combined=np.sort((original+signal-medians[:,None,None])/scales[:,None,None],axis=0)[1]
                yymax,xxmax=np.unravel_index(np.argmax(combined),combined.shape)
                added_snr=float(np.sort(signal[:,yymax,xxmax]/scales)[1])
                recovered=combined[yymax,xxmax]>=6 and added_snr>=6 and np.hypot(xxmax+ix-6-x,yymax+iy-6-y)*PIXEL<1
                result.append(dict(field_gaia_id=str(row.source_id),separation_arcsec=sep,angle_deg=angle,
                                   flat_ab_mag=mag,recovered=bool(recovered),added_snr_second_band=added_snr,
                                   psf_fwhm_grz_arcsec=fwhm.tolist()))
    return result


def make_plot(sid, image, residual, records, center):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    for band in range(3):
        for row, cube in enumerate([image,residual]):
            scale = robust_scale(cube[band])
            axes[row,band].imshow(np.arcsinh(cube[band]/(3*scale)), origin='lower', cmap='gray',vmin=-2,vmax=5)
            axes[row,band].set_title(f'{"grz"[band]} {"coadd" if row==0 else "model residual"}')
            axes[row,band].plot(*center,'+',color='cyan')
            for rec in records:
                color = 'lime' if rec['status'].startswith('catalogued') else 'orange'
                axes[row,band].plot(rec['x'],rec['y'],'o',mfc='none',mec=color,ms=8)
    fig.suptitle(f'Gaia {sid}: green=catalog match, orange=unmatched peak (not a star claim)')
    fig.savefig(PLOTS/f'{sid}.png',dpi=110)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--expanded',action='store_true',help='64 stars in four magnitude strata, ordered by SHA256(source_id)')
    parser.add_argument('--run-name',default='pilot12_v2')
    args=parser.parse_args()
    global OUT,PLOTS
    OUT=OUT/args.run_name;PLOTS=PLOTS/args.run_name
    for directory in [RAW,OUT,PLOTS]: directory.mkdir(parents=True,exist_ok=True)
    selection_query=EXPANDED_QUERY if args.expanded else PILOT_QUERY
    selection_pool,meta=query('expanded_pool_gaia.csv' if args.expanded else 'pilot_gaia.csv',TAP,selection_query)
    if args.expanded:
        pilot_parts=[]
        for lo,hi in [(10,11),(11,12),(12,13.5),(13.5,15.000001)]:
            subset=selection_pool[(selection_pool.phot_g_mean_mag>=lo)&(selection_pool.phot_g_mean_mag<hi)].copy()
            subset['selection_hash']=subset.source_id.map(lambda s:sha256(str(s).encode()).hexdigest())
            subset=subset.sort_values('selection_hash').head(16)
            if len(subset)<16: raise ValueError(f'Insufficient fixed-stratum targets {lo} to {hi}: {len(subset)}')
            pilot_parts.append(subset)
        pilot=pd.concat(pilot_parts,ignore_index=True)
    else: pilot=selection_pool
    pilot.to_csv(OUT/'pilot.csv',index=False)
    records, summaries, failures, injections = [], [], [], []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [(row.source_id,pool.submit(load_field,row)) for row in pilot.itertuples(index=False)]
        for sid,future in futures:
            try:
                rec,summ,inj = search_field(future.result())
                records.extend(rec);summaries.append(summ);injections.extend(inj)
                print(sid,'analyzed',len(rec),'detections',flush=True)
            except Exception as error:
                failures.append(dict(gaia_id=str(sid),error=repr(error)))
                print(sid,'FAILED',repr(error),flush=True)
    pd.DataFrame(records).to_csv(OUT/'detections.csv',index=False)
    pd.DataFrame(injections).to_csv(OUT/'injections.csv',index=False)
    pd.DataFrame([{k:v for k,v in s.items() if k!='files'} for s in summaries]).to_csv(OUT/'fields.csv',index=False)
    manifest=[]
    for path in sorted(RAW.rglob('*.json')):
        info=json.loads(path.read_text()); info['local_path']=str(path.with_suffix(''))
        manifest.append(info)
    (OUT/'raw_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    result=dict(method_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),pilot_query=selection_query,
                selection_pool_count=len(selection_pool),selection='SHA256(source_id),16 per magnitude bin' if args.expanded else 'brightest12',
                pilot_count=len(pilot),analyzed_fields=len(summaries),failures=failures,
                detection_count=len(records),counts=pd.Series([r['status'] for r in records]).value_counts().to_dict(),
                controls=dict(catalogued_psf_detections=sum(s['matched_psf_positive_controls'] for s in summaries)),
                fields=summaries)
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='fields'},indent=2),flush=True)


if __name__=='__main__': main()
