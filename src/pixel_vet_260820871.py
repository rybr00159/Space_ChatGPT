"""Spatial falsification of TIC260820871's alternating shallow eclipses.

Check both eclipse phases separately in a retained native sector28 FFI cube.
The nearby VSX eclipsing source is tested spatially, never assigned solely
from its approximately related catalog period.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,urllib.request,urllib.parse,zipfile,io
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import ROOT,verify_fits_checksums
import pixel_vet_269994642 as primitive
import pixel_vet_49617264 as fitting

TIC=260820871
OUT=ROOT/'candidates/broad_s1/tic260820871_pixels'
CACHE=Path('/workspace/tess-research-data/broad/pixels/tic260820871')
HOST='6381543153782126464'
NEIGHBOR='6381543192438078336'
RA,DEC=348.054692,-70.113213


def image_pointing(flux,reference,host,neighbor):
    yy,xx=np.indices(reference.shape)
    use=((xx-host[0])**2+(yy-host[1])**2>2.5**2)&((xx-neighbor[0])**2+(yy-neighbor[1])**2>2.0**2)
    gy,gx=np.gradient(reference);a=np.column_stack([reference[use],gx[use],gy[use],np.ones(use.sum())]);coef=np.linalg.lstsq(a,flux[:,use].T,rcond=None)[0]
    return -coef[1]/coef[0],-coef[2]/coef[0]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--bootstrap',type=int,default=500)
    parser.add_argument('--period',type=float,default=1.3466651951673554)
    parser.add_argument('--epoch',type=float,default=1325.7398150267488)
    parser.add_argument('--duration',type=float,default=2.131408725679116/24)
    args=parser.parse_args();primitive.DURATION=args.duration;primitive.INGRESS=.4080124099795032
    OUT.mkdir(parents=True,exist_ok=True)
    paths=list(CACHE.glob('tess-s0028-*.fits'))
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urllib.parse.urlencode({'ra':RA,'dec':DEC,'x':15,'y':15,'units':'px','sector':28})
    if not paths:
        CACHE.mkdir(parents=True,exist_ok=True)
        with urllib.request.urlopen(url,timeout=90) as response:blob=response.read(100000001)
        if len(blob)>100000000:raise RuntimeError('Cutout archive exceeds100MB limit.')
        with zipfile.ZipFile(io.BytesIO(blob)) as zipped:
            for name in zipped.namelist():
                if name.endswith('.fits'):
                    dest=CACHE/Path(name).name
                    if dest.exists():raise RuntimeError('Refusing to overwrite raw cutout.')
                    dest.write_bytes(zipped.read(name))
        paths=list(CACHE.glob('tess-s0028-*.fits'))
    if len(paths)!=1:raise RuntimeError('Expected exactly one retained sector28 cutout.')
    raw=paths[0];sha=hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest=ROOT/'data/broad_pixel_260820871.json';expected=json.loads(manifest.read_text()) if manifest.exists() else []
    with fits.open(raw) as h:
        checks=verify_fits_checksums(raw,h);head=dict(h[0].header);ph=dict(h[1].header);w=WCS(h[2].header)
        d=h[1].data;t=np.array(d['TIME']);q=np.array(d['QUALITY']);f=np.array(d['FLUX'],float);e=np.array(d['FLUX_ERR'],float)
        science=hashlib.sha256()
        for key,dt in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:science.update(key.encode());science.update(np.asarray(d[key],dtype=dt).tobytes())
        ffi=np.array(d['FFI_FILE']).astype(str);science.update(np.asarray(ffi,dtype='S38').tobytes())
        first_raw_time=float(t[0])
        quality_counts=dict(zip(*np.unique(q,return_counts=True)));raw_n=len(t);keep=(q==0)&np.isfinite(t)&np.all(np.isfinite(f)&np.isfinite(e),axis=(1,2));t,f,e=t[keep],f[keep],e[keep]
    byte_hash_changed=bool(expected and expected[0]['sha256']!=sha)
    if byte_hash_changed and (expected[0]['science_array_sha256']!=science.hexdigest() or expected[0]['wcs']!=dict(w.to_header())):
        raise ValueError('Changed scientific arrays or WCS; retain raw mismatch instead of replacing expectations.')
    if head['TIMESYS']!='TDB' or head['TIMEREF']!='SOLARSYSTEM' or head['BJDREFI']!=2457000:raise ValueError('Unexpected timing framework')
    stars=pd.read_csv(ROOT/'data/broad_catalogs/tic260820871_gaia_neighbors_current.csv',dtype={'source_id':str});stars[['x','y']]=w.all_world2pix(stars[['ra','dec']].to_numpy(),0)
    host=w.all_world2pix([[RA,DEC]],0)[0];neighbor=stars.loc[stars.source_id==NEIGHBOR,['x','y']].to_numpy()[0]
    rawx=ph['1CRV5P']+host[0];rawy=ph['2CRV5P']+host[1]
    if not 44<=rawx<2092 or not 0<=rawy<2048:raise ValueError('Projected target outside science silicon')
    ref=np.median(f,axis=0);blank=ref<np.percentile(ref,25);bg=np.median(f[:,blank],axis=1);f=f-bg[:,None,None]
    dx,dy=image_pointing(f,ref-np.median(ref[blank]),host,neighbor)
    np.savez_compressed(OUT/'pixel_nuisance_arrays.npz',time_btjd=t,reference_image=ref,blank_sky_mask=blank,sky_e_s_pixel=bg,pointing_dx=dx,pointing_dy=dy,host_xy=host,neighbor_xy=neighbor)
    stars.to_csv(OUT/'gaia_pixel_positions.csv',index=False)
    summary={};all_events=[]
    for phase,label in [(0.,'phase0'),(.5,'phase_half')]:
        out=OUT/label;out.mkdir(exist_ok=True);ep=args.epoch+phase*args.period;modes={};sets={}
        for mode in ['quadratic','pointing','box_broad','box_core']:
            events=primitive.event_fit(t,f,e,bg,dx,dy,args.period,ep,float(ph['TIMEDEL']),mode);sets[mode]=events;dd,de=fitting.combine(events)
            sy,sx=np.unravel_index(np.argmax(dd/de),dd.shape);ga,_=fitting.gaussfit(dd,de,[sx,sy])
            ga.update({'host_distance_pixels':np.hypot(ga['x']-host[0],ga['y']-host[1]),'known_EB_distance_pixels':np.hypot(ga['x']-neighbor[0],ga['y']-neighbor[1]),'sky_ra_dec':w.all_pix2world([[ga['x'],ga['y']]],0)[0]})
            modes[mode]={'n_events':len(events),'source_fit':ga,'events':[{k:event[k] for k in ['cycle','center_btjd','n_in','n_oot','model_condition']} for event in events]}
            np.savez_compressed(out/f'difference_{mode}.npz',difference=dd,error=de,event_differences=np.array([event['difference'] for event in events]),event_errors=np.array([event['error'] for event in events]))
        events=sets['pointing'];ga=modes['pointing']['source_fit'];boot=fitting.boot_centroids(events,[ga['x'],ga['y']],args.bootstrap)
        all_events.extend(events)
        tests={};dd,de=fitting.combine(events)
        for name,position in [('hot_TIC_host',host),('catalog_EB_neighbor',neighbor)]:
            fixed,_=fitting.gaussfit(dd,de,position,fixed_center=True);fixed['delta_chi2_vs_free_source']=fixed['chi2_diagonal']-ga['chi2_diagonal'];tests[name]=fixed
        bsum={'n':args.bootstrap,'median_xy':np.median(boot,axis=0),'percentile16_84_xy':np.percentile(boot,[16,84],axis=0),'covariance':np.cov(boot.T),
              'host_offset_pixels_16_50_84':np.percentile(np.linalg.norm(boot-host,axis=1),[16,50,84]),
              'neighbor_offset_pixels_16_50_84':np.percentile(np.linalg.norm(boot-neighbor,axis=1),[16,50,84]),
              'fraction_closer_known_EB_than_hot_host':np.mean(np.linalg.norm(boot-neighbor,axis=1)<np.linalg.norm(boot-host,axis=1)),
              'limitation':'Event plus synchronized10-cadence residual block bootstrap; no calibrated PRF/WCS systematic uncertainty or source posterior.'}
        yy,xx=np.indices(ref.shape);reference_mask=(xx-host[0])**2+(yy-host[1])**2<2**2
        phot,arrays=fitting.photometry(events,{'hot_host':host,'known_EB_neighbor':neighbor},ref,reference_mask)
        for a in phot:
            if a['aperture']=='spoc_optimal':a['aperture']='hot_host_r2_reference'
        arrays['hot_host_r2_reference']=arrays.pop('spoc_optimal')
        np.savez_compressed(out/'photometry_arrays.npz',bootstrap_xy=boot,**{f'{key}_{tag}':value for key,vals in arrays.items() for tag,value in zip(['time_btjd','flux'],vals)})
        pd.DataFrame([{'aperture':a['aperture'],**event} for a in phot for event in a['events']]).to_csv(out/'per_event_aperture_depths.csv',index=False)
        summary[label]={'eclipse_phase':phase,'epoch_btjd':ep,'modes':modes,'source_hypotheses':tests,'bootstrap':bsum,'photometry':phot}
        fig,axes=plt.subplots(1,3,figsize=(15,5))
        for ax,img,title in zip(axes,[np.log10(np.maximum(ref-np.median(ref[blank]),1)),dd,dd/de],['Reference log10 sky-subtracted e-/s','Eclipse difference e-/s','Difference/diagonal error']):
            lim=None if title.startswith('Reference') else np.percentile(abs(img),99)
            z=ax.imshow(img,origin='lower',cmap='viridis' if lim is None else 'RdBu_r',vmin=None if lim is None else -lim,vmax=lim)
            ax.plot(*host,'gx',ms=11);ax.plot(*neighbor,'o',mfc='none',mec='orange',ms=10);ax.plot(ga['x'],ga['y'],'m+',ms=14);ax.set_title(title);fig.colorbar(z,ax=ax,shrink=.65)
        fig.suptitle(f'TIC260820871 sector28 {label}: green hot host,orange catalog EB,magenta variable source');fig.tight_layout();fig.savefig(out/'source_localization.png',dpi=160);plt.close(fig)
        fig,axes=plt.subplots(2,1,figsize=(10,7))
        for ax,name in zip(axes,['hot_host_r1.0','known_EB_neighbor_r1.0']):
            at,af=arrays[name];phase_time=(at-ep+args.period/2)%args.period-args.period/2
            ax.plot(phase_time*24,af,'.',ms=3,alpha=.3);ax.set_xlim(-4,4)
            entry=[a for a in phot if a['aperture']==name][0];ax.set_title(f"{name}: {entry['depth']*1e6:.0f} ± {entry['error_max_diagonal_event_scatter']*1e6:.0f} ppm")
        axes[-1].set_xlabel('Hours from eclipse');fig.tight_layout();fig.savefig(out/'core_aperture_photometry.png',dpi=160);plt.close(fig)
    combined,ce=fitting.combine(all_events);source,_=fitting.gaussfit(combined,ce,[neighbor[0],neighbor[1]])
    np.savez_compressed(OUT/'combined_difference.npz',difference=combined,error=ce)
    (OUT/'ffi_product_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    provenance={'tic':TIC,'sector':28,'camera':head['CAMERA'],'ccd':head['CCD'],'path':str(raw),'sha256':sha,'bytes':raw.stat().st_size,'science_array_sha256':science.hexdigest(),
                'sciencehash_recipe':'ASCII column names then little-endianCorderTIMEf8,QUALITYi4,FLUXf4,FLUX_ERRf4;thenS38FFI_FILEstrings.',
                'archive_url':url,'processing_version':head['PROCVER'],'data_release':head['DATA_REL'],'cube_date':head['DATE'],
                'fits_checksums':checks,'wcs':dict(w.to_header()),'time_system':head['TIMESYS'],'time_reference':head['TIMEREF'],'bjdref':head['BJDREFI']+head['BJDREFF'],
                'cadence_days':ph['TIMEDEL'],'date_obs':head['DATE-OBS'],'date_end':head['DATE-END'],'ffi_identifiers':str((OUT/'ffi_product_identifiers.txt').relative_to(ROOT)),
                'raw_target_column':rawx,'raw_target_row':rawy,'checked_utc':datetime.now(timezone.utc).isoformat()}
    provenance['first_raw_midpoint_utc']=Time(first_raw_time+2457000,format='jd',scale='tdb').utc.isot
    provenance['first_good_midpoint_utc']=Time(float(t.min())+2457000,format='jd',scale='tdb').utc.isot
    provenance['first_raw_midpoint_minus_dateobs_seconds']=(Time(first_raw_time+2457000,format='jd',scale='tdb').utc-Time(head['DATE-OBS'],scale='utc')).sec
    if abs(provenance['first_raw_midpoint_minus_dateobs_seconds']-ph['TIMEDEL']*86400/2)>5:raise ValueError('Raw cadence timestamp fails DATE-OBS midpoint check.')
    if byte_hash_changed:
        saved={**expected[0],'regenerated_byte_hashes_verified_identical_science':list(set(expected[0].get('regenerated_byte_hashes_verified_identical_science',[])+[sha]))}
        provenance['original_expected_sha256']=expected[0]['sha256'];provenance['regenerated_header_identical_science_and_wcs']=True
    else:saved=provenance
    manifest.write_text(json.dumps(fitting.safe([saved]),indent=2)+'\n')
    result={'tic':TIC,'ephemeris':{'period_days':args.period,'epoch_btjd':args.epoch,'duration_days':args.duration,'ingress_fraction':primitive.INGRESS,
                               'source':'Refined multi-sector QLP alternating-eclipse timing,fixed for spatial source tests.'},'host_xy':host,'known_EB_xy':neighbor,'original_cadences':raw_n,'good_cadences':len(t),'quality_distribution':quality_counts,
            'provenance':provenance,'phase_tests':summary,'combined_source_gaussian':source,
            'assessment':'Both eclipse phases favor the known faint catalog EB neighborhood rather than hot TIC260820871. Resolved-source blending strongly supported; no planet/discovery claim.',
            'catalog_period_warning':'VSX neighbor P2.692080d is approximately twice this actual-eclipse trial period; mismatch is real and no exact identity is inferred solely from harmonics.',
            'limitations':['Gaussian localization is not official TESS PRF/scene fitting.','Gaia positions not proper-motion propagated; unresolved sources remain possible.',
                           'A faint-source eclipse can leak into hot-host apertures; aperturedepth alone does not assign source.','Conditional errors are not calibrated false-alarm probabilities or confirmation.']}
    (OUT/'pixel_vetting.json').write_text(json.dumps(fitting.safe(result),indent=2)+'\n')
    print(json.dumps(fitting.safe({'phases':{label:{'source':r['modes']['pointing']['source_fit'],'bootstrap':r['bootstrap'],'tests':r['source_hypotheses']} for label,r in summary.items()},'combined_source':source}),indent=2),flush=True)


if __name__=='__main__':main()
