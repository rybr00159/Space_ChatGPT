"""Native-pixel falsification of a selected isolated dip in TIC353933809.

No recurrence period, discovery claim, calibrated FAP or PRF precision assumed.
"""
from pathlib import Path
import sys,json,csv,hashlib,io,zipfile,argparse
from urllib.parse import urlencode
import numpy as np
import requests
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import multisector_followup as mf
from broad_vet import integrated_trapezoid
from pixel_vet_49617264 import gaussfit
import matplotlib.pyplot as plt
TIC=353933809;SECTOR=86;EPOCH=3637.2184787932
OUT=mf.OUT/f'tic{TIC}'/'isolated_s86'
CACHE=mf.CACHE/'followup'/f'tic{TIC}'/'pixels_isolated'
SOURCE=Path(__file__).read_bytes();SHA=hashlib.sha256(SOURCE).hexdigest()

def retrieve():
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urlencode(dict(ra=271.832661386,dec=54.5919743002,x=9,y=9,units='px',sector=SECTOR))
    archive=CACHE/'sector86_9x9_response.zip';manifest=OUT/'retrieval_manifest.json'
    if archive.exists():
        m=json.loads(manifest.read_text());assert m['url']==url and m['status']==200 and m['sha256']==mf.digest(archive)
    else:
        response=requests.get(url,timeout=90);archive.write_bytes(response.content)
        mf.save(manifest,dict(url=url,status=response.status_code,path=str(archive),sha256=mf.digest(archive),bytes=archive.stat().st_size,checked_utc=mf.utc()));response.raise_for_status()
    with zipfile.ZipFile(archive)as z:
        paths=[]
        for name in z.namelist():
            if name.endswith('.fits'):
                path=CACHE/Path(name).name;body=z.read(name)
                if path.exists():assert path.read_bytes()==body,'Cached raw FITS changed'
                else:path.write_bytes(body)
                paths.append(path)
    if len(paths)!=1:raise ValueError('Exactly one native sector expected')
    return paths[0],url

def shape_model(x,pars,degree,exposure):
    offset,duration,ingress,depth=pars[:4];basis=np.column_stack([x**k for k in range(degree+1)])
    q=integrated_trapezoid(x-offset,duration,ingress*duration,exposure)
    return basis@pars[4:]-depth*q

def fit_shape(t,f,error,half,degree,exposure,bootstrap=0):
    use=abs(t-EPOCH)<half;x=t[use]-EPOCH;flux=f[use];err=error[use];norm=float(np.median(flux));y=flux/norm;err=err/norm
    good=np.isfinite(y)&np.isfinite(err)&(err>0);x,y,err=x[good],y[good],err[good]
    lower=np.r_[-.3,.03,.025,-.03,[.5]+[-1.]*degree];upper=np.r_[.3,1.0,.499,.03,[1.5]+[1.]*degree]
    results=[]
    for duration in [.1,.3,.5,.8]:
        for offset in [-.04,.06]:
            start=np.r_[offset,duration,.3,.001,1.,[0.]*degree]
            fit=least_squares(lambda pars:(shape_model(x,pars,degree,exposure)-y)/err,start,bounds=(lower,upper),max_nfev=450)
            results.append(fit)
    best=min(results,key=lambda fit:sum(fit.fun**2));model=shape_model(x,best.x,degree,exposure);res=y-model
    dof=len(y)-len(best.x);chi2=float(sum(best.fun**2));cov=np.linalg.pinv(best.jac.T@best.jac)*max(1.,chi2/max(1,dof));se=np.sqrt(np.maximum(0,np.diag(cov)))
    poly=np.column_stack([x**k for k in range(degree+1)]);beta=np.linalg.lstsq(poly/err[:,None],y/err,rcond=None)[0];chi2flat=float(sum(((poly@beta-y)/err)**2))
    q=integrated_trapezoid(x-best.x[0],best.x[1],best.x[2]*best.x[1],exposure)
    out=dict(half_window_days=half,degree=degree,n_cadences=len(y),left_oot=int(((q==0)&(x<best.x[0])).sum()),right_oot=int(((q==0)&(x>best.x[0])).sum()),midtime_btjd=EPOCH+float(best.x[0]),midtime_error_days=float(se[0]),duration_days=float(best.x[1]),duration_error_days=float(se[1]),ingress_fraction=float(best.x[2]),depth=float(best.x[3]),depth_error=float(se[3]),depth_diagonal_snr=float(best.x[3]/se[3]),baseline_coefficients=best.x[4:],chi2_diagonal=chi2,dof=dof,delta_chi2_flat=chi2flat-chi2,fit_parameters=best.x,bound_hit=bool(np.any(best.x-lower<1e-4)|np.any(upper-best.x<1e-4)),warning='Selected single-event integrated trapezoid, no limb darkening. Covariance is conditional and incomplete for correlated noise. Both dip and brightening depth signs allowed. No planetary interpretation.')
    if bootstrap:
        rng=np.random.default_rng(TIC);draws=[];block=30
        for _ in range(bootstrap):
            starts=rng.integers(0,len(res),(len(res)+block-1)//block);indices=np.concatenate([(s+np.arange(block))%len(res)for s in starts])[:len(res)]
            fake=model+res[indices];fit=least_squares(lambda pars:(shape_model(x,pars,degree,exposure)-fake)/err,best.x,bounds=(lower,upper),max_nfev=250)
            draws.append(fit.x[:4])
        out['residual_block_bootstrap']=dict(n_draws=bootstrap,block_cadences=block,block_days=block*exposure,parameters=['midtime_offset_days','duration_days','ingress_fraction','depth'],q16_q50_q84=np.percentile(draws,[16,50,84],axis=0),warning='Conditional synchronized residual-block draws of one selected event; blocks cannot capture all longer correlations/systematics or survey trials.')
        np.savez_compressed(OUT/'shape_bootstrap.npz',draws=np.array(draws))
    return out,dict(time=x+EPOCH,flux=y,error=err,model=model,baseline=poly@best.x[4:])

def difference(t,cube,errors,half,degree,template,exposure,mask):
    use=abs(t-EPOCH)<half;x=t[use]-EPOCH;q=integrated_trapezoid(t[use]-template['midtime_btjd'],template['duration_days'],template['duration_days']*template['ingress_fraction'],exposure)
    a=np.column_stack([x**k for k in range(degree+1)]+[-q]);inv=np.linalg.pinv(a);f=cube[use].reshape(len(x),-1);coef=inv@f;res=f-a@coef
    scatter=np.maximum(np.std(res,axis=0,ddof=a.shape[1]),np.sqrt(np.mean(errors[use].reshape(len(x),-1)**2,axis=0)));projection=inv[-1];de=scatter*np.linalg.norm(projection)
    image=coef[-1].reshape(cube.shape[1:]);de=de.reshape(image.shape)
    aperture_residual=res[:,mask.ravel()].sum(axis=1);sigma=max(np.std(aperture_residual,ddof=a.shape[1]),np.sqrt(np.mean(np.sum(errors[use][:,mask]**2,axis=1))))
    amplitude=float(image[mask].sum());error=float(sigma*np.linalg.norm(projection));baseline=float(np.median((a[:,:-1]@coef[:-1])[:,mask.ravel()].sum(axis=1)))
    return dict(half_window_days=half,degree=degree,n_inside=int((q>.05).sum()),n_oot=int((q==0).sum()),condition=float(np.linalg.cond(a)),aperture_difference_e_s=amplitude,aperture_error_e_s=error,aperture_snr=amplitude/error,aperture_depth=amplitude/baseline,aperture_depth_error=error/baseline,positive_sum=float(np.maximum(image,0).sum()),negative_abs_sum=float(np.maximum(-image,0).sum())),image,de,res,projection

def redistribution(background):
    """Conditional empirical image-scale versus image-gradient diagnostic."""
    folder=OUT/background;data=np.load(folder/'diagnostic_arrays.npz');source=json.loads((folder/'isolated_pixel_validation.json').read_text());ref=data['reference'];target=data['target_xy'];yy,xx=np.indices(ref.shape);roi=(xx-target[0])**2+(yy-target[1])**2<3.5**2
    gy,gx=np.gradient(ref);basis=[ref,gx,gy,np.ones_like(ref),xx-target[0],yy-target[1]];design=np.column_stack([v[roi]for v in basis]);results=[];plot=[]
    for half in [.8,1.2]:
        image=data[f'w{half:g}_difference'];error=data[f'w{half:g}_error'];y=image[roi];e=error[roi];a=design/e[:,None];beta=np.linalg.lstsq(a,y/e,rcond=None)[0];res=y-design@beta;chi2=float(sum((res/e)**2));dof=len(y)-len(beta);cov=np.linalg.pinv(a.T@a)*max(1.,chi2/dof);se=np.sqrt(np.diag(cov));compare={}
        for name,columns in [('no_gradients',[0,3,4,5]),('no_flux_scaling',[1,2,3,4,5])]:
            reduced=design[:,columns];coef=np.linalg.lstsq(reduced/e[:,None],y/e,rcond=None)[0];compare[name]=float(sum(((y-reduced@coef)/e)**2))
        results.append(dict(half_window_days=half,parameters=['reference_fractional_decrement','x_gradient_coefficient_pixels','y_gradient_coefficient_pixels','constant_background_e_s','x_background_gradient','y_background_gradient'],coefficients=beta,errors_scaled_diagonal=se,chi2_diagonal=chi2,dof=dof,reduced_model_chi2=compare,warning='Selected difference map and empirical finite-difference gradients; these are explanatory fits, not source identification or a calibrated null test. Gradient templates approximate pointing, focus and PRF errors but may also absorb a blended real event. Pixel/background covariance omitted.'))
        model=sum(b*v for b,v in zip(beta,basis));plot.append((half,image,model,image-model));print('REDISTRIBUTION',background,half,'depthppm',beta[0]*1e6,'errorppm',se[0]*1e6,'grad',beta[1:3],'chi2',chi2,dof,compare,flush=True)
    mf.save(folder/'difference_redistribution_test.json',dict(driver_source_sha256=SHA,input_map=dict(path=str((folder/'diagnostic_arrays.npz').relative_to(ROOT)),sha256=mf.digest(folder/'diagnostic_arrays.npz')),reference_policy='Empirical out-of-event image from same selected native cube.',tests=results))
    snapshot=folder/'method_snapshots'/(SHA+'.py');snapshot.parent.mkdir(exist_ok=True)
    if not snapshot.exists():snapshot.write_bytes(SOURCE)
    fig,axes=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    for row,(half,image,model,residual)in zip(axes,plot):
        limit=np.percentile(abs(image),98)
        for ax,img,title in zip(row,[image,model,residual],['Selected difference','Reference scale + image gradients + plane','Residual']):
            im=ax.imshow(img,origin='lower',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.plot(*target,'r+',ms=10);ax.set_title(f'{title},w{half}d',fontsize=8);fig.colorbar(im,ax=ax,shrink=.7)
    fig.suptitle('Flux-versus-redistribution diagnostic; no calibrated significance or source claim');fig.savefig(folder/'difference_redistribution_test.png',dpi=140);plt.close(fig)

def main(draws=100,background='none'):
    global OUT
    path,url=retrieve();old=OUT/'native_manifest.json';oldraw=json.loads(old.read_text())if old.exists()else None
    if oldraw and oldraw['sha256']!=mf.digest(path):raise ValueError('Raw native pixel SHA mismatch')
    with fits.open(path,memmap=False)as hd:
        hd.verify('exception');h=hd[1].header;d=hd[1].data;checks=mf.verify_fits_checksums(path,hd);wcs=WCS(hd[2].header);exposure=float(h['TIMEDEL'])
        if h.get('TIMESYS')!='TDB' or abs(h.get('BJDREFI',0)+h.get('BJDREFF',0)-2457000)>1e-8:raise ValueError('Unrecognized native timing')
        good=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.all(np.isfinite(d['FLUX'])&np.isfinite(d['FLUX_ERR'])&(d['FLUX_ERR']>0),axis=(1,2))
        t=np.asarray(d['TIME'][good],float);cube=np.asarray(d['FLUX'][good],float);errors=np.asarray(d['FLUX_ERR'][good],float)
        original_time=np.asarray(d['TIME'],float);quality=np.asarray(d['QUALITY'],int);inside=abs(original_time-EPOCH)<.25
        flags={str(int(k)):int(v)for k,v in zip(*np.unique(quality,return_counts=True))};ffi=[str(v)for v in d['FFI_FILE']]
        science=hashlib.sha256()
        for key,typ in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:science.update(key.encode());science.update(np.asarray(d[key],dtype=typ).tobytes())
        creator=hd[0].header.get('CREATOR');col0=float(h['1CRV5P']);row0=float(h['2CRV5P'])
    manifest=dict(url=url,path=str(path),sha256=mf.digest(path),science_sha256=science.hexdigest(),bytes=path.stat().st_size,n_original=len(original_time),n_quality_finite=len(t),quality_counts=flags,near_selected_event_quality_counts={str(int(k)):int(v)for k,v in zip(*np.unique(quality[inside],return_counts=True))},exposure_days=exposure,creator=creator,first_detector_column=col0,first_detector_row=row0,fits_checksums=checks,checked_utc=mf.utc())
    mf.save(old,manifest);(OUT/'ffi_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    OUT=OUT/background;OUT.mkdir(parents=True,exist_ok=True)
    gp=ROOT/'data/multisector/catalogs/deep/20261016T054503162583Z/tic353933809_gaia_neighbors.csv';rows=list(csv.DictReader(gp.open()));year=float(Time(EPOCH+2457000,format='jd',scale='tdb').jyear)
    stars=[]
    for sid in ['2149774171239412352','2149774175535905792']:
        r=next(r for r in rows if r['source_id']==sid);ra,dec=float(r['ra']),float(r['dec']);dec+=float(r['pmdec'])*(year-2016.)/3600000.;ra+=float(r['pmra'])*(year-2016.)/(3600000.*np.cos(np.radians(float(r['dec']))));xy=wcs.all_world2pix([[ra,dec]],0)[0]
        stars.append(dict(source_id=sid,gaia2016_ra=float(r['ra']),gaia2016_dec=float(r['dec']),pmra_mas_yr=float(r['pmra']),pmdec_mas_yr=float(r['pmdec']),gmag=float(r['phot_g_mean_mag']),propagated_year=year,propagated_ra=ra,propagated_dec=dec,xy=xy))
    target=np.array(stars[0]['xy']);neighbor=np.array(stars[1]['xy']);yy,xx=np.indices(cube.shape[1:]);masks={f'target_r{r:g}':(xx-target[0])**2+(yy-target[1])**2<r**2 for r in [1.,1.5,2.5]};masks['neighbor_r1']=(xx-neighbor[0])**2+(yy-neighbor[1])**2<1.
    reference=np.median(cube[abs(t-EPOCH)>.4],axis=0);results={};arrays={};best=None
    blank=((xx-target[0])**2+(yy-target[1])**2>3.1**2)&(reference<np.percentile(reference,45))
    if blank.sum()<12:raise ValueError('Insufficient uncontaminated peripheral sky pixels')
    background_info=dict(policy=background,n_blank_pixels=int(blank.sum()),blank_mask=blank,selection='Radius>3.1pixels from Gaia primary and below45th percentile of native out-of-event reference image.',warning='Selected blank pixels may include constant faint sources and wings. Correlated sky uncertainty and large-scale spacecraft background remain incompletely modeled.')
    if background=='median':
        sky=np.median(cube[:,blank],axis=1);sky_variance=np.pi/2*np.mean(errors[:,blank]**2,axis=1)/blank.sum();cube=cube-sky[:,None,None];errors=np.sqrt(errors**2+sky_variance[:,None,None]);arrays.update(sky=sky,sky_variance=sky_variance)
    elif background=='rim_plane':
        basis=np.column_stack([np.ones(blank.sum()),xx[blank]-target[0],yy[blank]-target[1]]);inv=np.linalg.pinv(basis);allbasis=np.stack([np.ones_like(xx),xx-target[0],yy-target[1]],axis=-1);coeff=cube[:,blank]@inv.T;sky=np.einsum('tk,xyk->txy',coeff,allbasis);cube=cube-sky
        variances=errors[:,blank]**2;prediction=np.einsum('ki,xyk->xyi',inv,allbasis);sky_variance=np.einsum('ti,xyi->txy',variances,prediction**2);errors=np.sqrt(errors**2+sky_variance);arrays.update(sky_plane_coefficients=coeff,sky_at_primary=coeff[:,0],sky_plane_variance=sky_variance)
    elif background!='none':raise ValueError('Unknown background policy')
    reference=np.median(cube[abs(t-EPOCH)>.4],axis=0)
    for label,mask in masks.items():
        f=cube[:,mask].sum(axis=1);e=np.sqrt((errors[:,mask]**2).sum(axis=1));results[label]=[]
        for half in [.8,1.2]:
            fit,arr=fit_shape(t,f,e,half,2,exposure,bootstrap=draws if label=='target_r1.5'and half==1.2 else 0);results[label].append(fit)
            for key,v in arr.items():arrays[f'{label}_w{half:g}_{key}']=v
            if label=='target_r1.5'and half==1.2:best=fit
            print('NATIVE_SHAPE',label,half,{k:fit[k]for k in ['midtime_btjd','duration_days','depth','depth_error','depth_diagonal_snr','bound_hit']},flush=True)
    modes=[];maps={};roi=(xx-target[0])**2+(yy-target[1])**2<3.5**2
    for half in [.8,1.2]:
        stats,image,de,res,projection=difference(t,cube,errors,half,2,best,exposure,masks['target_r1.5']);label=f'w{half:g}';maps[label+'_difference']=image;maps[label+'_error']=de
        pos,pars=gaussfit(image,de,target,roi);neg,_=gaussfit(-image,de,target,roi)
        fixed={name:gaussfit(image,de,xy,roi,fixed_center=True)[0]for name,xy in [('target',target),('neighbor',neighbor)]}
        stats.update(positive_gaussian=pos,negative_gaussian=neg,fixed_positions=fixed,target_offset_pixels=float(np.hypot(pos['x']-target[0],pos['y']-target[1])),neighbor_offset_pixels=float(np.hypot(pos['x']-neighbor[0],pos['y']-neighbor[1])))
        stats['position_assessment']='diagnostic_only_no_reliable_localization_claim'
        if stats['aperture_snr']>=7 and not pos['width_bound_hit']:
            rng=np.random.default_rng(353933809);centers=[]
            for _ in range(draws):
                starts=rng.integers(0,len(res),(len(res)+29)//30);inds=np.concatenate([(s+np.arange(30))%len(res)for s in starts])[:len(res)];fake=image+(projection@res[inds]).reshape(image.shape)
                try:
                    fit,_=gaussfit(fake,de,[pos['x'],pos['y']],roi);centers.append([fit['x'],fit['y'],fit['sigma_x'],fit['sigma_y'],fit['chi2_diagonal']])
                except Exception:continue
            stats['diagnostic_block_bootstrap']=dict(n=len(centers),block_cadences=30,q16_q50_q84=np.percentile(centers,[16,50,84],axis=0),parameters=['x','y','sigma_x','sigma_y','chi2_diagonal'],warning='Not an official PRF; Gaussian mismatch, source blending, WCS and long instrumental covariance omitted. No confidence-level host assignment.')
            maps[label+'_bootstrap_xywidth']=np.array(centers)
        modes.append(stats);print('NATIVE_DIFFERENCE',half,stats,flush=True)
    output=dict(tic=TIC,sector=SECTOR,driver_source_sha256=SHA,checked_utc=mf.utc(),isolated_event_initial_epoch_btjd=EPOCH,source=manifest,background_treatment=background_info,gaia_input=dict(path=str(gp.relative_to(ROOT)),sha256=mf.digest(gp),proper_motion_policy='Linear proper motion from Gaia J2016 to event barycentric jyear; annual parallax and orbital photocentre displacement omitted, a few mas versus21arcsec native pixels.'),stars=stars,aperture_shapes=results,pixel_template=best,difference_modes=modes,assessment='isolated_selected_event_under_falsification_not_a_planet_candidate',limitations=['No recurrence period assumed.','Raw pixels share photons with QLP; not independent detection significance.','Single event, selected template, aperture, time, baseline and search trials are uncalibrated.','No stellar surface, eclipsing-binary or physical planetary transit model is established.','Pixel-diagonal errors and conditional block bootstraps do not capture full spacecraft, PRF or background covariance.','Nearby companion is separated by only about0.6native pixels; a smooth empirical Gaussian cannot by itself guarantee source attribution.'])
    mf.save(OUT/'isolated_pixel_validation.json',output);np.savez_compressed(OUT/'diagnostic_arrays.npz',**arrays,**maps,reference=reference,target_xy=target,neighbor_xy=neighbor)
    snap=OUT/'method_snapshots'/(SHA+'.py');snap.parent.mkdir(exist_ok=True)
    if not snap.exists():snap.write_bytes(SOURCE)
    fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    for ax,mode in zip(axes,modes):
        label=f'w{mode["half_window_days"]:g}'
        for a,img,title in zip(ax,[reference,maps[label+'_difference'],maps[label+'_difference']/maps[label+'_error']],['Native median flux','Isolated-event difference','Conditional diagonal S/N']):
            im=a.imshow(img,origin='lower');a.plot(*target,'r+',ms=12);a.plot(*neighbor,'cx',ms=10);a.set_title(f'{title}, halfwindow{mode["half_window_days"]}d',fontsize=9);fig.colorbar(im,ax=a,shrink=.7)
    fig.suptitle('TIC353933809 S86 isolated event: target(red), Gaia neighbor(cyan); no secure localization');fig.savefig(OUT/'native_difference_diagnostics.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,(label,tests)in zip(axes.ravel(),results.items()):
        for fit in tests:
            pre=f'{label}_w{fit["half_window_days"]:g}_';x=arrays[pre+'time']-EPOCH;y=arrays[pre+'flux'];b=arrays[pre+'baseline'];ax.plot(x,(y-b)*1000,'.',ms=1,alpha=.3,label=f'halfwindow{fit["half_window_days"]}d');ax.plot(x,(arrays[pre+'model']-b)*1000,lw=1.5)
        ax.set(xlabel='Days from initial selected epoch',ylabel='Local residual(ppt)',title=label);ax.legend(fontsize=8)
    fig.suptitle('Native aperture event shape, two quadratic baseline windows; selected single event');fig.savefig(OUT/'native_aperture_shapes.png',dpi=150);plt.close(fig)
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--bootstrap',type=int,default=100);parser.add_argument('--background',choices=['none','median','rim_plane'],default='none');parser.add_argument('--decompose-only',action='store_true');args=parser.parse_args()
    try:
        if args.decompose_only:redistribution(args.background)
        else:main(args.bootstrap,args.background)
    except Exception as exc:
        import traceback
        mf.save(OUT/'analysis_failure.json',dict(error=repr(exc),traceback=traceback.format_exc(),driver_source_sha256=SHA,checked_utc=mf.utc()));raise
