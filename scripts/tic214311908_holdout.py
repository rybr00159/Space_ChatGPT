"""Prospective fixed-ephemeris test on previously uninspected native FFIs.

QLP sectors103/104 were excluded for documented timing issues. Native FFI
TDB timestamps are used directly, with original FFI identifiers retained. The
frozen discovery ephemeris and analysis policy precede retrieval of holdout flux.
"""
from pathlib import Path
import os,sys,json,hashlib,io,zipfile
from urllib.parse import urlencode
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np,requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import multisector_followup as f
from archive import verify_fits_checksums
OUT=ROOT/'candidates/tic214311908/holdout';CACHE=Path('/workspace/tess-research-data/tic214311908/holdout')
RA=260.96969903043;DEC=-51.32383313824

def policy():
 OUT.mkdir(parents=True,exist_ok=True);path=OUT/'preregistration.json'
 if path.exists():
  d=json.loads(path.read_text())
  if d['script_sha256']!=f.digest(__file__):raise ValueError('Holdout driver changed after preregistration')
  return d
 pool=ROOT/'candidates/tic214311908/pooled_w0.75.json';p=json.loads(pool.read_text())['peaks'][0]
 d=dict(registered_utc=f.utc(),script_sha256=f.digest(__file__),discovery_source=dict(path=str(pool.relative_to(ROOT)),sha256=f.digest(pool)),discovery_sectors=[12,39,66,93],uninspected_sectors=[103,104],period_days=p['period'],epoch_btjd=p['t0_btjd'],duration_days=p['duration_days'],aperture_radii_pixels=[1.,1.5,2.],primary_radius_pixels=1.5,local_event_windows_days=[p['duration_days'],.07],background='median of lower25percentile reference-image pixels farther than3pixels from propagated target; constant per cadence, same mask across sector',time_policy='Native FFI TIME/BJDREFI/BJDREFF/TIMESYS verified; no QLP photometry or timestamp applied.',limitations='Same fixed ephemeris; apertures and two windows are correlated diagnostics, not independent tests. Real blend/stellar classification and calibrated red-noise FAP remain unresolved. All predicted coverage retained; a gap is not a null.')
 f.save(path,d);return d

def retrieve(sector):
 folder=CACHE/f'sector{sector}';folder.mkdir(parents=True,exist_ok=True);man=folder/'download_manifest.json'
 if man.exists():
  d=json.loads(man.read_text());p=Path(d['fits_path'])
  if f.digest(p)!=d['fits_sha256']:raise ValueError('Native holdout original changed')
  return p,d
 url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urlencode(dict(ra=RA,dec=DEC,x=9,y=9,units='px',sector=sector));r=requests.get(url,timeout=(20,180));path=folder/'original_response.zip';path.write_bytes(r.content);r.raise_for_status()
 with zipfile.ZipFile(io.BytesIO(r.content))as z:
  names=[n for n in z.namelist()if n.endswith('.fits')]
  if len(names)!=1:raise ValueError('Unexpected number of native cutouts')
  p=folder/Path(names[0]).name;p.write_bytes(z.read(names[0]))
 d=dict(url=url,http_status=r.status_code,zip_path=str(path),zip_sha256=f.digest(path),fits_path=str(p),fits_sha256=f.digest(p),retrieved_utc=f.utc());f.save(man,d);return p,d

def analyze(sector,prereg):
 path,source=retrieve(sector)
 with fits.open(path,memmap=False)as h:
  h.verify('exception');checks=verify_fits_checksums(path,h);head=h[1].header;d=h[1].data
  if head.get('TIMESYS')!='TDB'or head.get('BJDREFI',0)+head.get('BJDREFF',0)!=2457000:raise ValueError('Holdout timing is not explicit BTJD TDB')
  wcs=WCS(h[2].header);good=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.all(np.isfinite(d['FLUX']),axis=(1,2));t=np.asarray(d['TIME'][good],float);cube=np.asarray(d['FLUX'][good],float);error=np.asarray(d['FLUX_ERR'][good],float);ffi=[str(x)for x in d['FFI_FILE']]
  native_time=np.asarray(d['TIME'],float);native_quality=np.asarray(d['QUALITY'],int);exposure=float(head['TIMEDEL'])
 gaia=json.loads((ROOT/'candidates/tic214311908/catalogs/catalog_check.json').read_text())['target'];year=Time(float(np.median(t))+2457000,format='jd',scale='tdb').jyear;dt=year-2016.;ra=float(gaia['ra'])+float(gaia['pmra'])*dt/(3600000*np.cos(np.deg2rad(float(gaia['dec']))));dec=float(gaia['dec'])+float(gaia['pmdec'])*dt/3600000;xy=wcs.all_world2pix([[ra,dec]],0)[0];ref=np.median(cube,axis=0);yy,xx=np.indices(ref.shape);dist=np.hypot(xx-xy[0],yy-xy[1]);blank=(ref<np.percentile(ref,25))&(dist>3)
 if blank.sum()<3:raise ValueError('Too few predeclared sky pixels')
 sky=np.median(cube[:,blank],axis=1);skyerr=1.253*np.sqrt(np.sum(error[:,blank]**2,axis=1))/blank.sum();p,e0=prereg['period_days'],prereg['epoch_btjd'];rows=[];arrays=dict(t=t,sky=sky,reference=ref,target_xy=xy,native_time=native_time,native_quality=native_quality)
 for radius in prereg['aperture_radii_pixels']:
  aperture=dist<=radius;raw=np.sum(cube[:,aperture],axis=1)-aperture.sum()*sky;err=np.sqrt(np.sum(error[:,aperture]**2,axis=1)+(aperture.sum()*skyerr)**2);norm=np.median(raw);flux=raw/norm;err=err/norm;arrays.update({f'r{radius:g}_f':flux,f'r{radius:g}_e':err})
  for duration in prereg['local_event_windows_days']:
   stats,local=f.local_quadratic(t,flux,err,p,e0,duration);secondary=f.local_quadratic(t,flux,err,p,e0+p/2,duration)[0];rows.append(dict(aperture_radius_pixels=radius,n_aperture_pixels=int(aperture.sum()),duration_days=duration,primary=stats,secondary=secondary))
 coverage=[]
 for n in range(int(np.floor((t.min()-e0)/p)),int(np.ceil((t.max()-e0)/p))+1):
  center=e0+n*p
  if t.min()<=center<=t.max():coverage.append(dict(cycle=n,center_btjd=center,n_good_inside=int(np.sum(abs(t-center)<prereg['duration_days']/2)),nearest_good_cadence_days=float(np.min(abs(t-center)))))
 dest=OUT/f'sector{sector}';dest.mkdir(exist_ok=True);np.savez_compressed(dest/'native_aperture_arrays.npz',**arrays);(dest/'ffi_identifiers.txt').write_text('\n'.join(ffi)+'\n');f.save(dest/'assessment.json',dict(sector=sector,preregistration_sha256=f.digest(OUT/'preregistration.json'),source=source,fits_checksums=checks,timing={k:head.get(k)for k in ['TIMESYS','BJDREFI','BJDREFF','TIMEDEL','TSTART','TSTOP']},n_native=len(d),n_good=len(t),epoch_jyear=year,target_xy=xy,astrometric_policy='Linear Gaia2016 proper motion; parallax/orbital photocentre and WCS systematics omitted.',n_sky_pixels=int(blank.sum()),tests=rows,predicted_coverage=coverage,limitations=prereg['limitations']))
 for r in rows:print('HOLDOUT',sector,r['aperture_radius_pixels'],r['duration_days'],{k:r['primary'].get(k)for k in ['n_events','n_positive_events','depth','depth_err','snr']},flush=True)

if __name__=='__main__':
 pre=policy()
 for s in pre['uninspected_sectors']:analyze(s,pre)
