"""Checkpointed falsification of new multi-sector transit leads.

These conditional diagnostics are not independent observations or calibrated
false-alarm probabilities. Old broad-survey scientific helpers are imported
read-only; their hashes, raw FITS bytes, and every attempted retrieval are saved.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, hashlib, json, os, sqlite3, traceback
import numpy as np
import requests
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares
from archive import verify_fits_checksums
from fast_screen import FastScreenConfig, preprocess_curve, event_stats, phase
from broad_vet import local_quadratic, fit_shape, centroid_diagnostic, variability_diagnostics
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path(os.environ.get('TESS_DATA_DIR','/workspace/tess-research-data'))/'multisector'
OUT=ROOT/'candidates/multisector'
DB=CACHE/'survey.sqlite'
VERSION='multisector-raw-alternative-detrending-falsification-v1'
LOADED_SOURCE=Path(__file__).read_bytes()

def utc(): return datetime.now(timezone.utc).isoformat()
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def clean(x):
    if isinstance(x, dict): return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [clean(v) for v in x]
    if isinstance(x,np.ndarray): return clean(x.tolist())
    if isinstance(x,np.generic): return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x): return None
    return x

def save(path,payload):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(clean(payload),indent=2,allow_nan=False)+'\n');tmp.replace(path)

def method():
    source=LOADED_SOURCE;sha=hashlib.sha256(source).hexdigest();snapshot=OUT/'method_snapshots'/(sha+'.py');snapshot.parent.mkdir(parents=True,exist_ok=True)
    if not snapshot.exists():snapshot.write_bytes(source)
    return dict(version=VERSION,source_sha256=sha,source_snapshot=str(snapshot.relative_to(ROOT)),frozen_helpers={n:digest(ROOT/'src'/n)for n in ['fast_screen.py','broad_vet.py','archive.py']})

def read_source(record):
    path=Path(record['path'])
    if digest(path)!=record['sha256']: raise ValueError('Raw SHA mismatch: '+str(path))
    with fits.open(path,memmap=False)as hd:
        hd.verify('exception');integrity=verify_fits_checksums(path,hd)
        h=hd[1].header;d=hd[1].data
        if h.get('TIMESYS',hd[0].header.get('TIMESYS'))!='TDB' or abs(h.get('BJDREFI',0)+h.get('BJDREFF',0)-2457000)>1e-8: raise ValueError('Timing is not explicitly BTJD TDB')
        corr='DET_FLUX'if'DET_FLUX'in d.names else'KSPSAP_FLUX';err=corr+'_ERR'
        aliases={'TIME':'TIME','QUALITY':'QUALITY','corrected':corr,'error':err,'raw':'SAP_FLUX'}
        for key,col in [('small',corr+'_SML'),('large',corr+'_LAG'),('centroid_x','SAP_X'),('centroid_y','SAP_Y'),('background','SAP_BKG')]:
            if col in d.names:aliases[key]=col
        a={key:np.asarray(d[col],float)for key,col in aliases.items()}
        metadata={k:h.get(k,hd[0].header.get(k))for k in ['TICID','SECTOR','TESSMAG','RADIUS','TEFF','CAMERA','CCD','TIMEDEL','PROCVER','DATE','RA_OBJ','DEC_OBJ']}
    if int(metadata['TICID'])!=int(record['tic']) or int(metadata['SECTOR'])!=int(record['sector']):raise ValueError('FITS target/sector identity differs from retrieval record')
    return a,dict(columns=aliases,metadata=metadata,fits_checksums=integrity)

def flatten(a,column,width):
    return preprocess_curve(a['TIME'],a[column],a['error'],a['QUALITY'],config=FastScreenConfig(width_days=width,bin_minutes=10.,edge_cut_days=.06))

def refine(a,peak):
    flat=flatten(a,'corrected',1.5);p=float(peak.get('period',peak.get('period_days')));duration=float(peak.get('duration_days',float(peak.get('duration_hours',2.))/24))
    span=max(.002*p,.6*duration*p/np.ptp(flat['t']))
    periods=np.linspace(max(.1,p-span),p+span,401)
    durations=np.unique(np.clip(np.r_[duration*np.array([.5,.75,1.,1.25,1.5]),.04,.08,.12],.015,min(.25*p,.35)))
    res=BoxLeastSquares(flat['t'],flat['f'],flat['e']).power(periods,durations,oversample=12,objective='snr');j=int(np.argmax(res.power))
    return dict(period=float(res.period[j]),epoch_btjd=float(res.transit_time[j]),duration_days=float(res.duration[j]),box_depth=float(res.depth[j]),diagonal_snr=float(res.depth_snr[j]),period_window=[float(periods[0]),float(periods[-1])],period_grid_count=len(periods),warning='Targeted local optimization; not a global false-alarm probability.')

def quality_context(a,p,epoch,duration):
    t=a['TIME'];good=np.isfinite(t);t=t[good];q=a['QUALITY'][good];order=np.argsort(t);t=t[order];q=q[order]
    inside=np.abs(phase(t,p,epoch))<duration/2;bad=q!=0
    dt=np.median(np.diff(t));near_bad=np.zeros(len(t),bool)
    bt=t[bad]
    if len(bt):
        j=np.searchsorted(bt,t);dist=np.minimum(abs(t-bt[np.clip(j,0,len(bt)-1)]),abs(t-bt[np.clip(j-1,0,len(bt)-1)]));near_bad=dist<max(.05,3*dt)
    return dict(n_inside_all_quality=int(inside.sum()),inside_bad_quality_fraction=float(bad[inside].mean())if inside.any()else None,outside_bad_quality_fraction=float(bad[~inside].mean())if(~inside).any()else None,inside_near_bad_fraction=float(near_bad[inside].mean())if inside.any()else None,note='Descriptive association with spacecraft flags; not a calibrated p-value.')

def strongest_secondary(flat,p,epoch,duration):
    outside=abs(phase(flat['t'],p,epoch))>1.5*duration
    if outside.sum()<30:return None
    durations=np.unique(np.clip([duration*.75,duration,duration*1.25],.015,.2*p))
    result=BoxLeastSquares(flat['t'][outside],flat['f'][outside],flat['e'][outside]).power([p],durations,oversample=15,objective='snr')
    return dict(epoch_btjd=float(result.transit_time[0]),duration_days=float(result.duration[0]),phase=float((result.transit_time[0]-epoch)%p/p),diagonal_snr=float(result.depth_snr[0]),warning='Secondary epoch selected after a phase search; conditional score requires trials correction.')

def vet_observation(record,peak,output=True):
    a,provenance=read_source(record);ephem=refine(a,peak);p,epoch,duration=ephem['period'],ephem['epoch_btjd'],ephem['duration_days']
    good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['error'])&(a['error']>0)
    t,e=a['TIME'][good],a['error'][good];tests={};intermediate={};local_data={}
    additional=strongest_secondary(flatten(a,'corrected',1.5),p,epoch,duration)
    for column in ['corrected','raw','small','large']:
        if column not in a:continue
        primary,loc=local_quadratic(t,a[column][good],e,p,epoch,duration);secondary,_=local_quadratic(t,a[column][good],e,p,epoch+p/2,duration)
        tested=dict(primary=primary,half_phase_secondary=secondary,alternative_windows=[],variability=variability_diagnostics(t,a[column][good],e,p,epoch,duration))
        if additional:tested['additional_phase_secondary']=local_quadratic(t,a[column][good],e,p,additional['epoch_btjd'],additional['duration_days'])[0]
        local_data[column]=loc
        for width in [.75,1.5]:
            flat=flatten(a,column,width);stats=event_stats(flat['t'],flat['f'],flat['e'],p,epoch,duration)
            tested['alternative_windows'].append(dict(width_days=width,stats=stats))
            for key in ['t','f','e']:intermediate[f'{column}_w{width:g}_{key}']=flat[key]
        tests[column]=tested
    cadence=float(np.median(np.diff(np.sort(t))));shape=fit_shape(local_data['corrected'],p,epoch,duration,cadence)
    centroids={k:centroid_diagnostic(t,a[k][good],p,epoch,duration)for k in ['centroid_x','centroid_y']if k in a}
    doubled={}
    for column in ['corrected','raw']:
        d1=local_quadratic(t,a[column][good],e,2*p,epoch,duration)[0];d2=local_quadratic(t,a[column][good],e,2*p,epoch+p,duration)[0]
        doubled[column]=dict(first=d1,alternating=d2,difference_sigma=abs(d1['depth']-d2['depth'])/np.hypot(d1['depth_err'],d2['depth_err'])if d1.get('depth')is not None and d2.get('depth')is not None else None)
    flags=[];primary=tests['corrected']['primary'];raw=tests['raw']['primary']
    if primary['snr']<6 or raw['snr']<5:flags.append('weak_or_inconsistent_raw_local_events')
    if primary['n_events']<3 or primary.get('n_positive_events',0)<3:flags.append('fewer_than_three_positive_locally_bracketed_events')
    if any(v['stats']['local_snr']<5 for v in tests['corrected']['alternative_windows']):flags.append('detrending_sensitivity_or_weak_signal')
    if(primary.get('odd_even_sigma')or 0)>4:flags.append('alternating_depth_warning')
    if shape and shape['ingress_fraction']>.45:flags.append('v_shaped_or_cadence_unresolved')
    for key in ['half_phase_secondary','additional_phase_secondary']:
        sec=tests['corrected'].get(key,{});rsec=tests['raw'].get(key,{})
        if sec.get('snr',0)>5 and rsec.get('snr',0)>4 and(sec.get('depth')or 0)>.15*(primary.get('depth')or 1):flags.append(key+'_warning')
    if primary.get('depth',0)>0:
        variability=tests['raw']['variability']
        for name in ['orbital','half_period']:
            if(variability.get(name+'_semiamplitude')or 0)>.5*primary['depth']and(variability.get(name+'_snr')or 0)>5:flags.append('large_raw_'+name+'_variability')
        if'small'in tests and'large'in tests:
            sm,la=tests['small']['primary'],tests['large']['primary']
            if sm.get('depth')is not None and la.get('depth')is not None and abs(sm['depth']-la['depth'])>max(.3*primary['depth'],4*np.hypot(sm['depth_err'],la['depth_err'])):flags.append('aperture_depth_dependence')
    if any((v.get('sigma')or 0)>5 and abs(v.get('offset_pixels')or 0)>.001 for v in centroids.values()):flags.append('centroid_shift_warning')
    result=dict(tic=int(record['tic']),sector=int(record['sector']),checked_utc=utc(),method=method(),source=record,source_metadata=provenance,input_peak=peak,refined=ephem,column_tests=tests,doubled_period_alternating_test=doubled,integrated_trapezoid=shape,centroid_tests=centroids,additional_phase_search=additional,quality_context=quality_context(a,p,epoch,duration),warnings=flags,assessment='unexplained_transit_like_lead_requires_recognition_recurrence_and_pixels'if not flags else'fails_or_incomplete_initial_falsification',limitations=['Same light curve through multiple detrenders is not independent evidence.','Scores are conditional diagonal estimates; red noise and the search trials are not calibrated here.','Raw SAP and extra apertures use corrected-flux errors with local residual floors.','No transit source localization, dilution correction, companion mass or novelty claim.'])
    if output:
        directory=OUT/f'tic{record["tic"]}'/f'sector{record["sector"]}_p{float(peak.get("period",peak.get("period_days"))):.6f}';directory.mkdir(parents=True,exist_ok=True)
        save(directory/'vetting.json',result);np.savez_compressed(directory/'diagnostic_arrays.npz',**intermediate)
        fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
        for row,column in enumerate(['corrected','raw','small']):
            if column not in tests:continue
            x,y,z=[intermediate[f'{column}_w1.5_{key}']for key in ['t','f','e']]
            axes[row,0].plot(x,(y-1)*1000,'.',ms=1,color='.4');axes[row,0].set(xlabel='BTJD TDB',ylabel=column+' (ppt)')
            ph=phase(x,p,epoch)/p;axes[row,1].plot(ph,(y-1)*1000,'.',ms=1,color='.65')
            edges=np.linspace(-.5,.5,151);idx=np.digitize(ph,edges)
            xx,yy=[],[]
            for j in range(1,len(edges)):
                select=idx==j
                if select.any():xx.append(np.mean(ph[select]));yy.append(np.mean(y[select]-1)*1000)
            axes[row,1].plot(xx,yy,'o-',ms=2);axes[row,1].set(xlabel=f'Phase P={p:.6f} d',ylabel='1.5 d median normalization (ppt)')
        fig.suptitle(f'TIC {record["tic"]}, sector {record["sector"]}: '+('; '.join(flags)or'Unexplained lead; validation incomplete'),fontsize=9);fig.savefig(directory/'diagnostics.png',dpi=130);plt.close(fig)
        (directory/'README.md').write_text(f'# TIC {record["tic"]}, sector {record["sector"]}\n\nPeriod {p:.7f} d; duration {duration*24:.2f} h. Local corrected/raw S/N {primary["snr"]:.2f}/{raw["snr"]:.2f}; {primary["n_events"]} bracketed events.\n\nAssessment: **{result["assessment"]}**.\n\nWarnings: '+(', '.join(flags)or'None in these limited checks.')+'\n\nNo planet, source-location or novelty claim. See `vetting.json` for provenance and individual event measurements.\n')
        print(record['tic'],record['sector'],round(p,6),'local',round(primary['snr'],2),round(raw['snr'],2),'warnings',flags,flush=True)
    return result

def stellar_variability_check(record, ephemeris):
    """Check pulsation/rotation alternatives with an out-of-event baseline.

    An eighth-order Fourier series may approximate many waveforms. It is an
    explanatory diagnostic, not evidence that an overlapping planet is absent.
    """
    from astropy.timeseries import LombScargle
    arrays,metadata=read_source(record);p=ephemeris['period'];epoch=ephemeris['epoch_btjd'];duration=ephemeris['duration_days'];output=[];plotdata=[]
    for col in ['corrected','raw']:
        flat=flatten(arrays,col,1.5);t,f,e=flat['t'],flat['f'],flat['e'];ph=phase(t,p,epoch)/p
        oot=abs(phase(t,p,epoch))>1.5*duration;theta=2*np.pi*(t-epoch)/p
        design=np.column_stack([np.ones(len(t))]+[fun(k*theta)for k in range(1,9)for fun in [np.sin,np.cos]])
        beta=np.linalg.lstsq(design[oot]/e[oot,None],f[oot]/e[oot],rcond=None)[0];model=design@beta;residual=1+f-model
        original=event_stats(t,f,e,p,epoch,duration);remaining=event_stats(t,residual,e,p,epoch,duration)
        frequency=np.arange(.2,min(30.,.45/np.median(np.diff(t))),1/(5*np.ptp(t)))
        ls=LombScargle(t[oot],f[oot],e[oot]);power=ls.power(frequency);peaks=[]
        for j in np.argsort(power)[::-1]:
            if any(abs(frequency[j]-v['frequency_per_day'])<1/np.ptp(t)for v in peaks):continue
            peaks.append(dict(frequency_per_day=float(frequency[j]),period_days=float(1/frequency[j]),power=float(power[j]),ratio_to_trial_orbital_frequency=float(frequency[j]*p)))
            if len(peaks)==5:break
        dominant=peaks[0]['frequency_per_day'];theta_ls=2*np.pi*(t-t.min())*dominant;design_ls=np.column_stack([np.ones(len(t)),np.sin(theta_ls),np.cos(theta_ls)])
        beta_ls=np.linalg.lstsq(design_ls[oot]/e[oot,None],f[oot]/e[oot],rcond=None)[0];model_ls=design_ls@beta_ls
        single_residual=event_stats(t,1+f-model_ls,e,p,epoch,duration)
        amplitudes=[float(np.hypot(beta[2*k-1],beta[2*k]))for k in range(1,9)]
        output.append(dict(column=col,highest_out_of_event_ls_peaks=peaks,eight_harmonic_semiamplitudes=amplitudes,original_primary=original,harmonic_baseline_residual_primary=remaining,dominant_frequency_residual_primary=single_residual,dominant_frequency_model_nparameters=3,n_out_of_event_fit=int(oot.sum()),warning='Harmonic order and period are selected diagnostics. Model flexibility, detrending, correlated noise and search trials are not calibrated; no exclusion of a coexisting transit.'))
        plotdata.append((col,ph,f,model,model_ls,frequency,power))
    directory=OUT/f'tic{record["tic"]}';destination=directory/f'sector{record["sector"]}_stellar_variability_v2.json';save(destination,dict(tic=int(record['tic']),sector=int(record['sector']),checked_utc=utc(),source=record,source_metadata=metadata,method=method(),ephemeris=ephemeris,tests=output))
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,(col,ph,f,model,model_ls,frequency,power)in zip(axes,plotdata):
        order=np.argsort(ph);ax[0].plot(ph,(f-1)*1000,'.',ms=1,color='.6');ax[0].plot(ph[order],(model[order]-1)*1000,color='tab:orange',lw=1.3,label='8 harmonics');ax[0].plot(ph[order],(model_ls[order]-1)*1000,color='tab:green',lw=.8,label='dominant frequency');ax[0].legend(fontsize=7);ax[0].set(xlabel='Phase',ylabel=col+' (ppt)',title='8 harmonics fit outside candidate events')
        ax[1].plot(frequency,power,lw=.8);ax[1].set(xlabel='Frequency (cycles/day)',ylabel='Out-of-event LS power')
    fig.suptitle('Stellar-variability alternative; flexible baseline does not exclude a planet');fig.savefig(directory/f'sector{record["sector"]}_stellar_variability_v2.png',dpi=140);plt.close(fig)
    print('STELLAR_BASELINE',record['tic'],record['sector'],[(test['column'],test['highest_out_of_event_ls_peaks'][0],test['harmonic_baseline_residual_primary']['local_snr'])for test in output],flush=True)
    return output


def joint_recurrence(records, peak, name='joint'):
    """Targeted long-baseline refinement with one common epoch across sectors.

    The refinement interval is declared from the original sector baseline. This
    is follow-up optimization on selected data, not blind recovery or a FAP.
    """
    parsed=[];curves=[]
    for record in records:
        arrays,metadata=read_source(record);flat=flatten(arrays,'corrected',1.5)
        if len(flat['t'])<40:continue
        parsed.append((record,arrays,metadata));curves.append(flat)
    if len(parsed)<2:raise ValueError('At least two usable sectors required')
    t,f,e=[np.concatenate([flat[k]for flat in curves])for k in ['t','f','e']]
    order=np.argsort(t);t,f,e=t[order],f[order],e[order]
    p=float(peak.get('period',peak.get('period_days')));duration=float(peak.get('duration_days',.1))
    original_baseline=np.ptp(curves[0]['t']);baseline=np.ptp(t)
    span=max(.003*p,.7*duration*p/original_baseline)
    requested_step=.2*duration*p/baseline
    count=min(6001,max(401,int(np.ceil(2*span/requested_step))+1))
    periods=np.linspace(p-span,p+span,count)
    durations=np.unique(np.clip([.7*duration,duration,1.3*duration],.015,.2*p))
    result=BoxLeastSquares(t,f,e).power(periods,durations,oversample=15,objective='snr');j=int(np.argmax(result.power))
    period,epoch,dur=float(result.period[j]),float(result.transit_time[j]),float(result.duration[j]);per_sector=[];local_arrays={}
    for record,arrays,metadata in parsed:
        good=(arrays['QUALITY']==0)&np.isfinite(arrays['TIME'])&np.isfinite(arrays['error'])&(arrays['error']>0)
        tests={}
        for col in ['corrected','raw']:
            primary,loc=local_quadratic(arrays['TIME'][good],arrays[col][good],arrays['error'][good],period,epoch,dur)
            tests[col]=primary
            if loc is not None:
                for key in ['t','f','e']:local_arrays[f's{record["sector"]}_{col}_{key}']=loc[key]
        per_sector.append(dict(sector=int(record['sector']),tests=tests,source=record,metadata=metadata))
    pooled={}
    for col in ['corrected','raw']:
        events=[event for sector in per_sector for event in sector['tests'][col]['events']]
        if not events:
            pooled[col]=dict(depth=None,depth_err=None,snr=0.,n_events=0,n_positive_events=0);continue
        ds=np.array([v['depth']for v in events]);es=np.array([v['depth_err']for v in events]);w=1/es**2
        depth=float(np.sum(w*ds)/np.sum(w));diagonal=float(1/np.sqrt(np.sum(w)));scatter=float(np.std(ds,ddof=1)/np.sqrt(len(ds)))if len(ds)>1 else diagonal
        error=max(diagonal,scatter)
        pooled[col]=dict(depth=depth,depth_err=error,snr=depth/error,n_events=len(events),n_positive_events=int((ds>0).sum()),diagonal_error=diagonal,event_scatter_error=scatter)
    output=dict(checked_utc=utc(),method=method(),tic=int(records[0]['tic']),period_days=period,epoch_btjd=epoch,duration_days=dur,diagonal_bls_snr=float(result.depth_snr[j]),targeted_search=dict(period_window=[float(periods[0]),float(periods[-1])],period_count=count,requested_step=requested_step,actual_step=float(periods[1]-periods[0]),baseline_days=float(baseline),original_baseline_days=float(original_baseline)),per_sector=per_sector,pooled_local_events=pooled,warning='A common epoch was optimized on all selected sectors. Neither selection, red noise, nor multiple trials is calibrated; no independent confirmation or planet claim.')
    directory=OUT/f'tic{records[0]["tic"]}';save(directory/f'{name}.json',output);np.savez_compressed(directory/f'{name}_arrays.npz',**local_arrays)
    fig,axes=plt.subplots(len(parsed),2,figsize=(12,3*len(parsed)),squeeze=False,layout='constrained')
    for ax,(record,arrays,metadata),flat in zip(axes,parsed,curves):
        ax[0].plot(flat['t'],(flat['f']-1)*1000,'.',ms=1);ax[0].set(xlabel='BTJD TDB',ylabel='Corrected normalized flux (ppt)',title=f'Sector {record["sector"]}')
        ph=phase(flat['t'],period,epoch)/period;ax[1].plot(ph,(flat['f']-1)*1000,'.',ms=1,color='.6');ax[1].set(xlabel=f'Common phase P={period:.7f} d',ylabel='Corrected normalized flux (ppt)')
        ax[1].set_xlim(-min(.2,2*dur/period),min(.2,2*dur/period))
    fig.suptitle('Targeted common-epoch follow-up, not a significance or source-localization test');fig.savefig(directory/f'{name}.png',dpi=140);plt.close(fig)
    print('JOINT',records[0]['tic'],period,[(r['sector'],r['tests']['corrected']['snr'],r['tests']['raw']['snr'])for r in per_sector],flush=True)
    return output


def qlp_url(tic,sector):
    version='v02'if 74<=sector<=79 else'v01';s=str(tic).zfill(16);sub='/'.join(s[i:i+4]for i in range(0,16,4));name=f'hlsp_qlp_tess_ffi_s{sector:04}-{s}_tess_{version}_llc.fits'
    return f'https://archive.stsci.edu/hlsps/qlp/s{sector:04}/{sub}/{name}'

def retrieve_independent(record,max_sectors=4):
    """Coverage is geometric only. Every missing/failed public URL is retained."""
    tic=int(record['tic']);folder=CACHE/'followup'/f'tic{tic}';folder.mkdir(parents=True,exist_ok=True);manifest=OUT/f'tic{tic}'/'retrieval.json'
    old=json.loads(manifest.read_text())if manifest.exists()else{}
    query=folder/'tesscut_coverage.json';url=f'https://mast.stsci.edu/tesscut/api/v0.1/sector?ra={float(record["ra"]):.10f}&dec={float(record["dec"]):.10f}&units=d'
    if not query.exists():
        response=requests.get(url,timeout=45);query.write_bytes(response.content);save(query.with_suffix('.manifest.json'),dict(url=url,status=response.status_code,sha256=digest(query),checked_utc=utc()))
        response.raise_for_status()
    query_meta=json.loads(query.with_suffix('.manifest.json').read_text())
    if query_meta['status']!=200 or query_meta['sha256']!=digest(query):raise ValueError('Coverage query missing, failed or changed')
    coverage=json.loads(query.read_text());sectors=sorted({int(r['sector'])for r in coverage.get('results',[])})
    # QLP 99--104 carries a documented processing warning; do not use it as
    # independent confirmation until a corrected release is verified.
    eligible=[s for s in sectors if s!=int(record['sector'])and not 99<=s<=104]
    # Bracket the discovery epoch when both earlier/later products exist.
    earlier=sorted([s for s in eligible if s<int(record['sector'])],reverse=True);later=sorted([s for s in eligible if s>int(record['sector'])]);selected=[]
    while len(selected)<max_sectors and(earlier or later):
        if earlier:selected.append(earlier.pop(0))
        if later and len(selected)<max_sectors:selected.append(later.pop(0))
    rows={int(r['sector']):r for r in old.get('observations',[])}
    for sector in selected:
        if sector in rows and rows[sector].get('status')=='retained':
            if digest(rows[sector]['path'])!=rows[sector]['sha256']:raise ValueError('Follow-up raw SHA mismatch')
            continue
        product_url=qlp_url(tic,sector);path=folder/Path(product_url).name
        row=dict(tic=tic,sector=sector,url=product_url,checked_utc=utc())
        try:
            if not path.exists():
                response=requests.get(product_url,timeout=60);row['http_status']=response.status_code
                if response.status_code!=200 or not response.content.startswith(b'SIMPLE  '):
                    error=folder/f'sector{sector}_response.bin';error.write_bytes(response.content);row.update(status='failed_or_unavailable',response_path=str(error),response_sha256=digest(error));rows[sector]=row;continue
                tmp=path.with_suffix('.part');tmp.write_bytes(response.content);tmp.replace(path)
            row.update(path=str(path),sha256=digest(path),bytes=path.stat().st_size,status='retained',ra=record['ra'],dec=record['dec'])
            _,meta=read_source(row);row['verified']=meta
        except Exception as exc:row.update(status='retrieval_error',error=repr(exc))
        rows[sector]=row
        save(manifest,dict(tic=tic,coverage=dict(url=url,path=str(query),sha256=digest(query),geometric_sectors=sectors),excluded_qlp_warning_sectors=[s for s in sectors if 99<=s<=104],observations=list(rows.values())))
    result=dict(tic=tic,coverage=dict(url=url,path=str(query),sha256=digest(query),geometric_sectors=sectors),excluded_qlp_warning_sectors=[s for s in sectors if 99<=s<=104],observations=list(rows.values()))
    save(manifest,result);return result

def localize_modulation(record, frequency, other_ra, other_dec, size=9):
    """Native-pixel sinusoidal amplitude localization, with explicit caveats."""
    import io,zipfile
    from urllib.parse import urlencode
    from astropy.wcs import WCS
    from pixel_vet_49617264 import gaussfit
    tic=int(record['tic']);sector=int(record['sector']);folder=CACHE/'followup'/f'tic{tic}'/'pixels';folder.mkdir(parents=True,exist_ok=True);out=OUT/f'tic{tic}'/f'sector{sector}_pixel_modulation';out.mkdir(parents=True,exist_ok=True)
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urlencode(dict(ra=record['ra'],dec=record['dec'],x=size,y=size,units='px',sector=sector))
    manifest=out/'raw_manifest.json';stored=json.loads(manifest.read_text())if manifest.exists()else{}
    paths=list(folder.glob(f'tess-s{sector:04d}-*_{size}x{size}_astrocut.fits'))
    if not paths:
        response=requests.get(url,timeout=90);raw=folder/f'sector{sector}_{size}x{size}_response.zip';raw.write_bytes(response.content);save(out/'http_manifest.json',dict(url=url,status=response.status_code,path=str(raw),sha256=digest(raw),bytes=raw.stat().st_size,checked_utc=utc()));response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(response.content))as archive:
            for name in archive.namelist():
                if name.endswith('.fits'):
                    path=folder/Path(name).name
                    if path.exists():raise ValueError('Refusing raw pixel overwrite')
                    path.write_bytes(archive.read(name));paths.append(path)
    if len(paths)!=1:raise ValueError('Expected exactly one native cutout')
    path=paths[0];sha=digest(path)
    if stored and stored['sha256']!=sha:raise ValueError('Pixel raw SHA mismatch')
    with fits.open(path,memmap=False)as hd:
        hd.verify('exception');checks=verify_fits_checksums(path,hd);h=hd[1].header;d=hd[1].data
        if h.get('TIMESYS',hd[0].header.get('TIMESYS'))!='TDB' or h.get('BJDREFI',0)+h.get('BJDREFF',0)!=2457000:raise ValueError('Pixel time is not BTJD TDB')
        wcs=WCS(hd[2].header);positions=wcs.all_world2pix(np.array([[record['ra'],record['dec']],[other_ra,other_dec]],dtype=float),0)
        good=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.all(np.isfinite(d['FLUX']),axis=(1,2))
        t=np.asarray(d['TIME'][good],float);cube=np.asarray(d['FLUX'][good],float);errors=np.asarray(d['FLUX_ERR'][good],float);exposure=float(h['TIMEDEL'])
        ffi=[str(x)for x in d['FFI_FILE']]if'FFI_FILE'in d.names else[]
        science=hashlib.sha256()
        for key,typ in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:science.update(key.encode());science.update(np.asarray(d[key],dtype=typ).tobytes())
        column=h['1CRV5P']+(positions[0,0]+1-h['1CRP5P'])*h['1CDL5P'];row=h['2CRV5P']+(positions[0,1]+1-h['2CRP5P'])*h['2CDL5P']
    source=dict(url=url,path=str(path),sha256=sha,science_sha256=science.hexdigest(),bytes=path.stat().st_size,checksums=checks,n_native=len(d),n_quality_finite=len(t),exposure_days=exposure,detector_column=float(column),detector_row=float(row),on_science_silicon=bool(44<=column<2092 and 0<=row<2048),checked_utc=utc());save(manifest,source);(out/'ffi_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    if not source['on_science_silicon']:raise ValueError('Cutout target is outside science silicon')
    # Constant+linear baselines per contiguous spacecraft interval.
    splits=np.r_[0,np.flatnonzero(np.diff(t)>max(.08,3*exposure))+1,len(t)];segment=np.zeros(len(t),int)
    for j,(a,b)in enumerate(zip(splits[:-1],splits[1:])):segment[a:b]=j
    basis=[]
    for j in range(segment.max()+1):
        mask=segment==j;basis.extend([mask.astype(float),np.where(mask,t-np.mean(t[mask]),0)])
    theta=2*np.pi*(t-t.min())*frequency;design=np.column_stack(basis+[np.sin(theta),np.cos(theta)]);inverse=np.linalg.pinv(design);matrix=cube.reshape(len(t),-1);coeff=inverse@matrix;residual=matrix-design@coeff
    scatter=np.maximum(np.std(residual,axis=0,ddof=design.shape[1]),np.sqrt(np.mean(errors.reshape(len(t),-1)**2,axis=0)))
    yy,xx=np.indices(cube.shape[1:]);central=((xx-positions[0,0])**2+(yy-positions[0,1])**2)<2.5**2
    summed=coeff[-2:,central.ravel()].sum(axis=1);unit=summed/np.linalg.norm(summed);projection=unit@inverse[-2:];image=(unit@coeff[-2:]).reshape(cube.shape[1:]);image_error=(scatter*np.linalg.norm(projection)).reshape(image.shape)
    roi=(xx-positions[0,0])**2+(yy-positions[0,1])**2<3.8**2;free,_=gaussfit(image,image_error,positions[1],roi);fixed={}
    for label,point in zip(['catalog_target','known_ew_neighbor'],positions):
        fit,_=gaussfit(image,image_error,point,roi,fixed_center=True);fixed[label]=fit
    center=np.array([free['x'],free['y']]);sky=wcs.all_pix2world([center],0)[0]
    result=dict(tic=tic,sector=sector,method=method(),gaussian_helper_sha256=digest(ROOT/'src/pixel_vet_49617264.py'),source=source,frequency_per_day=frequency,period_days=1/frequency,target_xy=positions[0],known_ew_xy=positions[1],known_ew_radec=[other_ra,other_dec],free_gaussian=free,fixed_positions=fixed,free_radec=sky,offset_target_pixels=float(np.linalg.norm(center-positions[0])),offset_known_ew_pixels=float(np.linalg.norm(center-positions[1])),method_description='Project each native-pixel sine/cosine coefficient onto the summed-aperture modulation phase; fit an integrated empirical Gaussian plus background plane.',limitations=['Conditioned frequency/phase and pixel-diagonal errors; not a calibrated localization confidence interval.','Empirical Gaussian is not the TESS PRF; source blending, WCS and proper-motion systematics remain.','No pixel-level correction for all spacecraft trends; constant/linear baselines are fitted per contiguous interval.','Localizing the dominant periodic variable does not rule out a weaker overlapping planetary transit.'])
    save(out/'localization.json',result);np.savez_compressed(out/'maps.npz',modulation=image,error=image_error,reference=np.median(cube,axis=0),coefficients=coeff[-2:].reshape(2,*image.shape),target_xy=positions[0],known_ew_xy=positions[1])
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    for ax,img,title in zip(axes,[np.median(cube,axis=0),image,image/image_error],['Native median flux','Projected modulation (e-/s)','Conditional diagonal S/N']):
        im=ax.imshow(img,origin='lower');ax.plot(*positions[0],'r+',ms=12,label='catalog target');ax.plot(*positions[1],'cx',ms=10,label='known EW');ax.plot(*center,'wo',ms=5,mfc='none',label='Gaussian fit');ax.set_title(title,fontsize=9);fig.colorbar(im,ax=ax,shrink=.7)
    axes[0].legend(fontsize=7);fig.suptitle(f'TIC{tic}, S{sector}, modulation P={1/frequency:.7f}d; illustrative source localization');fig.savefig(out/'localization.png',dpi=150);plt.close(fig)
    print('LOCALIZATION',tic,sector,result['offset_target_pixels'],result['offset_known_ew_pixels'],free,fixed,flush=True);return result


def write_summary():
    rows=[]
    for path in sorted(OUT.glob('tic*/sector*/vetting.json')):
        result=json.loads(path.read_text());primary=result['column_tests']['corrected']['primary'];raw=result['column_tests']['raw']['primary']
        rows.append(dict(tic=result['tic'],sector=result['sector'],period_days=result['refined']['period'],duration_hours=result['refined']['duration_days']*24,local_depth_ppm=(primary['depth']or 0)*1e6,local_depth_error_ppm=(primary['depth_err']or 0)*1e6,local_corrected_snr=primary['snr'],local_raw_snr=raw['snr'],n_local_events=primary['n_events'],n_positive_local_events=primary.get('n_positive_events',0),warnings=';'.join(result['warnings']),assessment=result['assessment'],report=str(path.relative_to(ROOT))))
    if rows:
        destination=OUT/'followup_summary.csv'
        with destination.open('w',newline='')as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    return rows


def records_from_db(tic=None):
    con=sqlite3.connect(f'file:{DB}?mode=ro',uri=True);tables=[r[0]for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")];records=[]
    for table in tables:
        names=[r[1]for r in con.execute('PRAGMA table_info('+table+')')]
        if'result_json'not in names:continue
        sql='SELECT result_json FROM '+table+' WHERE result_json IS NOT NULL';args=[]
        if tic and'tic'in names:sql+=' AND tic=?';args=[tic]
        for raw, in con.execute(sql,args):
            record=json.loads(raw)
            if not tic or int(record['tic'])==tic:records.append(record)
    con.close();return records

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tic',type=int);parser.add_argument('--sector',type=int);parser.add_argument('--period',type=float);parser.add_argument('--limit',type=int,default=12);parser.add_argument('--followup',type=int,default=0);parser.add_argument('--min-snr',type=float,default=7.);parser.add_argument('--summarize',action='store_true');args=parser.parse_args()
    if args.summarize:
        print(json.dumps(dict(summarized_observations=len(write_summary()))));return
    queue=[]
    for record in records_from_db(args.tic):
        if args.sector and int(record['sector'])!=args.sector:continue
        for scan in record.get('scans',[record.get('screen',record.get('result',{}))]):
            for peak in scan.get('peaks',[]):
                if args.period and abs(float(peak['period'])-args.period)>.03*args.period:continue
                if float(peak.get('local_snr',peak.get('snr',0)))<args.min_snr:continue
                queue.append((float(peak.get('local_snr',peak.get('snr',0))),record,peak))
    queue.sort(key=lambda x:x[0],reverse=True);seen=set();results=[]
    for score,record,peak in queue:
        key=(int(record['tic']),int(record['sector']))
        if key in seen:continue
        seen.add(key)
        try:
            result=vet_observation(record,peak);results.append(result)
            if args.followup and result['column_tests']['raw']['primary']['snr']>=5 and len(result['warnings'])<3:
                extra=retrieve_independent(record,args.followup)
                for obs in extra['observations']:
                    if obs.get('status')=='retained':vet_observation(obs,dict(period=result['refined']['period'],duration_days=result['refined']['duration_days']))
        except Exception as exc:
            save(OUT/f'tic{record["tic"]}'/f'sector{record["sector"]}_error.json',dict(error=repr(exc),traceback=traceback.format_exc(),checked_utc=utc()));print('ERROR',key,repr(exc),flush=True)
        if len(seen)>=args.limit:break
    summary=[dict(tic=r['tic'],sector=r['sector'],period=r['refined']['period'],local_corrected_snr=r['column_tests']['corrected']['primary']['snr'],local_raw_snr=r['column_tests']['raw']['primary']['snr'],warnings=r['warnings'],assessment=r['assessment'])for r in results]
    save(OUT/'latest_queue_summary.json',dict(checked_utc=utc(),method=method(),evaluated=summary));write_summary()

if __name__=='__main__':main()
