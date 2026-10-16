"""Fixed pooled-ephemeris falsification of seven northern TESS search leads.

Imports frozen scientific helpers read-only. Extra sectors are genuine new
observations but the pooled search ephemeris has selection/period uncertainty.
"""
from pathlib import Path
import sys, json, csv, argparse, traceback, hashlib, sqlite3
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import multisector_followup as mf
from astropy.timeseries import LombScargle
from scipy.ndimage import median_filter
TARGETS={353807236:6.224717707910657,359679661:1.9134095518192817,353933809:18.978402288292596,366490373:2.810238855348628,287086419:12.613387995535525,356821129:3.871952352751222,353783909:6.397423709551385}
NEIGHBOR_TARGETS={280163751:.5743168974311555,199610815:.6554869878286899,1271331520:.8828000628488375}
DRIVER_SOURCE=Path(__file__).read_bytes()
SOURCE_SHA=hashlib.sha256(DRIVER_SOURCE).hexdigest()

def native_variability(record,ephem):
    """Out-of-event native sampling LS, no binning, conditional diagnostic."""
    a,meta=mf.read_source(record);result=[]
    for col in ['corrected','raw']:
        good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a[col])&np.isfinite(a['error'])&(a['error']>0)
        t,f,e=a['TIME'][good],a[col][good],a['error'][good]
        order=np.argsort(t);t,f,e=t[order],f[order],e[order];dt=np.median(np.diff(t));norm=np.median(f)
        f=f/norm;e=e/norm
        width=int(np.ceil(1.5/dt));width+=width%2==0
        trend=median_filter(f,size=width,mode='nearest');f=f/trend;e=e/trend
        mask=(abs(mf.phase(t,ephem['period'],ephem['epoch_btjd']))>1.5*ephem['duration_days'])&(t>t.min()+.06)&(t<t.max()-.06)
        t,f,e=t[mask],f[mask],e[mask]
        freq=np.arange(.2,min(180.,.45/dt),1/(5*np.ptp(t)))
        ls=LombScargle(t,f,e);power=ls.power(freq);peaks=[]
        for j in np.argsort(power)[::-1]:
            if any(abs(freq[j]-p['frequency_per_day'])<1/np.ptp(t)for p in peaks):continue
            amp=float(np.hypot(*ls.model_parameters(freq[j])[1:]))
            peaks.append(dict(frequency_per_day=float(freq[j]),period_days=float(1/freq[j]),power=float(power[j]),semiamplitude_ppm=amp*1e6,ratio_to_orbital_frequency=float(freq[j]*ephem['period'])))
            if len(peaks)==8:break
        result.append(dict(column=col,native_good_out_of_event=len(t),cadence_days=float(dt),frequency_range_per_day=[float(freq.min()),float(freq.max())],top_peaks=peaks,warning='Conditional exploratory LS; no correlated-noise/search FAP or rejection of coexisting planet.'))
    return result

def no_bracket_fallback(record,ephem):
    """Document zero conservatively bracketed events without mutating helper."""
    arrays,meta=mf.read_source(record);good=(arrays['QUALITY']==0)&np.isfinite(arrays['TIME'])&np.isfinite(arrays['error'])&(arrays['error']>0)
    tests={}
    for col in ['corrected','raw','small','large']:
        if col not in arrays:continue
        primary,_=mf.local_quadratic(arrays['TIME'][good],arrays[col][good],arrays['error'][good],ephem['period'],ephem['epoch_btjd'],ephem['duration_days'])
        alternative=[]
        for width in [.75,1.5]:
            flat=mf.flatten(arrays,col,width);alternative.append(dict(width_days=width,stats=mf.event_stats(flat['t'],flat['f'],flat['e'],ephem['period'],ephem['epoch_btjd'],ephem['duration_days'])))
        tests[col]=dict(primary=primary,alternative_windows=alternative)
    if tests['corrected']['primary']['n_events']!=0:raise ValueError('Fallback only applies to zero locally bracketed corrected events')
    return dict(tic=int(record['tic']),sector=int(record['sector']),refined=ephem,column_tests=tests,source=record,source_metadata=meta,warnings=['zero_contiguous_locally_bracketed_events','frozen_helper_None_depth_guard_unsupported'],assessment='insufficient_contiguous_local_event_coverage',limitations=['Frozen helper raises TypeError for depth=None. This wrapper records its zero-event result and alternative median-window measurements without changing helper source.','Zero local events is missing diagnostic coverage, not zero depth or excluded planet.'])

def run(tic):
    table=list(csv.DictReader(open(ROOT/'results/multisector/combined_signals.csv')))
    rows=[r for r in table if int(r['tic'])==tic and abs(float(r['period_days'])-TARGETS[tic])<1e-6]
    row=max(rows,key=lambda r:float(r['local_snr']))
    pooledpath=Path(row['evidence_path']);pooled=json.loads(pooledpath.read_text())
    ephem=dict(period=float(row['period_days']),epoch_btjd=float(row['epoch_btjd']),duration_days=float(row['duration_days']))
    peak=dict(period=ephem['period'],duration_days=ephem['duration_days'])
    records=mf.records_from_db(tic);records=[r for r in records if int(r['sector'])in[83,86]]
    if len(records)!=2:raise ValueError('Expected exactly two original sectors')
    known={int(s['sector']):s for s in pooled['identity']['sources']}
    for r in records:
        if r['sha256']!=known[int(r['sector'])]['sha256'] or r['path']!=known[int(r['sector'])]['path']:raise ValueError('Pooled and database raw identity mismatch')
    d=mf.OUT/f'tic{tic}';d.mkdir(parents=True,exist_ok=True)
    snapshot=d/'method_snapshots'/(SOURCE_SHA+'.py');snapshot.parent.mkdir(exist_ok=True)
    if not snapshot.exists():snapshot.write_bytes(DRIVER_SOURCE)
    output=dict(tic=tic,checked_utc=mf.utc(),driver_source_sha256=SOURCE_SHA,pooled_source=dict(path=str(pooledpath),sha256=mf.digest(pooledpath)),selected_blind_row=row,fixed_ephemeris=ephem,original_records=records,observations=[],errors=[],limitations=['All original sectors participated in the blind pooled selection.','Extra sectors use the exact pooled ephemeris without sector retuning. Long-baseline period uncertainty can cause phase drift; a nondetection is not in itself exclusion.','Conditional scores, raw/corrected variants and multiple detrenders are not independent false-alarm probabilities.'])
    mf.save(d/'pooled_north_input.json',output)
    retrieval=mf.retrieve_independent(records[0],max_sectors=2,exclude_sectors=[83,86])
    if sum(r.get('status')=='retained' for r in retrieval['observations'])<2:
        retrieval=mf.retrieve_independent(records[0],max_sectors=3,exclude_sectors=[83,86])
    output['retrieval']=retrieval
    extra=[r for r in retrieval['observations']if r.get('status')=='retained' and int(r['sector'])not in[83,86]]
    for record in records+extra:
        try:
            try:result=mf.vet_observation(record,peak,fixed_ephemeris=ephem)
            except TypeError as exc:
                if "'>' not supported between instances of 'NoneType' and 'int'" not in str(exc):raise
                result=no_bracket_fallback(record,ephem)
                folder=d/f'sector{record["sector"]}_fixed_p{ephem["period"]:.6f}';mf.save(folder/'vetting.json',result)
            variability=mf.stellar_variability_check(record,ephem)
            native=native_variability(record,ephem)
            mf.save(d/f'sector{record["sector"]}_native_variability.json',dict(tic=tic,sector=record['sector'],source=record,ephemeris=ephem,driver_source_sha256=SOURCE_SHA,tests=native))
            output['observations'].append(dict(sector=int(record['sector']),discovery_sector=int(record['sector'])in[83,86],vetting=result,stellar_variability=variability,native_variability=native))
        except Exception as exc:
            output['errors'].append(dict(sector=int(record['sector']),error=repr(exc),traceback=traceback.format_exc()));print('FAILED',tic,record['sector'],repr(exc),flush=True)
        mf.save(d/'pooled_north_validation.json',output)
    return output

def instrumental_event_check():
    """Check a single selected S86 event against all retained same-detector stars.

    Same-time sky controls diagnose shared artifacts; they are not a null
    experiment for the blind survey, and no binomial/normal FAP is inferred.
    """
    tic=353933809;d=mf.OUT/f'tic{tic}';src=json.loads((d/'pooled_north_validation.json').read_text())
    event=next(r for r in src['observations']if r['sector']==86)
    epoch=event['vetting']['column_tests']['corrected']['primary']['events'][0]['time']
    duration=src['fixed_ephemeris']['duration_days']
    con=sqlite3.connect(f'file:{mf.DB}?mode=ro',uri=True)
    records=[json.loads(r[0])for r in con.execute('SELECT result_json FROM observations WHERE sector=86 AND result_json IS NOT NULL')];con.close()
    selected=[r for r in records if r['metadata'].get('CAMERA')==4 and r['metadata'].get('CCD')==4]
    results=[];errors=[];local={}
    for record in selected:
        try:
            arrays,meta=mf.read_source(record)
            good=(arrays['QUALITY']==0)&np.isfinite(arrays['TIME'])&np.isfinite(arrays['error'])&(arrays['error']>0)
            tests={}
            for col in ['corrected','raw']:
                tests[col],loc=mf.local_quadratic(arrays['TIME'][good],arrays[col][good],arrays['error'][good],1000.,epoch,duration)
                if loc is not None:local[int(record['tic']),col]=loc
            results.append(dict(tic=int(record['tic']),ra=record['ra'],dec=record['dec'],tmag=record['tmag'],source=dict(path=record['path'],sha256=record['sha256'],url=record['url']),tests=tests))
        except Exception as exc:errors.append(dict(tic=int(record['tic']),error=repr(exc)))
    out=dict(tic=tic,sector=86,checked_utc=mf.utc(),driver_source_sha256=SOURCE_SHA,epoch_btjd=epoch,duration_days=duration,camera=4,ccd=4,selection='All 79 already-screened sector86 camera4 CCD4 observations, including target; selection frozen before viewing controls.',results=results,errors=errors,limitations=['Selected event and fixed time; not a survey false-alarm probability.','Stars can have independent real transits/variability at this time. Shared trends require assessment of depth, morphology and background, not simply a count.','Native local polynomial uses residual diagonal errors; correlated noise incomplete.'])
    snapshot=d/'method_snapshots'/(SOURCE_SHA+'.py');snapshot.parent.mkdir(exist_ok=True)
    if not snapshot.exists():snapshot.write_bytes(DRIVER_SOURCE)
    mf.save(d/'sector86_same_detector_event_check.json',out)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for ax,col in zip(axes,['corrected','raw']):
        vals=[r for r in results if r['tests'][col].get('depth')is not None]
        x=[float(r['tmag'])for r in vals];y=[r['tests'][col]['snr']for r in vals]
        ax.scatter(x,y,s=12,alpha=.7);target=next(r for r in vals if r['tic']==tic);ax.scatter(float(target['tmag']),target['tests'][col]['snr'],marker='*',s=120,color='r',label='Selected target');ax.axhline(0,color='.7');ax.set(xlabel='TESS magnitude',ylabel='Conditional single-event local S/N',title=col);ax.legend()
    fig.suptitle(f'BTJD {epoch:.5f}, S86 camera4 CCD4: same-time diagnostic, no calibrated FAP');fig.savefig(d/'sector86_same_detector_event_check.png',dpi=140);plt.close(fig)
    control_checks={};common_arrays={}
    fig,axes=plt.subplots(2,1,figsize=(12,7),layout='constrained')
    for ax,col in zip(axes,['corrected','raw']):
        target=local[tic,col];t=target['t'];controls=[];identities=[]
        for (other,c),loc in local.items():
            if c!=col or other==tic:continue
            residual=loc['f']-1;scale=1.4826*np.median(abs(residual-np.median(residual)))
            if not np.isfinite(scale)or scale<=0:continue
            sampled=np.interp(t,loc['t'],residual/scale,left=np.nan,right=np.nan)
            if np.isfinite(sampled).mean()<.9:continue
            controls.append(sampled);identities.append(other)
        median=np.nanmedian(controls,axis=0);tar=target['f']-1;scale=1.4826*np.median(abs(tar-np.median(tar)));inside=abs(t-epoch)<duration/2
        coef=np.polyfit(median[~inside],tar[~inside],1);res=tar-np.polyval(coef,median)
        err=target['e'];w=1/err[inside]**2;dep=-np.sum(w*res[inside])/w.sum();ee=1/np.sqrt(w.sum())
        control_checks[col]=dict(n_controls=len(controls),control_tics=identities,median_control_vs_target_correlation=float(np.corrcoef(median,tar)[0,1]),out_of_event_common_mode_slope=float(coef[0]),residual_single_event_depth=float(dep),residual_diagonal_error=float(ee),residual_diagonal_snr=float(dep/ee),warning='Same-detector median control is an exploratory explanatory trend, with target excluded. Its subtraction does not calibrate correlated-noise uncertainty or prove instrumentality.')
        common_arrays.update({col+'_time':t,col+'_target':tar,col+'_target_error':err,col+'_median_control_standardized':median,col+'_target_common_mode_residual':res})
        ax.plot(t-epoch,tar*1e6,'.',ms=2,label='Target local quadratic residual');ax.plot(t-epoch,np.polyval(coef,median)*1e6,lw=1,label='OOT-scaled median of other same-detector stars');ax.axvspan(-duration/2,duration/2,color='.7',alpha=.2);ax.set(xlabel='Days from selected event',ylabel=col+' (ppm)');ax.legend(fontsize=8)
    out['common_mode_experiment']=control_checks;mf.save(d/'sector86_same_detector_event_check.json',out);np.savez_compressed(d/'sector86_same_detector_common_mode.npz',**common_arrays)
    fig.suptitle('Selected single-event common-mode experiment; target excluded from comparison median');fig.savefig(d/'sector86_same_detector_common_mode.png',dpi=140);plt.close(fig)
    print('SAME_DETECTOR',len(results),len(errors),{col:sorted([(r['tic'],round(r['tests'][col]['snr'],2))for r in results],key=lambda x:x[1],reverse=True)[:12]for col in ['corrected','raw']},flush=True)
    return out

def neighbor_fast_checks():
    table=list(csv.DictReader(open(ROOT/'results/multisector/combined_signals.csv')))
    for tic,period in NEIGHBOR_TARGETS.items():
        rows=[r for r in table if int(r['tic'])==tic and abs(float(r['period_days'])-period)<1e-8]
        row=max(rows,key=lambda r:float(r['local_snr']))
        ephem=dict(period=float(row['period_days']),epoch_btjd=float(row['epoch_btjd']),duration_days=float(row['duration_days']))
        records=[r for r in mf.records_from_db(tic)if int(r['sector'])in[83,86]]
        d=mf.OUT/f'tic{tic}';d.mkdir(parents=True,exist_ok=True)
        out=dict(tic=tic,checked_utc=mf.utc(),driver_source_sha256=SOURCE_SHA,selected_blind_row=row,fixed_ephemeris=ephem,observations=[],limitations=['Fixed period/epoch from pooled selected data; no independent confirmation.','Known variable neighbors require pixels for source attribution; these are light-curve variability diagnostics.'])
        for record in records:
            v=mf.vet_observation(record,dict(period=period,duration_days=ephem['duration_days']),output=False,fixed_ephemeris=ephem)
            native=native_variability(record,ephem)
            out['observations'].append(dict(sector=int(record['sector']),vetting=v,native_variability=native))
            print('NEIGHBOR',tic,record['sector'],[(col,v['column_tests'][col]['primary']['snr'])for col in ['corrected','raw']],[(t['column'],t['top_peaks'][0])for t in native],flush=True)
        snapshot=d/'method_snapshots'/(SOURCE_SHA+'.py');snapshot.parent.mkdir(exist_ok=True)
        if not snapshot.exists():snapshot.write_bytes(DRIVER_SOURCE)
        mf.save(d/'pooled_north_neighbor.json',out)

def event_depth_summary(observations,column,original=None):
    events=[e for ob in observations if original is None or ob['discovery_sector']==original for e in ob['vetting']['column_tests'][column]['primary']['events']]
    if not events:return dict(n_events=0,depth=None,error=None,snr=None)
    ds=np.array([e['depth']for e in events]);es=np.array([e['depth_err']for e in events]);w=1/es**2
    mean=float(np.sum(w*ds)/sum(w));diagonal=float(1/np.sqrt(sum(w)));scatter=float(np.std(ds,ddof=1)/np.sqrt(len(ds)))if len(ds)>1 else diagonal;error=max(diagonal,scatter)
    return dict(n_events=len(events),n_positive_events=int(sum(ds>0)),depth=mean,error=error,snr=mean/error,diagonal_error=diagonal,event_scatter_error=scatter,constant_depth_diagonal_chi2=float(sum(((ds-mean)/es)**2)),constant_depth_dof=len(ds)-1,warning='Constant-depth chi2 and S/N use incomplete event-noise covariance and a selected ephemeris. No formal null probability; camera-dependent dilution/stellar/systematic variability remain.')

def summarize():
    rows=[];lines=['# Northern pooled-search follow-up','', 'These are fixed-ephemeris, conditional falsification tests, not calibrated significance, independent confirmation or planet claims. The exact period and epoch are selected once by the blind pooled search over sectors 83/86. Extra sectors retain that ephemeris, with period uncertainty and possible long-baseline drift.','', '| TIC | Period (d) | Sector | Original? | Corrected local depth (ppm) | Corrected/raw S/N | Positive/bracketed events | Warnings |','|---|---:|---:|---|---:|---|---|---|']
    for tic in TARGETS:
        path=mf.OUT/f'tic{tic}'/'pooled_north_validation.json'
        if not path.exists():continue
        obj=json.loads(path.read_text())
        for ob in obj['observations']:
            v=ob['vetting'];c=v['column_tests']['corrected']['primary'];r=v['column_tests']['raw']['primary']
            row=dict(tic=tic,period_days=obj['fixed_ephemeris']['period'],sector=ob['sector'],original_sector=ob['discovery_sector'],corrected_depth_ppm=(c.get('depth')or 0)*1e6,corrected_depth_error_ppm=(c.get('depth_err')or 0)*1e6,corrected_snr=c['snr'],raw_depth_ppm=(r.get('depth')or 0)*1e6,raw_snr=r['snr'],n_events=c['n_events'],n_positive_events=c.get('n_positive_events',0),warnings=';'.join(v['warnings']),evidence=str(path.relative_to(ROOT)))
            rows.append(row);lines.append(f'| {tic} | {row["period_days"]:.6f} | {ob["sector"]} | {ob["discovery_sector"]} | {row["corrected_depth_ppm"]:.0f} ± {row["corrected_depth_error_ppm"]:.0f} | {c["snr"]:.2f}/{r["snr"]:.2f} | {row["n_positive_events"]}/{row["n_events"]} | {row["warnings"]} |')
    lines.extend(['','## Event-depth consistency across sectors','','The following combines locally quadratic, individually bracketed event depths. The uncertainty is the larger of the diagonal error and event-depth standard deviation divided by sqrt(N). It remains conditional and does not calibrate red noise or search trials. Fresh-sector columns use only newly retrieved observations, with the original pooled ephemeris fixed.','', '| TIC | Original corrected/raw S/N | Fresh corrected/raw S/N | All corrected/raw S/N | Interpretation |','|---|---|---|---|---|'])
    notes={353807236:'Original corrected signal weakens in raw flux and fails to repeat in corrected fresh sectors; not a credible periodic transit.',359679661:'Residual weak lead: some positive fresh events, but very large event-depth variability and broad out-of-event harmonics. No validated planetary waveform.',353933809:'The 18.9784-day period fails to repeat. A stronger isolated S86 event remains unexplained; known Gaia binary and nearby star require source attribution.',366490373:'Residual weak lead: mostly positive raw event depths, but substantial fresh-sector scatter, broad phase structure and other-phase dips. No validated planetary waveform.',287086419:'Four original corrected events do not repeat in fresh S84. S81 lacks events satisfying the strict contiguous local-baseline requirement; median-window check is also weak/negative.',356821129:'Original positive events do not repeat in fresh S82/S84.',353783909:'Original positive events do not repeat in fresh S82/S84.'}
    for tic in TARGETS:
        path=mf.OUT/f'tic{tic}'/'pooled_north_validation.json'
        if not path.exists():continue
        obj=json.loads(path.read_text());depths={label:{col:event_depth_summary(obj['observations'],col,original)for col in ['corrected','raw']}for label,original in [('original',True),('fresh',False),('all',None)]}
        snapshot=path.parent/'method_snapshots'/(SOURCE_SHA+'.py');snapshot.parent.mkdir(exist_ok=True)
        if not snapshot.exists():snapshot.write_bytes(DRIVER_SOURCE)
        mf.save(path.parent/'pooled_north_event_depth_summary.json',dict(tic=tic,driver_source_sha256=SOURCE_SHA,fixed_ephemeris=obj['fixed_ephemeris'],tests=depths,assessment=notes[tic]))
        fmt=lambda key:'/'.join(f'{depths[key][col]["snr"]:.2f}'if depths[key][col]['snr']is not None else'NA'for col in ['corrected','raw'])
        lines.append(f'| {tic} | {fmt("original")} | {fmt("fresh")} | {fmt("all")} | {notes[tic]} |')
        (path.parent/'README.md').write_text(f'# TIC {tic}: northern pooled-search follow-up\n\n**{notes[tic]}**\n\nThe exact blind pooled-search ephemeris is period {obj["fixed_ephemeris"]["period"]:.9f} d, epoch BTJD {obj["fixed_ephemeris"]["epoch_btjd"]:.9f}, duration {24*obj["fixed_ephemeris"]["duration_days"]:.2f} h. It was held fixed across the original and extra observations.\n\nOriginal corrected/raw local-event S/N: {fmt("original")}; fresh: {fmt("fresh")}; all: {fmt("all")}. These use an event-scatter floor but remain conditional, without calibrated red-noise or survey-trial significance.\n\nSee [full measurements and raw provenance](pooled_north_validation.json), [event-depth consistency](pooled_north_event_depth_summary.json), per-sector fixed-ephemeris diagnostic plots, and native-cadence variability files. Missing/failed retrievals remain in [retrieval.json](retrieval.json).\n\nNo planet, source attribution or novelty claim. A nondetection at a fixed ephemeris is not an unconditional exclusion when period uncertainty, transit timing changes or dilution remain.\n')
    lines.extend(['','## Important limits','','- Original sectors 83 and 86 were both used to select the period and epoch. They are not independent recoveries. New sectors 81/82/84 are genuine additional data, evaluated without sector-specific retuning. Long-baseline period uncertainty can still cause phase drift.','- Generic variability flags are not automatic planet rejections. In particular, TIC359679661 and TIC366490373 retain some weak event support and are recorded as unresolved, variability-contaminated leads.','- The frozen eighth-order harmonic extrapolation can be unstable inside excluded event windows. Its occasionally huge residual scores are not used as evidence. Native-cadence out-of-event periodograms provide a separate descriptive variability check.','- TIC287086419 S82 returned HTTP404 and is preserved as unavailable. S81 was retrieved and verified but has zero strictly contiguous locally bracketed events. A wrapper preserves this result and median-window measurements because the frozen helper does not handle depth=None in its final flag guard. No helper source was changed.',''])
    common=mf.OUT/'tic353933809/sector86_same_detector_event_check.json'
    if common.exists():
        x=json.loads(common.read_text());tests=x.get('common_mode_experiment',{})
        if tests:
            lines.extend(['## Isolated event in TIC353933809','','The S86 dip at BTJD3637.21848 has local corrected/raw depths925/1007ppm, while the two S83 events are much shallower and fresh S82/S84 show no corresponding dips. Its proposed18.9784-day period therefore is not supported. A same-time check read and verified all79 retained S86 camera4/CCD4 products;77 comparison curves covered the target interval sufficiently for a target-excluded median control.','',f'After scaling that control outside the selected event, conditional corrected/raw event depths are {tests["corrected"]["residual_single_event_depth"]*1e6:.0f}/{tests["raw"]["residual_single_event_depth"]*1e6:.0f}ppm, S/N {tests["corrected"]["residual_diagonal_snr"]:.2f}/{tests["raw"]["residual_diagonal_snr"]:.2f}. This test does **not** establish a shared median instrumental cause. The event remains an unexplained single broad dip, not a planet candidate. Gaia recognition separately identifies a known stellar orbit and a nearby comoving source; source localization and the published binary geometry must be considered before any physical interpretation.',''])
    neighbor_lines=[]
    for tic in NEIGHBOR_TARGETS:
        path=mf.OUT/f'tic{tic}'/'pooled_north_neighbor.json'
        if not path.exists():continue
        x=json.loads(path.read_text())
        for ob in x['observations']:
            a=next(t for t in ob['native_variability']if t['column']=='raw')['top_peaks'][0]
            neighbor_lines.append(f'| {tic} | {ob["sector"]} | {x["fixed_ephemeris"]["period"]:.6f} | {a["period_days"]:.6f} | {a["ratio_to_orbital_frequency"]:.4f} | {a["semiamplitude_ppm"]:.0f} |')
    if neighbor_lines:lines.extend(['## Catalog-neighbor cases: native variability checks','','These use the fixed blind pooled period in both original sectors, plus native-cadence out-of-event LS. Source attribution remains unresolved without pixels. For the first two objects, the large coherent oscillation occurs at six times the selected transit frequency, explaining the dense BLS period aliases; a coexisting much weaker transit is not excluded.','','| TIC | Sector | Selected BLS P(d) | Dominant raw LS P(d) | LS frequency / selected orbital frequency | Raw semiamplitude(ppm) |','|---|---:|---:|---:|---:|---:|']+neighbor_lines+[''])
    out=ROOT/'reports/multisector';out.mkdir(parents=True,exist_ok=True)
    (out/'pooled_north.md').write_text('\n'.join(lines)+'\n')
    if rows:
        with (out/'pooled_north.csv').open('w',newline='')as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tic',type=int);p.add_argument('--summarize',action='store_true');p.add_argument('--instrumental-event',action='store_true');p.add_argument('--neighbor-cases',action='store_true');a=p.parse_args()
    if a.instrumental_event:instrumental_event_check()
    elif a.neighbor_cases:neighbor_fast_checks()
    elif not a.summarize:
        for tic in [a.tic]if a.tic else TARGETS:
            try:run(tic)
            except Exception as exc:
                mf.save(mf.OUT/f'tic{tic}'/'pooled_north_failure.json',dict(error=repr(exc),traceback=traceback.format_exc(),checked_utc=mf.utc()));print('TARGET_FAILED',tic,repr(exc),flush=True)
    summarize()
