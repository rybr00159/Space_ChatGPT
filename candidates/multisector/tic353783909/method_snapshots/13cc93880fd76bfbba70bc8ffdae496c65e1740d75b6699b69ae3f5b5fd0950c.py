"""Fixed pooled-ephemeris falsification of seven northern TESS search leads.

Imports frozen scientific helpers read-only. Extra sectors are genuine new
observations but the pooled search ephemeris has selection/period uncertainty.
"""
from pathlib import Path
import sys, json, csv, argparse, traceback, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import multisector_followup as mf
from astropy.timeseries import LombScargle
from scipy.ndimage import median_filter
TARGETS={353807236:6.224717707910657,359679661:1.9134095518192817,353933809:18.978402288292596,366490373:2.810238855348628,287086419:12.613387995535525,356821129:3.871952352751222,353783909:6.397423709551385}
SOURCE_SHA=mf.digest(__file__)

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
    output=dict(tic=tic,checked_utc=mf.utc(),driver_source_sha256=SOURCE_SHA,pooled_source=dict(path=str(pooledpath),sha256=mf.digest(pooledpath)),selected_blind_row=row,fixed_ephemeris=ephem,original_records=records,observations=[],errors=[],limitations=['All original sectors participated in the blind pooled selection.','Extra sectors use the exact pooled ephemeris without sector retuning. Long-baseline period uncertainty can cause phase drift; a nondetection is not in itself exclusion.','Conditional scores, raw/corrected variants and multiple detrenders are not independent false-alarm probabilities.'])
    mf.save(d/'pooled_north_input.json',output)
    retrieval=mf.retrieve_independent(records[0],max_sectors=2,exclude_sectors=[83,86]);output['retrieval']=retrieval
    extra=[r for r in retrieval['observations']if r.get('status')=='retained' and int(r['sector'])not in[83,86]]
    for record in records+extra:
        try:
            result=mf.vet_observation(record,peak,fixed_ephemeris=ephem)
            variability=mf.stellar_variability_check(record,ephem)
            native=native_variability(record,ephem)
            mf.save(d/f'sector{record["sector"]}_native_variability.json',dict(tic=tic,sector=record['sector'],source=record,ephemeris=ephem,driver_source_sha256=SOURCE_SHA,tests=native))
            output['observations'].append(dict(sector=int(record['sector']),discovery_sector=int(record['sector'])in[83,86],vetting=result,stellar_variability=variability,native_variability=native))
        except Exception as exc:
            output['errors'].append(dict(sector=int(record['sector']),error=repr(exc),traceback=traceback.format_exc()));print('FAILED',tic,record['sector'],repr(exc),flush=True)
        mf.save(d/'pooled_north_validation.json',output)
    return output

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
    out=ROOT/'reports/multisector';out.mkdir(parents=True,exist_ok=True)
    (out/'pooled_north.md').write_text('\n'.join(lines)+'\n')
    if rows:
        with (out/'pooled_north.csv').open('w',newline='')as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tic',type=int);p.add_argument('--summarize',action='store_true');a=p.parse_args()
    if not a.summarize:
        for tic in [a.tic]if a.tic else TARGETS:
            try:run(tic)
            except Exception as exc:
                mf.save(mf.OUT/f'tic{tic}'/'pooled_north_failure.json',dict(error=repr(exc),traceback=traceback.format_exc(),checked_utc=mf.utc()));print('TARGET_FAILED',tic,repr(exc),flush=True)
    summarize()
