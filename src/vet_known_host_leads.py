"""Independent gap/epoch/SAP vetting for the two strongest known-host pilot leads.

This reads existing FITS and screen JSON, and writes only its own candidate folder.
It never estimates a formal false-positive probability. Fixed-ephemeris tests are
retrospective because all sectors have already been screened.
"""
from pathlib import Path
import json
import math
import os
import sys
from datetime import datetime, timezone
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tess-known-host-vetting-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from search import ROOT, read_curve, flatten, known_for, mask_known, phase, scan

OUT = ROOT / 'candidates' / 'known_host_rejected'
TICS = [79748331, 354944123]
EDGE_BUFFER = .5


def safe(value):
    if isinstance(value, dict): return {str(k): safe(v) for k,v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)): return [safe(v) for v in value]
    if isinstance(value, np.generic): return safe(value.item())
    if isinstance(value, float) and not math.isfinite(value): return None
    return value


def segments(t):
    split = np.flatnonzero(np.diff(t) > .2)
    lo = np.r_[0, split + 1]; hi = np.r_[split, len(t)-1]
    return [(float(t[a]),float(t[b])) for a,b in zip(lo,hi)]


def edge_keep(t, bounds, buffer=EDGE_BUFFER):
    keep = np.zeros(len(t), bool)
    for lo,hi in bounds: keep |= (t >= lo+buffer) & (t <= hi-buffer)
    return keep


def fixed_depth(t,f,e,p,d,t0):
    if len(t)<10: return {'depth':None,'depth_err_diagonal':None,'snr_diagonal':None}
    inside=np.abs(phase(t,p,t0))<d/2
    if inside.sum()<2:
        return {'status':'fewer than two remaining in-transit bins','n_inside':int(inside.sum()),
                'depth':None,'depth_err_diagonal':None,'snr_diagonal':None}
    try:
        value=BoxLeastSquares(t,f,e).compute_stats(p,d,t0)['depth']
        return {'depth':float(value[0]),'depth_err_diagonal':float(value[1]),
                'snr_diagonal':float(value[0]/value[1]) if value[1]>0 else None}
    except Exception as exc:
        return {'error':str(exc),'depth':None,'depth_err_diagonal':None,'snr_diagonal':None}


def event_measure(t,f,e,center,duration,mask,bounds,known):
    inside=np.abs(t-center)<duration/2
    local=(np.abs(t-center)>duration)&(np.abs(t-center)<duration*4)
    selected=inside&mask; baseline=local&mask
    depth=error=None
    if selected.sum()>=2 and baseline.sum()>=3:
        wi=1/e[selected]**2; wo=1/e[baseline]**2
        depth=float(np.sum(wo*f[baseline])/wo.sum()-np.sum(wi*f[selected])/wi.sum())
        error=float(np.sqrt(1/wi.sum()+1/wo.sum()))
    edges=np.array(bounds).reshape(-1)
    offsets=[{'planet':k['name'],'phase_offset_days':float(phase(center,k['period'],k['epoch_bjd']-2457000)),
              'mask_halfwidth_days':max(k['duration_hours']/24*.85,.05)} for k in known]
    return {'center_btjd':float(center),'n_inside_before_mask':int(inside.sum()),'n_inside_after_mask':int(selected.sum()),
            'n_baseline':int(baseline.sum()),'depth_local':depth,'depth_err_diagonal':error,
            'snr_local_diagonal':depth/error if depth is not None and error else None,
            'nearest_usable_segment_edge_days':float(np.min(np.abs(edges-center))),
            'inside_edge_buffer':bool(np.min(np.abs(edges-center))<EDGE_BUFFER),
            'known_phase_offsets':offsets}


def predictions(t,f,e,period,duration,t0,mask,bounds,known):
    cycles=range(int(np.floor((t.min()-t0)/period))-1,int(np.ceil((t.max()-t0)/period))+2)
    return [event_measure(t,f,e,t0+cycle*period,duration,mask,bounds,known)
            for cycle in cycles if t.min()-duration/2 <= t0+cycle*period <= t.max()+duration/2]


def plot_sector(tic, sector, lead, curves,bounds,known,result):
    fig,axes=plt.subplots(4,1,figsize=(12,12),layout='constrained')
    color={'PDCSAP_FLUX':'tab:blue','SAP_FLUX':'tab:orange'}
    for column in ('PDCSAP_FLUX','SAP_FLUX'):
        t,f,e,k,trim=curves[(1.0,column)]
        axes[0].plot(t,f,'.',ms=1,label=column,color=color[column],alpha=.6)
        axes[2].plot(phase(t[k],lead['period'],lead['t0_btjd']),f[k],'.',ms=2,label=column,color=color[column],alpha=.6)
    t,f,e,k,trim=curves[(1.0,'PDCSAP_FLUX')]
    for lo,hi in bounds:
        axes[0].axvspan(lo,min(lo+EDGE_BUFFER,hi),color='tab:red',alpha=.14)
        axes[0].axvspan(max(hi-EDGE_BUFFER,lo),hi,color='tab:red',alpha=.14)
    strongest=max(lead['events'],key=lambda row:row['depth'])['time'] if lead['events'] else lead['t0_btjd']
    for column in ('PDCSAP_FLUX','SAP_FLUX'):
        tx,fx,ex,kx,tr=curves[(1.0,column)]
        close=np.abs(tx-strongest)<.8
        axes[1].plot(tx[close]-strongest,fx[close],'.',ms=2,color=color[column],label=column)
    for edge in np.array(bounds).ravel():
        if abs(edge-strongest)<.8: axes[1].axvline(edge-strongest,color='red',lw=1)
    axes[1].axvspan(-lead['duration_days']/2,lead['duration_days']/2,color='.8',alpha=.5)
    epochs=predictions(t,f,e,lead['period'],lead['duration_days'],lead['t0_btjd'],k,bounds,known)
    for event in epochs:
        if event['depth_local'] is not None:
            axes[3].errorbar(event['center_btjd'],event['depth_local']*1e6,yerr=event['depth_err_diagonal']*1e6,
                             fmt='o',color='tab:red' if event['inside_edge_buffer'] else 'black')
    axes[0].set(xlabel='BTJD (BJD_TDB−2457000)',ylabel='Flattened flux',title=f'TIC {tic}, sector {sector}; shaded usable segment edges ±0.5d')
    axes[0].legend(markerscale=5)
    axes[1].set(xlabel=f'Days from strongest screened event ({strongest:.6f} BTJD)',ylabel='Flattened flux',title='Raw SAP and PDCSAP independently flattened; red lines are usable-segment boundaries')
    axes[1].legend(markerscale=4)
    axes[2].set(xlabel='Phase time (days)',ylabel='Flattened flux',xlim=(-lead['duration_days']*3,lead['duration_days']*3),title=f"Original P={lead['period']:.6f}d; diagonal S/N={lead['snr']:.1f}")
    axes[3].set(xlabel='Predicted event BTJD',ylabel='Local depth (ppm)',title='Individual predicted events; red = within 0.5d of usable edge; errors assume independent bins')
    axes[3].axhline(0,color='.6',lw=1)
    path=OUT/f'tic{tic}_s{sector:04d}_boundary_diagnostics.png';fig.savefig(path,dpi=140);plt.close(fig)
    return str(path.relative_to(ROOT))


def check_trimmed_alias(prepared):
    """Follow the new sector-102 maximum and its two-event doubled-period alias."""
    source=prepared[(354944123,102)]['result']
    lead=next(row for row in source['trimmed_scans'] if row['width']==1.0)['peaks'][0]
    p,t0,d=lead['period'],lead['t0_btjd'],lead['duration_days']
    checks=[]
    for (tic,sector),prep in prepared.items():
        if tic!=354944123:continue
        for col in ('PDCSAP_FLUX','SAP_FLUX'):
            for width in (1.0,2.0):
                t,f,e,k,tr=prep['curves'][(width,col)]
                for period in (p,2*p):
                    events=predictions(t,f,e,period,d,t0,k,prep['bounds'],prep['known'])
                    # Explore a declared timing tolerance; maxima are screening
                    # statistics and this extra search increases the trial count.
                    for event in events:
                        grid=[]
                        for shift in np.linspace(-.3,.3,31):
                            test=event_measure(t,f,e,event['center_btjd']+shift,d,k,prep['bounds'],prep['known'])
                            if test['snr_local_diagonal'] is not None:
                                grid.append({'shift_days':float(shift),'depth':test['depth_local'],
                                             'snr_local_diagonal':test['snr_local_diagonal']})
                        event['timing_tolerance_search_days']=.3
                        event['timing_tolerance_maximum']=max(grid,key=lambda test:test['snr_local_diagonal']) if grid else None
                    checks.append({'sector':sector,'column':col,'width':width,'period':period,
                        't0_btjd':t0,'duration_days':d,
                        'fixed_edge_trimmed':fixed_depth(t[tr],f[tr],e[tr],period,d,t0),'events':events})
    record={'tic':354944123,'source_sector':102,'source_trimmed_peak':lead,'checks':checks,
        'assessment':'Reduction-dependent two-event alignment. Opposite-sign SAP events and absent fixed-ephemeris repetition prevent promotion.',
        'limitations':['Both 6.671d and 13.342d ephemerides share two observed events; the intervening short-period event is in a gap.',
            'The original two-event period has not received a full posterior fit; null fixed times alone do not exclude every alias.',
            'The +/-0.3d timing search is a sensitivity check, not a formal confidence interval or false-alarm probability.',
            'SAP could contain strong additive backgrounds; opposite signs justify independent pixel reduction if this lead is revisited.']}
    (OUT/'toi4342_trimmed_alias_check.json').write_text(json.dumps(safe(record),indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(3,2,figsize=(12,11),layout='constrained')
    for index,center in enumerate([t0,t0+2*p]):
        prep=prepared[(354944123,102)]
        for col,color in [('PDCSAP_FLUX','tab:blue'),('SAP_FLUX','tab:orange')]:
            t,f,e,k,tr=prep['curves'][(1.0,col)];close=np.abs(t-center)<.65
            axes[0,index].plot(t[close]-center,f[close],'.',ms=2,color=color,label=col)
            rawt,rawf,_,_=read_curve(prep['result']['source_fits'],col)
            local=np.abs(rawt-center)<.65
            normalization=np.median(rawf[local&(np.abs(rawt-center)>d)])
            axes[1,index].plot(rawt[local]-center,rawf[local]/normalization,'.',ms=1,color=color,label=col)
        axes[0,index].set(title=f'Sector 102, event BTJD {center:.5f}',ylabel='Independently flattened flux',xlabel='Days from predicted center')
        axes[1,index].set(ylabel='Raw local-normalized flux',xlabel='Days from predicted center')
        for ax in [axes[0,index],axes[1,index]]:ax.axvspan(-d/2,d/2,color='.8',alpha=.4);ax.legend(markerscale=4)
    for index,sector in enumerate([101,103]):
        prep=prepared[(354944123,sector)]
        t,f,e,k,tr=prep['curves'][(1.0,'PDCSAP_FLUX')]
        ph=phase(t,2*p,t0);close=np.abs(ph)<.6
        axes[2,index].plot(ph[close&k],f[close&k],'.',ms=2,color='tab:blue',label='PDCSAP')
        axes[2,index].axvspan(-d/2,d/2,color='.8',alpha=.4)
        axes[2,index].set(title=f'Sector {sector}: fixed 13.342d predicted phase',xlabel='Phase time (days)',ylabel='Flattened flux')
    fig.savefig(OUT/'toi4342_trimmed_alias_shapes.png',dpi=140);plt.close(fig)
    lines=['# TOI-4342: follow-up of the edge-trimmed two-event alignment','',
        f'After boundary rejection, sector 102 gives P={p:.8f} d, epoch BTJD {t0:.6f}, duration {d:.3f} d. The observed events are separated by {2*p:.6f} d. The intervening 6.671d prediction lies in the large quality gap, so 13.342d is equally compatible with these two event times. This specific signal remains unpromoted.', '',
        '| Sector | Flux | Width (d) | Fixed 13.342d diagonal S/N |','| --- | --- | --- | --- |']
    for row in checks:
        if row['period']>10:
            snr=row['fixed_edge_trimmed'].get('snr_diagonal')
            lines.append(f"| {row['sector']} | {row['column']} | {row['width']:g} | {snr:.2f} |")
    lines+=['', 'In the source sector, the 1d-detrended PDC local depths are approximately 1055 and 1433 ppm; corresponding SAP depths are approximately −228 and −401 ppm. The same ephemeris has no significant PDC repetition at fixed times in neighboring sectors. Exact event measurements and a ±0.3d timing-sensitivity search are saved in `toi4342_trimmed_alias_check.json`. Extra timing searches introduce extra trials; their maxima are not discovery significance.', '',
        '[Event-shape comparison](toi4342_trimmed_alias_shapes.png) includes both original raw flux columns, independently normalized, and both independently flattened columns. Opposite signs can be caused by additive instrumental backgrounds, so this is evidence of reduction dependence rather than a theorem excluding all transits. Pixel photometry and longer-baseline alias fitting would be required if future observations revive the signal. No current evidence justifies calling it a credible new candidate.', '',
        'The two-event period uncertainty was not fitted. Fixed-time null results alone do not exclude every nearby period or timing variation. No unpublished/discovered status, formal false-positive probability, or planetary radius is assigned.']
    (OUT/'TOI4342_TRIMMED_ALIAS.md').write_text('\n'.join(lines)+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records=json.loads((ROOT/'data'/'observations.json').read_text())
    all_results=[];prepared={}
    for rec in records:
        if rec['tic'] not in TICS: continue
        tic,sector=rec['tic'],rec['sector'];known=known_for(tic)
        with fits.open(rec['local_path']) as hdul:
            metadata={k:hdul[1].header.get(k) for k in ('TIMESYS','BJDREFI','BJDREFF','CROWDSAP','FLFRCSAP','TSTART','TSTOP')}
            metadata.update({k:hdul[0].header.get(k) for k in ('TICID','SECTOR','TEFF','RADIUS','TESSMAG','PROCVER')})
            assert metadata['TIMESYS']=='TDB' and metadata['BJDREFI']==2457000 and metadata['BJDREFF']==0
        original={}
        for width in (1.0,2.0):
            p=ROOT/'results'/'scans'/f'tic{tic}_s{sector:04}_w{width:g}_PDCSAP_FLUX_residual.json'
            if p.exists(): original[width]=json.loads(p.read_text())
        curves={};rawtimes=None
        for column in ('PDCSAP_FLUX','SAP_FLUX'):
            t,f,e,_=read_curve(rec['local_path'],column)
            bounds=segments(t)
            if column=='PDCSAP_FLUX': rawtimes=t.copy();reference_bounds=bounds
            for width in (1.0,2.0):
                tb,fb,eb=flatten(t,f,e,width)
                k=mask_known(tb,known);trim=k&edge_keep(tb,bounds)
                curves[(width,column)]=(tb,fb,eb,k,trim)
                np.savez_compressed(OUT/f'tic{tic}_s{sector:04d}_w{width:g}_{column}_vetting.npz',
                                    t=tb,f=fb,e=eb,known_mask=k,edge_trimmed_mask=trim)
        result={'tic':tic,'sector':sector,'source_fits':rec['local_path'],'source_sha256':rec['sha256'],
                'metadata':metadata,'known_planets':known,'usable_segments_btjd':reference_bounds,
                'edge_buffer_days':EDGE_BUFFER,'leads':[],'trimmed_scans':[],'known_ephemeris_checks':[]}
        for width,screen in original.items():
            tb,fb,eb,k,trim=curves[(width,'PDCSAP_FLUX')]
            peaks,_=scan(tb[trim],fb[trim],eb[trim],n_peaks=3)
            result['trimmed_scans'].append({'width':width,'n_before_trim':int(k.sum()),'n_after_trim':int(trim.sum()),'peaks':peaks})
            for lead in screen['peaks']:
                entry={'original_width':width,'screen':lead,'tests':[]}
                for w in (1.0,2.0):
                    for col in ('PDCSAP_FLUX','SAP_FLUX'):
                        tx,fx,ex,kx,tr=curves[(w,col)]
                        entry['tests'].append({'width':w,'column':col,
                            'fixed_ephemeris_untrimmed':fixed_depth(tx[kx],fx[kx],ex[kx],lead['period'],lead['duration_days'],lead['t0_btjd']),
                            'fixed_ephemeris_edge_trimmed':fixed_depth(tx[tr],fx[tr],ex[tr],lead['period'],lead['duration_days'],lead['t0_btjd']),
                            'events':predictions(tx,fx,ex,lead['period'],lead['duration_days'],lead['t0_btjd'],kx,reference_bounds,known)})
                result['leads'].append(entry)
        # Known ephemeris audit uses full flux after edge trimming, retaining known planets.
        tb,fb,eb,k,trim=curves[(1.0,'PDCSAP_FLUX')];onlyedge=edge_keep(tb,reference_bounds)
        for planet in known:
            d=planet['duration_hours']/24;epoch=planet['epoch_bjd']-2457000
            given=fixed_depth(tb[onlyedge],fb[onlyedge],eb[onlyedge],planet['period'],d,epoch)
            phasefit=BoxLeastSquares(tb[onlyedge],fb[onlyedge],eb[onlyedge]).power([planet['period']],[d],oversample=20,objective='snr')
            result['known_ephemeris_checks'].append({'name':planet['name'],'catalog_fixed':given,
                'free_phase_depth':float(phasefit.depth[0]),'free_phase_snr':float(phasefit.depth_snr[0]),
                'free_phase_t0_btjd':float(phasefit.transit_time[0]),
                'free_phase_offset_from_catalog_days':float(phase(float(phasefit.transit_time[0]),planet['period'],epoch)),
                'warning':'Free phase can select another signal or correlated noise; inspect consistency across sectors.'})
        if original and original.get(1.0,{}).get('peaks'):
            result['plot']=plot_sector(tic,sector,original[1.0]['peaks'][0],curves,reference_bounds,known,result)
        prepared[(tic,sector)]={'curves':curves,'bounds':reference_bounds,'known':known,'result':result}
        all_results.append(result)
        print('Vetted',tic,sector,'original',original.get(1.0,{}).get('peaks',[{}])[0].get('snr'),
              'trimmed',[(x['width'],[(round(y['period'],3),round(y['snr'],2)) for y in x['peaks']]) for x in result['trimmed_scans']],flush=True)
        (OUT/f'tic{tic}_s{sector:04d}.json').write_text(json.dumps(safe(result),indent=2,allow_nan=False)+'\n')
    # Freeze each original ephemeris and measure all other already-screened sectors.
    crosschecks=[]
    for source in all_results:
        for entry in source['leads']:
            lead=entry['screen'];tests=[]
            for (tic,sector),prep in prepared.items():
                if tic!=source['tic'] or sector==source['sector']:continue
                t,f,e,k,tr=prep['curves'][(1.0,'PDCSAP_FLUX')]
                events=predictions(t,f,e,lead['period'],lead['duration_days'],lead['t0_btjd'],k,prep['bounds'],prep['known'])
                tests.append({'sector':sector,'fixed_ephemeris_untrimmed':fixed_depth(t[k],f[k],e[k],lead['period'],lead['duration_days'],lead['t0_btjd']),
                              'fixed_ephemeris_trimmed':fixed_depth(t[tr],f[tr],e[tr],lead['period'],lead['duration_days'],lead['t0_btjd']),
                              'predicted_events':events})
            crosschecks.append({'tic':source['tic'],'source_sector':source['sector'],'width':entry['original_width'],
                                'period':lead['period'],'t0_btjd':lead['t0_btjd'],'duration_days':lead['duration_days'],'tests':tests})
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),'sectors':all_results,'retrospective_fixed_ephemeris_checks':crosschecks,
             'interpretation':'Screening scores require event-level falsification; edge trimming is a diagnostic, not a blanket proof that astrophysical transits cannot occur near gaps.',
             'statistics_limit':'All S/N values assume diagonal errors. No global false-alarm probability or confirmation is inferred.',
             'data_preservation':'Raw FITS remain unchanged; normalized variants and masks are retained beside these reports.'}
    (OUT/'known_host_vetting.json').write_text(json.dumps(safe(summary),indent=2,allow_nan=False)+'\n')
    check_trimmed_alias(prepared)
    lines=['# Known-host pilot: independent residual vetting','','Analysis date: '+summary['created_utc'], '',
      'The strongest residuals on TOI-1064 (TIC 79748331) and TOI-4342 (TIC 354944123) were re-examined in raw SAP and PDCSAP, both median-detrending widths, individual usable segments, and other observed sectors. The prominent long-period screen peaks are dominated by different quality-gap or sector-boundary excursions. These residual periods are not supported as new planets.', '',
      'Time units were verified directly from every FITS: `TIMESYS=TDB`, `BJDREFI=2457000`, `BJDREFF=0`. Catalogue BJD_TDB epochs were converted to BTJD once. `CROWDSAP` was read from the light-curve extension; high values here are aperture contamination estimates, not blend exclusion.', '',
      '| TIC | Sector | Original top period (d) | Original diagonal S/N | Top S/N after 0.5d gap/edge trim, width 1d | CROWDSAP |',
      '| --- | --- | --- | --- | --- | --- |']
    for row in all_results:
        leads=[entry for entry in row['leads'] if entry['original_width']==1.0]
        peak=leads[0]['screen'] if leads else None
        tr=next((entry for entry in row['trimmed_scans'] if entry['width']==1.0),{})
        after=tr.get('peaks',[{}])[0].get('snr') if tr.get('peaks') else None
        if peak:lines.append(f"| {row['tic']} | {row['sector']} | {peak['period']:.6f} | {peak['snr']:.2f} | {after:.2f} | {row['metadata']['CROWDSAP']:.4f} |")
    lines += ['', '## Event-level evidence', '']
    for row in all_results:
        top=next((entry for entry in row['leads'] if entry['original_width']==1.0),None)
        if not top:continue
        first=top['tests'][0];events=[e for e in first['events'] if e['depth_local'] is not None]
        strongest=max(events,key=lambda event:event['depth_local']) if events else None
        if strongest:
            lines.append(f"- TIC {row['tic']}, sector {row['sector']}: deepest predicted event at BTJD {strongest['center_btjd']:.6f}, local depth {strongest['depth_local']*1e6:.0f} ppm, nearest usable-segment boundary {strongest['nearest_usable_segment_edge_days']:.4f} d. Fixed original-ephemeris S/N changes from {first['fixed_ephemeris_untrimmed'].get('snr_diagonal')} to {first['fixed_ephemeris_edge_trimmed'].get('snr_diagonal')} after trimming. [Diagnostic plot]({Path(row['plot']).name}).")
    lines += ['', '## Scope and limitations', '',
      'Each original period/epoch is propagated without retuning to the other sectors in `known_host_vetting.json`. Individual predicted events, masks, local depths and uncertainties are retained, including missing predictions. This is a retrospective consistency check; the other sectors were already searched and do not constitute independent prospective validation. Per-sector periods do not justify an ephemeris uncertainty small enough for a distant-season null result to exclude every possible period alias.', '',
      'The 0.5d exclusion applies to every contiguous good-quality segment, including 0.217d gaps. It tests the hypothesis that instrumental edge excursions dominate the BLS score. True transits can occur near gaps; trimming sensitivity is reduced, and this diagnostic rejects the promoted residual interpretation rather than proving either star lacks additional planets.', '',
      'Odd/even statistics with only one supported event are undefined in scientific practice even when the screen reports zero. The TOI-4342 sector-102 ~9.379d interpretation also alternates deep events with an essentially absent middle event; inspect individual events rather than counting an available time window as a transit. Exact figures for both detrending widths are preserved in JSON.', '',
      'No formal transit fit, localization, or false-positive probability was justified for these rejected artifacts. Known ephemerides have also been audited with fixed-phase and freely shifted-phase depths; a free-phase detection is not automatically timing variation, and must be consistent across sectors. TOI-4342 has a previously reported ~47.5d RV candidate according to the parallel literature review; a future long-period transit search must treat that as prior recognition.', '',
      'Reproduce with `python src/vet_known_host_leads.py` after downloading the same FITS products listed by SHA256 in each sector JSON. The script writes only this folder and does not modify the archive manifest or screen implementation. Raw data are retained in the project cache.']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
