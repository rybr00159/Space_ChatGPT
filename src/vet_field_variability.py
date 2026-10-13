"""Gap-aware Fourier/box vetting of two recurrent field-search signals.

This is a falsification experiment, not a calibrated detection/absence test.
The first two sectors choose a sinusoidal frequency and harmonic complexity;
the final two sectors test that frozen waveform, and a separate per-sector
amplitude refit supplies an exploratory residual search.  Results are local.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import numpy as np
from astropy.timeseries import LombScargle, BoxLeastSquares
from scipy.optimize import minimize_scalar
from scipy.ndimage import median_filter
from search import ROOT, CACHE, read_curve, phase
import matplotlib.pyplot as plt

OUT = ROOT / 'candidates' / 'field_variability'
DURATIONS = np.array([.01, .02, .035, .05])


def prepare(record, width=5, column='PDCSAP_FLUX', injection=None):
    t, f, e, _ = read_curve(record['local_path'], column)
    if injection is not None:
        ip,it0,iduration,idepth=injection
        inside=np.abs(phase(t,ip,it0))<iduration/2
        f=f*(1-idepth*inside)
    edges = np.r_[0, np.flatnonzero(np.diff(t) > .2) + 1, len(t)]
    fn, en = np.empty_like(f), np.empty_like(e)
    for a, b in zip(edges[:-1], edges[1:]):
        cadence = np.median(np.diff(t[a:b])) if b-a>2 else .0014
        size=max(3,int(width/cadence)//2*2+1)
        trend=median_filter(f[a:b],size=size,mode='reflect')
        fn[a:b],en[a:b]=f[a:b]/trend,e[a:b]/trend
    # Keep photon uncertainties here. An empirical residual-noise floor is
    # applied after the variable-star model, not before its removal.
    idx=np.floor((t-t.min())/(5/1440)).astype(int)
    _,inv=np.unique(idx,return_inverse=True)
    weight=1/en**2; total=np.bincount(inv,weights=weight)
    tb=np.bincount(inv,weights=weight*t)/total
    fb=np.bincount(inv,weights=weight*fn)/total
    eb=1/np.sqrt(total)
    keep = np.zeros(len(tb), bool)
    for a, b in zip(edges[:-1], edges[1:]):
        keep |= (tb >= t[a] + .25) & (tb <= t[b-1] - .25)
    scatter=1.4826*np.median(np.abs(fb-np.median(fb)))
    keep &= fb < np.median(fb)+7*scatter
    return tb[keep], fb[keep], eb[keep]


def design(t, frequency, harmonics, origin):
    angle = 2*np.pi*frequency*(t-origin)
    return np.column_stack([np.ones(len(t))] +
        [v for k in range(1, harmonics+1) for v in (np.sin(k*angle), np.cos(k*angle))])


def fit(t, f, e, frequency, harmonics, origin):
    matrix = design(t, frequency, harmonics, origin)
    coeff = np.linalg.lstsq(matrix/e[:, None], f/e, rcond=None)[0]
    model = matrix @ coeff
    chi2 = np.sum(((f-model)/e)**2)
    return coeff, model, float(chi2), float(chi2 + (2*harmonics+1)*np.log(len(t)))


def spectrum(t, f, e):
    frequency, power = LombScargle(t, f, e).autopower(
        minimum_frequency=1/3, maximum_frequency=10, samples_per_peak=12)
    j = np.argmax(power)
    return frequency, power, int(j)


def scan_short(t, f, e, maximum=.7, minimum=.1, n=5):
    # A log grid limits transit-time drift to <= shortest duration/3 over the
    # full actual baseline. Unlike uniformly indexed FFTs, gaps remain gaps.
    baseline = t.max()-t.min()
    count = int(np.ceil(np.log(maximum/minimum)*3*baseline/DURATIONS.min())) + 1
    periods = np.geomspace(minimum, maximum, count)
    bls = BoxLeastSquares(t, f, e)
    result = bls.power(periods, DURATIONS, objective='snr', oversample=10)
    available = np.ones(len(periods), bool)
    peaks = []
    for _ in range(n):
        j = int(np.argmax(np.where(available, result.power, -np.inf)))
        p, t0, d = float(result.period[j]), float(result.transit_time[j]), float(result.duration[j])
        # Exact slow evaluation guards fast-method approximation errors,
        # particularly when a .05d trial reaches a large duty cycle.
        exact=bls.power([p],DURATIONS,objective='snr',oversample=10,method='slow')
        t0,d=float(exact.transit_time[0]),float(exact.duration[0])
        stats = bls.compute_stats(p, d, t0)
        odd, even = stats['depth_odd'], stats['depth_even']
        peaks.append(dict(period_days=p, epoch_btjd=t0, duration_days=d,
            depth_ppm=float(exact.depth[0]*1e6), depth_err_ppm=float(exact.depth_err[0]*1e6),
            diagonal_snr=float(exact.depth_snr[0]),fast_grid_snr=float(result.depth_snr[j]),
            odd_even_sigma=float(abs(odd[0]-even[0])/np.hypot(odd[1], even[1])),
            grid_count=count))
        available &= np.abs(periods/p-1)>.01
    return peaks, result


def fixed_box(t, f, e, peak):
    p, t0, d = (peak[k] for k in ['period_days','epoch_btjd','duration_days'])
    inside = np.abs(phase(t, p, t0)) < d/2
    if not inside.any() or inside.all(): return None
    w = 1/e**2
    fi, fo = np.average(f[inside], weights=w[inside]), np.average(f[~inside], weights=w[~inside])
    error = np.sqrt(1/w[inside].sum()+1/w[~inside].sum())
    return dict(depth_ppm=float((fo-fi)*1e6), depth_error_ppm=float(error*1e6),
                diagonal_snr=float((fo-fi)/error), n_in=int(inside.sum()))


def folded(t, f, p, origin, bins=80):
    phi = ((t-origin)/p) % 1
    edges=np.linspace(0,1,bins+1)
    values=[]; errors=[]
    for a,b in zip(edges[:-1],edges[1:]):
        x=f[(phi>=a)&(phi<b)]
        values.append(np.mean(x) if len(x) else np.nan)
        errors.append(np.std(x)/np.sqrt(len(x)) if len(x)>1 else np.nan)
    return .5*(edges[:-1]+edges[1:]), np.array(values), np.array(errors)


def run(tic, records, column='PDCSAP_FLUX', width=5):
    chunks=[prepare(r, width, column) for r in records]
    train=[c for r,c in zip(records,chunks) if r['sector'] in [95,96]]
    test=[c for r,c in zip(records,chunks) if r['sector'] in [97,98]]
    tt, tf, te = (np.concatenate([c[j] for c in train]) for j in range(3))
    origin=float(tt.min())
    freq,power,j=spectrum(tt,tf,te)
    # TIC348839788 has instrumental/slow variability in sector96; constrain
    # to the empirically recurring <0.7d frequency when selecting that family.
    if tic==348839788:
        j=np.argmax(np.where(freq>1/.7,power,-np.inf))
    grid_step=freq[1]-freq[0]
    opt=minimize_scalar(lambda x: fit(tt,tf,te,x,1,origin)[2],
                        bounds=(freq[j]-2*grid_step,freq[j]+2*grid_step),
                        method='bounded',options={'xatol':1e-11})
    frequency=float(opt.x); period=1/frequency
    choices=[]
    for h in range(1,5):
        coeff,model,chi2,bic=fit(tt,tf,te,frequency,h,origin)
        choices.append(dict(harmonics=h,chi2=chi2,bic_diagonal=bic,coefficients=coeff.tolist()))
    h=min(choices,key=lambda x:x['bic_diagonal'])['harmonics']
    coefficients=np.array(choices[h-1]['coefficients'])
    summary=dict(tic=tic,column=column,width_days=width,created_utc=datetime.now(timezone.utc).isoformat(),
        train_sectors=[95,96],test_sectors=[97,98],frequency_per_day=frequency,period_days=period,
        origin_btjd=origin,selected_harmonics=h,harmonic_choices=choices,
        cautions=['All S/N and BIC values use diagonal errors and are descriptive, not false-alarm probabilities.',
                  'Median detrending, 5-minute bins and harmonic removal reduce transit completeness.',
                  'Positive-only 7-MAD clipping before and after Fourier modeling limits flare leverage; negative excursions remain.',
                  'Per-sector amplitude refits are exploratory and may absorb coincident harmonic transits.',
                  'No claim of absolute absence, planetary origin, novelty, or stellar source localization.'],
        observations=[dict(tic=r['tic'],sector=r['sector'],data_uri=r['data_uri'],sha256=r.get('sha256'),
                           metadata=r.get('metadata')) for r in records],sectors=[])
    fig,axes=plt.subplots(4,3,figsize=(15,13),layout='constrained')
    residual_chunks=[]
    for row,c,ax in zip(records,chunks,axes):
        t,f,e=c; sf,sp,sj=spectrum(t,f,e)
        coeff,model,chi2,bic=fit(t,f,e,frequency,h,origin)
        one_coeff,one_model,one_chi2,_=fit(t,f,e,frequency,1,origin)
        frozen=design(t,frequency,h,origin)@coefficients
        # Robust empirical 5-min residual scatter caps ranking significance.
        residual=f-model+1
        sigma=1.4826*np.median(abs(residual-np.median(residual)))
        e_search=e*max(1,sigma/np.median(e))
        keep=residual<1+7*sigma
        peaks,bls=scan_short(t[keep],residual[keep],e_search[keep])
        frozen_residual=f-frozen+1
        frozen_sigma=1.4826*np.median(abs(frozen_residual-np.median(frozen_residual)))
        frozen_e=e*max(1,frozen_sigma/np.median(e))
        frozen_peaks,_=scan_short(t[keep],frozen_residual[keep],frozen_e[keep],n=3)
        # Full waveform: a sinusoid versus the best periodic constant box.
        box=BoxLeastSquares(t,f,e)
        box_result=box.power([period],[.01,.02,.035,.05,min(.15,period*.45),period*.45],oversample=10,method='slow')
        boxmodel=box.model(t,period,float(box_result.duration[0]),float(box_result.transit_time[0]))
        box_chi2=float(np.sum(((f-boxmodel)/e)**2))
        null_chi2=float(np.sum(((f-np.average(f,weights=1/e**2))/e)**2))
        sector=dict(sector=row['sector'],n_points=len(t),t_min=float(t.min()),t_max=float(t.max()),
            ls_best_period_days=float(1/sf[sj]),ls_peak_power=float(sp[sj]),
            fixed_period_sinusoid_explained_fraction=float(1-one_chi2/null_chi2),
            sinusoid_peak_to_peak_ppm=float(np.ptp(one_model)*1e6),
            harmonics_explained_fraction=float(1-chi2/null_chi2),
            best_box_explained_fraction=float(1-box_chi2/null_chi2),best_box_duration_days=float(box_result.duration[0]),
            coefficients=coeff.tolist(),residual_robust_scatter_ppm=float(sigma*1e6),
            frozen_waveform_residual_scatter_ppm=float(frozen_sigma*1e6),
            residual_peaks=peaks,frozen_model_residual_peaks=frozen_peaks,
            frozen_model_rms_ppm=float(np.std(f-frozen)*1e6))
        summary['sectors'].append(sector)
        residual_chunks.append((t[keep],residual[keep],e_search[keep]))
        ax[0].plot(1/sf,sp,lw=.5);ax[0].axvline(period,color='r',lw=.7)
        ax[0].set(xlabel='Period (d)',ylabel='Lomb–Scargle power',title=f'TIC{tic} sector{row["sector"]}',xlim=(.1,3))
        phi,y,se=folded(t,f,period,origin)
        ax[1].errorbar(phi,(y-1)*1e6,yerr=se*1e6,fmt='.',ms=3)
        tx=origin+phi*period
        ax[1].plot(phi,(design(tx,frequency,h,origin)@coeff-1)*1e6,color='r',label=f'{h} Fourier harmonics')
        ax[1].plot(phi,(box.model(tx,period,float(box_result.duration[0]),float(box_result.transit_time[0]))-1)*1e6,color='.4',label='best box')
        ax[1].set(xlabel=f'Phase at train P={period:.7f}d',ylabel='Flux (ppm)');ax[1].legend(fontsize=7)
        best=peaks[0];ph,y,se=folded(t[keep],residual[keep],best['period_days'],best['epoch_btjd'])
        ax[2].errorbar(ph,(y-1)*1e6,yerr=se*1e6,fmt='.',ms=3)
        ax[2].set(xlabel=f'Residual phase P={best["period_days"]:.6f}d',ylabel='Residual (ppm)',
                  title=f'Exploratory max short BLS S/N={best["diagonal_snr"]:.1f}')
        print(tic,row['sector'],'trainP',period,'sinfrac',sector['fixed_period_sinusoid_explained_fraction'],
              'boxfrac',sector['best_box_explained_fraction'],'resBLS',[(round(x['period_days'],5),round(x['diagonal_snr'],1)) for x in peaks],flush=True)
    # Fixed train ephemerides, evaluated without choosing a test-sector phase.
    train_peaks=[]
    for row,chunk,sector in zip(records,residual_chunks,summary['sectors']):
        if row['sector'] in [95,96]:train_peaks+=sector['residual_peaks'][:3]
    summary['train_peak_test_checks']=[]
    for peak in train_peaks:
        checks=[]
        for row,chunk in zip(records,residual_chunks):
            checks.append(dict(sector=row['sector'],measurement=fixed_box(*chunk,peak)))
        summary['train_peak_test_checks'].append(dict(peak=peak,sectors=checks))
    OUT.mkdir(parents=True,exist_ok=True)
    label=f'tic{tic}_{column}_w{width:g}'
    fig.savefig(OUT/f'{label}.png',dpi=135);plt.close(fig)
    (OUT/f'{label}.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    inter=CACHE/'intermediate'/f'{label}_variability_residuals.npz'
    np.savez_compressed(inter,**{f's{r["sector"]}_{v}':c[j] for r,c in zip(records,residual_chunks) for j,v in enumerate(['t','f','e'])})
    return summary


def probe_commensurability(records, period=.689773):
    """Falsify a historical high-ranked 2/3-waveform-period residual alias.

    Training sector96 and 2d detrending choose phase/duration. All alternative
    preprocessing and sectors evaluate that same frozen ephemeris.
    """
    from validate import local_depth, noise_beta
    tic=150359500
    variants=[('PDCSAP_FLUX',2),('PDCSAP_FLUX',5),('SAP_FLUX',5)]
    summary=dict(tic=tic,period_days=period,
                 motivation='Historical 2d-detrending residual peak near 2/3 of the dominant 1.033d waveform; selected before positive flare clipping correction.',
                 train_sector=96,tests=[],created_utc=datetime.now(timezone.utc).isoformat())
    ephemeris=None
    fig,axes=plt.subplots(3,4,figsize=(16,10),layout='constrained')
    for vi,(column,width) in enumerate(variants):
        source=json.loads((OUT/f'tic{tic}_{column}_w{width}.json').read_text())
        freq,h,origin=(source[k] for k in ['frequency_per_day','selected_harmonics','origin_btjd'])
        chunks=[]
        for row in records:
            t,f,e=prepare(row,width,column)
            _,model,_,_=fit(t,f,e,freq,h,origin)
            residual=f-model+1
            sigma=1.4826*np.median(np.abs(residual-np.median(residual)))
            e*=max(1,sigma/np.median(e));keep=residual<1+7*sigma
            chunks.append((t[keep],residual[keep],e[keep]))
        if ephemeris is None:
            k=next(i for i,r in enumerate(records) if r['sector']==96)
            t,f,e=chunks[k]
            z=BoxLeastSquares(t,f,e).power([period],DURATIONS,method='slow',objective='snr',oversample=20)
            ephemeris=dict(period_days=period,epoch_btjd=float(z.transit_time[0]),duration_days=float(z.duration[0]))
            summary['frozen_ephemeris']=ephemeris
        for row,chunk,ax in zip(records,chunks,axes[vi]):
            t,f,e=chunk
            direct=fixed_box(t,f,e,ephemeris)
            local=local_depth(t,f,e,period,ephemeris['epoch_btjd'],ephemeris['duration_days'])
            beta=noise_beta(t,f,e,60)
            summary['tests'].append(dict(column=column,width_days=width,sector=row['sector'],
                global_box=direct,local_linear_baseline=local,red_noise=beta))
            ph,y,err=folded(t,f,period,ephemeris['epoch_btjd'],100)
            ax.errorbar(ph,(y-1)*1e6,yerr=err*1e6,fmt='.',ms=3)
            ax.set(xlabel='Frozen phase P=.689773d',ylabel='Residual (ppm)',
                   title=f'{column} w{width} s{row["sector"]}: boxS/N={direct["diagonal_snr"]:.1f}')
    (OUT/'tic150359500_commensurability_probe.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    fig.savefig(OUT/'tic150359500_commensurability_probe.png',dpi=135);plt.close(fig)
    print('COMMENSURABILITY',[(x['column'],x['width_days'],x['sector'],round(x['global_box']['diagonal_snr'],1),round(x['local_linear_baseline'].get('event_scatter_snr',0),1)) for x in summary['tests']],flush=True)

    # A small, fixed-phase nonblind preservation control; not a completeness
    # curve and not a calibrated injection-recovery blind-search experiment.
    controls=[]
    for width in [2,5]:
        source=json.loads((OUT/f'tic{tic}_PDCSAP_FLUX_w{width}.json').read_text())
        freq,h,origin=(source[k] for k in ['frequency_per_day','selected_harmonics','origin_btjd'])
        for row in records:
            depths=[]
            for injected in [0,.001]:
                spec=(period,ephemeris['epoch_btjd'],ephemeris['duration_days'],injected) if injected else None
                t,f,e=prepare(row,width,'PDCSAP_FLUX',spec)
                _,model,_,_=fit(t,f,e,freq,h,origin)
                residual=f-model+1
                sigma=1.4826*np.median(np.abs(residual-np.median(residual)))
                keep=residual<1+7*sigma
                e*=max(1,sigma/np.median(e))
                depths.append(fixed_box(t[keep],residual[keep],e[keep],ephemeris)['depth_ppm'])
            controls.append(dict(width_days=width,sector=row['sector'],injected_depth_ppm=1000,
                                 original_depth_ppm=depths[0],injected_measurement_ppm=depths[1],
                                 recovered_difference_ppm=depths[1]-depths[0]))
    summary['nonblind_preservation_controls']=controls
    (OUT/'tic150359500_commensurability_probe.json').write_text(json.dumps(summary,indent=2,allow_nan=False))


def summary_tables():
    sectors,checks=[],[]
    for path in sorted(OUT.glob('tic*_FLUX_w*.json')):
        source=json.loads(path.read_text())
        for s in source['sectors']:
            sectors.append(dict(tic=source['tic'],column=source['column'],width_days=source['width_days'],
                sector=s['sector'],train_period_days=source['period_days'],harmonics=source['selected_harmonics'],
                n_points=s['n_points'],sinusoid_ptp_ppm=s['sinusoid_peak_to_peak_ppm'],
                sinusoid_explained_fraction=s['fixed_period_sinusoid_explained_fraction'],
                best_box_explained_fraction=s['best_box_explained_fraction'],
                residual_scatter_ppm=s['residual_robust_scatter_ppm'],
                strongest_short_fast_selected_period=s['residual_peaks'][0]['period_days'],
                strongest_short_exact_snr=s['residual_peaks'][0]['diagonal_snr']))
        for item in source['train_peak_test_checks']:
            peak=item['peak']; measurements={f"s{s['sector']}_snr":s['measurement']['diagonal_snr'] for s in item['sectors']}
            checks.append(dict(tic=source['tic'],column=source['column'],width_days=source['width_days'],
                period_days=peak['period_days'],epoch_btjd=peak['epoch_btjd'],duration_days=peak['duration_days'],
                **measurements,status=('follow-up required' if min(measurements.get('s97_snr',0),measurements.get('s98_snr',0))>=5
                    else 'no >=5 diagonal S/N recurrence in both held-out sectors; threshold not calibrated')))
    for name,rows in [('waveform_summary.csv',sectors),('frozen_ephemeris_checks.csv',checks)]:
        if rows:
            with (OUT/name).open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--tic',type=int,nargs='+',default=[348839788,150359500]);p.add_argument('--column',default='PDCSAP_FLUX');p.add_argument('--width',type=float,default=5);p.add_argument('--probe-only',action='store_true')
    a=p.parse_args();rows=json.loads((ROOT/'data'/'observations.json').read_text())
    OUT.mkdir(parents=True,exist_ok=True)
    if a.probe_only:
        selected=sorted([r for r in rows if r['tic']==150359500 and r['sector'] in [95,96,97,98]],key=lambda r:r['sector'])
        probe_commensurability(selected);summary_tables();return
    for tic in a.tic:
        selected=sorted([r for r in rows if r['tic']==tic and r['sector'] in [95,96,97,98]],key=lambda r:r['sector'])
        if len(selected)!=4:raise ValueError(f'Need sectors95-98 for TIC{tic}')
        run(tic,selected,a.column,a.width)
    summary_tables()


if __name__=='__main__':main()
