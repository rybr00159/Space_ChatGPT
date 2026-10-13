"""Falsification of a broad 4.59d TOI-1803 residual, without discovery claims."""
import csv
from datetime import datetime, timezone
import json
import numpy as np
from astropy.timeseries import BoxLeastSquares, LombScargle
from search import ROOT, CACHE, read_curve, flatten, phase, known_for, mask_known
from validate import noise_beta
import matplotlib.pyplot as plt

OUT=ROOT/'candidates'/'toi1803'
P=4.590202809336783
T0=1903.0191751754433
D=.2295
DURATIONS=np.array([.045,.065,.09,.13,.18,.23])


def local_poly(t,f,e,p,t0,d,degree=1,balanced=False):
    events=[]
    cycles=np.round((t-t0)/p).astype(int)
    for cycle in np.unique(cycles):
        center=t0+cycle*p;x=(t-center)/d
        inside=np.abs(x)<.5
        outside=(np.abs(x)>1)&(np.abs(x)<3)
        n_before=int((outside&(x<0)).sum());n_after=int((outside&(x>0)).sum())
        if balanced and min(n_before,n_after)<3:continue
        use=inside|outside
        if inside.sum()<2 or outside.sum()<6:continue
        design=np.column_stack([x[use]**k for k in range(degree+1)]+[-inside[use].astype(float)])
        w=1/e[use]**2;normal=design.T@(w[:,None]*design)
        if np.linalg.cond(normal)>1e10:continue
        covariance=np.linalg.inv(normal);coeff=covariance@(design.T@(w*f[use]))
        model=design@coeff
        chi2=np.sum(w*(f[use]-model)**2)/max(1,use.sum()-degree-2)
        events.append(dict(cycle=int(cycle),center_btjd=float(center),depth_ppm=float(coeff[-1]*1e6),
            formal_error_ppm=float(np.sqrt(covariance[-1,-1])*1e6),
            residual_inflated_error_ppm=float(np.sqrt(covariance[-1,-1]*max(1,chi2))*1e6),
            residual_chi2=float(chi2),n_in=int(inside.sum()),n_out=int(outside.sum())))
        events[-1].update(n_baseline_before=n_before,n_baseline_after=n_after,
                         bilateral_baseline_available=bool(min(n_before,n_after)>=3))
    if not events:return dict(n_events=0,events=[])
    dep=np.array([x['depth_ppm'] for x in events]);err=np.array([x['residual_inflated_error_ppm'] for x in events]);w=1/err**2
    mean=np.average(dep,weights=w);formal=1/np.sqrt(w.sum())
    scatter_chi2=np.sum(w*(dep-mean)**2)/max(1,len(dep)-1)
    inflated=formal*np.sqrt(max(1,scatter_chi2))
    return dict(depth_ppm=float(mean),residual_inflated_error_ppm=float(formal),
        event_scatter_error_ppm=float(inflated),event_scatter_snr=float(mean/inflated),
        event_scatter_chi2=float(scatter_chi2),n_events=len(events),events=events)


def known_mask(t,extra_half_width=0):
    keep=mask_known(t,known_for(144401492))
    if extra_half_width:
        for k in known_for(144401492):
            keep &= np.abs(phase(t,k['period'],k['epoch_bjd']-2457000))>extra_half_width
    return keep


def prepare(row,column,width,extra_half_width=0,injection=None):
    t,f,e,_=read_curve(row['local_path'],column)
    if injection is not None:
        ip,it0,iduration,depth=injection
        f=f*(1-depth*(np.abs(phase(t,ip,it0))<iduration/2))
    edges=np.r_[0,np.flatnonzero(np.diff(t)>.2)+1,len(t)]
    boundaries=np.array([v for a,b in zip(edges[:-1],edges[1:]) for v in (t[a],t[b-1])])
    if width==0:
        med=np.median(f);tb,fb,eb=t,f/med,e/med
    else:tb,fb,eb=flatten(t,f,e,width)
    keep=np.zeros(len(tb),bool)
    for a,b in zip(edges[:-1],edges[1:]):
        keep|=(tb>=t[a]+.25)&(tb<=t[b-1]-.25)
    keep &= known_mask(tb,extra_half_width)
    return tb[keep],fb[keep],eb[keep],boundaries


def fixed_box(t,f,e,p=P,t0=T0,d=D):
    stats=BoxLeastSquares(t,f,e).compute_stats(p,d,t0)
    dep,err=stats['depth']
    return dict(depth_ppm=float(dep*1e6),formal_error_ppm=float(err*1e6),diagonal_snr=float(dep/err))


def family_scan(t,f,e):
    periods=np.linspace(4.3,4.9,1201)
    bls=BoxLeastSquares(t,f,e);z=bls.power(periods,DURATIONS,objective='snr',oversample=10)
    j=np.argmax(z.power)
    exact=bls.power([float(z.period[j])],DURATIONS,objective='snr',oversample=10,method='slow')
    top=dict(period_days=float(exact.period[0]),epoch_btjd=float(exact.transit_time[0]),
             duration_days=float(exact.duration[0]),depth_ppm=float(exact.depth[0]*1e6),
             diagonal_snr=float(exact.depth_snr[0]))
    plateau=periods[z.power>=.95*z.power.max()]
    top['descriptive_95percent_peak_period_range_days']=[float(plateau.min()),float(plateau.max())]
    return top,z


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records=sorted([r for r in json.loads((ROOT/'data'/'observations.json').read_text()) if r['tic']==144401492],key=lambda r:r['sector'])
    summary=dict(tic=144401492,system='TOI-1803',created_utc=datetime.now(timezone.utc).isoformat(),
        frozen_ephemeris=dict(period_days=P,epoch_btjd=T0,duration_days=D),known_planets=known_for(144401492),
        observations=records,measurements=[],warning='Diagonal statistics and period-profile plateau are descriptive; not calibrated FAP or period credible intervals. Current known-planet masks do not include measured transit timing variations.')
    prepared={};csvrows=[]
    for row in records:
        for col in ['PDCSAP_FLUX','SAP_FLUX']:
            for width in [0,1,2,5]:
                t,f,e,b=prepare(row,col,width)
                box=fixed_box(t,f,e)
                linear=local_poly(t,f,e,P,T0,D,1);quad=local_poly(t,f,e,P,T0,D,2)
                for event in linear['events']:event['nearest_good_segment_boundary_days']=float(np.min(np.abs(b-event['center_btjd'])))
                measurement=dict(sector=row['sector'],column=col,width_days=width,n_searched=len(t),
                    fixed_box=box,linear=linear,quadratic=quad,
                    bilateral_quadratic=local_poly(t,f,e,P,T0,D,2,balanced=True),
                    narrow_quadratic=[dict(duration_days=d,measurement=local_poly(t,f,e,P,T0,d,2)) for d in [.08,.11,.16]],
                    red_noise=[noise_beta(t,f,e,minutes) for minutes in [60,120]],good_segment_boundaries_btjd=b.tolist())
                if width in [1,2]:
                    fam,z=family_scan(t,f,e);measurement['family_scan']=fam
                    prepared[(row['sector'],col,width)]=(t,f,e)
                summary['measurements'].append(measurement)
                csvrows.append(dict(sector=row['sector'],column=col,width_days=width,
                    fixed_depth_ppm=box['depth_ppm'],fixed_diagonal_snr=box['diagonal_snr'],
                    local_linear_depth_ppm=linear.get('depth_ppm'),local_linear_eventscatter_snr=linear.get('event_scatter_snr'),
                    local_quadratic_depth_ppm=quad.get('depth_ppm'),local_quadratic_eventscatter_snr=quad.get('event_scatter_snr')))
                print(row['sector'],col,width,'box',round(box['diagonal_snr'],2),'linear',round(linear.get('event_scatter_snr',0),2),
                      'quadratic',round(quad.get('event_scatter_snr',0),2),'family',measurement.get('family_scan'),flush=True)
    # Broader known-planet masks test leakage from timing-varying known planets.
    summary['known_mask_sensitivity']=[]
    for width in [1,2]:
        for halfwidth in [.25,.5]:
            row=next(r for r in records if r['sector']==22)
            t,f,e,_=prepare(row,'PDCSAP_FLUX',width,halfwidth)
            summary['known_mask_sensitivity'].append(dict(width_days=width,known_half_mask_days=halfwidth,
                n_searched=len(t),fixed_box=fixed_box(t,f,e),quadratic=local_poly(t,f,e,P,T0,D,2)))
    # The long separation permits many cycle aliases; a narrow family search
    # is a diagnostic consistency test, not independent confirmation.
    summary['joint_family_scans']=[]
    for col in ['PDCSAP_FLUX','SAP_FLUX']:
        for width in [1,2]:
            t,f,e=(np.concatenate([prepared[(s,col,width)][k] for s in [22,49]]) for k in range(3))
            z=BoxLeastSquares(t,f,e).power(np.linspace(4.3,4.9,4801),DURATIONS,objective='snr',oversample=10)
            j=np.argmax(z.power);p,t0,d=(float(z[k][j]) for k in ['period','transit_time','duration'])
            exact=BoxLeastSquares(t,f,e).power([p],DURATIONS,objective='snr',oversample=10,method='slow')
            t0,d=float(exact.transit_time[0]),float(exact.duration[0])
            checks=[]
            for sector in [22,49]:
                ts,fs,es=prepared[(sector,col,width)]
                checks.append(dict(sector=sector,global_box=fixed_box(ts,fs,es,p,t0,d),
                    linear=local_poly(ts,fs,es,p,t0,d,1),quadratic=local_poly(ts,fs,es,p,t0,d,2)))
            summary['joint_family_scans'].append(dict(column=col,width_days=width,period_days=p,
                epoch_btjd=t0,duration_days=d,diagonal_snr=float(exact.depth_snr[0]),sector_checks=checks))
    # A broad BLS epoch can hide a narrower dip offset within its box. Search
    # modest training-only refinements with a quadratic local baseline rather
    # than interpreting the fixed narrow-box null as an absence statement.
    trial_results=[]
    t,f,e=prepared[(22,'PDCSAP_FLUX',2)]
    for p in np.linspace(4.53,4.6265,15):
        for t0 in T0+np.linspace(-.15,.15,21):
            for d in [.065,.09,.13]:
                result=local_poly(t,f,e,p,t0,d,2)
                trial_results.append(dict(period_days=float(p),epoch_btjd=float(t0),duration_days=d,
                                         local_event_scatter_snr=result.get('event_scatter_snr',0),measurement=result))
    trial_results.sort(key=lambda x:x['local_event_scatter_snr'],reverse=True)
    summary['narrow_quadratic_refinement']=dict(n_trials=len(trial_results),
        caveat='Training-only exploratory ranking; many correlated trials; no global FAP or independent validation.',
        top_trials=trial_results[:5],best_trial_checks=[])
    best=trial_results[0]
    for (sector,col,width),(t,f,e) in prepared.items():
        summary['narrow_quadratic_refinement']['best_trial_checks'].append(dict(sector=sector,column=col,width_days=width,
            measurement=local_poly(t,f,e,best['period_days'],best['epoch_btjd'],best['duration_days'],2)))
    print('NARROW REFINEMENT',best['period_days'],best['epoch_btjd'],best['duration_days'],best['local_event_scatter_snr'],flush=True)
    # Nonblind preservation controls determine whether local curvature removal
    # itself automatically erases a box at the investigated ephemeris.
    summary['nonblind_preservation_controls']=[]
    row=next(r for r in records if r['sector']==22)
    for col in ['PDCSAP_FLUX','SAP_FLUX']:
        for width in [0,1,2]:
            fits=[]
            for depth in [0,.001]:
                spec=(P,T0,D,depth) if depth else None
                t,f,e,_=prepare(row,col,width,injection=spec)
                fits.append(local_poly(t,f,e,P,T0,D,2))
            old={r['cycle']:r['depth_ppm'] for r in fits[0]['events']}
            changes=[dict(cycle=r['cycle'],recovered_increment_ppm=r['depth_ppm']-old[r['cycle']]) for r in fits[1]['events'] if r['cycle'] in old]
            summary['nonblind_preservation_controls'].append(dict(column=col,width_days=width,injected_depth_ppm=1000,
                median_recovered_increment_ppm=float(np.median([r['recovered_increment_ppm'] for r in changes])),events=changes))
    # Individual event plots retain surrounding waveform and known-mask points.
    fig,axes=plt.subplots(6,3,figsize=(13,17),layout='constrained')
    row=next(r for r in records if r['sector']==22)
    for c in range(6):
        center=T0+c*P
        for col in ['PDCSAP_FLUX','SAP_FLUX']:
            tr,fr,er,_=read_curve(row['local_path'],col);med=np.median(fr);fr/=med;er/=med
            use=np.abs(tr-center)<.85
            baseline=np.median(fr[use])
            axes[c,0].plot(tr[use]-center,(fr[use]-baseline)*1e6,'.',ms=1,label=col)
            for a,w in zip(axes[c,1:],[1,2]):
                t,f,e=prepared[(22,col,w)];near=np.abs(t-center)<.85
                a.plot(t[near]-center,(f[near]-1)*1e6,'.',ms=2,label=col)
        for a in axes[c]:
            a.axvspan(-D/2,D/2,color='r',alpha=.1);a.axvline(0,color='r',lw=.5)
            a.set(xlabel='Time from fixed event center (d)',ylabel='Flux offset (ppm)')
        axes[c,0].set_title(f'Native, center BTJD{center:.3f}');axes[c,0].legend(fontsize=7)
        axes[c,1].set_title('1d normalization; standard known mask');axes[c,2].set_title('2d normalization; standard known mask')
    fig.savefig(OUT/'sector22_individual_events.png',dpi=130);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for sector,axrow in zip([22,49],axes):
        for width,ax in zip([1,2],axrow):
            for col in ['PDCSAP_FLUX','SAP_FLUX']:
                t,f,e=prepared[(sector,col,width)];phi=phase(t,P,T0)
                idx=np.floor((phi+P/2)/.025).astype(int)
                x,y,err=[],[],[]
                for k in np.unique(idx):
                    use=idx==k;x.append(np.mean(phi[use]));y.append(np.mean(f[use]));err.append(np.std(f[use])/np.sqrt(use.sum()))
                ax.errorbar(x,(np.array(y)-1)*1e6,yerr=np.array(err)*1e6,fmt='.',ms=3,label=col)
            ax.axvspan(-D/2,D/2,color='r',alpha=.1);ax.legend();ax.set(xlabel='Frozen phase time (d)',ylabel='Flux (ppm)',
                title=f'Sector{sector}, width{width}d, fixedP4.5902d',xlim=(-1,1))
    fig.savefig(OUT/'frozen_phase_comparison.png',dpi=130);plt.close(fig)
    (OUT/'vetting.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    with (OUT/'measurements.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(csvrows[0]));writer.writeheader();writer.writerows(csvrows)
    np.savez_compressed(CACHE/'intermediate'/'toi1803_vetting.npz',**{f's{s}_{col}_w{w}_{name}':values[j] for (s,col,w),values in prepared.items() for j,name in enumerate(['t','f','e'])})


if __name__=='__main__':main()
