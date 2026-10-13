"""Vetting of TIC 150299576's sector-13 two-event ~11.27d alignment.

Write only candidates/tic150299576. Use original FITS and existing screen
results. Regional BLS and bootstrap maxima are screening diagnostics, not a
full-search false-alarm probability or a validated transit model.
"""
from pathlib import Path
import argparse
import gzip
import json
import math
import os
from datetime import datetime, timezone
os.environ.setdefault('MPLCONFIGDIR','/tmp/tic150299576-matplotlib')
import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from search import ROOT, read_curve, flatten, phase
from vet_known_host_leads import safe, fixed_depth, predictions, segments

TIC=150299576
OUT=ROOT/'candidates'/f'tic{TIC}'
PERIODS=np.linspace(10.5,12,1201)
DURATIONS=np.array([.065,.09,.13,.18,.23])


def local_raw_box(t,f,e,center,duration,degree=1,halfwindow=1.0):
    selected=np.abs(t-center)<halfwindow
    tx=t[selected]-center;fx=f[selected];ex=e[selected]
    inside=np.abs(tx)<duration/2
    if inside.sum()<5 or (~inside).sum()<30:
        return {'status':'insufficient data','n_in':int(inside.sum()),'n_out':int((~inside).sum())}
    normalization=np.median(fx)
    fx,ex=fx/normalization,ex/normalization
    design=np.column_stack([tx**i for i in range(degree+1)]+[-inside.astype(float)])
    w=1/ex**2
    normal=design.T@(w[:,None]*design)
    covariance=np.linalg.inv(normal)
    coefficient=covariance@(design.T@(w*fx))
    residual=fx-design@coefficient
    inflation=max(1,float(np.sum(w*residual**2)/(len(tx)-design.shape[1])))
    covariance*=inflation
    depth=float(coefficient[-1]);error=float(np.sqrt(covariance[-1,-1]))
    return {'depth':depth,'depth_err_diagonal_residual_scaled':error,'snr_diagonal':depth/error,
            'degree':degree,'halfwindow_days':halfwindow,'n_in':int(inside.sum()),'n_out':int((~inside).sum()),
            'residual_rms':float(np.std(residual)),'chi2_error_inflation':inflation,
            'limitation':'Polynomial baseline and independent residuals; transit shape is a fixed box, not a physical planet fit.'}


def rednoise_beta(t,f,period,t0,duration):
    # Nonoverlapping sample-count blocks only within contiguous usable intervals.
    oot=np.abs(phase(t,period,t0))>duration
    sigma=float(np.std(f[oot],ddof=1))
    step=float(np.median(np.diff(t)))
    n=max(2,round(duration/step))
    means=[]
    for lo,hi in segments(t):
        idx=np.flatnonzero((t>=lo)&(t<=hi)&oot)
        # Split again where masking has removed a candidate event.
        if not len(idx):continue
        boundaries=np.r_[0,np.flatnonzero(np.diff(t[idx])>step*1.8)+1,len(idx)]
        for first,last in zip(boundaries[:-1],boundaries[1:]):
            use=f[idx[first:last]]
            for offset in range(0,len(use)-n+1,n):means.append(float(np.mean(use[offset:offset+n])))
    if len(means)<3:return {'beta':None,'n_blocks':len(means)}
    expected=sigma/np.sqrt(n)*np.sqrt(len(means)/(len(means)-1))
    beta=float(np.std(means,ddof=1)/expected) if expected>0 else None
    return {'beta':max(1,beta) if beta is not None else None,'raw_beta':beta,'n_blocks':len(means),
            'block_points':n,'duration_days':duration,
            'limitation':'Binning diagnostic for correlated variability, not a likelihood model.'}


def region(t,f,e):
    result=BoxLeastSquares(t,f,e).power(PERIODS,DURATIONS,oversample=10,objective='snr')
    j=int(np.argmax(result.power))
    stats=BoxLeastSquares(t,f,e).compute_stats(float(result.period[j]),float(result.duration[j]),float(result.transit_time[j]))
    return {'period':float(result.period[j]),'t0_btjd':float(result.transit_time[j]),'duration_days':float(result.duration[j]),
            'depth':float(result.depth[j]),'depth_err_diagonal':float(result.depth_err[j]),'snr_diagonal':float(result.depth_snr[j]),
            'stats':safe(dict(stats))},result


def resample_blocks(t,f,rng,block_days):
    result=np.empty_like(f)
    for lo,hi in segments(t):
        idx=np.flatnonzero((t>=lo)&(t<=hi))
        step=np.median(np.diff(t[idx])) if len(idx)>1 else 10/1440
        n=max(2,round(block_days/step));values=f[idx]
        pieces=[];length=0
        while length<len(values):
            start=int(rng.integers(0,max(1,len(values)-n+1)))
            piece=values[start:start+n];pieces.append(piece);length+=len(piece)
        result[idx]=np.concatenate(pieces)[:len(values)]
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--null-draws',type=int,default=64)
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    recs=[rec for rec in json.loads((ROOT/'data'/'observations.json').read_text()) if rec['tic']==TIC]
    screen=json.loads((ROOT/'results'/'scans'/f'tic{TIC}_s0013_w1_PDCSAP_FLUX_residual.json').read_text())
    lead=screen['peaks'][0];p,t0,d=lead['period'],lead['t0_btjd'],lead['duration_days']
    curves={};raw={};results=[]
    for rec in recs:
        sector=rec['sector'];bounds=None
        with fits.open(rec['local_path']) as hdul:
            metadata={key:hdul[1].header.get(key) for key in ['TIMESYS','BJDREFI','BJDREFF','CROWDSAP','FLFRCSAP']}
            metadata.update({key:hdul[0].header.get(key) for key in ['TICID','SECTOR','TEFF','RADIUS','TESSMAG','PROCVER']})
            assert metadata['TIMESYS']=='TDB' and metadata['BJDREFI']==2457000 and metadata['BJDREFF']==0
        row={'sector':sector,'data_uri':rec['data_uri'],'raw_sha256':rec['sha256'],'metadata':metadata,'reductions':[]}
        for column in ['PDCSAP_FLUX','SAP_FLUX']:
            t,f,e,h=read_curve(rec['local_path'],column);raw[(sector,column)]=(t,f,e);bounds=segments(t)
            for width in [1.,2.]:
                tb,fb,eb=flatten(t,f,e,width);curves[(sector,column,width)]=(tb,fb,eb)
                events=predictions(tb,fb,eb,p,d,t0,np.ones(len(tb),bool),bounds,[])
                for event in events:
                    event['raw_baseline_models']=[local_raw_box(t,f,e,event['center_btjd'],d,degree,halfwindow)
                        for degree in [1,2] for halfwindow in [.65,1.25]]
                row['reductions'].append({'column':column,'width_days':width,
                    'fixed_ephemeris':fixed_depth(tb,fb,eb,p,d,t0),'events':events,'correlated_noise':rednoise_beta(tb,fb,p,t0,d)})
                path=OUT/f'sector{sector:02d}_w{width:g}_{column}.csv.gz'
                with gzip.open(path,'wt') as stream:
                    np.savetxt(stream,np.column_stack([tb,fb,eb]),delimiter=',',header='btjd,normalized_flux,error_diagonal',comments='',fmt='%.11g')
        results.append(row)
        print('Prepared sector',sector,flush=True)
    joint=[];periodograms={}
    for column in ['PDCSAP_FLUX','SAP_FLUX']:
        for width in [1.,2.]:
            merged=[np.concatenate([curves[(rec['sector'],column,width)][i] for rec in recs]) for i in range(3)]
            best,pg=region(*merged);periodograms[(column,width)]=pg
            best.update({'column':column,'width_days':width,'period_grid_min':10.5,'period_grid_max':12.,'period_grid_size':len(PERIODS)})
            best['original_ephemeris_fixed']=fixed_depth(*merged,p,d,t0)
            joint.append(best)
            np.savez_compressed(OUT/f'joint_w{width:g}_{column}.npz',t=merged[0],f=merged[1],e=merged[2],period=pg.period,snr=pg.power)
            print('Joint',column,width,best['period'],best['snr_diagonal'],flush=True)
    # Empirical check for the selected star and regional search. It does not
    # include the earlier many-star/full-period selection and cannot validate.
    rng=np.random.default_rng(20261013);null=[]
    observed=max(entry['snr_diagonal'] for entry in joint if entry['column']=='PDCSAP_FLUX')
    for draw in range(args.null_draws):
        block_days=.3 if draw<args.null_draws//2 else .75
        maxima=[]
        for width in [1.,2.]:
            resampled=[]
            for rec in recs:
                tx,fx,ex=curves[(rec['sector'],'PDCSAP_FLUX',width)]
                resampled.append((tx,resample_blocks(tx,fx,rng,block_days),ex))
            merged=[np.concatenate([arr[i] for arr in resampled]) for i in range(3)]
            best,_=region(*merged);maxima.append(best['snr_diagonal'])
        null.append({'draw':draw,'block_days':block_days,'max_snr':max(maxima)})
        if (draw+1)%16==0:print('Null draws',draw+1,'/',args.null_draws,flush=True)
    exceed=sum(row['max_snr']>=observed for row in null)
    null_summary={'seed':20261013,'draws':args.null_draws,'observed_max_joint_pdc_snr':observed,
        'exceedances':exceed,'empirical_tail_plus_one':(exceed+1)/(len(null)+1),
        'calibration_scope':'One already selected star, 10.5-12d, five durations, maxima across two PDC detrenders; original full-period/many-star trials are absent.',
        'limitations':['Moving blocks preserve short-range variability only, introduce boundaries, and do not model all instrumental covariance.',
                       'Finite ensemble; tail estimate is diagnostic, not a global discovery false-alarm probability.'],'samples':null}
    summary={'tic':TIC,'created_utc':datetime.now(timezone.utc).isoformat(),'original_sector13_peak':lead,
             'sectors':results,'joint_regional_scans':joint,'conditional_null':null_summary,
             'assessment':'Not a credible recurring planetary candidate at the sector13 ephemeris; reduction dependence, correlated variability and nonrepetition require further data before promotion.',
             'time_convention':'BTJD=BJD_TDB-2457000, verified in each FITS',
             'novelty_status':'Catalog/literature assessment is maintained separately; flux tests do not establish novelty.'}
    (OUT/'vetting.json').write_text(json.dumps(safe(summary),indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(3,2,figsize=(12,11),layout='constrained')
    sourcecenters=[t0,t0+p]
    for index,center in enumerate(sourcecenters):
        for col,color in [('PDCSAP_FLUX','tab:blue'),('SAP_FLUX','tab:orange')]:
            tx,fx,ex=curves[(13,col,1.)];sel=np.abs(tx-center)<.8
            axes[0,index].plot(tx[sel]-center,fx[sel],'.',ms=2,color=color,label=col)
            tx,fx,ex=raw[(13,col)];sel=np.abs(tx-center)<.8
            axes[1,index].plot(tx[sel]-center,fx[sel]/np.median(fx[sel]),'.',ms=1,color=color,label=col)
        for ax in [axes[0,index],axes[1,index]]:
            ax.axvspan(-d/2,d/2,color='.8',alpha=.4);ax.legend(markerscale=4);ax.set(xlabel='Days from predicted event center',ylabel='Normalized flux')
        axes[0,index].set_title(f'Sector13 event {center:.5f}: independent1d median detrending')
        axes[1,index].set_title('Raw locally normalized photometry')
    for (col,width),pg in periodograms.items():
        axes[2,0].plot(pg.period,pg.power,label=f'{col} {width:g}d',lw=.7)
    axes[2,0].axvline(p,color='black',ls='--',lw=1);axes[2,0].set(xlabel='Joint period (d)',ylabel='BLS diagonal S/N',title='Four-sector regional search');axes[2,0].legend(fontsize=8)
    for row in results:
        red=next(red for red in row['reductions'] if red['column']=='PDCSAP_FLUX' and red['width_days']==1)
        for event in red['events']:
            if event['depth_local'] is not None:
                axes[2,1].errorbar(event['center_btjd'],event['depth_local']*1e6,yerr=event['depth_err_diagonal']*1e6,
                                   fmt='o',color='tab:red' if row['sector']==13 else 'black')
    axes[2,1].axhline(0,color='.6');axes[2,1].set(xlabel='Predicted event BTJD',ylabel='Local PDC depth(ppm)',title='Original sector13 ephemeris; red=source sector')
    fig.savefig(OUT/'diagnostics.png',dpi=140);plt.close(fig)
    lines=['# TIC 150299576: sector 13 weak long-period alignment','',
        'The two central sector 13 events are investigated as a weak lead. They do not currently support a credible recurring planet. Original and independently normalized scientific curves, exact FITS hashes, sector tests and conditional null samples are retained here.', '',
        f'Original screen: P = {p:.8f} d, epoch BTJD = {t0:.6f}, box duration = {d*24:.2f} h. The maximum duration grid point was selected, making baseline choice especially relevant. Both SAP and PDC were analyzed, with 1 d and 2 d rolling-median detrending and independent local linear/quadratic baseline fits to raw 2 min data.', '',
        '| Sector | PDC fixed S/N, 1 d | PDC fixed S/N, 2 d | SAP fixed S/N, 1 d | SAP fixed S/N, 2 d |', '| --- | --- | --- | --- | --- |']
    for row in results:
        vals=[next(red['fixed_ephemeris']['snr_diagonal'] for red in row['reductions'] if red['column']==col and red['width_days']==width) for col in ['PDCSAP_FLUX','SAP_FLUX'] for width in [1.,2.]]
        lines.append(f"| {row['sector']} | "+' | '.join(f'{value:.2f}' for value in vals)+' |')
    lines += ['', 'The joint 10.5–12 d search uses 1201 periods and five durations. Its best PDC scores are '+', '.join(f"{row['snr_diagonal']:.2f} at {row['period']:.5f} d (width {row['width_days']:g} d)" for row in joint if row['column']=='PDCSAP_FLUX')+'. They shift away from the sector 13 winning period. This does not establish a coherent extension of the original two events.', '',
        f'An empirical moving-block screen used {args.null_draws} realizations (0.3 d and 0.75 d block sizes) and maxima over the regional grid and both PDC detrenders. {exceed} null maxima equaled/exceeded the observed joint maximum; the plus-one empirical tail is {null_summary["empirical_tail_plus_one"]:.3f}. It omits the earlier many-star/full-period selection and is not a global false-alarm probability. The small ensemble does not support precise tail claims.', '',
        'Independent local raw SAP baselines give only 52–332 ppm for the first source event (0.2–1.9 diagonal sigma), compared with 568–1366 ppm in PDC. The second event has 544–655 ppm in SAP versus 1012–1106 ppm in PDC. These values are conditional on the fixed box and baseline window; their discrepancy argues against promoting the PDC interpretation. Sector 13 PDC red-noise binning factors at the proposed duration are 1.31 and 1.57 for the 1 d and 2 d detrenders.', '',
        'Each original ephemeris prediction, local raw-baseline depth/error, and red-noise binning factor is in `vetting.json`. The observations in sectors 10–12 were already available for analysis, so this is retrospective robustness. Two-event period uncertainty, period aliases and transit timing variations have not been eliminated by a formal posterior; fixed-time null results alone cannot exclude every nearby ephemeris.', '',
        'The inferred circular central-transit duration for a roughly 0.62 solar-radius/0.60 solar-mass TIC dwarf at 11.27 d is about 3 h. The selected 5.5 h box is long compared with that scale, but eccentricity and stellar-parameter uncertainty can lengthen it. No radius or physical planet parameters are claimed for a reduction-dependent, nonrepeating lead. No confirmation or formal validation is inferred.', '',
        '[Diagnostic figure](diagnostics.png). All normalized 10 min scientific curves are retained as CSV.GZ; raw FITS retrieval products and SHA256 values remain in `vetting.json`. Reproduce with `python src/vet_tic150299576.py --null-draws 64`. Investigate pixels or additional sectors only if a transit-preserving independent reduction produces consistent events at a jointly constrained ephemeris.']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('Saved',OUT,'conditional tail',null_summary['empirical_tail_plus_one'],flush=True)

if __name__=='__main__':main()
