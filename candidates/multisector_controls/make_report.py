"""Regenerate exact-retention summaries and figures from retained control scans.

The experiment JSONL's ``queued`` field is the documented SNR+screening_pass
proxy. This audit applies the new search's actual additional positive-event
guard, without changing the original experiment records.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from multisector_screen import seed_reasons
from multisector_controls import read_observation, combine_fixed
from fast_screen import preprocess_curve, event_stats
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = ROOT / 'results/multisector_controls'
OUT = Path(__file__).resolve().parent


def phase_offset_control():
    """Fixed-period wrong-phase checks; no blind period or global null claim."""
    source=json.loads((DATA/'actual_controls_new.json').read_text())
    candidate=next(c for c in source['controls'] if c['tic']==4206066)
    observations=[r for r in source['observations'] if r['tic']==4206066]
    prepared={}
    for record in observations:
        arrays,_=read_observation(record)
        for width in [.75,1.5]:
            prepared[record['sector'],width]=preprocess_curve(*arrays[:3],quality=arrays[3],config=dict(width_days=width))
    rows=[]
    for width in [.75,1.5]:
        for k in range(17):
            epoch=candidate['epoch_btjd']+candidate['period']*k/17
            sectors=[]
            for record in observations:
                d=prepared[record['sector'],width]
                s=event_stats(d['t'],d['f'],d['e'],candidate['period'],epoch,candidate['duration_days'])
                sectors.append(dict(sector=record['sector'],**s))
            rows.append(dict(width_days=width,phase_offset_cycles=k/17,epoch_btjd=epoch,
                             sectors=sectors,combined=combine_fixed(sectors)))
    (DATA/'tic4206066_phase_offsets.json').write_text(json.dumps(dict(
        target=4206066,period_days=candidate['period'],rows=rows,
        scope='Published period held fixed, 16 displaced epochs plus published epoch. Illustrates phase specificity only; not a blind search or trial-corrected significance.',
        limitations=['Phase sidebands can overlap real transits. The candidate was chosen because its signal is known.',
                     'Widths are correlated transformations of the same flux; do not multiply their evidence.',
                     'Period and epoch uncertainties not propagated.']),indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(8,4))
    for width in [.75,1.5]:
        rr=[r for r in rows if r['width_days']==width]
        ax.plot([r['phase_offset_cycles'] for r in rr],[r['combined']['diagonal_snr'] for r in rr],'o-',ms=3,label=f'{width:g} d detrend')
    ax.set(xlabel='Epoch offset / published period',ylabel='Pooled diagonal local SNR',
           title='TIC 4206066: published-period phase-specificity control\nSelected candidate; 17 fixed epochs; no global false-alarm interpretation')
    ax.legend();fig.tight_layout();fig.savefig(OUT/'phase_offset_control.png',dpi=160);plt.close(fig)


def exact_retained(row, method):
    selected = []
    for m in row['matches']:
        p = next(p for p in row['peaks'] if p['variant']==m['variant'] and p['period']==m['period'])
        if method == 'new':
            reasons = seed_reasons(p)
        else:
            reasons = list(p['guard_failures'])
            if p['local_snr'] < 8:
                reasons.append('local_snr_below_8_control_proxy')
        selected.append(dict(variant=m['variant'], period=p['period'], reasons=reasons))
    return any(not v['reasons'] for v in selected), selected


def retained_at_snr8(row, method):
    for m in row['matches']:
        p=next(p for p in row['peaks'] if p['variant']==m['variant'] and p['period']==m['period'])
        reasons=seed_reasons(p) if method=='new' else p['guard_failures']
        if not reasons and p['local_snr']>=8:
            return True
    return False


def main():
    rows, audits = [], []
    for filename in sorted(DATA.glob('injections_*.jsonl')):
        method = filename.stem.removeprefix('injections_')
        cases = [json.loads(s) for s in filename.read_text().splitlines()]
        for r in cases:
            retained, matches = exact_retained(r, method)
            audits.append(dict(method=method,case_id=r['case_id'],retained=retained,
                               retained_at_snr8=retained_at_snr8(r,method),matches=matches))
        groups=sorted({(r['tic'],r['sector'],r['period_days'],r['depth_ppm']) for r in cases})
        for tic,sector,period,depth in groups:
            g=[r for r in cases if (r['tic'],r['sector'],r['period_days'],r['depth_ppm'])==(tic,sector,period,depth)]
            rows.append(dict(method=method,tic=tic,sector=sector,period_days=period,depth_ppm=depth,
                             n=len(g),found=sum(r['found_period_and_epoch'] for r in g),
                             retained=sum(exact_retained(r,method)[0] for r in g),
                             retained_at_snr8=sum(retained_at_snr8(r,method) for r in g),
                             fewer_than_two_sampled_events=sum(r['fixed_ephemeris_n_events']<2 for r in g),
                             median_fixed_ephemeris_diagonal_snr=float(np.median([r['fixed_ephemeris_snr'] for r in g]))))
    with (DATA/'exact_retention_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    actual=[]
    for filename in sorted(DATA.glob('actual_controls_*.json')):
        method=filename.stem.removeprefix('actual_controls_')
        for control in json.loads(filename.read_text())['controls']:
            for sector in control['sectors']:
                peaks=[dict(variant=i,**p) for i,s in enumerate(sector['scans']) for p in s['peaks']]
                retained,matches=exact_retained(dict(peaks=peaks,matches=sector['recovery']['matches']),method)
                actual.append(dict(method=method,tic=control['tic'],name=control['name'],sector=sector['sector'],
                                   found=sector['recovery']['found_period_and_epoch'],retained=retained,matches=matches))
    source=['src/multisector_screen.py','src/multisector_controls.py','candidates/multisector_controls/make_report.py']
    (DATA/'exact_retention_audit.json').write_text(json.dumps(dict(
        methods_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source},
        injections=audits,actual_controls=actual,
        old_retention_scope='Period+epoch, screening_pass and localSNR>=8; old catalog/radius/literature classification is not reproduced.',
        new_retention_scope='Exact src.multisector_screen.seed_reasons; no literature/classification validation.',
        limitation='Retention statistics are candidate-generation tests, not detection false-alarm probabilities.'),indent=2)+'\n')
    backgrounds=sorted({(r['tic'],r['sector']) for r in rows})
    groups=sorted({(r['period_days'],r['depth_ppm']) for r in rows})
    fig,axes=plt.subplots(len(backgrounds),2,figsize=(14,6),squeeze=False)
    for i,(tic,sector) in enumerate(backgrounds):
        for j,(metric,label) in enumerate([('found','Period and epoch recovered'),('retained','Retained by SNR and event guards')]):
            matrix=np.full((2,len(groups)),np.nan)
            labels={}
            for m,method in enumerate(['old','new']):
                for k,(period,depth) in enumerate(groups):
                    matches=[r for r in rows if (r['tic'],r['sector'],r['method'],r['period_days'],r['depth_ppm'])==(tic,sector,method,period,depth)]
                    if matches:
                        r=matches[0];matrix[m,k]=r[metric]/r['n'];labels[m,k]=f'{r[metric]}/{r["n"]}'
            ax=axes[i,j]
            ax.imshow(matrix,vmin=0,vmax=1,cmap='Blues',aspect='auto')
            for (m,k),labeltext in labels.items():
                ax.text(k,m,labeltext,ha='center',va='center',color='white' if matrix[m,k]>.5 else 'black',fontsize=8)
            ax.set(yticks=[0,1],yticklabels=['Old','Expanded'],xticks=range(len(groups)),
                   xticklabels=[f'{p:g} d\n{d} ppm' for p,d in groups],
                   title=f'TIC {tic}, sector {sector}: {label}')
            ax.tick_params(axis='x',labelsize=8)
    fig.suptitle('Two real-noise backgrounds, 16 fixed phases per cell\nOld retention SNR≥8; expanded SNR≥5.5 + ≥60% positive events. Conditional tests, not TESS completeness.',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.89));fig.savefig(OUT/'injection_recovery.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,4,figsize=(14,4),sharey=False)
    controls=json.loads((DATA/'actual_controls_new.json').read_text())['controls']
    for ax,c in zip(axes,controls):
        ss=c['sectors'];x=list(range(len(ss)))
        for i,width in enumerate([.75,1.5]):
            y=[s['fixed_ephemeris'][i]['local_depth']*1e6 for s in ss]
            e=[s['fixed_ephemeris'][i]['local_depth_err']*1e6 for s in ss]
            ax.errorbar(np.array(x)+(i-.5)*.1,y,yerr=e,fmt='o',ms=4,label=f'{width:g} d detrend')
        pooled=c['combined_fixed'][0]
        ax.axhline(pooled['depth_ppm'],color='gray',lw=1)
        ax.axhspan(pooled['depth_ppm']-pooled['depth_error_ppm'],pooled['depth_ppm']+pooled['depth_error_ppm'],color='gray',alpha=.15)
        ax.set(xticks=x,xticklabels=[str(s['sector']) for s in ss],xlabel='Sector',ylabel='Local depth (ppm)',
               title=c['name'].replace(' reported candidate','')+f'\nPooled diagonal SNR {pooled["diagonal_snr"]:.1f}')
    axes[0].legend(fontsize=8)
    fig.suptitle('Published fixed-ephemeris consistency controls; correlated detrenders are not independent evidence',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.86));fig.savefig(OUT/'fixed_ephemeris_controls.png',dpi=160);plt.close(fig)
    phase_offset_control()
    inputs=['data/observations.json','data/known_planets.csv','data/catalogs/confirmed_transiting.csv',
            'candidates/tic4206066/recovery.json']
    inputs.extend(str(p.relative_to(ROOT)) for p in (ROOT/'candidates/tic4206066').glob('*manifest.json'))
    products=[p for directory in [DATA,OUT] for p in directory.rglob('*')
              if p.is_file() and p.name!='checkpoint.json' and '__pycache__' not in p.parts]
    (DATA/'checkpoint.json').write_text(json.dumps(dict(
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source},
        original_input_manifests={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        artifacts=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                        bytes=p.stat().st_size) for p in sorted(products)],
        verification='Every original observation is SHA256 and FITS-structure checked before analysis; existing FITS CHECKSUM/DATASUM checked where present. Original QLP can lack FITS checksum cards.',
        known_limitations='Public archive and existing raw-cache persistence are inherited from the project; original retrieval URLs and exact SHA256 retained in observation outputs.'),indent=2)+'\n')
    print(json.dumps(dict(injections=len(audits),actual_sector_cases=len(actual),groups=len(rows))))


if __name__=='__main__': main()
