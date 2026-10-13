"""Search contiguous-sector residuals; deliberately avoid huge multi-year grids."""
from pathlib import Path
import argparse,json
import numpy as np
from search import ROOT, CACHE, flatten, read_curve, known_for, mask_known, scan, phase
import matplotlib.pyplot as plt

def grouped(records):
    groups=[]
    for r in sorted(records,key=lambda r:r['sector']):
        if groups and r['sector']<=groups[-1][-1]['sector']+1: groups[-1].append(r)
        else: groups.append([r])
    return [g for g in groups if len(g)>1]

def process_group(group,width,minimum=.7,maximum=25):
    tic=group[0]['tic'];sectors=[r['sector'] for r in group]
    chunks=[]
    for r in group:
        t,f,e,h=read_curve(r['local_path']);tb,fb,eb=flatten(t,f,e,width)
        keep=mask_known(tb,known_for(tic));chunks.append((tb[keep],fb[keep],eb[keep]))
    t,f,e=(np.concatenate([c[j] for c in chunks]) for j in range(3))
    peaks,r=scan(t,f,e,minimum,maximum,n_peaks=5)
    label=f"tic{tic}_s{'-'.join(map(str,sectors))}_w{width:g}"
    out=ROOT/'results'/'joint';out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(CACHE/'intermediate'/f'{label}_joint.npz',t=t,f=f,e=e)
    result={'tic':tic,'sectors':sectors,'width':width,'n_searched':len(t),'t_min':float(t.min()),'t_max':float(t.max()),'peaks':peaks,'intermediate':str(CACHE/'intermediate'/f'{label}_joint.npz')}
    (out/f'{label}.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    if peaks:
        best=peaks[0];fig,ax=plt.subplots(3,1,figsize=(10,9),layout='constrained')
        ax[0].plot(t,f,'.',ms=1);ax[0].set(xlabel='BTJD',ylabel='Normalized flux',title=label)
        ax[1].plot(r.period,r.power,lw=.5);ax[1].set(xlabel='Period (d)',ylabel='BLS S/N')
        ph=phase(t,best['period'],best['t0_btjd']);ax[2].plot(ph,f,'.',ms=1)
        ax[2].set(xlabel='Phase time (d)',ylabel='Flux',xlim=(-.5,.5),title=f"P{best['period']:.5f} SNR{best['snr']:.1f} events{best['n_events']}")
        fig.savefig(out/f'{label}.png',dpi=120);plt.close(fig)
    print(label,[(round(p['period'],4),round(p['snr'],1),p['n_events']) for p in peaks],flush=True)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--tic',type=int,nargs='+');p.add_argument('--width',type=float,nargs='+',default=[1.,2.]);p.add_argument('--max-period',type=float,default=25)
    a=p.parse_args();rows=json.loads((ROOT/'data'/'observations.json').read_text())
    for tic in sorted({r['tic'] for r in rows}):
        if a.tic and tic not in a.tic:continue
        for group in grouped([r for r in rows if r['tic']==tic]):
            for width in a.width:
                label=f"tic{tic}_s{'-'.join(str(r['sector']) for r in group)}_w{width:g}"
                if (ROOT/'results'/'joint'/f'{label}.json').exists():continue
                process_group(group,width,maximum=a.max_period)
if __name__=='__main__':main()
