"""Independent local-baseline check of final gap-trimmed TOI-4342 joint maximum."""
import json
import numpy as np
from search import ROOT,read_curve,flatten,known_for,mask_known,phase
from validate import local_depth,noise_beta
from vet_tic150360085 import quadratic_depth,gap_keep
import matplotlib.pyplot as plt

def main():
    tic=354944123
    selected=json.loads((ROOT/'results'/'robust'/'tic354944123_s101-102-103_w2_edge0.25.json').read_text())['peaks'][0]
    p,t0,d=selected['period'],selected['t0_btjd'],selected['duration_days']
    records=[r for r in json.loads((ROOT/'data'/'observations.json').read_text()) if r['tic']==tic]
    output={'tic':tic,'selected_peak':selected,'interpretation':'unpromoted; local baseline/raw flux diagnostics below','measurements':[]}
    fig,axes=plt.subplots(3,1,figsize=(10,9),layout='constrained')
    for col in ('PDCSAP_FLUX','SAP_FLUX'):
        for width in (1,2):
            chunks=[]
            for rec in records:
                t,f,e,_=read_curve(rec['local_path'],col);tb,fb,eb=flatten(t,f,e,width)
                keep=gap_keep(t,tb)&mask_known(tb,known_for(tic))
                chunks.append((tb[keep],fb[keep],eb[keep]))
            t,f,e=(np.concatenate([c[j] for c in chunks]) for j in range(3))
            linear=local_depth(t,f,e,p,t0,d);quadratic=quadratic_depth(t,f,e,p,t0,d)
            output['measurements'].append({'column':col,'width':width,'linear':linear,'quadratic':quadratic,'noise':noise_beta(t,f,e,minutes=round(d*1440))})
            print(col,width,'linear',linear.get('depth',0)*1e6,linear.get('formal_snr'),'quadratic',quadratic.get('depth',0)*1e6,quadratic.get('formal_snr'),flush=True)
            if col=='PDCSAP_FLUX' and width==2:
                axes[0].plot(t,f,'.',ms=1);axes[1].plot(phase(t,p,t0),f,'.',ms=1);axes[1].set_xlim(-d*3,d*3)
                events=linear['events'];axes[2].errorbar([k['center_btjd'] for k in events],[k['depth']*1e6 for k in events],yerr=[k['formal_error']*1e6 for k in events],fmt='o')
    native=[]
    for rec in records:
        tt,ff,ee,_=read_curve(rec['local_path'],'SAP_FLUX')
        scale=np.median(ff);ff=ff/scale;ee=ee/scale
        noise=1.4826*np.median(np.abs(np.diff(ff)-np.median(np.diff(ff))))/np.sqrt(2)
        ee*=max(1,noise/np.median(ee))
        keep=gap_keep(tt,tt)&mask_known(tt,known_for(tic))
        native.append((tt[keep],ff[keep],ee[keep]))
    nt,nf,ne=(np.concatenate([c[j] for c in native]) for j in range(3))
    output['native_sap']={'linear':local_depth(nt,nf,ne,p,t0,d),'quadratic':quadratic_depth(nt,nf,ne,p,t0,d)}
    axes[0].set(xlabel='BTJD',ylabel='PDC normalized flux',title=f'TOI-4342 final joint gap-trim maximum P={p:.5f} d')
    axes[1].set(xlabel='Phase time (days)',ylabel='Flux');axes[2].axhline(0,color='.5');axes[2].set(xlabel='Event center (BTJD)',ylabel='Local depth (ppm)')
    out=ROOT/'candidates'/'known_host_rejected';out.mkdir(parents=True,exist_ok=True)
    (out/'toi4342_joint_local_check.json').write_text(json.dumps(output,indent=2,allow_nan=False));fig.savefig(out/'toi4342_joint_local_check.png',dpi=140);plt.close(fig)
    lines=['# TOI-4342 final joint gap-trimmed peak','',f'P={p:.6f} d, epoch BTJD {t0:.6f}, duration {d*24:.2f} h; screening S/N {selected["snr"]:.2f}, supported events {selected["n_events"]}. Selected using all three sectors, so no independent heldout claim applies.','', '|Reduction|Width (d)|Local baseline|Depth (ppm)|Diagonal S/N|','|---|---|---|---|---|']
    for row in output['measurements']:
        for baseline in ('linear','quadratic'):
            m=row[baseline];lines.append(f'|{row["column"]}|{row["width"]}|{baseline}|{m.get("depth",0)*1e6:.0f}|{m.get("formal_snr",0):.2f}|')
    for baseline in ('linear','quadratic'):
        m=output['native_sap'][baseline];lines.append(f'|Native SAP|—|{baseline}|{m.get("depth",0)*1e6:.0f}|{m.get("formal_snr",0):.2f}|')
    lines+=['','Assessment: the corrected-flux dips lack a positive native-SAP counterpart (linear S/N about -1.37; quadratic about 0.62). No credible candidate is promoted. This diagnoses reduction dependence; it does not prove that the star lacks further planets.','The two supported cycles have only 2 versus 31 ten-minute in-box samples and the same parity; odd/even agreement cannot be evaluated. Partial coverage is not evidence of two well-observed transits.','These are local boxes with separate nuisance baselines, not calibrated false-alarm probabilities. Selection spans the full joint data and multiple grids. SAP/PDC share photons and local fits can still depend on baseline window. No pixel localization or physical transit inference is justified by weak or reduction-dependent support. Original boundary/alias diagnostics remain preserved separately.','', 'Reproduce: `python src/vet_toi4342_joint.py`.']
    (out/'TOI4342_JOINT.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()
