"""Regenerate compact moving-block diagnostic outputs, keeping core controls frozen."""
from pathlib import Path
import csv,hashlib,json,os
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
DATA=ROOT/'results/multisector_controls/nulls'
rows=[json.loads(x) for x in (DATA/'iterations.jsonl').read_text().splitlines()]
keys=['case_id','tic','sector','iteration','seed','max_local_snr','max_retained_local_snr',
      'any_retained','retained_peak_records','n_input','surrogate_path','surrogate_sha256']
with (DATA/'iterations.csv').open('w',newline='')as stream:
    writer=csv.DictWriter(stream,fieldnames=keys,lineterminator='\n');writer.writeheader()
    writer.writerows({k:r[k] for k in keys} for r in rows)
tics=list(dict.fromkeys(r['tic']for r in rows))
fig,ax=plt.subplots(figsize=(10,4.5))
for i,tic in enumerate(tics):
    rr=[r for r in rows if r['tic']==tic]
    y=[r['max_local_snr'] for r in rr]
    colors=['tab:red'if r['any_retained']else'tab:blue'for r in rr]
    ax.scatter(i+np.linspace(-.15,.15,len(rr)),y,c=colors,s=25)
ax.axhline(5.5,color='gray',ls='--',label='Retention threshold 5.5')
ax.axhline(8,color='black',ls=':',label='SNR 8 reference')
ax.set(xticks=range(len(tics)),xticklabels=[str(t)for t in tics],ylabel='Maximum reported local SNR across both widths',
       xlabel='First eight pre-flux SHA-ordered sector-83 targets',
       title='32 moving-block surrogate searches; red = at least one retained seed\nAstrophysical fragments and block artifacts remain: not a formal false-alarm calibration')
ax.tick_params(axis='x',labelsize=8);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(OUT/'surrogate_diagnostic.png',dpi=160);plt.close(fig)
sources=['src/multisector_nulls.py','candidates/multisector_controls/summarize_nulls.py']
artifacts=[p for p in DATA.rglob('*')if p.is_file() and p.name!='checkpoint.json']
artifacts.extend([OUT/'surrogate_diagnostic.png',OUT/'NULLS.md'])
manifest=dict(methods=json.loads((DATA/'config.json').read_text())['methods_sha256'],
              sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},
              artifacts=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size)for p in artifacts],
              surrogate_arrays=[dict(path=r['surrogate_path'],sha256=r['surrogate_sha256'])for r in rows],
              core_controls_checkpoint_sha256=hashlib.sha256((DATA.parent/'checkpoint.json').read_bytes()).hexdigest())
(DATA/'checkpoint.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Saved',len(rows),'surrogate rows and checked-source manifest')
