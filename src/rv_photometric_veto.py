"""Independent TESS variability checks for the two sparse RV leads.

This targets pulsation/activity alternatives, not planet discovery. A peak in
one TESS sector alone cannot identify the RV mechanism or exclude a companion.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.ndimage import median_filter
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import verify_fits_checksums

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/new_stars/rv/tess'
OUT=ROOT/'candidates/new_stars/rv'


def main():
    results=[];curves=[]
    for manifest in sorted(DATA.glob('*_qlp_observations.json')):
        for record in json.loads(manifest.read_text()):
            if record['http_status']!=200:continue
            path=Path(record['path'])
            if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('Raw SHA mismatch')
            with fits.open(path)as hd:
                integrity=verify_fits_checksums(path,hd);d=hd[1].data
                key='DET_FLUX'if'DET_FLUX'in d.names else'KSPSAP_FLUX'
                good=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d['SAP_FLUX'])&np.isfinite(d[key])
                t=np.asarray(d['TIME'][good],float);raw=np.asarray(d['SAP_FLUX'][good],float);corr=np.asarray(d[key][good],float)
                exposure=float(hd[1].header['TIMEDEL'])
                if hd[0].header.get('TIMESYS',hd[1].header.get('TIMESYS'))!='TDB':raise ValueError('Unexpected time standard')
            raw/=np.median(raw);corr/=np.median(corr)
            n=max(5,int(round(1/exposure))//2*2+1)
            high=raw/median_filter(raw,size=n,mode='nearest')-1
            mad=1.4826*np.median(abs(high-np.median(high)));keep=abs(high-np.median(high))<8*mad
            t=t[keep];high=high[keep];corr=corr[keep];raw=raw[keep]
            max_frequency=min(70.,.45/exposure)
            frequency=np.arange(.3,max_frequency,1/(np.ptp(t)*8))
            powers={};measurements=[]
            for name,y in [('raw_highpass',high),('corrected',corr-1)]:
                ls=LombScargle(t,y);power=ls.power(frequency);powers[name]=power
                order=np.argsort(power)[::-1];selected=[]
                for j in order:
                    if all(abs(frequency[j]-frequency[k])>1/np.ptp(t)for k in selected):selected.append(j)
                    if len(selected)==5:break
                peaks=[]
                for j in selected:
                    params=ls.model_parameters(frequency[j]);amp=float(np.hypot(params[1],params[2]))
                    peaks.append(dict(frequency_per_day=float(frequency[j]),period_days=float(1/frequency[j]),power=float(power[j]),sinusoid_semiamplitude_ppm=amp*1e6))
                measurements.append(dict(method=name,peaks=peaks))
            result=dict(tic=record['tic'],sector=record['sector'],raw=record,checksums=integrity,n_points=len(t),
                exposure_days=exposure,baseline_days=float(np.ptp(t)),highpass_width_days=1.,
                measurements=measurements,limitations=['One-day raw median detrending attenuates slow rotation and can distort gap edges.',
                'Corrected aperture flux and raw flux are the same observations; neither localizes an astrophysical variable.',
                'Peak powers are targeted descriptive diagnostics, not calibrated survey false-positive probabilities.'])
            results.append(result);curves.append((result,t,high,frequency,powers))
            print(record['tic'],record['sector'],[(m['method'],m['peaks'][0])for m in measurements],flush=True)
    fig,axes=plt.subplots(len(curves),2,figsize=(12,3.1*len(curves)),squeeze=False)
    for ax,(r,t,y,f,pow)in zip(axes,curves):
        ax[0].scatter(t,y*1000,s=1);ax[0].set(xlabel='BTJD TDB',ylabel='Raw high-pass flux (ppt)',title=f'TIC{r["tic"]},sector{r["sector"]}')
        for name,p in pow.items():ax[1].plot(f,p,label=name,lw=.8)
        ax[1].set(xlabel='Frequency (cycles/day)',ylabel='Lomb–Scargle power');ax[1].legend(fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'photometric_veto.png',dpi=160);plt.close(fig)
    summary=dict(checked_utc=datetime.now(timezone.utc).isoformat(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),lightcurves=results)
    (OUT/'photometric_veto.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
