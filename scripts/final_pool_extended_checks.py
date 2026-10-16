"""Additional falsification of two final pooled-search leads; no discovery claim.

TIC307085108 has broad recurring troughs at two rationally related trial periods;
TIC183373155 has one deep event per discovery sector and no secure unique period.
Both are investigated rather than vetoed just for warning flags or sparse events.
"""
from pathlib import Path
import json,os,sys,traceback
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import numpy as np
from astropy.timeseries import LombScargle
import multisector_followup as f
from pooled_south_vet import explicit_followup_error_policy,incomplete_observation


def main():
    for tic in [307085108,183373155]:
        output=f.OUT/f'tic{tic}'/'final_pool_extended.json'
        if output.exists():continue
        base=ROOT/'results/multisector_final_pool_followup'/f'tic{tic}.json';d=json.loads(base.read_text());ephem=d['ephemeris'];peak=d['selection'];records=[r['source']for r in d['observations']]
        result=dict(tic=tic,checked_utc=f.utc(),input_path=str(base.relative_to(ROOT)),input_sha256=f.digest(base),script_sha256=f.digest(__file__),ephemeris=ephem,discovery_spectra=[],external_tests=[],errors=[],limitations=['Exploratory follow-up on selected signals; spectra and phase checks are conditional.','Variable-star harmonics and timing uncertainty must be examined before classifying a non-recovery.'])
        for record in records:
            a,meta=f.read_source(record);record.update(ra=meta['metadata']['RA_OBJ'],dec=meta['metadata']['DEC_OBJ'])
            for col in ['corrected','raw']:
                flat=f.flatten(a,col,1.5);t,y,e=[flat[k]for k in ['t','f','e']];oot=abs(f.phase(t,ephem['period'],ephem['epoch_btjd']))>min(ephem['duration_days'],.20*ephem['period']);freq=np.arange(.05,20,1/(10*np.ptp(t)))
                ls=LombScargle(t[oot],y[oot],e[oot]);power=ls.power(freq);j=int(np.argmax(power));fr=float(freq[j]);theta=2*np.pi*(t-t.min())*fr
                design=np.column_stack([np.ones(len(t)),np.sin(theta),np.cos(theta),np.sin(2*theta),np.cos(2*theta)]);weight=design[oot]/e[oot,None];coef=np.linalg.lstsq(weight,y[oot]/e[oot],rcond=None)[0];model=design@coef
                stats=f.event_stats(t,1+y-model,e,ephem['period'],ephem['epoch_btjd'],ephem['duration_days'],include_details=False)
                result['discovery_spectra'].append(dict(sector=record['sector'],column=col,dominant_period_days=1/fr,bls_period_over_dominant_period=ephem['period']*fr,condition=float(np.linalg.cond(weight)),fitted_peak_to_peak=float(np.ptp(model)),residual_primary=stats,warning='Two-harmonic fit to out-of-event data; a flexible baseline can absorb real transits.'))
        try:
            retrieval=f.retrieve_independent(records[0],max_sectors=2,exclude_sectors=[r['sector']for r in records]);result['retrieval']=retrieval
            for record in retrieval['observations']:
                if record.get('status')!='retained':continue
                record=explicit_followup_error_policy(record)
                try:result['external_tests'].append(f.vet_observation(record,peak,output=False,fixed_ephemeris=ephem))
                except Exception as exc:
                    result['errors'].append(dict(sector=record['sector'],error=repr(exc),traceback=traceback.format_exc()));result['external_tests'].append(incomplete_observation(record,ephem,exc))
        except Exception as exc:result['errors'].append(dict(stage='retrieval',error=repr(exc),traceback=traceback.format_exc()))
        f.save(output,result);print('EXTENDED',tic,'spectra',[(r['sector'],r['column'],round(r['dominant_period_days'],6),round(r['bls_period_over_dominant_period'],3),round(r['residual_primary']['local_snr'],2))for r in result['discovery_spectra']],'external',[(r['sector'],r.get('warnings'))for r in result['external_tests']],'errors',result['errors'],flush=True)

if __name__=='__main__':main()
