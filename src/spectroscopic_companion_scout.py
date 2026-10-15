"""Fresh companion search in repeat GALAH DR4 radial velocities.

These are catalogued visible stars, not uncatalogued sky sources. A new
radial-velocity variability assessment is only a lead for an unseen companion;
it does not measure a secondary mass or establish first recognition.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import io
import json
import sys
import requests
import numpy as np
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/new_stars/rv'
OUT=ROOT/'candidates/new_stars/rv'
ENDPOINT='https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
TABLE='"J/other/PASA/42.51/allspec"'
QUALITY='flagSp = 0 and e_RVcomp1 > 0 and e_RVcomp1 < 2 and SNRpxCCD2 > 30 and logg > 3.7 and Teff between 4000 and 7000'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def query(name,adql,download=False):
    DATA.mkdir(parents=True,exist_ok=True)
    path=DATA/(name+'.csv');manifest=path.with_suffix('.csv.manifest.json')
    if path.exists():
        meta=json.loads(manifest.read_text())
        if digest(path)!=meta['sha256']or meta['query']!=adql:raise ValueError('Cache identity mismatch '+name)
    elif download:
        params=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=adql)
        response=requests.get(ENDPOINT,params=params,timeout=150)
        path.write_bytes(response.content)
        manifest.write_text(json.dumps(dict(url=response.url,query=adql,http_status=response.status_code,
             bytes=len(response.content),sha256=digest(path),checked_utc=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
        response.raise_for_status()
    else:raise FileNotFoundError(str(path)+'; run --download to retrieve missing public input')
    text=path.read_text()
    if text.lstrip().startswith('<'):raise ValueError('Failed TAP response is not an empty table: '+name)
    reader=csv.DictReader(io.StringIO(text))
    if not reader.fieldnames:raise ValueError('Missing CSV header: '+name)
    return list(reader)


def retrieve(download=False):
    jobs={
        'quality_positive_sample':f'SELECT TOP 10 GALAH, GaiaDR3, MJD, RVcomp1, e_RVcomp1, RVcompNr, flagSp, flagRed, Teff, logg FROM {TABLE} WHERE {QUALITY}',
        'quality_repeat_groups':f'SELECT GaiaDR3, COUNT(*) AS n, MIN(MJD) AS tmin, MAX(MJD) AS tmax, MIN(RVcomp1) AS rvmin, MAX(RVcomp1) AS rvmax FROM {TABLE} WHERE {QUALITY} GROUP BY GaiaDR3 HAVING COUNT(*) >= 3',
        'quality_count':f'SELECT COUNT(*) AS n FROM {TABLE} WHERE {QUALITY}',
    }
    # Sequential public queries avoid sending several full-table aggregates
    # to the archive simultaneously. Cached successful products are reused.
    result={}
    for name,q in jobs.items():
        result[name]=query(name,q,download)
        print(name,'rows',len(result[name]),flush=True)
    return result


def summarize(download=False,max_leads=40):
    OUT.mkdir(parents=True,exist_ok=True)
    initial=retrieve(download)
    groups=[]
    for raw in initial['quality_repeat_groups']:
        r=dict(raw);r.update(n=int(raw['n']),tmin=float(raw['tmin']),tmax=float(raw['tmax']),
                            rvmin=float(raw['rvmin']),rvmax=float(raw['rvmax']))
        r['span_days']=r['tmax']-r['tmin'];r['amplitude_kms']=r['rvmax']-r['rvmin'];groups.append(r)
    groups.sort(key=lambda r:(r['amplitude_kms'],r['n']),reverse=True)
    selected=[r for r in groups if r['amplitude_kms']>=20 and r['span_days']>=2][:max_leads]
    print('Quality repeat stars',len(groups),'selected amplitude leads',len(selected),flush=True)
    if not selected:raise ValueError('No selected leads; preserve aggregate then adjust strategy explicitly')
    ids=','.join(r['GaiaDR3']for r in selected)
    q=f'SELECT GALAH, GaiaDR3, "2MASS", RAJ2000, DEJ2000, MJD, RVcomp1, e_RVcomp1, RVcompNr, RVcomp2, RVGDR3, e_RVGDR3, flagSp, flagRed, Teff, logg, vsini, SNRpxCCD2, VbaryEff FROM {TABLE} WHERE GaiaDR3 IN ({ids})'
    visits=query('selected_repeat_visits',q,download)
    assessments=[]
    for group in selected:
        rows=[r for r in visits if r['GaiaDR3']==group['GaiaDR3']]
        good=[]
        for row in rows:
            try:
                good_quality=int(row['flagSp'])==0 and 0<float(row['e_RVcomp1'])<2 and float(row['SNRpxCCD2'])>30 and float(row['logg'])>3.7 and 4000<=float(row['Teff'])<=7000
                if good_quality:good.append(row)
            except(ValueError,TypeError):pass
        rv=np.array([float(r['RVcomp1'])for r in good]);err=np.array([float(r['e_RVcomp1'])for r in good])
        time=np.array([float(r['MJD'])for r in good])
        # A deliberately conservative >=1 km/s visit floor guards against
        # underestimated formal errors; it is not an empirical full noise model.
        total_error=np.hypot(err,1.)
        mean=float(np.average(rv,weights=total_error**-2));stat=float(np.sum(((rv-mean)/total_error)**2))
        nights=len(np.unique(np.floor(time)))
        first=good[0]
        assessments.append(dict(**group,ra=float(first['RAJ2000']),dec=float(first['DEJ2000']),tmass=first['2MASS'],
             usable_visits=len(good),distinct_mjd_nights=nights,all_visits=len(rows),
             visits_with_nonzero_spectral_flag=sum(int(r['flagSp'])!=0 for r in rows),
             weighted_constant_rv=mean,constant_chi2=stat,constant_dof=len(rv)-1,
             constant_pvalue_conditional=float(chi2.sf(stat,len(rv)-1)),
             velocity_floor_kms=1.,minimum_pair_difference=float(np.max(rv)-np.min(rv)),
             median_teff=float(np.median([float(r['Teff'])for r in good])),
             median_logg=float(np.median([float(r['logg'])for r in good])),
             status='rv_variable_screening_lead_not_new_star_claim',
             limits='Binarity, mass, alias period, calibration, source identity and prior recognition remain to test.'))
    def write_csv(path,rows):
        with path.open('w',newline='')as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    write_csv(OUT/'all_repeat_groups.csv',groups);write_csv(OUT/'rv_variable_leads.csv',assessments)
    summary=dict(checked_utc=datetime.now(timezone.utc).isoformat(),source_sha256=digest(Path(__file__)),
        fixed_quality_predicate=QUALITY,selection='All quality repeat groups>=3, then amplitude>=20km/s and span>=2d; amplitude-ranked maximum40 for detailed follow-up.',
        quality_visits=int(initial['quality_count'][0]['n']),repeat_stars=len(groups),
        repeat_visits=sum(r['n']for r in groups),selected_stars=len(selected),all_selected_visit_rows=len(visits),
        measurements=assessments,limitations=['Visible sources are already catalogued; only unseen multiplicity could be newly recognized.',
        'Quality flags and formal errors do not eliminate all RV failures; raw spectra and repeated epochs are needed.',
        'Conditional constant-velocity p-values do not classify companion mass or establish discovery.',
        'Selection is outcome-ranked; small uncorrected p-values are not novelty evidence.'])
    (OUT/'screen_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');parser.add_argument('--max-leads',type=int,default=40)
    args=parser.parse_args();result=summarize(args.download,args.max_leads)
    print({k:v for k,v in result.items()if k not in ['measurements','limitations']},flush=True)


if __name__=='__main__':main()
