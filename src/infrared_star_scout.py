#!/usr/bin/env python3
"""Reproducible independent infrared proper-motion reconnaissance.

No object is called an uncatalogued star merely for lacking a Gaia/SIMBAD match.
CatWISE is itself a catalogue; its apparent motions require independent checks.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, hashlib, io, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/new_stars/infrared'
CAND = ROOT / 'candidates/new_stars/infrared'
CACHE = Path('/workspace/tess-research-data/new_stars/infrared')
TAP = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
for directory in (OUT, CAND, CACHE): directory.mkdir(parents=True, exist_ok=True)

def retrieve(label, url, params, timeout=120):
    path = CACHE / (label + '.txt')
    meta = CACHE / (label + '.json')
    if path.exists() and meta.exists():
        record=json.loads(meta.read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('Retained response hash mismatch: '+str(path))
        return path.read_text(), record
    started = time.time()
    try:
        response = requests.get(url, params=params, timeout=timeout)
        content, code = response.text, response.status_code
        final_url = response.url
    except requests.RequestException as exc:
        content, code, final_url = repr(exc), None, url
    record = dict(label=label, url=final_url, params=params, http_status=code,
                  utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                  elapsed_seconds=time.time()-started,
                  sha256=hashlib.sha256(content.encode()).hexdigest(),
                  raw_path=str(path))
    path.write_text(content); meta.write_text(json.dumps(record, indent=2))
    return content, record

def tap(label, query, timeout=120):
    content, meta = retrieve(label, TAP, dict(REQUEST='doQuery', LANG='ADQL',
                                             FORMAT='csv', QUERY=query), timeout)
    if meta['http_status'] != 200 or content.lstrip().startswith('<'):
        raise RuntimeError(f'{label}: invalid TAP response {meta["http_status"]}')
    # Gaia source IDs are identifiers, never floating-point measurements.
    frame = pd.read_csv(io.StringIO(content), dtype={'Source': 'string'})
    return frame, meta

def fixed_fields():
    # Galactic poles-adjacent fields avoid Galactic-plane crowding. These are
    # fixed BEFORE looking at results, not cones chosen around intriguing stars.
    rows = []
    for latitude in (-60, 60):
        for longitude in range(0, 360, 60):
            c = SkyCoord(l=longitude*u.deg, b=latitude*u.deg, frame='galactic').icrs
            rows.append(dict(label=f'l{longitude:03d}_b{latitude:+03d}',
                             ra=c.ra.deg, dec=c.dec.deg, radius_deg=5.,
                             galactic_longitude=longitude, galactic_latitude=latitude))
    return rows

def field_query(field):
    return f'''select top 2000 objID,Name,RA_ICRS,DE_ICRS,RAPMdeg,DEPMdeg,
pmRA,pmDE,e_pmRA,e_pmDE,W1mproPM,W2mproPM,snrW1pm,snrW2pm,
chi2W1pm,chi2W2pm,ccf,abf,pmQual,MJD from "II/365/catwise"
where 1=contains(point('ICRS',RA_ICRS,DE_ICRS),circle('ICRS',
{field['ra']:.9f},{field['dec']:.9f},{field['radius_deg']}))
and pmRA*pmRA+pmDE*pmDE>0.0225
and pmRA*pmRA/(e_pmRA*e_pmRA)+pmDE*pmDE/(e_pmDE*e_pmDE)>100
and snrW1pm>30 and snrW2pm>15
and chi2W1pm<1.5 and chi2W2pm<1.5
and abf='00' and (ccf='0000' or ccf is null)
and W1mproPM-W2mproPM between -0.1 and 1.0 and W1mproPM>12'''

def fetch():
    fields = fixed_fields(); outputs=[]; records=[]
    def one(field):
        try:
            table, meta = tap('catwise_'+field['label'], field_query(field))
            table['field'] = field['label']
            meta['row_count'] = len(table); meta['field'] = field
            return table, meta
        except Exception as exc:
            return pd.DataFrame(), dict(field=field, error=repr(exc))
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for table, meta in pool.map(one, fields):
            outputs.append(table); records.append(meta)
            print(meta.get('field'), meta.get('row_count', meta.get('error')), flush=True)
    good=[t for t in outputs if len(t)]
    allrows=pd.concat(good, ignore_index=True) if good else pd.DataFrame()
    if len(allrows):
        allrows=allrows.drop_duplicates('objID').reset_index(drop=True)
        allrows['motion_arcsec_per_year']=np.hypot(allrows.pmRA, allrows.pmDE)
        allrows['motion_snr_diagonal']=np.sqrt((allrows.pmRA/allrows.e_pmRA)**2+
                                             (allrows.pmDE/allrows.e_pmDE)**2)
        allrows['W1_W2']=allrows.W1mproPM-allrows.W2mproPM
    allrows.to_csv(OUT/'catwise_screen.csv', index=False)
    (OUT/'field_manifest.json').write_text(json.dumps(records, indent=2))
    (OUT/'selection.json').write_text(json.dumps(dict(
        fields=fields,
        planned_solid_angle_deg2=12*2*np.pi*(1-np.cos(np.deg2rad(5)))*(180/np.pi)**2,
        successfully_queried_solid_angle_deg2=sum('row_count' in r for r in records)*2*np.pi*(1-np.cos(np.deg2rad(5)))*(180/np.pi)**2,
        unique_selected_rows=len(allrows),
        caveat='Count is returned quality/motion-selected rows, not all underlying CatWISE detections. '
               'Motion errors treated diagonally; fitted motion and significance are not independent proof of real motion.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2))

def gaia_query(ra, dec, radius_arcsec=5):
    return f'''select Source,RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,e_pmRA,pmDE,e_pmDE,
Gmag,"BP-RP",RUWE from "I/355/gaiadr3" where
1=contains(point('ICRS',RA_ICRS,DE_ICRS),circle('ICRS',{ra},{dec},{radius_arcsec/3600}))'''

def gaia():
    table=pd.read_csv(OUT/'catwise_screen.csv'); records=[]; matches=[]
    def one(pair):
        index,row=pair
        name=f'gaia_{index:04d}'
        try:
            g,meta=tap(name,gaia_query(row.RAPMdeg,row.DEPMdeg))
            meta['catwise_objID']=row.objID; meta['row_count']=len(g)
            result=dict(catwise_objID=row.objID,gaia_count=len(g))
            if len(g):
                c=SkyCoord(row.RAPMdeg*u.deg,row.DEPMdeg*u.deg)
                sep=c.separation(SkyCoord(g.RA_ICRS.to_numpy()*u.deg,g.DE_ICRS.to_numpy()*u.deg)).arcsec
                k=int(np.argmin(sep)); nearest=g.iloc[k]
                result.update({f'gaia_{key}':value for key,value in nearest.items()})
                result['gaia_separation_arcsec']=float(sep[k])
                result['gaia_motion_disagreement_masyr']=float(np.hypot(
                    1000*row.pmRA-nearest.pmRA,1000*row.pmDE-nearest.pmDE))
            return result,meta
        except Exception as exc:
            return dict(catwise_objID=row.objID,gaia_count=-1),dict(error=repr(exc),catwise_objID=row.objID)
    with cf.ThreadPoolExecutor(max_workers=5) as pool:
        for n,(result,meta) in enumerate(pool.map(one,table.iterrows())):
            matches.append(result); records.append(meta)
            if n%20==0: print('Gaia checked',n+1,'of',len(table),flush=True)
    joined=table.merge(pd.DataFrame(matches),left_on='objID',right_on='catwise_objID')
    joined.to_csv(OUT/'gaia_crossmatch.csv',index=False)
    (OUT/'gaia_manifest.json').write_text(json.dumps(records,indent=2,default=str))
    print(joined.groupby('gaia_count').size(),flush=True)

def two_mass_query(ra, dec):
    return f'''select RAJ2000,DEJ2000,"2MASS",Jmag,e_Jmag,Hmag,e_Hmag,
Kmag,e_Kmag,Qflg,Cflg,Xflg,Aflg,JD,Jpsfchi,Hpsfchi,Kpsfchi from "II/246/out"
where 1=contains(point('ICRS',RAJ2000,DEJ2000),circle('ICRS',{ra},{dec},{10/3600}))'''

def counterparts():
    table=pd.read_csv(OUT/'gaia_crossmatch.csv',dtype={'gaia_Source':'string'})
    # All successful no-match cases, plus ten Gaia-matched controls spanning
    # motion/color. Failed Gaia queries are retained but cannot be absences.
    sample=table[table.gaia_count==0].copy()
    controls=table[(table.gaia_count>0)&(table.gaia_motion_disagreement_masyr<50)].sort_values('W1_W2').iloc[::max(1,int((table.gaia_count>0).sum()/10))].head(10)
    sample=pd.concat([sample.assign(role='gaia_unmatched'),controls.assign(role='positive_control')])
    results=[];manifest=[]
    def one(pair):
        index,row=pair
        # CatWISE motion solution positions are explicitly at MJD 57170,
        # not at the per-source mean static-photometry epoch.
        c=SkyCoord(ra=row.RAPMdeg*u.deg,dec=row.DEPMdeg*u.deg,
                   pm_ra_cosdec=row.pmRA*u.arcsec/u.yr,pm_dec=row.pmDE*u.arcsec/u.yr,
                   obstime=Time(57170,format='mjd'))
        pred=c.apply_space_motion(new_obstime=Time(2000.,format='jyear'))
        try:
            t,meta=tap(f'tmass_{index:04d}',two_mass_query(pred.ra.deg,pred.dec.deg))
            meta['catwise_objID']=row.objID;meta['row_count']=len(t)
            result=dict(objID=row.objID,role=row.role,tmass_count=len(t),
                        predicted_ra_2000=pred.ra.deg,predicted_dec_2000=pred.dec.deg)
            if len(t):
                best=None
                for j,star in t.iterrows():
                    epoch=Time(star.JD,format='jd')
                    predicted=c.apply_space_motion(new_obstime=epoch)
                    observed=SkyCoord(star.RAJ2000*u.deg,star.DEJ2000*u.deg)
                    sep=predicted.separation(observed).arcsec
                    if best is None or sep<best[0]:best=(sep,j,predicted)
                sep,j,predicted=best;star=t.iloc[j]
                result.update({f'tmass_{key}':value for key,value in star.items()})
                dt=(57170+2400000.5-float(star.JD))/365.25
                dra=(row.RAPMdeg-star.RAJ2000)*3600*np.cos(np.deg2rad(row.DEPMdeg))
                ddec=(row.DEPMdeg-star.DEJ2000)*3600
                result.update(tmass_prediction_residual_arcsec=sep,
                              baseline_years=dt,independent_pmra=dra/dt,
                              independent_pmdec=ddec/dt,
                              raw_epoch_displacement_arcsec=np.hypot(dra,ddec),
                              J_W2=float(star.Jmag-row.W2mproPM))
            return result,meta
        except Exception as exc:
            return dict(objID=row.objID,role=row.role,tmass_count=-1),dict(error=repr(exc),objID=row.objID)
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for result,meta in pool.map(one,sample.iterrows()):results.append(result);manifest.append(meta)
    final=sample.merge(pd.DataFrame(results),on=['objID','role'])
    final.to_csv(OUT/'independent_epoch_checks.csv',index=False)
    (OUT/'independent_epoch_manifest.json').write_text(json.dumps(manifest,indent=2,default=str))
    print(final[['Name','role','tmass_count','tmass_prediction_residual_arcsec','J_W2']].to_string(index=False),flush=True)

def novelty():
    table=pd.read_csv(OUT/'independent_epoch_checks.csv',dtype={'gaia_Source':'string'})
    # These numerical cuts only prioritize follow-up. The response of every
    # tested source is retained, including galaxies and known brown dwarfs.
    sample=table[(table.tmass_count>0)&
                 (table.tmass_prediction_residual_arcsec<1.)&
                 (table.raw_epoch_displacement_arcsec>1.5)].copy()
    manifest=[];results=[]
    def one(pair):
        index,row=pair
        source_id=str(row.objID)
        label=hashlib.sha256(source_id.encode()).hexdigest()[:12]
        query=f'''select top 20 TIC,RAJ2000,DEJ2000,"S/G",Teff,Rad,Mass,Tmag,"2MASS"
from "IV/39/tic82" where 1=contains(point('ICRS',RAJ2000,DEJ2000),
circle('ICRS',{row.tmass_RAJ2000},{row.tmass_DEJ2000},{3/3600}))'''
        result=dict(objID=source_id)
        try:
            tic,tm=tap('novelty_tic_'+label,query,30)
            exact=tic[tic['2MASS'].astype(str).str.strip()==str(row['tmass_2MASS']).strip()]
            result['tic_exact_2mass_matches']=len(exact)
            if len(exact):
                result.update(tic_id=str(exact.TIC.iloc[0]),
                              tic_classification=str(exact['S/G'].iloc[0]).strip(),
                              tic_teff=exact.Teff.iloc[0])
        except Exception as exc:tm=dict(error=repr(exc));result['tic_exact_2mass_matches']=-1
        query=f'''select main_id,otype,ra,dec,pmra,pmdec,plx_value,sp_type,nbref
from basic where contains(point('ICRS',ra,dec),circle('ICRS',
{row.tmass_RAJ2000},{row.tmass_DEJ2000},{20/3600}))=1'''
        content,sm=retrieve('novelty_simbad_'+label,
            'https://simbad.cds.unistra.fr/simbad/sim-tap/sync',
            dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query),30)
        try:
            if sm['http_status']!=200 or content.lstrip().startswith('<'):raise ValueError('invalid response')
            sim=pd.read_csv(io.StringIO(content));result['simbad_cone_count']=len(sim)
            if len(sim):
                c=SkyCoord(row.tmass_RAJ2000*u.deg,row.tmass_DEJ2000*u.deg)
                sep=c.separation(SkyCoord(sim.ra.to_numpy()*u.deg,sim.dec.to_numpy()*u.deg)).arcsec
                j=int(np.argmin(sep));near=sim.iloc[j]
                result.update(simbad_nearest_separation_arcsec=float(sep[j]),
                    simbad_main_id=near.main_id,simbad_otype=near.otype,
                    simbad_spectral_type=near.sp_type,simbad_nbref=near.nbref)
        except Exception as exc:sm['parse_error']=repr(exc);result['simbad_cone_count']=-1
        return result,dict(objID=source_id,tic=tm,simbad=sm)
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for result,meta in pool.map(one,sample.iterrows()):results.append(result);manifest.append(meta)
    final=sample.merge(pd.DataFrame(results),on='objID')
    final.to_csv(OUT/'novelty_checks.csv',index=False)
    (OUT/'novelty_manifest.json').write_text(json.dumps(manifest,indent=2,default=str))
    print('Independently consistent selected rows:',len(final),flush=True)
    print(final[['Name','role','tic_id','tic_classification','simbad_main_id','simbad_otype','simbad_spectral_type']].to_string(index=False),flush=True)

def optical():
    table=pd.read_csv(OUT/'independent_epoch_checks.csv')
    # A predeclared numerical follow-up subset of the unmatched sample. Its
    # thresholds were chosen after seeing the pilot, so this is NOT a blind
    # occurrence-rate experiment.
    sample=table[(table.role=='gaia_unmatched')&(table.tmass_count==0)&
                 (table.motion_arcsec_per_year>.38)&(table.W1_W2<.4)].copy()
    results=[];records=[]
    def one(pair):
        index,row=pair;label=hashlib.sha256(str(row.objID).encode()).hexdigest()[:12]
        if row.DEPMdeg>=-30:
            columns=['objID','raMean','decMean','raMeanErr','decMeanErr','epochMean',
                     'nDetections','gMeanPSFMag','rMeanPSFMag','iMeanPSFMag',
                     'zMeanPSFMag','yMeanPSFMag','rMeanKronMag']
            url='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv'
            params=dict(ra=row.RAPMdeg,dec=row.DEPMdeg,radius=5/3600,
                        **{'nDetections.gte':3,'columns':json.dumps(columns)})
            rac,decc,epochc='raMean','decMean','epochMean'
            survey='PS1_DR2'
        else:
            url='https://api.skymapper.nci.org.au/public/tap/sync'
            params=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',
                QUERY=f"select top 20 * from dr4.master where 1=contains(point('ICRS',raj2000,dej2000),circle('ICRS',{row.RAPMdeg},{row.DEPMdeg},{5/3600}))")
            rac,decc,epochc='raj2000','dej2000','mean_epoch';survey='SkyMapper_DR4'
        raw,meta=retrieve('optical_'+label,url,params,40)
        result=dict(objID=row.objID,survey=survey)
        try:
            if meta['http_status']!=200 or raw.lstrip().startswith('<'):raise ValueError('invalid response')
            t=pd.read_csv(io.StringIO(raw));result['optical_count']=len(t)
            if len(t):
                current=SkyCoord(row.RAPMdeg*u.deg,row.DEPMdeg*u.deg)
                pos=SkyCoord(t[rac].to_numpy()*u.deg,t[decc].to_numpy()*u.deg)
                sep=current.separation(pos).arcsec;j=int(np.argmin(sep));source=t.iloc[j]
                dt=(source[epochc]-57170)/365.25
                dra=(source[rac]-row.RAPMdeg)*3600*np.cos(np.deg2rad(row.DEPMdeg))
                ddec=(source[decc]-row.DEPMdeg)*3600
                predicted=np.hypot(row.pmRA*dt,row.pmDE*dt)
                residual=np.hypot(dra-row.pmRA*dt,ddec-row.pmDE*dt)
                result.update(optical_ra=source[rac],optical_dec=source[decc],
                    optical_epoch_mjd=source[epochc],baseline_from_catwise_years=dt,
                    static_position_residual_arcsec=sep[j],
                    moving_prediction_residual_arcsec=residual,
                    predicted_displacement_arcsec=predicted,
                    dx_observed_arcsec=dra,dy_observed_arcsec=ddec,
                    class_star=source.get('class_star',None),
                    optical_motion_interpretation=(
                        'stationary_position_preferred' if sep[j]<.4 and residual>.8 else
                        'motion_prediction_supported' if residual<.4 and sep[j]>.8 else
                        'ambiguous_or_blended'))
        except Exception as exc:result['optical_count']=-1;meta['error']=repr(exc)
        return result,dict(objID=row.objID,**meta)
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for result,meta in pool.map(one,sample.iterrows()):results.append(result);records.append(meta)
    final=sample.merge(pd.DataFrame(results),on='objID')
    final.to_csv(OUT/'optical_epoch_checks.csv',index=False)
    (OUT/'optical_manifest.json').write_text(json.dumps(records,indent=2,default=str))
    print(final[['Name','survey','optical_count','static_position_residual_arcsec',
                 'moving_prediction_residual_arcsec','optical_motion_interpretation']].to_string(index=False),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['fetch','gaia','counterparts','novelty','optical'])
    args=p.parse_args();globals()[args.stage]()
if __name__=='__main__': main()
