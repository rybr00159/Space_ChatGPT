#!/usr/bin/env python3
"""Derive transparent counts and ranked leads without promoting discoveries."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/new_stars'


def read(path):return json.loads((ROOT/path).read_text())


def field_controls():
    table=pd.read_csv(ROOT/'data/new_stars/rv/same_field_epoch_controls.csv',dtype={'GaiaDR3':str,'GALAH':str})
    good=table[(table.flagSp==0)&(table.flagRed==0)&table.Teff.between(4000,7000)&(table.logg>3.7)&
               (table.SNRpxCCD2>30)&(table.e_RVcomp1>0)&(table.e_RVcomp1<2)].copy()
    good['prefix']=good.GALAH.str[:12];results={}
    for name,old,new in [('field543',['140312002701','140314003201'],['230113003201']),
                         ('field469',['140806004701','140823002701'],['230801003301','230803003801'])]:
        a=good[good.prefix.isin(old)].groupby('GaiaDR3').RVcomp1.median()
        b=good[good.prefix.isin(new)].groupby('GaiaDR3').RVcomp1.median()
        values=(b-a).dropna();median=float(values.median());scatter=float(1.4826*np.median(abs(values-median)))
        results[name]=dict(paired_stars=len(values),median_new_minus_old_kms=median,robust_scatter_kms=scatter,
                           caveat='Control includes real binaries and selection effects; scatter is not a calibrated survey false-alarm model.')
    return results


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    im1=read('data/new_stars/imaging/pilot12_v2/summary.json');im2=read('data/new_stars/imaging/expanded64_v4/summary.json')
    selection=read('data/new_stars/infrared/selection.json')
    gaia=pd.read_csv(ROOT/'data/new_stars/infrared/gaia_crossmatch.csv')
    epoch=pd.read_csv(ROOT/'data/new_stars/infrared/independent_epoch_checks.csv')
    novelty=pd.read_csv(ROOT/'data/new_stars/infrared/novelty_checks.csv')
    optical=pd.read_csv(ROOT/'data/new_stars/infrared/optical_epoch_checks.csv')
    pairs=read('candidates/new_stars/companions/recognition_and_robustness.json')
    quality=read('candidates/new_stars/companions/quality_and_catalog_scope.json')
    rv=read('candidates/new_stars/rv/screen_summary.json')
    strict=pd.read_csv(ROOT/'data/new_stars/rv/strict_unflagged_repeat_groups.csv')
    stars=pd.read_csv(ROOT/'data/new_stars/rv/selected_allstar_classification.csv')
    flags=stars.flagSp.to_numpy(int);recognized=int(np.sum((flags&(4|8|128))!=0))
    imaging=dict(distinct_primary_fields=im1['analyzed_fields']+im2['analyzed_fields'],
       pilot_fields=im1['analyzed_fields'],expanded_fields=im2['analyzed_fields'],
       pilot_peaks=im1['detection_count'],expanded_peaks=im2['detection_count'],
       gaia_missing_but_deeper_catalogued_peaks=im1['counts']['catalogued_deeper_survey']+im2['counts']['catalogued_deeper_survey'],
       validated_uncatalogued_stars=0,one_faint_peak_unverified=True,expanded_last_peak_saturated_in_grz=True)
    infrared=dict(selected_catwise_rows=len(gaia),successful_fields=10,failed_fixed_fields=2,
       successful_solid_angle_deg2=selection['successfully_queried_solid_angle_deg2'],
       gaia_unmatched=int((gaia.gaia_count==0).sum()),gaia_query_failures=int((gaia.gaia_count<0).sum()),
       epoch_tested=len(epoch),gaia_unmatched_epoch_consistent=int((novelty.role=='gaia_unmatched').sum()),
       epoch_consistent_unmatched_already_tic_star=int(((novelty.role=='gaia_unmatched')&(novelty.tic_classification=='STAR')).sum()),
       optical_tested=len(optical),apparent_optical_motion_matches=int((optical.optical_motion_interpretation=='motion_prediction_supported').sum()),
       last_optical_match_falsified_as_stationary_sdss_galaxy=True,validated_new_stars=0)
    rvcounts=dict(repeat_stars=rv['repeat_stars'],repeat_quality_visits=rv['repeat_visits'],
       strong_amplitude_screened=rv['selected_stars'],already_binary_flagged=recognized,
       remaining_strong_lead_rejected_for_reduction_quality=True,strict_unflagged_repeat_stars=len(strict),
       strict_amplitude_above3kms=int(((strict.rvmax-strict.rvmin)>3).sum()),
       documented_gaia_nss_solutions=3,additional_gaia_rv_variability_cases=2,
       provisional_sparse_rv_leads=2,component_masses_measured=0)
    summary=dict(method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),imaging=imaging,
       infrared=infrared,astrometric_companions=read('candidates/new_stars/companions/summary.json'),
       spectroscopy=rvcounts,same_field_rv_controls=field_controls(),
       conclusion='No previously unknown/uncatalogued star established. Two cleaner possible unreported associations and two provisional RV leads remain among already catalogued stars.',
       counting_caveat='Tracks overlap and have different selection units; do not add their star/peak/catalog-query counts into a claimed independent survey total.')
    (OUT/'search_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    ranked=[]
    for rank,number in enumerate([1,3,2,5,9,11],1):
        p=next(x for x in pairs if x['pair_number']==number);q=next(x for x in quality if x['pair_number']==number)
        ranked.append(dict(track='Gaia_association',rank_within_track=rank,source=p['primary'],other_source=p['neighbor'],
            assessment=q['assessment'],previously_uncatalogued=False,priority_proven=False,
            key_measurement=f"separation={p['separation_arcsec']:.6f}arcsec; projected={p['projected_separation_au']:.1f}AU",
            evidence='candidates/new_stars/companions/README.md'))
    for rank,sid in enumerate(['4690783166942834816','5437472031746145536'],1):
        ranked.append(dict(track='sparse_RV',rank_within_track=rank,source=sid,other_source='',
            assessment='Unexplained long-term velocities; region-dependent spectral shifts, sparse cadence and stellar variability prevent stellar-companion classification.',
            previously_uncatalogued=False,priority_proven=False,key_measurement='~7.75km/s or ~5.80km/s catalogue range; no orbital period or mass',evidence='reports/new_stars/spectroscopy.md'))
    with (OUT/'ranked_leads.csv').open('w',newline='')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(ranked[0]),lineterminator='\n');writer.writeheader();writer.writerows(ranked)
    print(json.dumps({k:v for k,v in summary.items()if k!='astrometric_companions'},indent=2))


if __name__=='__main__':main()
