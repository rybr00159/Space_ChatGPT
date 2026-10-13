"""Vetting of two revised TOI-1064 peaks after a strict gap-edge search.

This does not inherit rejection of the earlier boundary maxima. It reads the
exact revised ephemerides and checks each event across reductions and sectors.
Writes only its JSON evidence and TOI1064_TRIMMED.md in known_host_rejected.
"""
from pathlib import Path
import json
from datetime import datetime, timezone
import numpy as np
from search import ROOT, read_curve, flatten, phase, known_for, mask_known
from vet_known_host_leads import safe, fixed_depth, predictions, segments, edge_keep
from vet_tic150299576 import local_raw_box, rednoise_beta

OUT=ROOT/'candidates'/'known_host_rejected'
TIC=79748331


def main():
    leads=[]
    for filename in ['tic79748331_s94_w2_edge0.25.json','tic79748331_s104-105_w2_edge0.25.json']:
        screen=json.loads((ROOT/'results'/'robust'/filename).read_text())
        leads.append({'source_screen':str((ROOT/'results'/'robust'/filename).relative_to(ROOT)),
                      'source_sectors':screen['sectors'],'peak':screen['peaks'][0]})
    records=[rec for rec in json.loads((ROOT/'data'/'observations.json').read_text()) if rec['tic']==TIC]
    known=known_for(TIC);prepared={}
    for rec in records:
        sector=rec['sector']
        for column in ['PDCSAP_FLUX','SAP_FLUX']:
            t,f,e,h=read_curve(rec['local_path'],column);bounds=segments(t)
            rawkeep=mask_known(t,known)
            for width in [1.,2.]:
                tb,fb,eb=flatten(t,f,e,width)
                keep=mask_known(tb,known)
                prepared[(sector,column,width)]={'curve':(tb,fb,eb),'known_mask':keep,'bounds':bounds,
                    'raw_curve_known_masked':(t[rawkeep],f[rawkeep],e[rawkeep]),'raw_source':rec}
    assessments=[]
    for signal in leads:
        lead=signal['peak'];p,t0,d=lead['period'],lead['t0_btjd'],lead['duration_days']
        checks=[];joint=[]
        for (sector,column,width),prep in prepared.items():
            t,f,e=prep['curve'];keep=prep['known_mask'];bounds=prep['bounds']
            events=predictions(t,f,e,p,d,t0,keep,bounds,known)
            rt,rf,re=prep['raw_curve_known_masked']
            for event in events:
                center=event['center_btjd']
                event['raw_known_masked_baselines']=[local_raw_box(rt,rf,re,center,d,degree,halfwindow)
                    for degree in [1,2] for halfwindow in [.65,1.25]]
                # A broad duration response is a variability warning, not a
                # physical ingress/egress fit or a formal model comparison.
                if sector in signal['source_sectors']:
                    event['raw_known_masked_duration_response']=[
                        {'duration_days':duration,**local_raw_box(rt,rf,re,center,duration,2,1.25)}
                        for duration in [.06,.10,.16,.23,.32,.45]]
                    grid=[]
                    for shift in np.linspace(-.3,.3,31):
                        test=local_raw_box(rt,rf,re,center+shift,d,2,1.25)
                        if test.get('snr_diagonal') is not None:
                            grid.append({'shift_days':float(shift),**test})
                    event['local_timing_scan_halfwidth_days']=.3
                    event['local_timing_maximum']=max(grid,key=lambda row:row['snr_diagonal']) if grid else None
            fixed=[]
            for cut in [.25,.5]:
                selected=keep&edge_keep(t,bounds,cut)
                fixed.append({'edge_buffer_days':cut,'n_searched':int(selected.sum()),
                              **fixed_depth(t[selected],f[selected],e[selected],p,d,t0)})
            checks.append({'sector':sector,'column':column,'width_days':width,'fixed_ephemeris':fixed,
                           'individual_events':events,'correlated_noise':rednoise_beta(t[keep],f[keep],p,t0,d)})
        for column in ['PDCSAP_FLUX','SAP_FLUX']:
            for width in [1.,2.]:
                for cut in [.25,.5]:
                    chunks=[]
                    for sector in signal['source_sectors']:
                        prep=prepared[(sector,column,width)];t,f,e=prep['curve']
                        selected=prep['known_mask']&edge_keep(t,prep['bounds'],cut)
                        chunks.append((t[selected],f[selected],e[selected]))
                    combined=[np.concatenate([arr[i] for arr in chunks]) for i in range(3)]
                    joint.append({'column':column,'width_days':width,'edge_buffer_days':cut,
                                  **fixed_depth(*combined,p,d,t0)})
        assessments.append({**signal,'checks':checks,'source_group_fixed_ephemeris':joint,
            'status':'unpromoted: fails agreement across reduction/detrending choices',
            'period_ambiguity':'No full transit timing posterior was fitted. Source-event timing searches are conditional on a box and baseline. A distant-season fixed-ephemeris null cannot exclude every nearby period alias.'})
    report={'tic':TIC,'created_utc':datetime.now(timezone.utc).isoformat(),'known_planets':known,
        'observations':[{'sector':rec['sector'],'data_uri':rec['data_uri'],'sha256':rec['sha256']} for rec in records],
        'signals':assessments,
        'conclusion':'Neither revised peak warrants a credible new planetary candidate. This conclusion is specific to the tested signals; it does not exclude additional planets around TOI-1064.',
        'limits':['Diagonal uncertainties understate correlated errors; red-noise binning is a diagnostic rather than a complete noise likelihood.',
                  'No pixel-level independent aperture or contamination localization was attempted for these unstable low-level residuals.',
                  'SAP can retain additive background effects, so a reduction disagreement is a reason to withhold promotion, not standalone proof of astrophysical absence.',
                  'Timing scans add trials. No global false-alarm or novelty claim is assigned.',
                  'Duration scans use boxes and polynomial baselines. No physical radius/eccentricity inference is warranted.']}
    (OUT/'toi1064_revised_trimmed_vetting.json').write_text(json.dumps(safe(report),indent=2,allow_nan=False)+'\n')
    lines=['# TOI-1064: separate vetting of revised gap-trimmed peaks','',
        'The 8.87447 d sector-94 and 14.85263 d joint sector-104/105 peaks were re-examined independently after the original boundary peaks were rejected. Neither revised peak is currently a credible new planetary candidate. This finding applies to these signals and does not exclude all additional planets in the system.', '',
        'Every test includes the known b/c transit masks. Raw SAP and PDC local models omit masked cadences before fitting. Both 1 d and 2 d rolling-median reductions, 0.25 d and 0.5 d good-segment edge exclusions, individual event depths, local polynomial baselines, duration response, and timing-tolerance searches are recorded in `toi1064_revised_trimmed_vetting.json`.', '',
        '| Source period (d) | Flux | Median width (d) | Source-group fixed S/N, edge 0.25 d | Source-group fixed S/N, edge 0.5 d |',
        '| --- | --- | --- | --- | --- |']
    for signal in assessments:
        for column in ['PDCSAP_FLUX','SAP_FLUX']:
            for width in [1.,2.]:
                values=[next(test.get('snr_diagonal') for test in signal['source_group_fixed_ephemeris']
                    if test['column']==column and test['width_days']==width and test['edge_buffer_days']==cut)
                    for cut in [.25,.5]]
                lines.append(f"| {signal['peak']['period']:.6f} | {column} | {width:g} | "+' | '.join('unavailable' if v is None else f'{v:.2f}' for v in values)+' |')
    lines += ['', '## 8.87447 d signal', '',
        'The source sector has three nominal supported windows, but detrending changes the PDC ranking statistic from 10.08 (2 d) to 4.49 (1 d). Under 1 d detrending, local PDC depths are approximately 450, 115, and 308 ppm, while SAP gives approximately −38, +186, and −516 ppm. The third window is only 0.36 d before a quality gap. The original formal odd/even discrepancy is 4.13 sigma, with three events and correlated-error limitations. These tests establish instability rather than a robust shallow transit sequence.', '',
        'At the exact source ephemeris, other seasons give no coherent depth: the 2 d PDC fixed scores are −2.12 in sector 104 and −0.53 in sector 105 with the 0.25 d edge cut. That is a retrospective exact-ephemeris check. The several-hundred-day separation amplifies an unfit period uncertainty, so it does not eliminate every neighboring alias. Source-event ±0.3 d raw baseline scans are saved rather than interpreting fixed predictions as a complete posterior.', '',
        '## 14.85263 d signal', '',
        'The first nominal event at BTJD 4193.86913 occurs only 0.27 d after a large quality gap; a 0.5 d exclusion removes its usable transit window. The remaining two sector-105 event windows are away from the edges. Their local 1 d PDC depths are approximately 357 and 349 ppm, compared with SAP depths of +149 and −579 ppm. Sector 105 alone gives PDC fixed scores of 5.08 (1 d) or 8.99 (2 d), while SAP gives −2.81 or −3.65. The joint score thus relies on broad, reduction-dependent variability, with one source event adjacent to a gap.', '',
        'The strongest reduction disagreement persists with direct fits to raw 2 min flux after known-transit masking: the event at BTJD 4223.574 has SAP depths from −458 to −1095 ppm across all four local linear/quadratic baselines (−3.1 to −8.5 diagonal sigma), while PDC gives +307 to +1300 ppm. A polynomial-baseline change does not restore a shared downward event. At the first event, raw PDC significance is only 0.1–1.6 diagonal sigma across these models. Source PDC red-noise binning factors at the proposed duration are approximately 1.7–2.1 even for the 1 d detrending.', '',
        'Both revised searches choose the longest original duration grid point, about 5.5 h. An eccentric planet or uncertain host density could give a long transit; duration alone is not a rejection. The raw polynomial baseline and duration-response measurements test whether a compact event survives alternative local trend choices. They provide no justified physical planet fit for these unstable residuals.', '',
        '## Scientific scope', '',
        'The earlier rejection of single boundary excursions was not reused as a blanket decision. The revised peaks fail their own SAP/PDC and detrending-consistency checks. If a future independent pixel reduction restores coherent, similarly shaped events, revisit them with a multi-season ephemeris posterior, injection recovery, and complete false-positive diagnostics. No global significance, planetary radius, discovery, or confirmation is claimed here. All raw FITS remain unchanged and retained.', '',
        'Reproduce with `python src/vet_toi1064_trimmed.py` after retaining the revised screen JSON and FITS products named by SHA256 in the evidence file.']
    (OUT/'TOI1064_TRIMMED.md').write_text('\n'.join(lines)+'\n')
    print('Saved revised TOI-1064 vetting',flush=True)

if __name__=='__main__':main()
