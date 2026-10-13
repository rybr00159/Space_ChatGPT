"""Falsification checks for a weak additional-period residual in TOI-904.

The period selected from sectors 93/94 is frozen when evaluating sector 88.
No additional period search, downloads, paid requests or physical fits occur.
"""
import json

import numpy as np
import pandas as pd
from astropy import constants as const, units as u
from astropy.io import fits

from search import ROOT, flatten, known_for, mask_known, read_curve
from validate import local_depth, noise_beta
from vet_tic150360085 import aggregate, gap_keep, quadratic_depth
import matplotlib.pyplot as plt

TIC=261257684
OUT=ROOT/"candidates"/"toi904"


def duration_check(period, observed_duration):
    table=pd.read_csv(ROOT/"data"/"catalogs"/"pilot_known_planets.csv")
    row=table[table["hostname"]=="TOI-904"].iloc[0]
    mass,radius=float(row["st_mass"]),float(row["st_rad"])
    p=period*u.day
    a=(const.G*mass*const.M_sun*p**2/(4*np.pi**2))**(1/3)
    a_over_r=(a/(radius*const.R_sun)).decompose().value
    central_hours=(p/np.pi*np.arcsin(1/a_over_r)).to_value(u.hour)
    ratio=observed_duration*24/central_hours
    minimum_eccentricity=(ratio**2-1)/(ratio**2+1) if ratio>1 else 0
    return {"stellar_mass_msun":mass,"stellar_radius_rsun":radius,
            "stellar_density_solar":mass/radius**3,"a_over_r_circular":float(a_over_r),
            "maximum_small_planet_circular_duration_hours":float(central_hours),
            "selected_bls_duration_hours":observed_duration*24,
            "duration_ratio":float(ratio),"illustrative_min_e_at_apastron":float(minimum_eccentricity),
            "source":"NASA confirmed snapshot pilot_known_planets.csv, Harris et al. 2023 stellar parameters",
            "warning":"Consistency calculation, not a transit fit or eccentricity measurement. Stellar uncertainties, finite planet radius, impact parameter, and BLS duration bias omitted."}


def write_report(result):
    ephem=result["ephemeris"]
    physical=result["duration_consistency"]
    lines=["# TOI-904 / TIC 261257684: additional-period residual fails screening", "",
           "A weak 13.46-day residual appears after known-planet masking in sectors 93/94. Its depth depends strongly on detrending width and local baseline curvature. The frozen-ephemeris comparison in sector 88 provides no significant additional support. These observations do not justify a credible new exoplanet candidate or physical transit fit.", "",
           f"Selection: P={ephem['period_days']:.8f} d, T0={ephem['epoch_btjd']:.8f} BTJD, duration={ephem['duration_days']*24:.3f} h, formal BLS S/N={ephem['selected_bls_snr']:.2f}. The selected duration sits at the longest searched box. The alternate 1-day detrending prefers P={ephem['alternate_width1_best_period_days']:.5f} d, demonstrating a broad/inconsistent period preference. Sectors 93/94 supplied the selected ephemeris; sector 88 is observationally separate.", "",
           "## Independent local baseline measurements", "",
           "The same selected period/epoch/duration are frozen for all comparisons. Each event has a separate local linear or quadratic baseline using points 1–3 durations from center. The candidate box uses points within half a duration. Known TOI-904 b/c transits and 0.25-day neighborhoods around raw-data boundaries/gaps are excluded. Formal uncertainties are diagonal-noise estimates and are not false-alarm probabilities.", "",
           "|Sectors|Reduction|Width (d)|Baseline|Depth (ppm)|Formal error (ppm)|Formal S/N|Events|",
           "|---|---|---|---|---|---|---|---|"]
    for row in result["measurements"]:
        for baseline in ("linear","quadratic"):
            m=row[baseline]
            lines.append(f"|{row['sectors']}|{row['column']}|{row['width_days']}|{baseline}|{m['depth']*1e6:.0f}|{m['formal_error']*1e6:.0f}|{m['formal_snr']:.2f}|{m['n_events']}|")
    for row in result["raw_sap"]:
        for baseline in ("linear","quadratic"):
            m=row[baseline]
            lines.append(f"|{row['sectors']}|Raw SAP|—|{baseline}|{m['depth']*1e6:.0f}|{m['formal_error']*1e6:.0f}|{m['formal_snr']:.2f}|{m['n_events']}|")
    lines += ["", "## Duration consistency", "",
              f"The saved confirmed-planet catalog gives Mstar={physical['stellar_mass_msun']:.3f} Msun, Rstar={physical['stellar_radius_rsun']:.3f} Rsun, density={physical['stellar_density_solar']:.2f} solar. At the selected period, the maximum central, circular, small-planet transit duration is approximately {physical['maximum_small_planet_circular_duration_hours']:.2f} hours. The BLS box duration is {physical['duration_ratio']:.2f} times this value. A central transit near apastron would need eccentricity approximately {physical['illustrative_min_e_at_apastron']:.2f} under these assumptions; other possibilities include incorrect stellar density or a biased box duration. This is a consistency warning, not a measured eccentricity or proof against a planet. A larger impact parameter shortens the circular duration.", "",
              "## Assessment and limitations", "",
              "- The local linear depth falls from roughly 285 ppm at 2-day detrending to 33 ppm at 1-day detrending; a quadratic baseline reduces the 2-day value further. This sensitivity is consistent with slow stellar/instrumental curvature entering the long box statistic. No repeated isolated transit-shaped events are established.",
              "- The 13.46-day trial is distinct from catalog periods 10.8772 and 83.9997 days, but that alone establishes neither novelty nor astrophysical origin. Current known-transit masks are used, with ephemeris uncertainties not fully propagated.",
              "- Sector 88 tests the frozen nominal ephemeris without a new search. Its weak positive measurements do not independently validate the signal. The selected broad peak has no rigorous timing/period posterior, so failure at the nominal phase would not alone exclude a true planet over the extrapolation interval.",
              "- SAP and PDCSAP are reductions of the same photons. Their comparison tests processing sensitivity, not independent astrophysical confirmation. Raw SAP baseline fits use empirical adjacent-point noise floors, while correlated stellar/instrumental variability remains incompletely modeled.",
              "- BLS explored 110,061 period values in the selected two-sector screen, plus epochs/durations and multiple stars, sectors and detrenders. An attractive formal maximum is insufficient. Local off-phase controls and binned noise ratios in the JSON are descriptive and cannot calibrate this global trial factor.",
              "- The source is not localized using pixels; background binaries and other blends remain untested. No exhaustive additional-candidate literature novelty claim is made. A physical fit, planet radius, impact parameter or false-positive probability would overinterpret the weak evidence.",
              "- This screening rejection is not evidence that TOI-904 lacks further planets. Longer baseline data and variability models with injection-tested short-transit preservation could revisit the system without presuming the selected period.", "",
              "Reproduce with `python src/vet_toi904.py`; inputs and product checksums are preserved in `data/observations.json`. The script does not download or rerun the period search.", ""]
    (OUT/"ASSESSMENT.md").write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    selected=json.loads((ROOT/"results"/"robust"/"tic261257684_s93-94_w2_edge0.25.json").read_text())
    alternate=json.loads((ROOT/"results"/"robust"/"tic261257684_s93-94_w1_edge0.25.json").read_text())
    best=selected["peaks"][0]
    p,t0,d=best["period"],best["t0_btjd"],best["duration_days"]
    records=[r for r in json.loads((ROOT/"data"/"observations.json").read_text()) if r["tic"]==TIC and r["sector"] in [88,93,94]]
    result={"tic":TIC,"status":"detrending_sensitive_residual_not_credible_candidate",
            "ephemeris":{"period_days":p,"epoch_btjd":t0,"duration_days":d,"selection_sectors":[93,94],
                         "held_out_sector":88,"selected_bls_snr":best["snr"],"period_grid_count":best["period_grid_count"],
                         "alternate_width1_best_period_days":alternate["peaks"][0]["period"]},
            "known_planets":known_for(TIC),"observations":records,"measurements":[],"raw_sap":[],
            "duration_consistency":duration_check(p,d)}
    fig,axes=plt.subplots(4,1,figsize=(10,11),layout="constrained")
    rng=np.random.default_rng(261257684)
    primary_plot=[]
    for column in ("PDCSAP_FLUX","SAP_FLUX"):
        for width in (1.,2.):
            chunks={}
            for record in records:
                t,f,e,_=read_curve(record["local_path"],column)
                tb,fb,eb=flatten(t,f,e,width)
                keep=gap_keep(t,tb)&mask_known(tb,known_for(TIC))
                chunks[record["sector"]]=(tb[keep],fb[keep],eb[keep])
            for sectors in ([93,94],[88]):
                t,f,e=(np.concatenate([chunks[sector][i]for sector in sectors])for i in range(3))
                linear=local_depth(t,f,e,p,t0,d)
                quadratic=quadratic_depth(t,f,e,p,t0,d)
                null=[]
                for shift in rng.uniform(.15*p,.85*p,48):
                    m=local_depth(t,f,e,p,t0+shift,d)
                    if m["n_events"]:null.append(m["formal_snr"])
                result["measurements"].append({"sectors":sectors,"column":column,"width_days":width,
                    "linear":linear,"quadratic":quadratic,"noise":noise_beta(t,f,e,d*1440),
                    "off_phase_controls":{"n":len(null),"formal_snr_std":float(np.std(null,ddof=1)),
                                          "max_formal_snr":float(max(null)),"warning":"Correlated descriptive phase controls; not search FAP"}})
                if column=="PDCSAP_FLUX" and sectors==[93,94]:
                    primary_plot.append((width,linear))
                    for event in linear["events"]:
                        if event["cycle"] not in range(4):continue
                        x=t-event["center_btjd"]
                        use=abs(x)<3*d
                        axes[event["cycle"]].plot(x[use],f[use],".",ms=2,alpha=.5,label=f"{width:g}-day median")
            print(column,width,[(row['sectors'],round(row['linear']['formal_snr'],2),round(row['quadratic']['formal_snr'],2))for row in result['measurements']if row['column']==column and row['width_days']==width],flush=True)
    for i,ax in enumerate(axes):
        ax.axvspan(-d/2,d/2,color="red",alpha=.1)
        ax.set(xlabel="Days from frozen event center",ylabel="Normalized residual flux",title=f"TOI-904 predicted cycle {i}; two correlated detrenders")
        ax.legend()
    fig.savefig(OUT/"event_windows.png",dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4),layout="constrained")
    for width,measurement in primary_plot:
        events=measurement["events"]
        ax.errorbar([event["cycle"]+.07*(width-1)for event in events],[event["depth"]*1e6 for event in events],
                    yerr=[event["formal_error"]*1e6 for event in events],fmt="o",label=f"{width:g}-day median")
    ax.axhline(0,color="black",lw=.7);ax.legend()
    ax.set(xlabel="Predicted cycle in sectors 93/94",ylabel="Local-linear box depth (ppm)",title="TOI-904: strong detrending sensitivity at frozen period")
    fig.savefig(OUT/"event_depths.png",dpi=140);plt.close(fig)
    raw_chunks={}
    for record in records:
        t,f,e,_=read_curve(record["local_path"],"SAP_FLUX")
        median=np.median(f);f=f/median;e=e/median
        adjacency=np.diff(t)<.01
        delta=np.diff(f)[adjacency]
        noise=1.4826*np.median(np.abs(delta-np.median(delta)))/np.sqrt(2)
        e*=max(1,noise/np.median(e))
        keep=gap_keep(t,t)&mask_known(t,known_for(TIC))
        raw_chunks[record["sector"]]=(t[keep],f[keep],e[keep])
    for sectors in ([93,94],[88]):
        t,f,e=(np.concatenate([raw_chunks[sector][i]for sector in sectors])for i in range(3))
        result["raw_sap"].append({"sectors":sectors,"linear":local_depth(t,f,e,p,t0,d),
                                  "quadratic":quadratic_depth(t,f,e,p,t0,d)})
    result["interpretation"]={"physical_fit_performed":False,"reason":"Signal lacks local transit significance and is strongly sensitive to detrending/local curvature.","novelty_claim":False}
    (OUT/"measurements.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    write_report(result)


if __name__=="__main__":
    main()
