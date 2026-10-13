"""Independent fixed-ephemeris vetting of the weak TIC 150360085 residual.

All sectors 93--96 selected the ephemeris; these are not held-out data. This
script performs no period search, uses no paid API, and preserves input FITS.
"""
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from search import ROOT, flatten, phase, read_curve
from validate import local_depth, noise_beta
import matplotlib.pyplot as plt

TIC = 150360085
OUT = ROOT / "candidates" / "tic150360085"


def aggregate(events):
    if not events:
        return {"n_events": 0, "events": []}
    depths = np.array([event["depth"] for event in events])
    errors = np.array([event["formal_error"] for event in events])
    weights = 1 / errors**2
    mean = float(np.average(depths, weights=weights))
    formal = float(np.sqrt(1 / weights.sum()))
    chi2 = float(np.sum(weights * (depths - mean)**2) / max(1, len(events)-1))
    error = formal * np.sqrt(max(1, chi2))
    return {"depth": mean, "formal_error": formal, "formal_snr": mean/formal,
            "event_reduced_chi2": chi2, "event_scatter_error": error,
            "event_scatter_snr": mean/error, "n_events": len(events), "events": events}


def quadratic_depth(t, f, e, period, epoch, duration):
    events = []
    for cycle in np.unique(np.round((t-epoch)/period).astype(int)):
        center = epoch + cycle*period
        x = (t-center)/duration
        inside = np.abs(x)<.5
        outside = (np.abs(x)>1)&(np.abs(x)<3)
        use = inside|outside
        if inside.sum()<2 or outside.sum()<6:
            continue
        design = np.column_stack([np.ones(use.sum()),x[use],x[use]**2,-inside[use].astype(float)])
        weight = 1/e[use]**2
        normal = design.T@(weight[:,None]*design)
        if np.linalg.cond(normal)>1e10:
            continue
        covariance = np.linalg.inv(normal)
        coef = covariance@(design.T@(weight*f[use]))
        residual = f[use]-design@coef
        events.append({"cycle":int(cycle),"center_btjd":float(center),"depth":float(coef[-1]),
                       "formal_error":float(np.sqrt(covariance[-1,-1])),
                       "n_in":int(inside.sum()),"n_baseline":int(outside.sum()),
                       "reduced_chi2":float(np.sum(weight*residual**2)/(use.sum()-4))})
    return aggregate(events)


def gap_keep(raw_t, binned_t, cut=.25):
    boundaries = np.r_[0,np.flatnonzero(np.diff(raw_t)>.2)+1,len(raw_t)]
    keep = np.zeros(len(binned_t),bool)
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        keep |= (binned_t>=raw_t[lo]+cut)&(binned_t<=raw_t[hi-1]-cut)
    return keep


def event_differences(measurement):
    events = measurement["events"]
    odd = aggregate([event for event in events if event["cycle"]%2])
    even = aggregate([event for event in events if not event["cycle"]%2])
    if not odd["n_events"] or not even["n_events"]:
        return None
    error = np.hypot(odd["formal_error"],even["formal_error"])
    return {"odd_depth_ppm":odd["depth"]*1e6,"even_depth_ppm":even["depth"]*1e6,
            "formal_difference_sigma":float(abs(odd["depth"]-even["depth"])/error)}


def write_report(result):
    lines = ["# TIC 150360085: weak residual rejected as a credible candidate", "",
             "A joint four-sector BLS screen selected a roughly 13.123-day signal. Independent local baseline measurements weaken it, and early predicted events lack significant individual depth. Their uncertainties are large enough that a shallow planet cannot be ruled out. This remains an exploratory residual, with no evidence sufficient to designate a credible new exoplanet candidate.", "",
             f"Frozen selection ephemeris: P={result['ephemeris']['period_days']:.8f} d, T0={result['ephemeris']['epoch_btjd']:.8f} BTJD, duration={result['ephemeris']['duration_days']*24:.3f} h. The 5.5-hour duration lies near the longest searched duration. Both detrenders select this period, but share observations and are correlated analyses. All sectors 93–96 contributed to selection; none are independent confirmation.", "",
             "## Independent local-baseline measurements", "",
             "Each predicted event uses a separate linear baseline from points 1–3 durations away; a box depth measures the central half-duration interval. The table uses the same frozen ephemeris for every reduction. Quadratic baselines are a diagnostic for local curvature. Formal errors assume independent points; event scatter is an additional warning, not a calibrated red-noise model.", "",
             "|Reduction|Baseline|Depth (ppm)|Formal error (ppm)|Formal S/N|Event-scatter S/N|Events|",
             "|---|---|---|---|---|---|---|"]
    for row in result["joint_measurements"]:
        for baseline in ("linear","quadratic"):
            m=row[baseline]
            lines.append(f"|{row['column']} width {row['width_days']} d|{baseline}|{m['depth']*1e6:.0f}|{m['formal_error']*1e6:.0f}|{m['formal_snr']:.2f}|{m['event_scatter_snr']:.2f}|{m['n_events']}|")
    for baseline,m in result["raw_sap_joint"].items():
        lines.append(f"|Raw SAP, no median detrending|{baseline}|{m['depth']*1e6:.0f}|{m['formal_error']*1e6:.0f}|{m['formal_snr']:.2f}|{m['event_scatter_snr']:.2f}|{m['n_events']}|")
    lines += ["", "## Per-sector support", "",
              "|Sector|Reduction|Depth (ppm)|Formal error (ppm)|Formal S/N|Events|",
              "|---|---|---|---|---|---|"]
    for row in result["sector_measurements"]:
        if row["width_days"]!=1:
            continue
        m=row["linear"]
        lines.append(f"|{row['sector']}|{row['column']} width 1 d|{m['depth']*1e6:.0f}|{m['formal_error']*1e6:.0f}|{m['formal_snr']:.2f}|{m['n_events']}|")
    lines += ["", "## Predicted event consistency", "",
              "|Cycle|Center (BTJD)|Depth (ppm)|Formal error (ppm)|Formal S/N|",
              "|---|---|---|---|---|"]
    primary=next(row for row in result["joint_measurements"] if row["column"]=="PDCSAP_FLUX" and row["width_days"]==1)
    for event in primary["linear"]["events"]:
        lines.append(f"|{event['cycle']}|{event['center_btjd']:.5f}|{event['depth']*1e6:.0f}|{event['formal_error']*1e6:.0f}|{event['depth']/event['formal_error']:.2f}|")
    lines += ["", "## Interpretation and limitations", "",
              "- Counted events mean sufficient sampling, not significant repeated transits. Local depths include zero/negative events. A constant-depth interpretation has excess event scatter; sector-by-sector support varies strongly.",
              "- Alternative detrending widths and local baseline curvature change the result. SAP comparison shares the same photons and cannot independently verify astrophysical origin. Raw SAP local-linear and quadratic measurements are included in the JSON; their diagonal errors can understate stellar/instrumental covariance.",
              "- Noise ratios on duration-scale bins and off-phase fixed-ephemeris control statistics are descriptive. They are neither a period-search false-alarm probability nor a proven error inflation rule. The ephemeris was selected among hundreds of thousands of period trials and multiple stars/sectors/reductions.",
              "- Odd/even or secondary-eclipse nondetections have little diagnostic power at this signal strength. Favorable aperture crowding fractions alone cannot rule out background blends, nor establish a transit on the target.",
              "- A physical transit/trapezoid fit is not justified: consistent transit-shaped events have not survived screening. No radius, impact parameter, or planetary interpretation is claimed. Pixel-level source localization and exhaustive novelty checks would be required if stronger independent evidence emerges.",
              "- No archived sectors 87–90 were used in this evaluation. A fixed-ephemeris comparison with independently retrieved older data could falsify or revive the hypothesis, provided ephemeris error from the selected broad peak is propagated. Lack of a nominal predicted event alone would be inconclusive over a long extrapolation.",
              "- This negative screening judgment does not prove the star lacks planets; it means the selected residual lacks sufficient evidence to survive this pilot's vetting.", "",
              "All inputs, observation identifiers, processing versions, and retrieval checksums are linked through `data/observations.json`. `src/vet_tic150360085.py` regenerates this assessment and diagnostic plots without a new period search.", ""]
    (OUT/"ASSESSMENT.md").write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    selection=json.loads((ROOT/"results"/"robust"/"tic150360085_s93-94-95-96_w1_edge0.25.json").read_text())
    peak=selection["peaks"][0]
    period,epoch,duration=peak["period"],peak["t0_btjd"],peak["duration_days"]
    records=[r for r in json.loads((ROOT/"data"/"observations.json").read_text()) if r["tic"]==TIC and r["sector"] in [93,94,95,96]]
    result={"tic":TIC,"status":"weak_selected_residual_not_credible_candidate",
            "ephemeris":{"period_days":period,"epoch_btjd":epoch,"duration_days":duration,
                         "selection_sectors":[93,94,95,96],"held_out_sectors":[],"selection_bls_snr":peak["snr"],
                         "period_grid_count":peak["period_grid_count"]},
            "observations":records,"joint_measurements":[],"sector_measurements":[],"raw_sap_measurements":[]}
    fig,axes=plt.subplots(4,2,figsize=(12,12),layout="constrained")
    rng=np.random.default_rng(150360085)
    for column in ("PDCSAP_FLUX","SAP_FLUX"):
        for width in (1.,2.):
            chunks=[]
            for record in records:
                t,f,e,_=read_curve(record["local_path"],column)
                tb,fb,eb=flatten(t,f,e,width)
                keep=gap_keep(t,tb)
                tt,ff,ee=tb[keep],fb[keep],eb[keep]
                chunks.append((tt,ff,ee))
                row={"sector":record["sector"],"column":column,"width_days":width,
                     "linear":local_depth(tt,ff,ee,period,epoch,duration),
                     "quadratic":quadratic_depth(tt,ff,ee,period,epoch,duration),
                     "noise":noise_beta(tt,ff,ee,duration*1440)}
                result["sector_measurements"].append(row)
                if width==1:
                    ax=axes[record["sector"]-93,0 if column=="PDCSAP_FLUX" else 1]
                    ph=phase(tt,period,epoch)
                    ax.plot(ph,ff,".",ms=1,alpha=.35)
                    bins=np.arange(-duration*3,duration*3+.02,.04)
                    centers,means,errors=[],[],[]
                    for a,b in zip(bins[:-1],bins[1:]):
                        use=(ph>=a)&(ph<b)
                        if use.sum()<2:continue
                        w=1/ee[use]**2
                        centers.append((a+b)/2);means.append(np.average(ff[use],weights=w));errors.append(np.sqrt(1/w.sum()))
                    ax.errorbar(centers,means,yerr=errors,fmt="o",ms=3,color="black")
                    ax.axvspan(-duration/2,duration/2,color="red",alpha=.1)
                    ax.set(xlim=(-duration*3,duration*3),xlabel="Time from frozen epoch modulo P (d)",ylabel="Flattened normalized flux",title=f"Sector {record['sector']} {column}")
            tt,ff,ee=(np.concatenate([chunk[i] for chunk in chunks]) for i in range(3))
            linear=local_depth(tt,ff,ee,period,epoch,duration)
            quadratic=quadratic_depth(tt,ff,ee,period,epoch,duration)
            null=[]
            for shift in rng.uniform(.15*period,.85*period,96):
                check=local_depth(tt,ff,ee,period,epoch+shift,duration)
                if check["n_events"]:null.append(check["formal_snr"])
            result["joint_measurements"].append({"column":column,"width_days":width,"linear":linear,"quadratic":quadratic,
                "noise":[noise_beta(tt,ff,ee,minutes) for minutes in (180,duration*1440,480)],
                "odd_even":event_differences(linear),"secondary":local_depth(tt,ff,ee,period,epoch+period/2,duration),
                "off_phase_controls":{"seed":150360085,"n":len(null),"formal_snr_std":float(np.std(null,ddof=1)),
                                      "maximum_formal_snr":float(max(null)),"n_at_least_observed":int(np.sum(np.array(null)>=linear["formal_snr"])),
                                      "warning":"Descriptive correlated phase controls; not search FAP"}})
            print(column,width,linear["depth"]*1e6,linear["formal_error"]*1e6,linear["formal_snr"],linear["event_scatter_snr"],flush=True)
    fig.savefig(OUT/"per_sector_frozen_fold.png",dpi=140);plt.close(fig)
    # Raw SAP local baselines, without a global median detrending or binning.
    raw_events={"linear":[],"quadratic":[]}
    for record in records:
        t,f,e,_=read_curve(record["local_path"],"SAP_FLUX")
        norm=np.median(f);f=f/norm;e=e/norm
        adjacent=np.diff(t)<.01
        scatter=1.4826*np.median(np.abs(np.diff(f)[adjacent]-np.median(np.diff(f)[adjacent])))/np.sqrt(2)
        e*=max(1,scatter/np.median(e))
        keep=gap_keep(t,t)
        row={"sector":record["sector"],"empirical_point_noise_ppm":float(scatter*1e6),
             "linear":local_depth(t[keep],f[keep],e[keep],period,epoch,duration),
             "quadratic":quadratic_depth(t[keep],f[keep],e[keep],period,epoch,duration)}
        result["raw_sap_measurements"].append(row)
        for baseline in raw_events:raw_events[baseline].extend(row[baseline]["events"])
        with fits.open(record["local_path"]) as hdul:
            row["crowdsap"]=hdul[1].header.get("CROWDSAP")
    result["raw_sap_joint"]={baseline:aggregate(events) for baseline,events in raw_events.items()}
    result["interpretation"]={"physical_fit_performed":False,"reason":"Inconsistent low-significance local events fail initial screening; no credible transit model is justified.",
                             "novelty_claim":False,"independent_observation_confirmation":False}
    (OUT/"measurements.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    # Event depths expose missing/variable events without folding them together.
    primary=next(row for row in result["joint_measurements"] if row["column"]=="PDCSAP_FLUX" and row["width_days"]==1)
    fig,ax=plt.subplots(figsize=(9,4),layout="constrained")
    for row in result["joint_measurements"]:
        if row["column"]!="PDCSAP_FLUX":continue
        events=row["linear"]["events"]
        ax.errorbar([event["cycle"]+.06*(row["width_days"]-1) for event in events],
                    [event["depth"]*1e6 for event in events],yerr=[event["formal_error"]*1e6 for event in events],
                    fmt="o",label=f"Median width {row['width_days']:g} d")
    ax.axhline(0,color="black",lw=.7)
    ax.set(xlabel="Predicted cycle",ylabel="Local-baseline box depth (ppm)",title="TIC 150360085: fixed selected ephemeris; formal errors")
    ax.legend();fig.savefig(OUT/"event_depths.png",dpi=140);plt.close(fig)
    write_report(result)


if __name__=="__main__":
    main()
