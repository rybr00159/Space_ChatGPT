"""Conditional S66/S93 source comparisons, independent of held-out sectors.

Uses existing native-difference maps and their exact raw cutout WCS. Propagates
Gaia DR3 proper motion from J2016 to each sector's median observing epoch.
An empirical Gaussian is not a calibrated PRF; chi-square differences are not
confidence levels, and overlapping reduction variants are not independent.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/workspace/tess-research-data/matplotlib")
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
from scipy.special import erf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pixel_vet_49617264 import gaussfit, safe

OUT = ROOT / "candidates/tic214311908/source_localization"
DATA = ROOT / "data/tic214311908/source_localization"
TARGET = "5925373272850784768"
NAMES = {
    TARGET: "target TIC214311908",
    "5925373272850782208": "comoving TIC214311910",
    "5925373268528684288": "Mira TIC214311906",
    "5925373272827260544": "3.15arcsec TIC1515709327",
    "5925373272827258752": "4.27arcsec TIC1515709326",
    "5925373272827278208": "known EB P0.316083d",
    "5925373204107739648": "known EB P0.260926d",
    "5925373272827253632": "11.29arcsec TIC1515709322",
}
IMPORTANT = list(NAMES)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, obj):
    Path(path).write_text(json.dumps(safe(obj), indent=2) + "\n")


def projected(row, epoch, wcs):
    c = SkyCoord(float(row["ra"])*u.deg, float(row["dec"])*u.deg)
    years = epoch.jyear - 2016.
    # Spherical angular propagation; perspective acceleration, parallax and
    # missing motion are not fitted. These corrections are tiny here compared
    # with the subpixel model/WCS uncertainty. Missing PM is explicitly flagged.
    pm_available = bool(row["pmra"] and row["pmdec"])
    c = c.spherical_offsets_by(float(row["pmra"] or 0)*years*u.mas,
                               float(row["pmdec"] or 0)*years*u.mas)
    xy = wcs.all_world2pix([[c.ra.deg, c.dec.deg]], 0)[0]
    return dict(source_id=row["source_id"], label=NAMES.get(row["source_id"], row["source_id"]),
                ra_epoch=c.ra.deg, dec_epoch=c.dec.deg, x=xy[0], y=xy[1],
                original_gaia_ra=float(row["ra"]), original_gaia_dec=float(row["dec"]),
                pmra_masyr=float(row["pmra"]) if row["pmra"] else None,
                pmdec_masyr=float(row["pmdec"]) if row["pmdec"] else None,
                proper_motion_applied=pm_available, reference_epoch_jyear=2016.,
                observing_jyear=epoch.jyear, gmag=float(row["phot_g_mean_mag"]) if row["phot_g_mean_mag"] else None,
                separation_arcsec_at_catalog_epoch=float(row["separation_arcsec"]))


def fixed_shape(image, error, pos, widths):
    """Linear positive source amplitude plus plane, with prescribed PSF widths."""
    yy,xx = np.indices(image.shape)
    gx=.5*(erf((xx+.5-pos[0])/(np.sqrt(2)*widths[0]))-erf((xx-.5-pos[0])/(np.sqrt(2)*widths[0])))
    gy=.5*(erf((yy+.5-pos[1])/(np.sqrt(2)*widths[1]))-erf((yy-.5-pos[1])/(np.sqrt(2)*widths[1])))
    matrix=np.column_stack([(gx*gy).ravel(),np.ones(image.size),(xx-xx.mean()).ravel(),(yy-yy.mean()).ravel()])
    scale=np.maximum(error.ravel(),.03)
    weighted=matrix/scale[:,None]; response=image.ravel()/scale
    coef=np.linalg.lstsq(weighted,response,rcond=None)[0]
    if coef[0]<0:
        coef[0]=0.;coef[1:]=np.linalg.lstsq(weighted[:,1:],response,rcond=None)[0]
    residual=(matrix@coef-image.ravel())/scale
    return dict(widths_pixels=list(widths),source_flux_e_s=coef[0],chi2_diagonal=np.sum(residual**2),dof=image.size-4,
                note="Fixed Gaussian shape plus background plane, not calibrated TESS PRF.")


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source_table=ROOT/"candidates/tic214311908/catalogs/gaia_neighbors_derived.csv"
    rows=list(csv.DictReader(source_table.open()))
    selection=[r for r in rows if float(r["separation_arcsec"])<=12. or r["source_id"] in IMPORTANT]
    summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),script_sha256=sha(__file__),
                 gaussian_helper_sha256=sha(ROOT/"src/pixel_vet_49617264.py"),
                 gaia_table=dict(path=str(source_table.relative_to(ROOT)),sha256=sha(source_table)),
                 selection="All Gaia sources within12arcsec plus named Mira and two known neighboring EBs, chosen from catalogues independently of difference-image preferred centers.",
                 no_holdout_data_accessed=True,sectors=[],limitations=[
                     "Two sectors and eight overlapping reductions reuse the same photons; variants are not independent evidence.",
                     "Gaussian model is not a calibrated TESS PRF; unresolved blending and unknown WCS residuals remain.",
                     "Diagonal errors/event-scatter floor do not model all pixel/event covariance. Do not translate delta-chi-square to sigma or FAP.",
                     "Proper motion is propagated on the sphere; annual parallax, perspective acceleration and missing-PM motion are not modeled.",
                     "Catalog associations and reference-image peaks do not assign the transit host."])
    flat=[]
    for sector in [66,93]:
        base=ROOT/f"candidates/tic214311908/sector{sector}_transit_difference"
        previous=json.loads((base/"localization.json").read_text())
        if previous["gaussian_helper_sha256"]!=summary["gaussian_helper_sha256"]:
            raise ValueError("Gaussian helper changed since input difference-map fit")
        raw=Path(previous["source"]["path"])
        if sha(raw)!=previous["source"]["sha256"]:raise ValueError("Raw pixel checksum mismatch")
        arrays=np.load(base/"maps.npz")
        with fits.open(raw,memmap=False) as hd:
            wcs=WCS(hd[2].header)
            epoch=Time(float(np.nanmedian(hd[1].data["TIME"]))+2457000,format="jd",scale="tdb")
        positions=[projected(r,epoch,wcs) for r in selection]
        position_map={p["source_id"]:p for p in positions}
        old_xy=arrays["target_xy"]
        new_xy=np.array([position_map[TARGET][k] for k in ["x","y"]])
        record=dict(sector=sector,epoch_btjd=epoch.jd-2457000,epoch_tdb=epoch.isot,
                    source=previous["source"],maps_sha256=sha(base/"maps.npz"),localization_sha256=sha(base/"localization.json"),
                    old_raw_ra_obj_target_xy=old_xy,proper_motion_corrected_target_xy=new_xy,
                    old_to_corrected_target_offset_pixels=np.linalg.norm(new_xy-old_xy),positions=positions,tests=[])
        for test in previous["tests"]:
            key=("qlp_quality_" if test["common_qlp_quality"] else "")+("sky_" if test["sky_subtracted"] else "")+f"degree{test['baseline_degree']}"
            image,error=arrays[key+"_difference"],arrays[key+"_error"]
            free=test["free_gaussian"]
            free_xy=np.array([free["x"],free["y"]])
            results=[]
            for pos in positions:
                xy=np.array([pos["x"],pos["y"]])
                fixed,_=gaussfit(image,error,xy,fixed_center=True)
                fixed["dof_fixed_position"]=image.size-6
                alt=[fixed_shape(image,error,xy,widths) for widths in [(.65,.65),(.85,.65)]]
                row=dict(source_id=pos["source_id"],label=pos["label"],source_xy=xy,free_center_distance_pixels=np.linalg.norm(free_xy-xy),
                         optimized_width_gaussian=fixed,prescribed_width_gaussians=alt,
                         delta_chi2_vs_free=fixed["chi2_diagonal"]-free["chi2_diagonal"])
                results.append(row)
            target_fit=next(r for r in results if r["source_id"]==TARGET)
            for r in results:
                r["delta_chi2_vs_target"]=r["optimized_width_gaussian"]["chi2_diagonal"]-target_fit["optimized_width_gaussian"]["chi2_diagonal"]
                flat.append(dict(sector=sector,reduction=key,source_id=r["source_id"],label=r["label"],
                    x=r["source_xy"][0],y=r["source_xy"][1],free_center_distance_pixels=r["free_center_distance_pixels"],
                    chi2=r["optimized_width_gaussian"]["chi2_diagonal"],delta_chi2_vs_target=r["delta_chi2_vs_target"],
                    flux_e_s=r["optimized_width_gaussian"]["source_flux_e_s"],sigma_x=r["optimized_width_gaussian"]["sigma_x"],sigma_y=r["optimized_width_gaussian"]["sigma_y"],width_bound_hit=r["optimized_width_gaussian"]["width_bound_hit"]))
            record["tests"].append(dict(reduction=key,n_events=test["n_events"],image_peak_snr=test["max_difference_snr"],
                previous_threshold_flag=test["localization_valid"],free_gaussian=free,
                corrected_target_offset_pixels=np.linalg.norm(free_xy-new_xy),fixed_source_tests=results))
        write(OUT/f"sector{sector}_source_comparison.json",record)
        summary["sectors"].append(dict(sector=sector,file=f"sector{sector}_source_comparison.json",sha256=sha(OUT/f"sector{sector}_source_comparison.json"),
            old_to_corrected_target_offset_pixels=record["old_to_corrected_target_offset_pixels"],
            target_offset_range_pixels=[min(t["corrected_target_offset_pixels"] for t in record["tests"]),max(t["corrected_target_offset_pixels"] for t in record["tests"])]))
        fig,axes=plt.subplots(1,3,figsize=(16,5),layout="constrained")
        for ax,arr,title in zip(axes,[np.log10(np.maximum(arrays["reference"],1)),arrays["degree1_difference"],arrays["sky_degree2_difference"]],
                               ["Native median: log10 flux","Linear local baseline: difference","Quadratic baseline + sky: difference"]):
            im=ax.imshow(arr,origin="lower",cmap="viridis");fig.colorbar(im,ax=ax,shrink=.7)
            for k,sid in enumerate(IMPORTANT):
                pp=position_map[sid];color={TARGET:"red","5925373272850782208":"orange","5925373268528684288":"magenta"}.get(sid,"cyan")
                ax.plot(pp["x"],pp["y"],"+",color=color,ms=9,mew=1.5)
                ax.annotate(str(k),(pp["x"],pp["y"]),xytext=(4,7+(k%2)*7),textcoords="offset points",color=color,fontsize=8)
            if ax is not axes[0]:
                for tt in record["tests"][:4]:ax.plot(tt["free_gaussian"]["x"],tt["free_gaussian"]["y"],"o",color="white",ms=4,mfc="none")
            ax.set_title(title,fontsize=9);ax.set_xlabel("Zero-based cutout x");ax.set_ylabel("Zero-based cutout y")
        fig.suptitle(f"TIC214311908 S{sector}; Gaia positions propagated to {epoch.isot[:10]}\n0 target;1 comoving;2 Mira;3/4 close faint sources;5/6 known EBs;7 possible companion-neighbor",fontsize=10)
        fig.savefig(OUT/f"sector{sector}_gaia_overlays.png",dpi=160);plt.close(fig)
    with (OUT/"fixed_source_comparisons.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    write(OUT/"summary.json",summary)
    ground=DATA/"decaps2_grz.fits"
    if ground.exists():
        with fits.open(ground,memmap=False) as hd:
            arr=hd[0].data;header=hd[0].header;wcs=WCS(header,naxis=2)
        plane=arr[-1] if arr.ndim==3 else arr
        fig,ax=plt.subplots(figsize=(8,8),layout="constrained")
        pos=np.maximum(plane,0);ax.imshow(np.arcsinh(pos/(np.nanpercentile(pos,80)+1e-8)),origin="lower",cmap="gray")
        for k,sid in enumerate(IMPORTANT):
            row=next(r for r in rows if r["source_id"]==sid)
            x,y=wcs.all_world2pix([[float(row["ra"]),float(row["dec"])]],0)[0]
            if 0<=x<plane.shape[1] and 0<=y<plane.shape[0]:
                ax.plot(x,y,"+",color="red" if sid==TARGET else "cyan",ms=8)
                ax.annotate(str(k),(x,y),xytext=(6,5),textcoords="offset points",color="red" if sid==TARGET else "cyan")
        ax.set_title("Public DECaPS2 z-band context; Gaia J2016 source labels\n0 target;1 comoving;2 bright Mira;3/4 close faint stars;7 neighbor of companion",fontsize=10)
        ax.set_xlabel("Zero-based archive cutout x");ax.set_ylabel("Zero-based archive cutout y")
        fig.savefig(OUT/"decaps2_resolved_context.png",dpi=180);plt.close(fig)
        write(OUT/"decaps2_context.json",dict(path=str(ground.relative_to(ROOT)),sha256=sha(ground),header_fits_text=header.tostring(sep="\n",padding=False),
            limitations="Resolved static field context only. Viewer coadd has no unique observation epoch; Gaia labels use J2016 without a false exact epoch claim. No transient, new-source or photometric-host conclusion."))
    print(json.dumps(safe(summary),indent=2))


if __name__=="__main__":
    main()
