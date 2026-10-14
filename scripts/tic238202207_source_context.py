"""Retain source context for TIC238202207 without changing its scientific ranking."""
import csv, io, json, hashlib, concurrent.futures
from pathlib import Path
from datetime import datetime,timezone
from astropy.coordinates import SkyCoord
import astropy.units as u
import requests
from scripts.secondary_leads_catalog_check import write_snapshot,tap,fetch
ROOT=Path("/workspace/Space_ChatGPT")
tic=238202207
j=json.loads((ROOT/f"data/broad_catalogs/tic{tic}_exofop_current.json").read_text())
ra=float(j["coordinates"]["ra"]);dec=float(j["coordinates"]["dec"]);ct=SkyCoord(ra*u.deg,dec*u.deg)
target=list(csv.DictReader((ROOT/f"data/broad_catalogs/tic{tic}_gaia_target_current.csv").open()))[0]
rows=list(csv.DictReader((ROOT/f"data/broad_catalogs/tic{tic}_gaia_neighbors_current.csv").open()))
for r in rows:
    r["separation_arcsec"]=float(ct.separation(SkyCoord(float(r["ra"])*u.deg,float(r["dec"])*u.deg)).arcsec)
    r["g_flux_ratio_to_target"]=10**(-.4*(float(r["phot_g_mean_mag"])-float(target["phot_g_mean_mag"]))) if r["phot_g_mean_mag"] else None
    r["is_gaia_target"]=r["source_id"]==target["source_id"]
rows.sort(key=lambda r:r["separation_arcsec"])
out=io.StringIO();w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
sources=[ROOT/f"data/broad_catalogs/tic{tic}_{x}"for x in ["exofop_current.json","gaia_target_current.csv","gaia_neighbors_current.csv"]]
write_snapshot(f"tic{tic}_gaia_neighbors_derived_current.csv",out.getvalue().encode(),{"derived_utc":datetime.now(timezone.utc).isoformat(),"scope":"120arcsec cone. Full G-band catalog flux ratios are not TESS aperture coupling.","sources":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources}})
for r in rows:
    if not r["is_gaia_target"] and (r["separation_arcsec"]<40 or r["g_flux_ratio_to_target"]>.1):
        print(r,flush=True)
def tic_cone():
    endpoint="https://mast.stsci.edu/api/v0/invoke"
    req={"service":"Mast.Catalogs.Tic.Cone","params":{"ra":ra,"dec":dec,"radius":1/30},"format":"json","pagesize":1000,"page":1}
    res=requests.post(endpoint,data={"request":json.dumps(req)},headers={"User-Agent":"Space_ChatGPT-public-research/0.3"},timeout=90)
    meta={"retrieved_utc":datetime.now(timezone.utc).isoformat(),"url":endpoint,"method":"POST","request":req,"http_status":res.status_code,"content_type":res.headers.get("Content-Type")}
    p=write_snapshot(f"tic{tic}_tic_cone_current.json",res.content,meta);res.raise_for_status()
    data=res.json()
    print("TIC cone status",data.get("status"),"rows",len(data.get("data",[])),flush=True)
    for r in data.get("data",[]):
        if str(r.get("GAIA"))in["4629517902303061632","4629516149956404480"]:
            print("BRIGHT TIC",r["ID"],r.get("GAIA"),r.get("Tmag"),r.get("rad"),r.get("Teff"),flush=True)
    return p
q="select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.bp_rp,g.ruwe,g.non_single_star,g.phot_variable_flag,g.parallax,g.parallax_error,g.ipd_frac_multi_peak,g.astrometric_excess_noise,g.astrometric_excess_noise_sig,g.radial_velocity,g.radial_velocity_error,g.rv_nb_transits,g.rv_amplitude_robust,g.rv_chisq_pvalue,g.rv_renormalised_gof,a.teff_gspphot,a.logg_gspphot,a.radius_gspphot,a.radius_flame,a.mass_flame from gaiadr3.gaia_source g left outer join gaiadr3.astrophysical_parameters a on g.source_id=a.source_id where g.source_id in (4629517902303061632,4629516149956404480)"
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    fs=[pool.submit(tic_cone),pool.submit(tap,tic,"gaia_bright_neighbors_current",q)]
    for f in fs:f.result()

