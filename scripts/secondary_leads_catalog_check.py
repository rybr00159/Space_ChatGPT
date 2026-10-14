"""Reproducible current public stellar lookup for four sector-1 secondary leads.

Preserves response bodies; never overwrites existing snapshots. This script only
writes candidate-prefixed sources and cached-match results, not source catalogs.
"""
import concurrent.futures, csv, hashlib, json, math, sqlite3, sys, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
import requests
ROOT=Path("/workspace/Space_ChatGPT")
sys.path.insert(0,str(ROOT/"src"))
from broad_crossmatch import CatalogIndex
OUT=ROOT/"data/broad_catalogs"
TICS=[79354493,231825824,259901286,308600649]
DB="file:/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite?mode=ro"
S=requests.Session(); S.headers["User-Agent"]="Space_ChatGPT-public-research/0.3"
STAMP=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
def write_snapshot(name,body,meta):
    p=OUT/name
    if p.exists(): p=p.with_name(p.stem+"_"+STAMP+p.suffix)
    if p.exists(): raise RuntimeError(f"Refuse overwrite {p}")
    p.write_bytes(body)
    meta.update(sha256=hashlib.sha256(body).hexdigest(),bytes=len(body))
    p.with_suffix(p.suffix+".manifest.json").write_text(json.dumps(meta,indent=2)+"\n")
    return p
def fetch(tic,label,url,query=None):
    meta=dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),url=url)
    if query: meta["query"]=query
    try:
        r=S.get(url,timeout=80)
        meta.update(http_status=r.status_code,content_type=r.headers.get("Content-Type"))
        p=write_snapshot(f"tic{tic}_{label}",r.content,meta)
        r.raise_for_status()
        print(tic,label,"OK",r.status_code,len(r.content),flush=True)
        return p
    except Exception as e:
        print(tic,label,"ERROR",type(e).__name__,str(e)[:180],flush=True)
        return None
def tap(tic,label,q):
    url="https://gea.esac.esa.int/tap-server/tap/sync?"+urlencode(dict(request="doQuery",lang="ADQL",format="csv",query=q))
    return fetch(tic,label+".csv",url,q)
def handle(tic):
    c=sqlite3.connect(DB,uri=True)
    raw=json.loads(c.execute("select result_json from stars where tic=?",(tic,)).fetchone()[0]); c.close()
    periods=[p["period"] for s in raw["scans"] for p in s.get("peaks",[])]
    ix=CatalogIndex(position_radius_arcsec=120)
    match={"checked_utc":datetime.now(timezone.utc).isoformat(),"tic":tic,"ra":raw["ra"],"dec":raw["dec"],
           "periods_days":periods,"scope":"Cached multi-catalog all-method target, positional and 2-5x integer harmonic/subharmonic matches within 1%; current live sources checked separately.",
           "catalog_metadata":ix.metadata(),"signals":[ix.match_signal(tic,p,raw["ra"],raw["dec"]) for p in periods]}
    write_snapshot(f"tic{tic}_cached_matches_current.json",(json.dumps(match,indent=2)+"\n").encode(),{"derived_utc":datetime.now(timezone.utc).isoformat(),"sources":ix.metadata()})
    p=fetch(tic,"exofop_current.json",f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json")
    if not p: return
    j=json.loads(p.read_text()); aliases=j["basic_info"].get("star_names","")
    fetch(tic,"exofop_nearby_current.csv",f"https://exofop.ipac.caltech.edu/tess/download_nearbytarget.php?id={tic}&output=csv")
    import re
    m=re.search(r"Gaia DR3 (\d+)",aliases)
    if not m: print(tic,"NO GAIA DR3 ALIAS",aliases,flush=True); return
    sid=m.group(1); ra=float(j["coordinates"]["ra"]); dec=float(j["coordinates"]["dec"])
    cols="source_id,ra,dec,phot_g_mean_mag,bp_rp,ruwe,non_single_star,phot_variable_flag,astrometric_params_solved,parallax,parallax_error,astrometric_excess_noise,astrometric_excess_noise_sig,duplicated_source,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,phot_bp_rp_excess_factor,phot_g_mean_flux_over_error,phot_g_n_obs,radial_velocity,radial_velocity_error,rv_nb_transits,rv_amplitude_robust,rv_chisq_pvalue,rv_renormalised_gof,rv_template_teff,rv_expected_sig_to_noise"
    tasks=[
      ("gaia_target_current",f"select {cols} from gaiadr3.gaia_source where source_id={sid}"),
      ("gaia_stellar_current",f"select source_id,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,radius_gspphot,radius_gspphot_lower,radius_gspphot_upper,radius_flame,mass_flame from gaiadr3.astrophysical_parameters where source_id={sid}"),
      ("gaia_class_current",f"select source_id,best_class_name,best_class_score from gaiadr3.vari_classifier_result where source_id={sid}"),
      ("gaia_eb_current",f"select * from gaiadr3.vari_eclipsing_binary where source_id={sid}"),
      ("gaia_neighbors_current",f"select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.bp_rp,g.ruwe,g.non_single_star,g.phot_variable_flag,v.in_vari_eclipsing_binary from gaiadr3.gaia_source g left outer join gaiadr3.vari_summary v on g.source_id=v.source_id where 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))")
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
      for f in [ex.submit(tap,tic,label,q) for label,q in tasks]: f.result()
if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        for f in [ex.submit(handle,t) for t in TICS]: f.result()

