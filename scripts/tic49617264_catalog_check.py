"""Public Gaia, ExoFOP, and cached catalog context for the public SPOC TIC49617264 TCE."""
import json, concurrent.futures
from scripts.secondary_leads_catalog_check import fetch,tap,write_snapshot,CatalogIndex,datetime,timezone
tic=49617264
ra=37.8931844608584; dec=-40.0058723004587; period=4.8319416046
ix=CatalogIndex(position_radius_arcsec=120)
check={"checked_utc":datetime.now(timezone.utc).isoformat(),"tic":tic,"ra":ra,"dec":dec,
       "scope":"Already public SPOC TCE; cached catalog absence is not novelty.","period_days":period,
       "signal":ix.match_signal(tic,period,ra,dec),"catalog_metadata":ix.metadata()}
write_snapshot(f"tic{tic}_cached_matches_current.json",(json.dumps(check,indent=2)+"\n").encode(),
               {"derived_utc":datetime.now(timezone.utc).isoformat(),"sources":ix.metadata()})
p=fetch(tic,"exofop_current.json",f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json")
j=json.loads(p.read_text())
import re
sid=re.search(r"Gaia DR3 (\d+)",j["basic_info"]["star_names"]).group(1)
ra=float(j["coordinates"]["ra"]);dec=float(j["coordinates"]["dec"])
cols="source_id,ra,dec,phot_g_mean_mag,bp_rp,ruwe,non_single_star,phot_variable_flag,astrometric_params_solved,parallax,parallax_error,astrometric_excess_noise,astrometric_excess_noise_sig,duplicated_source,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,phot_bp_rp_excess_factor,phot_g_mean_flux_over_error,phot_g_n_obs,radial_velocity,radial_velocity_error,rv_nb_transits,rv_amplitude_robust,rv_chisq_pvalue,rv_renormalised_gof,rv_template_teff,rv_expected_sig_to_noise"
jobs=[
("gaia_target_current",f"select {cols} from gaiadr3.gaia_source where source_id={sid}"),
("gaia_stellar_current",f"select source_id,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,radius_gspphot,radius_gspphot_lower,radius_gspphot_upper,radius_flame,mass_flame from gaiadr3.astrophysical_parameters where source_id={sid}"),
("gaia_class_current",f"select source_id,best_class_name,best_class_score from gaiadr3.vari_classifier_result where source_id={sid}"),
("gaia_eb_current",f"select * from gaiadr3.vari_eclipsing_binary where source_id={sid}"),
("gaia_neighbors_current",f"select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.bp_rp,g.ruwe,g.non_single_star,g.phot_variable_flag,v.in_vari_eclipsing_binary from gaiadr3.gaia_source g left outer join gaiadr3.vari_summary v on g.source_id=v.source_id where 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))")]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(tap,tic,label,q) for label,q in jobs]
  futures.append(pool.submit(fetch,tic,"exofop_nearby_current.csv",f"https://exofop.ipac.caltech.edu/tess/download_nearbytarget.php?id={tic}&output=csv"))
  for f in futures:f.result()

