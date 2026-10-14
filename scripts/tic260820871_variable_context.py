"""Finite public SIMBAD/alias/arXiv context for TIC260820871."""
import csv,io,json,concurrent.futures
from urllib.parse import urlencode
from astropy.coordinates import SkyCoord
import astropy.units as u
from scripts.secondary_leads_catalog_check import fetch
tic=260820871;ra=348.054595469913;dec=-70.1131990097929
def sim(label,q):
  url="https://simbad.cds.unistra.fr/simbad/sim-tap/sync?"+urlencode({"request":"doQuery","lang":"ADQL","format":"csv","query":q})
  return fetch(tic,label+".csv",url,q)
q=f"select oid,main_id,ra,dec,otype,nbref from basic where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))"
p=sim("simbad_current",q);rows=list(csv.DictReader(p.open()))
ct=SkyCoord(ra*u.deg,dec*u.deg);associated=[]
for r in rows:
  sep=float(ct.separation(SkyCoord(float(r["ra"])*u.deg,float(r["dec"])*u.deg)).arcsec)
  print("SIMBAD",sep,r,flush=True)
  if sep<=2:associated.append(r)
tokens=["TIC 260820871","260820871","23121312-7006475","6381543153782126464","9341-1858-1","100-104680","34531911"]
if associated:
  r=associated[0];oid=r["oid"];tokens.append(r["main_id"])
  with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:
    fs=[pool.submit(sim,"simbad_current_aliases",f"select id from ident where oidref={oid}"),
        pool.submit(sim,"simbad_current_references",f"select r.bibcode,r.doi,r.title,r.abstract from has_ref h join ref r on h.oidbibref=r.oidbib where h.oidref={oid}")]
    for f in fs:f.result()
query=" OR ".join('all:"'+x+'"'for x in tokens)
fetch(tic,"arxiv_current.xml","https://export.arxiv.org/api/query?"+urlencode({"search_query":query,"max_results":100}))

