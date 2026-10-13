"""Public metadata and eclipsing-binary checks for field pilot; no full novelty claim."""
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib,re,shutil,xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from check_field_catalogs import sep_arcsec
ROOT=Path(__file__).resolve().parent

def fetch(name,url):
 p=ROOT/name
 manifest=p.with_suffix(p.suffix+'.manifest.json')
 if p.exists() and manifest.exists():
  meta=json.loads(manifest.read_text());body=p.read_bytes()
  if hashlib.sha256(body).hexdigest()!=meta['sha256']:
   raise ValueError(f'checksum mismatch in {p}; preserve and investigate')
  if meta['url']==url:return body
 if p.exists():
  # Selection/query changes must not silently reuse a different scientific request.
  stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
  archival=p.with_name(p.name+'.'+stamp+'.previous')
  shutil.copy2(p,archival)
  if manifest.exists():shutil.copy2(manifest,archival.with_suffix(archival.suffix+'.manifest.json'))
 with urlopen(Request(url,headers={'User-Agent':'Space_ChatGPT-public-research/0.1'}),timeout=60) as r:b=r.read();status=r.status
 p.write_bytes(b)
 p.with_suffix(p.suffix+'.manifest.json').write_text(json.dumps({'retrieved_utc':datetime.now(timezone.utc).isoformat(),'url':url,'http_status':status,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()},indent=2)+'\n')
 return b

if __name__=='__main__':
 STARS=json.loads((ROOT.parent/'field_sample.json').read_text())['selected']
 conds=[f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r['ra']},{r['dec']},0.008333333333))" for r in STARS]
 q='select main_id,ra,dec,otype from basic where '+' OR '.join(conds)
 url='https://simbad.cds.unistra.fr/simbad/sim-tap/sync?'+urlencode({'request':'doQuery','lang':'adql','format':'csv','query':q})
 b=fetch('field_simbad.csv',url)
 source=list(csv.DictReader(b.decode().splitlines()));matched=[]
 for s in STARS:
  for r in source:
   d=sep_arcsec(s['ra'],s['dec'],float(r['ra']),float(r['dec']))
   if d<=30:matched.append({'tic':s['tic'],'main_id':r['main_id'],'otype':r['otype'],'separation_arcsec':d,'ra':r['ra'],'dec':r['dec']})
 with (ROOT/'field_simbad_matches.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['tic','main_id','otype','separation_arcsec','ra','dec']);w.writeheader();w.writerows(matched)
 print('SIMBAD locally filtered rows:',len(matched))
 def eb(star):
  tid=str(star['tic']);url='https://tessebs.villanova.edu/search_results?'+urlencode({'tic':tid})
  b=fetch('field_tesseb_'+tid+'.html',url);s=b.decode(errors='replace')
  # A detail-link matching the TIC, rather than its mere appearance in the search box, identifies a result.
  found=bool(re.search(r'href=[\"\'][^\"\']*0*'+tid+r'[\"\']',s))
  count=re.search(r'(\d+)\s*(?:result|object)',re.sub('<[^>]+>',' ',s),re.I)
  inclusion=None;signals=[]
  if found:
   detail_url='https://tessebs.villanova.edu/'+tid.zfill(10)
   detail=fetch('field_tesseb_detail_'+tid+'.html',detail_url).decode(errors='replace')
   if 'In catalog?' in detail:
    section=detail.split('In catalog?',1)[1].split('</table>',1)[0]
    for row in re.findall(r'<tr[^>]*>(.*?)</tr>',section,re.S):
     cells=[' '.join(re.sub('<[^>]+>',' ',cell).split()) for cell in re.findall(r'<td[^>]*>(.*?)</td>',row,re.S)]
     if len(cells)>=7:signals.append({'signal':cells[0],'in_catalog':cells[1]=='True','epoch':cells[2],'period':cells[4],'morphology':cells[6]})
    if signals:inclusion=any(row['in_catalog'] for row in signals)
  return {'tic':tid,'result_link_found':found,'eb_catalog_inclusion':inclusion,'signal_rows':signals,'count_text':count.group(0) if count else '', 'url':url,'note':'Search returns light-curve-only entries as well as EBs; only explicit In catalog=True counts as a catalog EB. Exact TIC query; absence does not exclude unlisted EBs or contamination.'}
 with ThreadPoolExecutor(max_workers=2) as pool:result=list(pool.map(eb,STARS))
 (ROOT/'field_tesseb_checks.json').write_text(json.dumps(result,indent=2)+'\n');print('TEB explicit inclusion checks:',[(r['tic'],r['eb_catalog_inclusion']) for r in result])
 query=' OR '.join('all:"'+str(s['tic'])+'"' for s in STARS)
 url='https://export.arxiv.org/api/query?'+urlencode({'search_query':query,'max_results':100})
 b=fetch('field_tic_arxiv.xml',url)
 n={'a':'http://www.w3.org/2005/Atom','op':'http://a9.com/-/spec/opensearch/1.1/'};root=ET.fromstring(b)
 total=root.findtext('op:totalResults','?',n);entries=root.findall('a:entry',n)
 (ROOT/'field_arxiv_check.json').write_text(json.dumps({'query':query,'total_results':total,'returned':len(entries),'scope':'arXiv searchable metadata for numeric TIC identifiers; not a full-text or all-alias search.'},indent=2)+'\n')
 print('ARXIV TIC numeric search',total,len(entries))

 ids=','.join(str(int(star['gaia'])) for star in STARS)
 q='select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.ruwe,g.non_single_star,g.phot_variable_flag,v.in_vari_eclipsing_binary,v.in_vari_rotation_modulation,v.in_vari_short_timescale from gaiadr3.gaia_source g left outer join gaiadr3.vari_summary v on g.source_id=v.source_id where g.source_id in ('+ids+')'
 url='https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode({'request':'doQuery','lang':'ADQL','format':'csv','query':q})
 try:
  b=fetch('field_gaia_dr3_flags.csv',url)
  print('GAIA DR3 flag rows:',len(list(csv.DictReader(b.decode().splitlines()))))
 except Exception as exc:
  print('GAIA check not complete',str(exc))
 names=','.join("'"+r['main_id'].replace("'","''")+"'" for r in matched)
 q='select b.main_id,t.otypes,b.nbref,f.G,f.V,f.J,f.H,f.K from basic b left outer join alltypes t on b.oid=t.oidref left outer join allfluxes f on b.oid=f.oidref where b.main_id in ('+names+')'
 url='https://simbad.cds.unistra.fr/simbad/sim-tap/sync?'+urlencode({'request':'doQuery','lang':'adql','format':'csv','query':q})
 b=fetch('field_simbad_types_fluxes.csv',url)
 print('SIMBAD all object types and magnitude rows',len(list(csv.DictReader(b.decode().splitlines()))))

 q='select source_id,best_class_name,best_class_score from gaiadr3.vari_classifier_result where source_id in ('+ids+')'
 url='https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode({'request':'doQuery','lang':'ADQL','format':'csv','query':q})
 b=fetch('field_gaia_variable_classifier.csv',url)
 print('GAIA variable classifications:',len(list(csv.DictReader(b.decode().splitlines()))))
