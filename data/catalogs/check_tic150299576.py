"""Focused prior recognition check for an unpromoted weak field signal."""
from check_field_external import fetch,ROOT
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
import csv,xml.etree.ElementTree as ET,json

if __name__=='__main__':
 aliases=[r['id'] for r in csv.DictReader((ROOT/'tic150299576_simbad_identifiers.csv').open())]
 query=' OR '.join('all:"'+a+'"' for a in aliases)
 requests=[
 ('tic150299576_alias_arxiv.xml','https://export.arxiv.org/api/query?'+urlencode({'search_query':query,'max_results':100})),
 ('tic150299576_exofop.json','https://exofop.ipac.caltech.edu/tess/target.php?id=150299576&json'),
 ('tic150299576_nearby_tic.csv','https://exofop.ipac.caltech.edu/tess/download_nearbytarget.php?id=150299576&output=csv'),
 ]
 def get(v):
  name,u=v
  try:
   b=fetch(name,u);print(name,len(b),'bytes cached with provenance')
   return name,b
  except Exception as e:print(name,type(e).__name__,str(e));return name,None
 with ThreadPoolExecutor(max_workers=2) as p:result=dict(p.map(get,requests))
 b=result['tic150299576_alias_arxiv.xml']
 if b:
  ns={'a':'http://www.w3.org/2005/Atom','op':'http://a9.com/-/spec/opensearch/1.1/'};root=ET.fromstring(b)
  out={'target':'TIC150299576','aliases':aliases,'query':query,'total_results':root.findtext('op:totalResults','?',ns),'returned':len(root.findall('a:entry',ns)),'scope':'arXiv metadata alias search only, not full text/journal-complete; no novelty inference from null search.'}
  (ROOT/'tic150299576_alias_arxiv_summary.json').write_text(json.dumps(out,indent=2)+'\n');print('Alias metadata matches:',out['total_results'])
 # Inspect authoritative SIMBAD reference table schemas before constructing the join.
 for table in ['has_ref','ref']:
  q="select column_name from TAP_SCHEMA.columns where table_name='"+table+"'"
  u='https://simbad.cds.unistra.fr/simbad/sim-tap/sync?'+urlencode({'request':'doQuery','lang':'adql','format':'csv','query':q})
  b=fetch('simbad_'+table+'_schema.csv',u);print(table,'schema columns:',len(list(csv.DictReader(b.decode().splitlines()))))

 q="select r.bibcode,r.doi,r.title,r.abstract from basic b join has_ref h on b.oid=h.oidref join ref r on h.oidbibref=r.oidbib where b.main_id='UCAC4 143-007406'"
 u='https://simbad.cds.unistra.fr/simbad/sim-tap/sync?'+urlencode({'request':'doQuery','lang':'adql','format':'csv','query':q})
 b=fetch('tic150299576_simbad_references.csv',u)
 print('SIMBAD reference rows',len(list(csv.DictReader(b.decode().splitlines()))))
