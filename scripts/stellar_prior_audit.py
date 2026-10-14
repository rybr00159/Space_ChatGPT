"""Finite recognition audit of six retained eclipsing-source leads.

Preserves source bytes and query/hash manifests. Missing catalog membership does
not establish novelty; matching a source does not attribute its TESS eclipse.
Run with --download for public queries, without it for retained-state derivation.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import requests

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/broad_catalogs'
OUT = ROOT / 'data/stellar_systems/prior_audit'
REPORT = ROOT / 'reports/stellar_systems'
TARGETS = {
    300814564: ('5273775087036017024',116.92191756103277,-68.62140624817651,1.860044776),
    273879358: ('4692911684016155008',32.30851048028379,-71.32001533152403,7.586894227),
    270312735: ('4617279547612110208',32.62890012752531,-83.14812279702222,8.370105786),
    370247818: ('4625697305555206912',55.45564700638868,-79.11323002464255,3.522194945921),
    238202207: ('4629516356114834944',58.37129206255529,-73.78568182519616,4.643651161748),
    269994642: ('6806240691223695744',311.4615173986208,-24.650049533874235,3.492744247184),
}
EXTRA_IDS = ['5273775091332140800','4617279616331585408','4625697168116255360']
SUPPLEMENT_IDS = ['4692910966756908288','5496090260794096128','5496090260793489024','5496090256498964096']

def utc():
    return datetime.now(timezone.utc).isoformat()


def fetch(name,url,params=None):
    """Reuse byte-validated response snapshot, including preserved failures."""
    OUT.mkdir(parents=True,exist_ok=True)
    p=OUT/name
    m=p.with_suffix(p.suffix+'.manifest.json')
    if p.exists() and m.exists():
        meta=json.loads(m.read_text())
        if hashlib.sha256(p.read_bytes()).hexdigest()!=meta['sha256']:
            raise RuntimeError('Hash mismatch: '+str(p))
        print('CACHED',name,meta.get('http_status'),flush=True)
        return p
    if p.exists() or m.exists():
        raise RuntimeError('Incomplete existing snapshot '+str(p))
    try:
        r=requests.get(url,params=params,timeout=(15,85),headers={'User-Agent':'Space_ChatGPT-public-stellar-research/1.0'})
        body=r.content
        meta=dict(retrieved_utc=utc(),requested_url=url,query_parameters=params,url=r.url,http_status=r.status_code,content_type=r.headers.get('Content-Type'),bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
        p.write_bytes(body);m.write_text(json.dumps(meta,indent=2)+'\n')
        print('FETCH',name,r.status_code,len(body),flush=True)
        return p
    except Exception as e:
        fail=OUT/(name+'.error.json')
        if not fail.exists():fail.write_text(json.dumps(dict(attempted_utc=utc(),url=url,query_parameters=params,error=type(e).__name__+': '+str(e)),indent=2)+'\n')
        print('ERROR',name,type(e).__name__,str(e)[:100],flush=True)
        return None


def tap(name,query,service='gaia'):
    bases={'gaia':'https://gea.esac.esa.int/tap-server/tap/sync','simbad':'https://simbad.cds.unistra.fr/simbad/sim-tap/sync','vizier':'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'}
    return fetch(name+'.csv',bases[service],dict(request='doQuery',lang='ADQL',format='csv',query=query))


def mast_tic(name,ids):
    p=OUT/(name+'.json');m=p.with_suffix('.json.manifest.json')
    if p.exists() and m.exists():
        if hashlib.sha256(p.read_bytes()).hexdigest()!=json.loads(m.read_text())['sha256']:
            raise RuntimeError('Hash mismatch: '+str(p))
        return p
    request={'service':'Mast.Catalogs.Filtered.Tic','params':{'columns':'*','filters':[{'paramName':'GAIA','values':ids}]},'format':'json','pagesize':100,'page':1}
    url='https://mast.stsci.edu/api/v0/invoke'
    response=requests.post(url,data={'request':json.dumps(request)},timeout=85)
    p.write_bytes(response.content)
    m.write_text(json.dumps(dict(retrieved_utc=utc(),url=url,method='POST',request=request,http_status=response.status_code,sha256=hashlib.sha256(response.content).hexdigest(),bytes=len(response.content)),indent=2)+'\n')
    response.raise_for_status()
    if response.json().get('status')!='COMPLETE':raise ValueError('Incomplete TIC query')
    return p


def download():
    sidlist=','.join(v[0] for v in TARGETS.values())+','+','.join(EXTRA_IDS)
    jobs=[]
    for table in ['nss_two_body_orbit','nss_acceleration_astro','nss_non_linear_spectro','nss_vim_fl','vari_eclipsing_binary','vari_classifier_result','vari_summary']:
        jobs.append((tap,('gaia_'+table,f'select * from gaiadr3.{table} where source_id in ({sidlist})')))
    jobs += [
        (tap,('vizier_stellar_binary_catalog_metadata',"select table_name,description from TAP_SCHEMA.tables where (description like '%APOGEE%' and (description like '%binary%' or description like '%Binary%' or description like '%Joker%')) or (description like '%GALAH%' and (description like '%binary%' or description like '%Binary%')) or description like '%ATLAS catalog of candidate variable%'",'vizier')),
        (tap,('gaia_nss_metadata',"select table_name,description from TAP_SCHEMA.tables where table_name like 'gaiadr3.nss%'")),
    ]
    ticlist=','.join(str(t) for t in TARGETS)
    for table in ['J/ApJS/259/50/table1','J/AJ/170/356/table1','J/MNRAS/511/4285/table2']:
        jobs.append((tap,(table.replace('/','_')+'_exact_tic',f'select * from "{table}" where TIC in ({ticlist})','vizier')))
    for label,table in [('kostov2025_published_table3_exact_tic','J/ApJS/279/50/table3'),('kostov2025_known_table4_exact_tic','J/ApJS/279/50/table4')]:
        jobs.append((tap,(label,f'select * from "{table}" where TIC in ({ticlist},273879360,260352037)','vizier')))
    jobs.append((tap,('vizier_apogee_galah_binary_metadata_v2',"select table_name,description from TAP_SCHEMA.tables where (description like '%APOGEE%' or description like '%GALAH%' or description like '%Joker%') and (description like '%binar%' or description like '%Binar%' or description like '%companion%' or description like '%spectroscop%')",'vizier')))
    for tic,(sid,ra,dec,period) in TARGETS.items():
        jobs.append((fetch,(f'tic{tic}_asassn_live.html','https://asas-sn.osu.edu/variables',{'ra':ra,'dec':dec,'radius':2,'show_non_periodic':'true','show_without_class':'true','asassn_discov_only':'false'})))
        for label,source in [('asassn','II/366'),('atlas2018','J/AJ/156/241/table4'),
                             ('shi2022','J/ApJS/259/50/table1'),
                             ('tang2025','J/AJ/170/356/table1'),
                             ('cdips2022','J/MNRAS/511/4285/table2'),
                             ('tian2020_wide','J/ApJS/246/4/catalog'),
                             ('10k2025','J/ApJS/279/50'),
                             ('wise2021','J/MNRAS/503/3975/catalog'),
                             ('rowan2022','J/MNRAS/517/2190/table'),
                             ('sb9','B/sb9/main'),('wds','B/wds/wds')]:
            jobs.append((fetch,(f'tic{tic}_{label}.tsv','https://vizier.cds.unistra.fr/viz-bin/asu-tsv',{'-source':source,'-c':f'{ra} {dec}','-c.rs':120,'-out.max':'unlimited','-out.all':'','-out.form':'tsv'})))
    with ThreadPoolExecutor(max_workers=5) as pool:
        for f in [pool.submit(fn,*args) for fn,args in jobs]:f.result()


def download_supplement():
    sidlist=','.join(SUPPLEMENT_IDS)
    jobs=[]
    jobs.append((mast_tic,('supplement_tic_sources',SUPPLEMENT_IDS)))
    for table in ['nss_two_body_orbit','nss_acceleration_astro','nss_non_linear_spectro','nss_vim_fl','vari_eclipsing_binary','vari_classifier_result','vari_summary']:
        jobs.append((tap,('supplement_gaia_'+table,f'select * from gaiadr3.{table} where source_id in ({sidlist})')))
    jobs.append((tap,('supplement_gaia_sources',f'select * from gaiadr3.gaia_source where source_id in ({sidlist})')))
    for tic in [273879360,260352037,737466840,260352041]:
        jobs.append((fetch,(f'tic{tic}_exofop.json',f'https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json')))
        jobs.append((fetch,(f'tic{tic}_tesseb.html',f'https://tessebs.villanova.edu/{tic:010d}')))
    for table in ['J/ApJS/259/50/table1','J/AJ/170/356/table1','J/MNRAS/511/4285/table2','J/ApJS/279/50/table3','J/ApJS/279/50/table4']:
        jobs.append((tap,('supplement_'+table.replace('/','_')+'_exact_tic',f'select * from "{table}" where TIC in (273879360,260352037,737466840,260352041)','vizier')))
    for name,ra,dec in [('background273',32.25246282511974,-71.32317034350537),('tic260352037',94.78791015296669,-57.07821822202367)]:
        jobs.append((fetch,(name+'_asassn_live.html','https://asas-sn.osu.edu/variables',{'ra':ra,'dec':dec,'radius':2,'show_non_periodic':'true','show_without_class':'true','asassn_discov_only':'false'})))
        jobs.append((tap,(name+'_simbad_cone',f"select oid,main_id,ra,dec,otype,nbref from basic where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))",'simbad')))
        for label,source in [('vsx','B/vsx/vsx'),('asassn','II/366'),('shi2022','J/ApJS/259/50/table1'),('tang2025','J/AJ/170/356/table1'),('cdips2022','J/MNRAS/511/4285/table2'),('sb9','B/sb9/main'),('wds','B/wds/wds'),('tian2020_wide','J/ApJS/246/4/catalog')]:
            jobs.append((fetch,(name+'_'+label+'.tsv','https://vizier.cds.unistra.fr/viz-bin/asu-tsv',{'-source':source,'-c':f'{ra} {dec}','-c.rs':120,'-out.max':'unlimited','-out.all':'','-out.form':'tsv'})))
        tname='tic273879360' if name=='background273' else name
        for label,source in [('10k2025','J/ApJS/279/50'),('wise2021','J/MNRAS/503/3975/catalog'),('rowan2022','J/MNRAS/517/2190/table')]:
            jobs.append((fetch,(tname+'_'+label+'.tsv','https://vizier.cds.unistra.fr/viz-bin/asu-tsv',{'-source':source,'-c':f'{ra} {dec}','-c.rs':120,'-out.max':'unlimited','-out.all':'','-out.form':'tsv'})))
    with ThreadPoolExecutor(max_workers=5) as pool:
        for f in [pool.submit(fn,*args) for fn,args in jobs]:f.result()


def verify_sources():
    records=[]
    for base in [BASE,OUT]:
        paths=sorted(base.glob('*.manifest.json'))
        for m in paths:
            if base==BASE and not any(m.name.startswith('tic'+str(t)+'_') for t in TARGETS):continue
            p=Path(str(m)[:-len('.manifest.json')]);meta=json.loads(m.read_text())
            if not p.exists() or 'sha256' not in meta:continue
            digest=hashlib.sha256(p.read_bytes()).hexdigest()
            if digest!=meta['sha256']:raise RuntimeError('Hash mismatch '+str(p))
            records.append(dict(file=str(p.relative_to(ROOT)),manifest=str(m.relative_to(ROOT)),sha256=digest,bytes=p.stat().st_size,retrieved_utc=meta.get('retrieved_utc',meta.get('derived_utc'))))
    REPORT.mkdir(parents=True,exist_ok=True)
    (REPORT/'prior_audit_sources.json').write_text(json.dumps(dict(verified_utc=utc(),source_count=len(records),sources=records),indent=2)+'\n')
    print('VERIFIED',len(records),'source snapshots',flush=True)


def parse_tsv(path, with_table=False):
    """Read VizieR TSV data without treating failed responses as empty tables."""
    text=path.read_text()
    if '#INFO\tError=' in text or 'QUERY_STATUS" value="ERROR' in text:
        raise ValueError('Query failed: '+str(path))
    if '#Name:' not in text:
        raise ValueError('Missing VizieR table metadata: '+str(path))
    blocks=[];block=[];table=None
    for line in text.splitlines():
        if line.startswith('#Table'):
            if block:blocks.append((table,block))
            block=[];table=None
        elif line.startswith('#Name:'):
            table=line.split(':',1)[1].strip()
        elif not line.startswith('#') and line.strip():
            block.append(line)
    if block:blocks.append((table,block))
    rows=[]
    for table,block in blocks:
        # Empty cones carry column metadata but no data header, units or rows.
        if len(block)<3:continue
        if not all(set(v.strip()) <= {'-'} for v in block[2].split('\t')):
            raise ValueError('Unexpected TSV layout: '+str(path))
        fields=block[0].split('\t')
        for line in block[3:]:
            row=dict(zip(fields,[v.strip() for v in line.split('\t')]))
            if with_table:row['vizier_table']=table
            rows.append(row)
    return rows


def csv_rows(path):
    if not path.exists():raise FileNotFoundError(path)
    meta=json.loads(path.with_suffix(path.suffix+'.manifest.json').read_text())
    if meta.get('http_status',200)!=200:raise ValueError('Unsuccessful response: '+str(path))
    return list(csv.DictReader(path.open()))


def retained_gaia(tic,sid):
    for ending in ['gaia_target_current.csv','gaia.csv','gaia_nearby.csv']:
        p=BASE/f'tic{tic}_{ending}'
        if p.exists():
            rows=[r for r in csv.DictReader(p.open()) if r['source_id']==sid]
            if rows:return rows[0]
    raise ValueError('No retained Gaia target '+str(tic))


def derive():
    REPORT.mkdir(parents=True,exist_ok=True)
    nss={table:csv_rows(OUT/('gaia_'+table+'.csv')) for table in ['nss_two_body_orbit','nss_acceleration_astro','nss_non_linear_spectro','nss_vim_fl','vari_eclipsing_binary','vari_classifier_result','vari_summary']}
    rows=[]
    scopes=['asassn','atlas2018','shi2022','tang2025','cdips2022','tian2020_wide','sb9','wds','wise2021','rowan2022']
    priority={370247818:1,273879358:2,300814564:3,238202207:4,269994642:5,270312735:6}
    for tic,(sid,ra,dec,period) in TARGETS.items():
        g=retained_gaia(tic,sid)
        row=dict(tic=tic,gaia_dr3_source_id=sid,ra_dr3=ra,dec_dr3=dec,photometric_or_binary_period_days=period,followup_priority=priority[tic],ruwe=g.get('ruwe',''),nss_flag=g.get('non_single_star',''))
        for scope in scopes:row[scope+'_cone_rows']=len(parse_tsv(OUT/f'tic{tic}_{scope}.tsv'))
        row['atlas2018_within_declination_coverage']=dec>=-30
        for scope,data in nss.items():row[scope+'_exact_source_rows']=sum(r.get('source_id')==sid for r in data)
        html=(OUT/f'tic{tic}_asassn_live.html').read_text()
        m=re.search(r'Total Light Curves Found:\s*<strong>(.*?)</strong>',html,re.S)
        if not m:raise ValueError('ASAS-SN result count not found for '+str(tic))
        row['asassn_live_120arcsec_rows']=int(m.group(1).replace(',',''))
        tenk=OUT/f'tic{tic}_10k2025.tsv'
        if tenk.exists():
            data=parse_tsv(tenk,with_table=True)
            row['kostov2025_published_eb_cone_rows']=sum(r['vizier_table'] in ['J/ApJS/279/50/table3','J/ApJS/279/50/table4'] for r in data)
            row['kostov2025_unvetted_inventory_cone_rows']=sum(r['vizier_table']=='J/ApJS/279/50/table2' for r in data)
        row['rv_amplitude_robust_kms']=g.get('rv_amplitude_robust','')
        row['rv_chisq_pvalue']=g.get('rv_chisq_pvalue','')
        if tic==370247818:
            rv=next(r for r in csv.DictReader((BASE/f'tic{tic}_gaia_rv.csv').open()) if r['source_id']==sid)
            row['rv_amplitude_robust_kms']=rv.get('rv_amplitude_robust','');row['rv_chisq_pvalue']=rv.get('rv_chisq_pvalue','')
        row['recognition_assessment']='No prior EB identification found in checked finite scope; novelty unestablished'
        rows.append(row)
    rows.sort(key=lambda r:r['followup_priority'])
    p=REPORT/'prior_audit_results.csv'
    with p.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    method=dict(generated_utc=utc(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),target_count=len(rows),cone_radius_arcsec=120,exact_tic_queries=['J_ApJS_259_50_table1_exact_tic.csv','J_AJ_170_356_table1_exact_tic.csv','J_MNRAS_511_4285_table2_exact_tic.csv'],limitations='Negative finite membership checks cannot establish novelty; catalog positions are not eclipse source localization. ATLAS outside −30degree south coverage is not an informative nonmatch. Gaia published variability statistics are prior evidence even when NSS/EB catalogs have no solution.')
    (REPORT/'prior_audit_results.manifest.json').write_text(json.dumps(method,indent=2)+'\n')
    print('DERIVED',len(rows),'recognition audit rows',flush=True)


def derive_supplement():
    import xml.etree.ElementTree as ET
    p=OUT/'supplement_tic_sources.json'
    if not p.exists():return
    targets={r['GAIA']:r for r in json.loads(p.read_text())['data']}
    gaia={r['source_id']:r for r in csv_rows(OUT/'supplement_gaia_sources.csv')}
    ntab={table:csv_rows(OUT/('supplement_gaia_'+table+'.csv'))for table in ['nss_two_body_orbit','nss_acceleration_astro','nss_non_linear_spectro','nss_vim_fl','vari_eclipsing_binary','vari_classifier_result','vari_summary']}
    rows=[]
    for sid in SUPPLEMENT_IDS:
        g=gaia[sid];t=targets[sid];tic=int(t['ID'])
        r=dict(tic=tic,gaia_dr3_source_id=sid,ra_dr3=g['ra'],dec_dr3=g['dec'],gmag=g['phot_g_mean_mag'],tmag=t.get('Tmag'),parallax_mas=g['parallax'],parallax_error_mas=g['parallax_error'],pmra_masyr=g['pmra'],pmdec_masyr=g['pmdec'],ruwe=g['ruwe'])
        for label,data in ntab.items():r[label+'_exact_source_rows']=sum(x.get('source_id')==sid for x in data)
        prefix='background273'if tic==273879360 else'tic260352037'
        for label in ['vsx','asassn','shi2022','tang2025','cdips2022','sb9','wds','tian2020_wide']:
            r[label+'_field120arcsec_rows']=len(parse_tsv(OUT/(prefix+'_'+label+'.tsv')))
        tprefix='tic273879360'if tic==273879360 else'tic260352037'
        for label in ['wise2021','rowan2022']:
            r[label+'_field120arcsec_rows']=len(parse_tsv(OUT/(tprefix+'_'+label+'.tsv')))
        data=parse_tsv(OUT/(tprefix+'_10k2025.tsv'),with_table=True)
        r['kostov2025_published_eb_field120arcsec_rows']=sum(x['vizier_table']in['J/ApJS/279/50/table3','J/ApJS/279/50/table4']for x in data)
        r['kostov2025_unvetted_inventory_field120arcsec_rows']=sum(x['vizier_table']=='J/ApJS/279/50/table2'for x in data)
        html=(OUT/(prefix+'_asassn_live.html')).read_text()
        m=re.search(r'Total Light Curves Found:\s*<strong>(.*?)</strong>',html,re.S)
        if not m:raise ValueError('Missing ASAS-SN count')
        r['asassn_live_field120arcsec_rows']=int(m.group(1).replace(',',''))
        j=json.loads((OUT/f'tic{tic}_exofop.json').read_text())
        r['aliases']=j['basic_info']['star_names']
        for label in ['tois','ctois','stellar_companions','files']:r['exofop_'+label+'_rows']=len(j.get(label,[]))
        ap=OUT/f'tic{tic}_arxiv.xml'
        r['arxiv_indexed_alias_results']='not_checked'
        if ap.exists():
            root=ET.fromstring(ap.read_bytes());value=root.findtext('{http://a9.com/-/spec/opensearch/1.1/}totalResults')
            if value is None:raise ValueError('No arXiv result count '+str(ap))
            r['arxiv_indexed_alias_results']=int(value)
        rows.append(r)
    p=REPORT/'prior_audit_supplement_results.csv'
    with p.open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('DERIVED',len(rows),'supplementary recognition rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--download',action='store_true');parser.add_argument('--supplement',action='store_true');args=parser.parse_args()
    if args.download:download()
    if args.supplement:download_supplement()
    verify_sources()
    derive()
    derive_supplement()
