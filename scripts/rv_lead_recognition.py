#!/usr/bin/env python3
"""Bounded prior-companion recognition audit of two GALAH RV leads.

Public GET requests only; immutable response snapshots include failed responses.
This checks recognition, never claims catalogue absence establishes novelty.
"""
from __future__ import annotations
import argparse, concurrent.futures, csv, hashlib, io, json, math
from datetime import datetime, timezone
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/new_stars/rv_recognition'
CAND=ROOT/'candidates/new_stars/rv_recognition'
TARGETS=[dict(source_id='4690783166942834816',ra=14.857739808017081,dec=-71.49272729328251,
              aliases=['UCAC4 093-001240','2MASS J00592586-7129339']),
         dict(source_id='5437472031746145536',ra=141.85594809943,dec=-35.12455248935,
              aliases=['TYC 7162-1658-1','2MASS J09272542-3507283'])]
CATALOGS={'sb9':'B/sb9/main','vsx':'B/vsx/vsx','wds':'B/wds/wds',
          'traven2020':'J/A+A/638/A145/catalog','asassn_v2':'II/366/catv2021'}
BASES={'gaia':'https://gea.esac.esa.int/tap-server/tap/sync',
       'vizier':'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'}

def fetch(name,url,params=None):
    OUT.mkdir(parents=True,exist_ok=True)
    p=OUT/name;m=p.with_name(p.name+'.manifest.json')
    if p.exists() and m.exists():
        meta=json.loads(m.read_text())
        if hashlib.sha256(p.read_bytes()).hexdigest()!=meta['sha256']:raise RuntimeError(p)
        return p
    if p.exists() or m.exists():raise RuntimeError('Incomplete cache '+str(p))
    try:
        r=requests.get(url,params=params,timeout=(15,90),headers={'User-Agent':'Space_ChatGPT-public-stellar-research/1.1'})
        p.write_bytes(r.content)
        m.write_text(json.dumps(dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),
            requested_url=url,parameters=params,url=r.url,http_status=r.status_code,
            content_type=r.headers.get('Content-Type'),bytes=len(r.content),
            sha256=hashlib.sha256(r.content).hexdigest()),indent=2)+'\n')
        print(name,r.status_code,len(r.content),flush=True)
        return p
    except requests.RequestException as exc:
        p.with_name(p.name+'.error.json').write_text(json.dumps(dict(url=url,parameters=params,
            attempted_utc=datetime.now(timezone.utc).isoformat(),error=str(exc)),indent=2)+'\n')
        print('FAILED',name,str(exc),flush=True)
        return None

def tap(name,query,service):
    return fetch(name+'.csv',BASES[service],dict(request='doQuery',lang='ADQL',format='csv',query=query))

def cone(name,table,ra,dec,radius=120):
    return fetch(name+'.tsv','https://vizier.cds.unistra.fr/viz-bin/asu-tsv',
        {'-source':table,'-c':f'{ra} {dec}','-c.rs':radius,'-out.all':'','-out.max':'unlimited','-out.form':'tsv'})

def csvrows(path):
    if not path.exists() or path.read_text().lstrip().startswith('<'):
        raise ValueError('Missing/failed CSV is not catalog absence: '+str(path))
    return list(csv.DictReader(io.StringIO(path.read_text())))

def tsvrows(path):
    """Parse ASU sections after their ruler, including zero-row catalogue responses."""
    if not path.exists():raise FileNotFoundError(path)
    text=path.read_text()
    if text.lstrip().startswith('<')or '#INFO\tError'in text or '#ERROR'in text:
        raise ValueError('Failed ASU response is not catalog absence: '+str(path))
    rows=[];header=None;after=False;table=None
    for line in text.splitlines():
        if line.startswith('#Name:'):table=line.split(':',1)[1].strip()
        if not line.strip() or line.startswith('#'):continue
        values=[x.strip() for x in line.split('\t')]
        if header is None:header=values;continue
        if all(not v or set(v)<=set('-') for v in values):after=True;continue
        if after and len(values)==len(header):rows.append(dict(zip(header,values),_table=table))
    return rows

def sep(a,d,b,e):
    a,d,b,e=map(math.radians,(a,d,b,e))
    h=math.sin((d-e)/2)**2+math.cos(d)*math.cos(e)*math.sin((a-b)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1,h))))*3600

def download():
    ids=','.join(t['source_id'] for t in TARGETS)
    jobs=[]
    for t in TARGETS:
        for key,cat in CATALOGS.items():jobs.append((cone,(t['source_id']+'_'+key,cat,t['ra'],t['dec'])))
    jobs += [(cone,('positive_algol_sb9',CATALOGS['sb9'],47.04221855,40.95564667)),
             (cone,('positive_algol_vsx',CATALOGS['vsx'],47.04221855,40.95564667)),
             (cone,('positive_algol_wds',CATALOGS['wds'],47.04221855,40.95564667)),
             (cone,('positive_traven',CATALOGS['traven2020'],52.928333,-69.362167)),
             (cone,('positive_asassn_validated',CATALOGS['asassn_v2'],45.09515,0.47689,30)),
             (tap,('gaia_dr2_neighbourhood',f'select * from gaiadr3.dr2_neighbourhood where dr3_source_id in ({ids})','gaia')),
             (tap,('wide_elbadry2021',f'select * from "J/MNRAS/506/2269/catalog" where Source1 in ({ids}) or Source2 in ({ids})','vizier')),
             (tap,('wide_positive','select * from "J/MNRAS/506/2269/catalog" where Source1=3770759051395978240 and Source2=3770759051395977472','vizier')),
             (tap,('galah_binary_metadata',"select table_name,description from TAP_SCHEMA.tables where description like '%GALAH%' or description like '%Traven%' or description like '%spectroscopic binar%' or description like '%Spectroscopic binar%'",'vizier'))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        for result in pool.map(lambda job:job[0](*job[1]),jobs):pass

def derive():
    CAND.mkdir(parents=True,exist_ok=True)
    manifests=[]
    for p in OUT.glob('*.manifest.json'):
        m=json.loads(p.read_text());raw=p.with_name(p.name.removesuffix('.manifest.json'))
        if hashlib.sha256(raw.read_bytes()).hexdigest()!=m['sha256']:raise RuntimeError(raw)
        manifests.append(dict(file=str(raw.relative_to(ROOT)),sha256=m['sha256'],http_status=m['http_status']))
    traven=ROOT/'data/stellar_systems/scout_novelty/traven2020_full.tsv'
    tm=json.loads(traven.with_name(traven.name+'.manifest.json').read_text())
    if hashlib.sha256(traven.read_bytes()).hexdigest()!=tm['sha256']:raise RuntimeError(traven)
    tr=tsvrows(traven)
    if len(tr)!=12760:raise RuntimeError(f'Full Traven catalogue incomplete {len(tr)}')
    dr2=csvrows(OUT/'gaia_dr2_neighbourhood.csv')
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'scope':'Bounded prior-companion recognition; null matches do not establish novelty.',
             'traven_full_rows':len(tr),'traven_full_sha256':tm['sha256'],'targets':[],
             'positive_controls':{},'manifests':manifests}
    summary['failed_positive_controls_preserved']={}
    for f in OUT.glob('positive*.tsv'):
        try:
            positive=tsvrows(f)
            summary['positive_controls'][f.name]=dict(rows=len(positive),sample=positive[:2])
        except ValueError as error:
            summary['failed_positive_controls_preserved'][f.name]=str(error)
    summary['positive_controls']['wide_positive.csv']=dict(rows=len(csvrows(OUT/'wide_positive.csv')))
    for name in ['positive_algol_wds.tsv','positive_algol_vsx.tsv','positive_algol_sb9.tsv',
                 'positive_traven.tsv','positive_asassn_validated.tsv','wide_positive.csv']:
        assert summary['positive_controls'][name]['rows']>0,'Positive retrieval control failed: '+name
    for t in TARGETS:
        t=dict(t);mapping=[r for r in dr2 if r.get('dr3_source_id')==t['source_id']]
        dids={r['dr2_source_id'] for r in mapping}
        t['dr2_mapping']=mapping
        t['traven_exact_dr2_matches']=[r for r in tr if r.get('GaiaDR2') in dids]
        t['traven_full_120arcsec_matches']=[dict(r,separation_arcsec=sep(t['ra'],t['dec'],float(r['RAJ2000']),float(r['DEJ2000']))) for r in tr if sep(t['ra'],t['dec'],float(r['RAJ2000']),float(r['DEJ2000']))<=120]
        t['catalog_cones']={key:tsvrows(OUT/(t['source_id']+'_'+key+'.tsv')) for key in CATALOGS}
        t['wide_exact_matches']=[r for r in csvrows(OUT/'wide_elbadry2021.csv') if t['source_id'] in (r.get('Source1'),r.get('Source2'))]
        summary['targets'].append(t)
    (CAND/'recognition.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('manifests','targets','positive_controls')},indent=2))
    for t in summary['targets']:
        print(t['source_id'],'DR2',t['dr2_mapping'],'traven_exact',len(t['traven_exact_dr2_matches']),
              'traven_cone',len(t['traven_full_120arcsec_matches']),
              'cones',{k:len(v) for k,v in t['catalog_cones'].items()},'wide',len(t['wide_exact_matches']))
    print('positive',{k:v['rows'] for k,v in summary['positive_controls'].items()})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--download',action='store_true');args=p.parse_args()
    if args.download:download()
    derive()
