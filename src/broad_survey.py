"""Checkpointed actual-flux QLP survey; catalog inventory is a separate count.

Retain every downloaded raw product. HTTP sessions reuse connections. A bounded
I/O queue feeds CPU workers; SQLite commits completed stars for interruption
recovery. No credential or paid-service access is used.
"""
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime,timezone
from pathlib import Path
import argparse,csv,gzip,hashlib,json,os,shutil,sqlite3,threading,time
import numpy as np
import requests
from astropy.io import fits
from archive import verify_fits_checksums
from fast_screen import screen_curve
from broad_crossmatch import CatalogIndex

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path(os.environ.get('TESS_DATA_DIR','/workspace/tess-research-data'))/'broad'
SURVEY=CACHE/'qlp_s0001'
_thread=threading.local(); _index=None; _method_hash=None
CONFIG={'n_refine':0,'n_peaks':2,'compute_full_stats':False}

def utc():return datetime.now(timezone.utc).isoformat()
def url_for(tic):
    s=str(tic).zfill(16)
    name=f'hlsp_qlp_tess_ffi_s0001-{s}_tess_v01_llc.fits'
    sub='/'.join(s[i:i+4] for i in range(0,16,4))
    return f'https://archive.stsci.edu/hlsps/qlp/s0001/{sub}/{name}'

def init_plan(size=160000):
    p=ROOT/'data/broad_inventory/qlp_s0001_selection.csv.gz'
    meta=p.with_suffix('.manifest.json')
    if p.exists():
        m=json.loads(meta.read_text())
        if m['selected_stars']!=size:raise ValueError('Existing selection has a different size; preserve rather than replace it')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=m['selection_sha256']:raise ValueError('Selection hash mismatch')
        with gzip.open(p,'rt') as f:return list(csv.DictReader(f)),m
    source=CACHE.parent/'broad_probe/qlp_s0001.csv'
    if not source.exists():
        source.parent.mkdir(parents=True,exist_ok=True)
        response=requests.get('https://archive.stsci.edu/hlsps/qlp/target_lists/s0001.csv',timeout=120)
        response.raise_for_status();source.write_bytes(response.content)
    with source.open() as f:rows=list(csv.DictReader(f))
    if len({r['#TIC_ID'] for r in rows})!=len(rows):raise ValueError('Duplicate TIC in inventory')
    rng=np.random.default_rng(20261013)
    ids=rng.permutation(len(rows))[:size]
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as z:
            import io
            text=io.TextIOWrapper(z,newline='')
            writer=csv.DictWriter(text,fieldnames=['order','tic','ra','dec','url'],lineterminator='\n');writer.writeheader()
            for order,i in enumerate(ids):
                r=rows[int(i)];writer.writerow(dict(order=order,tic=r['#TIC_ID'],ra=r['RA'],dec=r['DEC'],url=url_for(int(r['#TIC_ID']))))
            text.flush();text.detach()
    inv=ROOT/'data/broad_inventory/qlp_s0001_targets.csv.gz'
    with inv.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as z:z.write(source.read_bytes())
    m={'created_utc':utc(),'inventory_rows':len(rows),'selected_stars':len(ids),'seed':20261013,
       'selection_rule':'Fixed random permutation of all official sector-1 QLP targets; selected before flux screening',
       'source_url':'https://archive.stsci.edu/hlsps/qlp/target_lists/s0001.csv','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
       'selection_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'inventory_gzip_sha256':hashlib.sha256(inv.read_bytes()).hexdigest(),
       'sector':1,'product_version':'v01','cadence_minutes':30,'timing':'BTJD, TDB',
       'planned_raw_bytes_approx':size*80640,'warning':'Selection count is not the number of successfully analyzed flux series'}
    meta.write_text(json.dumps(m,indent=2))
    with gzip.open(p,'rt') as f:return list(csv.DictReader(f)),m

def download(row):
    tic=int(row['tic']);directory=SURVEY/'raw'/f'{tic%1000:03}'
    path=directory/Path(row['url']).name;directory.mkdir(parents=True,exist_ok=True)
    begin=time.monotonic()
    if path.exists():return {**row,'path':str(path),'download_seconds':0,'download_cached':True}
    if shutil.disk_usage(CACHE).free<5*1024**3:raise RuntimeError('Storage reserve reached: 5GiB required for candidate followup')
    if not hasattr(_thread,'session'):
        _thread.session=requests.Session()
        _thread.session.headers.update({'User-Agent':'Space_ChatGPT-public-TESS-research/0.3'})
    for attempt in range(4):
        try:
            response=_thread.session.get(row['url'],timeout=(15,60))
            if response.status_code in [429,500,502,503,504]:
                retry=float(response.headers.get('Retry-After',min(16,2**attempt)))
                response.close();time.sleep(min(30,max(1,retry)));continue
            response.raise_for_status();body=response.content;response.close()
            if not body.startswith(b'SIMPLE  '):raise ValueError('Not FITS bytes; do not cache an error document')
            tmp=path.with_suffix('.part');tmp.write_bytes(body);tmp.replace(path)
            return {**row,'path':str(path),'download_seconds':time.monotonic()-begin,'download_cached':False}
        except requests.RequestException:
            if attempt==3:raise
            time.sleep(2**attempt)
    raise RuntimeError('Archive retry limit reached; retry later after diagnosis')

def initialize_worker():
    global _index,_method_hash
    _index=CatalogIndex()
    _method_hash=hashlib.sha256((ROOT/'src/fast_screen.py').read_bytes()).hexdigest()

def finite(v):
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v

def analyze(record):
    tic=int(record['tic']);ra,dec=float(record['ra']),float(record['dec']);path=Path(record['path'])
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    with fits.open(path,memmap=False) as hd:
        hd.verify('exception');integrity=verify_fits_checksums(path,hd)
        header=hd[1].header;d=hd[1].data
        if header.get('TIMESYS',hd[0].header.get('TIMESYS'))!='TDB' or abs(header.get('BJDREFI',0)+header.get('BJDREFF',0)-2457000)>1e-8:
            raise ValueError('Invalid/undefined TESS timing; cannot assume BTJD')
        metadata={k:finite(header.get(k,hd[0].header.get(k))) for k in ['TICID','SECTOR','TESSMAG','TEFF','RADIUS','PROCVER','DATE','CAMERA','CCD']}
        t=np.asarray(d['TIME'],float);f=np.asarray(d['KSPSAP_FLUX'],float);e=np.asarray(d['KSPSAP_FLUX_ERR'],float);q=np.asarray(d['QUALITY'])
    target=_index.match_target(tic,ra,dec)
    result=screen_curve(t,f,e,q,config=CONFIG)
    result['label']='unmasked_broad'
    scans=[result]
    known=[]
    for r in target['matches']:
        if r['host_match'] and r['transiting'] and r['period_days'] and r['epoch_btjd'] is not None:
            duration=r['duration_days'] or 3/24
            if duration/r['period_days']>=.3:continue
            known.append({'period':r['period_days'],'epoch_bjd':r['epoch_btjd']+2457000,'duration_hours':duration*24})
    if known and not target['is_eb_host']:
        residual=screen_curve(t,f,e,q,known=known,config=CONFIG);residual['label']='known_masked_residual';scans.append(residual)
    for scan in scans:
        for peak in scan.get('peaks',[]):
            match=_index.match_signal(tic,peak['period'],ra,dec,tolerance=.01)
            peak['catalog']={k:match[k] for k in ['known_period_match','known_harmonic_match','known_rv_period_match','possible_nearby_period_contamination','period_status']}
            peak['matched_names']=[r['name'] for r in match['known_matches']]
    flags={k:v for k,v in target.items() if k.startswith('is_') or k in ['nearby_catalog_object','novelty_status']}
    return {**record,'tic':tic,'ra':ra,'dec':dec,'sha256':digest,'raw_bytes':path.stat().st_size,
            'fits_checksums':integrity,'method_sha256':_method_hash,'metadata':metadata,'catalog_flags':flags,'scans':scans,'checked_utc':utc()}

def connect():
    SURVEY.mkdir(parents=True,exist_ok=True)
    db=sqlite3.connect(SURVEY/'survey.sqlite');db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=NORMAL')
    db.execute('CREATE TABLE IF NOT EXISTS survey_meta(key TEXT PRIMARY KEY,value TEXT)')
    db.execute('CREATE TABLE IF NOT EXISTS stars(tic INTEGER PRIMARY KEY, status TEXT, result_json TEXT, error TEXT, updated_utc TEXT)')
    db.commit();return db

def counts(db):return dict(db.execute('SELECT status,count(*) FROM stars GROUP BY status'))

def export(db,plan_meta):
    out=ROOT/'results/broad';out.mkdir(parents=True,exist_ok=True)
    target_cols=['tic','sector','ra','dec','status','tmag','teff','radius','n_input','n_clean','n_searched','n_segments','period_grid_count','max_unvetted_snr','single_events','is_confirmed_host','is_toi_host','is_ctoi_host','is_eb_host','is_diamante_host','is_published_host','nearby_catalog_object','sha256','raw_bytes']
    signal_cols=['tic','sector','label','rank','period_days','epoch_btjd','duration_hours','depth_ppm','depth_error_ppm','snr_diagonal','local_snr','n_events','n_positive_events','screening_pass','guard_failures','known_period_match','known_harmonic_match','possible_nearby_period_contamination','matched_names']
    # Atomic replace; SQLite is the authoritative incremental checkpoint.
    with gzip.open(out/'targets.csv.gz.tmp','wt',newline='') as tf,gzip.open(out/'signals.csv.gz.tmp','wt',newline='') as sf,gzip.open(out/'observations.jsonl.gz.tmp','wt') as mf:
        tw=csv.DictWriter(tf,fieldnames=target_cols,lineterminator='\n');sw=csv.DictWriter(sf,fieldnames=signal_cols,lineterminator='\n');tw.writeheader();sw.writeheader()
        for tic,status,raw,error in db.execute('SELECT tic,status,result_json,error FROM stars ORDER BY tic'):
            if not raw:tw.writerow({'tic':tic,'sector':1,'status':status});continue
            d=json.loads(raw);scan=d['scans'][0];meta=d['metadata'];flags=d['catalog_flags']
            peaks=[p for s in d['scans'] for p in s['peaks']]
            row={'tic':tic,'sector':1,'ra':d['ra'],'dec':d['dec'],'status':status,'tmag':meta.get('TESSMAG'),'teff':meta.get('TEFF'),'radius':meta.get('RADIUS'),'n_input':scan['n_input'],'n_clean':scan['n_clean'],'n_searched':scan['n_searched'],'n_segments':scan['n_segments'],'period_grid_count':scan.get('search_metadata',{}).get('period_grid_count'),'max_unvetted_snr':max((p['snr'] for p in peaks),default=None),'single_events':len(scan.get('single_events',[])),'sha256':d['sha256'],'raw_bytes':d['raw_bytes']}
            row.update({k:flags.get(k) for k in target_cols if k.startswith('is_') or k=='nearby_catalog_object'});tw.writerow(row)
            for s in d['scans']:
                for rank,p in enumerate(s['peaks'],1):
                    sw.writerow({'tic':tic,'sector':1,'label':s['label'],'rank':rank,'period_days':p['period'],'epoch_btjd':p['t0_btjd'],'duration_hours':p['duration_days']*24,'depth_ppm':p['depth']*1e6,'depth_error_ppm':p['depth_err']*1e6,'snr_diagonal':p['snr'],'local_snr':p.get('local_snr'),'n_events':p['n_events'],'n_positive_events':p['n_positive_events'],'screening_pass':p['screening_pass'],'guard_failures':';'.join(p['guard_failures']),**{k:p['catalog'].get(k) for k in ['known_period_match','known_harmonic_match','possible_nearby_period_contamination']},'matched_names':';'.join(p['matched_names'])})
            obs={k:d[k] for k in ['tic','ra','dec','url','path','sha256','raw_bytes','fits_checksums','metadata','checked_utc']};mf.write(json.dumps(obs,allow_nan=False)+'\n')
    for name in ['targets.csv.gz','signals.csv.gz','observations.jsonl.gz']:(out/(name+'.tmp')).replace(out/name)
    stats=counts(db);summary={'updated_utc':utc(),'inventory_stars':plan_meta['inventory_rows'],'planned_stars':plan_meta['selected_stars'],'status_counts':stats,'actually_period_searched':stats.get('screened',0),'downloaded_but_insufficient':stats.get('insufficient_data',0),'errors':stats.get('error',0),'source':'QLP sector1 v01,30minute','source_selection':plan_meta,'config':CONFIG,'method_sha256':hashlib.sha256((ROOT/'src/fast_screen.py').read_bytes()).hexdigest(),'paid_calls':0,'discovery_claim':False,'counts_exclude_original_19_star_pilot':True}
    (out/'summary.json').write_text(json.dumps(summary,indent=2));return summary

def run(args):
    rows,meta=init_plan(args.planned);db=connect()
    method_hash=hashlib.sha256((ROOT/'src/fast_screen.py').read_bytes()).hexdigest()
    prior=db.execute("SELECT value FROM survey_meta WHERE key='method_sha256'").fetchone()
    if prior and prior[0]!=method_hash:raise ValueError('Fast screen implementation changed; preserve this database and start a separately versioned experiment')
    if not prior:
        db.execute("UPDATE stars SET status='needs_version_annotation' WHERE status IN ('screened','insufficient_data')")
        db.execute('INSERT INTO survey_meta VALUES(?,?)',('method_sha256',method_hash));db.commit()
    done={r[0] for r in db.execute("SELECT tic FROM stars WHERE status IN ('screened','insufficient_data')")}
    jobs=[r for r in rows[:args.limit or len(rows)] if int(r['tic']) not in done]
    print('PLAN',meta['inventory_rows'],'inventory;',meta['selected_stars'],'selected;',len(jobs),'pending in this run',flush=True)
    started=time.monotonic();completed=0;last_report=0;download_futures={};analysis_futures={};iterator=iter(jobs);exhausted=False
    with ThreadPoolExecutor(args.download_workers) as io,ProcessPoolExecutor(args.analysis_workers,initializer=initialize_worker) as cpu:
        while not exhausted or download_futures or analysis_futures:
            while not exhausted and len(download_futures)+len(analysis_futures)<args.queue:
                try:r=next(iterator)
                except StopIteration:exhausted=True;break
                download_futures[io.submit(download,r)]=r
            if not download_futures and not analysis_futures:break
            ready,_=wait(list(download_futures)+list(analysis_futures),timeout=10,return_when=FIRST_COMPLETED)
            for future in ready:
                if future in download_futures:
                    r=download_futures.pop(future)
                    try:rec=future.result();analysis_futures[cpu.submit(analyze,rec)]=r
                    except Exception as e:
                        db.execute('INSERT OR REPLACE INTO stars VALUES(?,?,?,?,?)',(int(r['tic']),'error',None,type(e).__name__+': '+str(e)[:180],utc()));completed+=1
                        if 'Storage reserve' in str(e):exhausted=True
                else:
                    r=analysis_futures.pop(future)
                    try:
                        result=future.result();status=result['scans'][0]['status']
                        db.execute('INSERT OR REPLACE INTO stars VALUES(?,?,?,?,?)',(int(r['tic']),status,json.dumps(result,allow_nan=False),None,utc()))
                    except Exception as e:db.execute('INSERT OR REPLACE INTO stars VALUES(?,?,?,?,?)',(int(r['tic']),'error',None,type(e).__name__+': '+str(e)[:180],utc()))
                    completed+=1
            db.commit()
            if completed-last_report>=args.report_every:
                last_report=completed
                c=counts(db);print('CHECKPOINT',c,'elapsed_s',round(time.monotonic()-started),'run_rate',round(completed/(time.monotonic()-started),2),flush=True)
                if completed>=args.report_every:export(db,meta)
            if args.max_run_stars and completed>=args.max_run_stars:exhausted=True
    summary=export(db,meta);db.close();print('FINAL',json.dumps(summary['status_counts']),round(time.monotonic()-started),'seconds',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--planned',type=int,default=160000);p.add_argument('--limit',type=int,default=0);p.add_argument('--download-workers',type=int,default=16);p.add_argument('--analysis-workers',type=int,default=4);p.add_argument('--queue',type=int,default=64);p.add_argument('--report-every',type=int,default=1000);p.add_argument('--max-run-stars',type=int,default=0)
    a=p.parse_args();run(a)
