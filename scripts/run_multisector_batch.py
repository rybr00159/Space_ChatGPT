#!/usr/bin/env python3
"""Named follow-on batches using the unchanged screening and recovery methods.

Only data/output paths are configured here; scientific parameters are imported
from the already frozen modules. A tighter metadata-only stellar selection is
recorded before any flux is inspected in this batch.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,hashlib,json,sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import multisector_survey as survey
import multisector_missing_errors as recovery


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def configure(name):
    if not name or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_'for c in name):
        raise ValueError('Use a simple lowercase batch name')
    cache=Path('/workspace/tess-research-data/multisector')/name
    data=ROOT/'data/multisector'/name
    output=ROOT/'results/multisector_batches'/name
    survey.CACHE=cache;survey.DATA=data;survey.OUT=output;survey.DB=cache/'survey.sqlite'
    recovery.CACHE=cache;recovery.MAIN_DB=survey.DB;recovery.RECOVERY_DB=cache/'error_recovery.sqlite';recovery.OUT=output/'error_recovery'
    return cache,data,output


def plan(name,sectors,per_sector):
    cache,data,output=configure(name);data.mkdir(parents=True,exist_ok=True)
    path=data/'bright_cool_metadata.csv';manifest=data/'batch.json'
    if path.exists():
        meta=json.loads(manifest.read_text())
        if meta['sectors']!=sectors or meta['metadata_sha256']!=sha(path):raise ValueError('Existing batch differs')
        return survey.plan(path,per_sector)
    sources=[];rows=[]
    for sector in sectors:
        source=ROOT/'data/multisector'/f'inventory_s{sector:04d}_metadata.csv'
        sources.append(dict(path=str(source.relative_to(ROOT)),sha256=sha(source)))
        for row in csv.DictReader(source.open()):
            try:
                good=0<float(row['tmag'])<=11.5 and 3000<=float(row['teff'])<=6000 and .15<=float(row['radius'])<=1.2
            except (ValueError,TypeError):good=False
            if good:rows.append(row)
    if not rows:raise ValueError('No metadata-eligible stars')
    with path.open('w',newline='')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    meta=dict(created_utc=datetime.now(timezone.utc).isoformat(),name=name,sectors=sectors,
              metadata_sha256=sha(path),source_metadata=sources,driver_sha256=sha(__file__),
              preflux_rule='Tmag<=11.5,Teff3000..6000K,radius.15..1.2Rsun,all finite; no flux or candidate score in selection.',
              requested_max_per_sector=per_sector,eligible_rows=len(rows),
              note='New dataset/path configuration only; original scientific screen and explicit empirical-error alternative are unchanged.')
    manifest.write_text(json.dumps(meta,indent=2)+'\n')
    return survey.plan(path,per_sector)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['plan','run','recover','export'])
    parser.add_argument('--name',required=True);parser.add_argument('--sectors',nargs='+',type=int)
    parser.add_argument('--per-sector',type=int,default=1024);parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--downloads',type=int,default=12);parser.add_argument('--limit',type=int,default=0)
    args=parser.parse_args()
    if args.mode=='plan':print('BATCH SELECTED',len(plan(args.name,args.sectors,args.per_sector)))
    else:
        configure(args.name)
        if args.mode=='run':survey.run(args.limit,args.workers,args.downloads,False)
        elif args.mode=='recover':recovery.run(limit=args.limit)
        else:
            conn=survey.connect();survey.export(conn);conn.close()
