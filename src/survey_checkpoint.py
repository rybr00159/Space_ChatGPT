"""Export/restore the full scientific SQLite state in checksum-verified chunks.

Original result_json strings are preserved exactly. Restoring does not download
observations or repeat searches. Existing databases/checkpoints are never replaced.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import uuid

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = Path('/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite')
DEFAULT_OUT = ROOT / 'data/broad_checkpoint'
COLUMNS = ['tic', 'status', 'result_json', 'error', 'updated_utc']


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def validate(directory):
    m = json.loads((directory / 'manifest.json').read_text())
    if m['format'] != 'Space_ChatGPT-broad-sqlite-chunks-v1' or m['columns'] != COLUMNS:
        raise ValueError('Unsupported checkpoint schema')
    for p in m['parts']:
        name = p['file']
        if Path(name).name != name:
            raise ValueError('Checkpoint part must be a local filename')
        f = directory / name
        if f.stat().st_size != p['bytes'] or digest(f) != p['sha256']:
            raise ValueError(f'Checkpoint integrity mismatch: {name}; preserve and investigate')
    return m


def export(db_path, directory, rows_per_part):
    if directory.exists():
        m = validate(directory)
        print('Existing verified checkpoint preserved:', m['status_counts'])
        return
    directory.mkdir(parents=True)
    c = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    c.execute('BEGIN')
    m = dict(format='Space_ChatGPT-broad-sqlite-chunks-v1', columns=COLUMNS,
             created_utc=datetime.now(timezone.utc).isoformat(), source_database=str(db_path),
             survey_meta=dict(c.execute('SELECT key,value FROM survey_meta')),
             status_counts=dict(c.execute('SELECT status,count(*) FROM stars GROUP BY status')),
             preservation='Original result_json strings; raw observations separately recorded in the observation manifest.',
             parts=[])
    records = iter(c.execute('SELECT tic,status,result_json,error,updated_utc FROM stars ORDER BY tic'))
    count = 0
    while True:
        first = next(records, None)
        if first is None:
            break
        p = directory / f'part{len(m["parts"]):03d}.jsonl.gz'
        n, last = 0, first
        with p.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as z:
                row = first
                while row is not None and n < rows_per_part:
                    z.write((json.dumps(list(row), separators=(',', ':')) + '\n').encode())
                    last = row
                    n += 1
                    if n < rows_per_part:
                        row = next(records, None)
        m['parts'].append(dict(file=p.name, rows=n, first_tic=first[0], last_tic=last[0], bytes=p.stat().st_size, sha256=digest(p)))
        count += n
        print('Preserved', count, 'rows', flush=True)
    c.close()
    if count != sum(m['status_counts'].values()):
        raise ValueError('Export row count mismatch; preserve incomplete checkpoint')
    m['rows'] = count
    m['compressed_bytes'] = sum(p['bytes'] for p in m['parts'])
    (directory / 'manifest.json').write_text(json.dumps(m, indent=2) + '\n')
    validate(directory)
    print('Export verified:', count, 'rows;', m['compressed_bytes'], 'bytes')


def restore(directory, db_path):
    m = validate(directory)
    if db_path.exists():
        raise FileExistsError('Existing database preserved; choose an unused --db path')
    db_path.parent.mkdir(parents=True, exist_ok=True)
    pending = db_path.with_name(db_path.name + '.' + uuid.uuid4().hex + '.partial')
    c = sqlite3.connect(pending)
    c.execute('CREATE TABLE survey_meta(key TEXT PRIMARY KEY,value TEXT)')
    c.execute('CREATE TABLE stars(tic INTEGER PRIMARY KEY,status TEXT,result_json TEXT,error TEXT,updated_utc TEXT)')
    c.executemany('INSERT INTO survey_meta VALUES(?,?)', m['survey_meta'].items())
    counts = Counter()
    try:
        for part in m['parts']:
            n = 0
            with gzip.open(directory / part['file'], 'rt') as f:
                for line in f:
                    row = json.loads(line)
                    if not isinstance(row, list) or len(row) != len(COLUMNS):
                        raise ValueError('Malformed checkpoint row')
                    if row[2] is not None and int(json.loads(row[2])['tic']) != row[0]:
                        raise ValueError('Result TIC differs from SQLite key')
                    c.execute('INSERT INTO stars VALUES(?,?,?,?,?)', row)
                    counts[row[1]] += 1
                    n += 1
            if n != part['rows']:
                raise ValueError('Part row count mismatch')
            c.commit()
        if dict(counts) != m['status_counts'] or sum(counts.values()) != m['rows']:
            raise ValueError('Restored checkpoint status/count mismatch')
        if c.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
            raise ValueError('Restored SQLite integrity failed')
        c.close()
        os.link(pending, db_path)  # Atomic non-overwriting publication.
        pending.unlink()  # Same verified bytes remain at the final path.
        print('Restored verified checkpoint:', dict(counts))
    except Exception:
        c.close()
        raise  # Partial database remains available for diagnosis.


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['export', 'restore', 'verify'])
    parser.add_argument('--directory', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--db', type=Path, default=DEFAULT_DB)
    parser.add_argument('--rows-per-part', type=int, default=40000)
    args = parser.parse_args()
    if args.rows_per_part < 1:
        parser.error('--rows-per-part must be positive')
    if args.action == 'export':
        export(args.db, args.directory, args.rows_per_part)
    elif args.action == 'restore':
        restore(args.directory, args.db)
    else:
        m = validate(args.directory)
        print('Verified checkpoint:', m['rows'], 'rows;', m['status_counts'])
