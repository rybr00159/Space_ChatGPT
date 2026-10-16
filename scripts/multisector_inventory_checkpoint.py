#!/usr/bin/env python3
"""Verify/restore exact public TIC snapshots and retrieve original QLP inventories.

No existing scientific file is overwritten. Changed live responses are retained
as unexpected files and are not accepted as the original research input.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import tarfile

import requests

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = Path('/workspace/tess-research-data/multisector')
MANIFEST = ROOT / 'data/multisector/archive_probe_query_archive.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_archive():
    meta = json.loads(MANIFEST.read_text())
    archive = ROOT / meta['archive_path']
    if sha(archive) != meta['sha256']:
        raise ValueError('Query archive checksum mismatch')
    for record in meta.get('source_files', []):
        if sha(ROOT / record['path']) != record['sha256']:
            raise ValueError('Source/artifact changed: ' + record['path'])
    return meta, archive


def verify(cache):
    meta, archive = validate_archive()
    expected = {r['path']:r for r in meta['files']}
    count, retained = 0, 0
    with tarfile.open(archive, 'r:xz') as stream:
        for member in stream:
            if not member.isfile() or member.name not in expected:
                raise ValueError('Unexpected archive member')
            record = expected.pop(member.name)
            body = stream.extractfile(member).read()
            if len(body) != record['bytes'] or hashlib.sha256(body).hexdigest() != record['sha256']:
                raise ValueError('Archived response changed: ' + member.name)
            target = cache / member.name
            if target.exists():
                if sha(target) != record['sha256']:
                    raise ValueError('Retained cache file changed: ' + str(target))
                retained += 1
            count += 1
    if expected:
        raise ValueError('Archive has missing members')
    print('VERIFIED exact archive members', count, 'retained cache members', retained,
          'source/artifact hashes', len(meta.get('source_files', [])))


def restore_queries(cache):
    meta, archive = validate_archive()
    expected = {r['path']:r for r in meta['files']}
    restored = 0
    with tarfile.open(archive, 'r:xz') as stream:
        for member in stream:
            if not member.isfile() or member.name not in expected:
                raise ValueError('Unexpected query archive member')
            target = (cache / member.name).resolve()
            if not target.is_relative_to(cache.resolve()):
                raise ValueError('Unsafe archived path')
            record = expected.pop(member.name)
            body = stream.extractfile(member).read()
            if hashlib.sha256(body).hexdigest() != record['sha256']:
                raise ValueError('Invalid archived response')
            if target.exists():
                if sha(target) != record['sha256']:
                    raise ValueError('Refusing to overwrite changed input: ' + str(target))
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(body); restored += 1
    if expected:
        raise ValueError('Archive has missing members')
    print('RESTORED exact query snapshots', restored)
    verify(cache)


def download_inventories(cache):
    probe = json.loads((ROOT / 'data/multisector/archive_probe.json').read_text())
    for record in probe['inventories']:
        target = cache / 'inventories' / Path(record['path']).name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if sha(target) != record['sha256']:
                raise ValueError('Refusing changed existing inventory: ' + str(target))
            print('CACHED inventory', record['sector']); continue
        response = requests.get(record['url'], timeout=(15, 150))
        response.raise_for_status(); body = response.content
        if hashlib.sha256(body).hexdigest() != record['sha256']:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
            unexpected = target.with_suffix('.unexpected_' + stamp + '.csv')
            unexpected.write_bytes(body)
            raise ValueError('Live inventory changed; saved separately: ' + str(unexpected))
        target.write_bytes(body)
        target.with_suffix('.manifest.json').write_text(json.dumps(record, indent=2) + '\n')
        print('RETRIEVED original inventory', record['sector'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['verify', 'restore-queries', 'download-inventories'])
    parser.add_argument('--cache-root', type=Path, default=DEFAULT_CACHE)
    args = parser.parse_args()
    {'verify': verify, 'restore-queries': restore_queries,
     'download-inventories': download_inventories}[args.mode](args.cache_root)
