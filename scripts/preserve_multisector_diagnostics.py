"""Retain bulk PNG/NPZ diagnostics outside ordinary Git, without deleting evidence.

All original bytes remain in this cache or the curated candidate folders. The
manifest maps original repository-relative paths to the retained physical file.
Reproduction uses the recorded raw observations and versioned source snapshots;
external cloud filesystem persistence has not been independently established.
"""
from pathlib import Path
import hashlib,json,shutil
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'candidates/multisector';CACHE=Path('/workspace/tess-research-data/multisector/diagnostics_archive');OUT=ROOT/'data/multisector/diagnostic_preservation.json'
KEEP_ALL={'tic229146531','tic368248545','tic183373155','tic71274048','tic4206066'}
KEEP_NAMES={'tic182517223':{'joint.png','joint_arrays.npz','sector95_empirical_p2.207684/diagnostics.png','sector95_empirical_p2.207684/diagnostic_arrays.npz'},'tic353933809':{'isolated_s86/rim_plane/difference_redistribution_test.png','isolated_s86/rim_plane/native_difference_diagnostics.png','sector86_same_detector_event_check.png'}}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 old=json.loads(OUT.read_text())if OUT.exists()else{};records={r['original_path']:r for r in old.get('files',[])}
 for path in sorted(BASE.rglob('*')):
  if not path.is_file()or path.suffix not in {'.png','.npz'}:continue
  relative=path.relative_to(BASE);tic=relative.parts[0];key=str(path.relative_to(ROOT));digest=sha(path);size=path.stat().st_size;keep=tic in KEEP_ALL or str(Path(*relative.parts[1:]))in KEEP_NAMES.get(tic,set())
  destination=path if keep else CACHE/key
  if not keep:
   destination.parent.mkdir(parents=True,exist_ok=True)
   if destination.exists():
    if sha(destination)!=digest:raise ValueError('Conflicting existing preserved evidence: '+str(destination))
    # Existing destination verified before original is removed; no byte loss.
    path.unlink()
   else:shutil.move(str(path),str(destination))
  if destination.stat().st_size!=size or sha(destination)!=digest:raise ValueError('Preservation identity failure')
  if key in records and records[key]['sha256']!=digest:raise ValueError('Previously preserved scientific artifact changed')
  records[key]=dict(original_path=key,retained_path=str(destination),sha256=digest,bytes=size,tracked_curated=keep)
 for r in records.values():
  if sha(Path(r['retained_path']))!=r['sha256']:raise ValueError('Prior preserved evidence differs')
 result=dict(updated_utc=datetime.now(timezone.utc).isoformat(),script_sha256=sha(Path(__file__)),files=list(records.values()),all_bytes_verified=True,external_persistence_verified=False,scope=__doc__,moved_files=sum(not r['tracked_curated']for r in records.values()),moved_bytes=sum(r['bytes']for r in records.values()if not r['tracked_curated']),curated_bytes=sum(r['bytes']for r in records.values()if r['tracked_curated']))
 OUT.write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items()if k not in ['files','scope']})

if __name__=='__main__':main()
