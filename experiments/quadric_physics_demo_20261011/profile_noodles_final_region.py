"""Read-only diagnostic of a FAILED boundary, never a preprocessing certificate."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time


def main(source,out):
    out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    stream=(out/'events.jsonl').open('x')
    def event(stage,**data):
        row=dict(stage=stage,seconds=time.monotonic()-start,**data)
        stream.write(json.dumps(row)+'\n');stream.flush();os.fsync(stream.fileno())
        print(json.dumps(row),flush=True)
    event('READ_FAILED_ARTIFACT_NOT_ADMITTED')
    raw=source.read_bytes();digest=hashlib.sha256(raw).hexdigest();packet=json.loads(raw)
    event('JSON_PARSED',sha256=digest,vertices=len(packet['vertices']),faces=len(packet['faces']))
    from make_it_stand.verification.spin_problem import MeshSnapshot
    def profile(frame,kind,arg):
        name=frame.f_code.co_name if kind in ('call','return') else getattr(arg,'__name__','')
        if name in ('_canonical_rational_vertices','make_exact_material_region_from_rationals'):
            event(kind+':'+name)
    sys.setprofile(profile)
    try:
        event('ORIGINAL_SNAPSHOT_CONSTRUCTOR_BEGIN')
        region=MeshSnapshot.from_rationals(packet['vertices'],packet['faces']).exact_region()
        event('ORIGINAL_SNAPSHOT_CONSTRUCTOR_END')
    finally:sys.setprofile(None)
    (out/'diagnostic.json').write_text(json.dumps(dict(scope='FAILED_ARTIFACT_COST_DIAGNOSTIC_ONLY',
        source=str(source),source_sha256=digest,facts=region.facts,
        geometry_admission=False,predecessor_admission=False,optimization_admission=False),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.source,a.output)
