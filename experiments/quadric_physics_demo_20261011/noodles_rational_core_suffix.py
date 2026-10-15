"""Retain rational operands after the witnessed 73rd Boolean; no float feedback."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def main(case, capture, predecessor, output):
    import gmpy2
    import numpy as np
    import trimesh
    import _make_it_stand as native
    from core_boolean_mesh import bounded_rational_core_mesh
    output.mkdir(parents=True,exist_ok=False)
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    manifest=json.loads((case/'manifest.json').read_text())
    prefix=json.loads(predecessor.read_text())
    witness=json.loads((capture/'diagnosis.json').read_text())
    if prefix['manifest_sha256']!=sha(case/'manifest.json') or len(prefix['excluded_parts'])!=72 or \
        prefix['excluded_parts']!=manifest['construction']['fixed_parts'][:72] or \
        witness['predecessor_sha256']!=sha(predecessor) or \
        witness['part']!=manifest['construction']['fixed_parts'][72]:
        raise ValueError('FROZEN_RATIONAL_SUFFIX_CONTEXT_MISMATCH')
    def load(name):
        entry=manifest['parts'][name];path=case/entry['model_obj']
        if sha(path)!=entry['sha256']['model']:raise ValueError('FIXED_INPUT_CHANGED:'+name)
        return trimesh.load(path,force='mesh',process=False)
    inputs=[load('outer-ellipsoid')]+[load(name) for name in manifest['construction']['fixed_parts']]
    lo=np.min([m.bounds[0] for m in inputs],axis=0);hi=np.max([m.bounds[1] for m in inputs],axis=0)
    length=float(np.nextafter(np.linalg.norm(hi-lo),0))
    started=time.monotonic()
    with (output/'events.jsonl').open('x') as stream:
        def event(status,**values):
            stream.write(json.dumps(dict(status=status,seconds=time.monotonic()-started,**values))+'\n')
            stream.flush();os.fsync(stream.fileno())
        def check_indices(packet):
            faces=np.asarray(packet['faces'],dtype=np.int64)
            if not packet['succeeded'] or not len(faces) or faces.min()<0 or faces.max()>=len(packet['vertices']):
                raise ValueError('INVALID_RATIONAL_BOOLEAN_RESULT')
            edges=np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
            keys=np.sort(edges,axis=1)
            _,inv,counts=np.unique(keys,axis=0,return_inverse=True,return_counts=True)
            directions=np.where(edges[:,0]<edges[:,1],1,-1)
            sums=np.bincount(inv,weights=directions)
            if np.any(counts!=2) or np.any(sums!=0):
                event('RATIONAL_EDGE_TOPOLOGY_REJECTED',abnormal_edges=int(np.count_nonzero(counts!=2)))
                raise ValueError('RATIONAL_BOUNDARY_NOT_CLOSED_ORIENTED_EDGE_MANIFOLD')
            return faces
        packet=json.loads((capture/'original-rational-output.json').read_text())
        check_indices(packet)
        excluded=manifest['construction']['fixed_parts'][:73]
        event('RATIONAL_SUFFIX_STARTED',excluded=73,capture_sha256=sha(capture/'original-rational-output.json'),
              profile='CGAL_CORE_RATIONAL_SUFFIX_V1',length_scale=length)
        for name in manifest['construction']['fixed_parts'][73:]:
            event('RATIONAL_EXCLUSION_STARTED',part=name,excluded=len(excluded))
            part=load(name)
            cutter=[[str(gmpy2.mpq(float(x))) for x in row] for row in part.vertices]
            packet=native.libigl_exact_boolean_construction_from_rationals(
                packet['vertices'],np.asarray(packet['faces'],dtype=np.int32),
                cutter,part.faces,'difference')
            with (output/f'domain-{len(excluded)+1:03d}-rational.json').open('x') as saved:
                json.dump(dict(packet),saved,default=lambda value:value.tolist(),allow_nan=False)
            check_indices(packet)
            excluded.append(name)
            event('RATIONAL_EXCLUSION_FINISHED',part=name,excluded=len(excluded),triangles=len(packet['faces']))
        # Only at this export boundary create a binary64 derived representation.
        mesh,representation=bounded_rational_core_mesh(packet['vertices'],packet['faces'],length_scale=length)
        with (output/'representation.json').open('x') as saved:json.dump(representation,saved,allow_nan=False)
        with (output/'domain-081.json').open('x') as saved:
            json.dump(dict(schema='noodles-core-binary64-domain-v1',profile='CGAL_CORE_RATIONAL_SUFFIX_V1',
                vertices=mesh.vertices.tolist(),faces=mesh.faces.tolist(),excluded_parts=excluded,
                manifest_sha256=sha(case/'manifest.json'),geometry_admission=False,optimization_admission=False),saved)
        mesh.export(output/'design-core.obj')
        with (output/'complete.json').open('x') as saved:
            json.dump(dict(core_domain_exclusions_complete=True,excluded_parts=excluded,
                seconds=time.monotonic()-started,geometry_admission=False,optimization_admission=False),saved)
        event('RATIONAL_DOMAIN_COMPLETE_PENDING_INDEPENDENT_ADMISSION',excluded=len(excluded))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('case','capture','predecessor','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.case,a.capture,a.predecessor,a.output)
