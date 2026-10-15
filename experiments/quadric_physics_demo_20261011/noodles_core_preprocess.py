"""Core CGAL design-domain construction, without repeated Nef certifications.

Each actual Boolean uses represented binary64 operands and exact CGAL
construction, then stores the binary64 mesh used by the next core operation.
This is a new numerical profile, not an upgrade of historical rational stages.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(case,output,previous=None,bounded_alias_chain=False):
    import gmpy2
    import numpy as np
    import trimesh
    import _make_it_stand as native
    from core_boolean_mesh import canonical_core_mesh, bounded_rational_core_mesh
    case,output=Path(case),Path(output);output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((case/'manifest.json').read_text())
    started=time.monotonic()
    with (output/'events.jsonl').open('x') as events:
        def event(status,**values):
            row=dict(status=status,seconds=time.monotonic()-started,**values)
            events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
            print(json.dumps(row),flush=True)
        def basic(vertices,faces,label):
            mesh=trimesh.Trimesh(vertices,faces,process=False)
            facts=dict(finite=bool(np.isfinite(mesh.vertices).all()),
                watertight=bool(mesh.is_watertight),winding_consistent=bool(mesh.is_winding_consistent),
                positive_volume=bool(mesh.volume>0),zero_area_faces=int(np.count_nonzero(mesh.area_faces==0)))
            if not all(facts[key] for key in ('finite','watertight','winding_consistent','positive_volume')) or facts['zero_area_faces']:
                event('CORE_GEOMETRY_REJECTED',label=label,facts=facts)
                raise ValueError('CORE_INVALID_BINARY64_BOUNDARY:'+label)
            return mesh
        def load(name):
            entry=manifest['parts'][name];path=case/entry['model_obj']
            if sha(path)!=entry['sha256']['model']:raise ValueError('SOURCE_CHANGED:'+name)
            mesh=trimesh.load(path,force='mesh',process=False)
            return basic(mesh.vertices,mesh.faces,name)
        excluded=[]
        domain=load('inner-ellipsoid-candidate')
        source_meshes=[load('outer-ellipsoid')]+[load(name) for name in manifest['construction']['fixed_parts']]
        lo=np.min([mesh.bounds[0] for mesh in source_meshes],axis=0)
        hi=np.max([mesh.bounds[1] for mesh in source_meshes],axis=0)
        frozen_length_scale=float(np.nextafter(np.linalg.norm(hi-lo),0.0))
        event('FROZEN_ASSEMBLY_LENGTH_SCALE',length_scale=frozen_length_scale,
              bounded_alias_chain=bounded_alias_chain)
        if previous is not None:
            previous=Path(previous);packet=json.loads(previous.read_text())
            excluded=list(packet['excluded_parts'])
            if packet['schema']!='noodles-core-binary64-domain-v1' or packet['manifest_sha256']!=sha(case/'manifest.json') or \
                    manifest['construction']['fixed_parts'][:len(excluded)]!=excluded:
                raise ValueError('CORE_PREDECESSOR_IDENTITY_OR_PREFIX_MISMATCH')
            domain=basic(packet['vertices'],packet['faces'],'core-predecessor')
            event('RESUME_CORE_PREFIX_NOT_EXACT_CERTIFICATE',excluded=len(excluded),predecessor_sha256=sha(previous))
        for name in manifest['construction']['fixed_parts'][len(excluded):]:
            event('CORE_CGAL_EXCLUSION_STARTED',part=name,excluded=len(excluded))
            part=load(name)
            result=native.libigl_exact_boolean_construction(domain.vertices,domain.faces,
                part.vertices,part.faces,'difference')
            if not result['succeeded']:raise ValueError('CORE_BOOLEAN_FAILED:'+name)
            if bounded_alias_chain:
                # Preserve the actual authority before deriving the next operand.
                with (output/f'boolean-{len(excluded)+1:03d}-rational.json').open('x') as stream:
                    json.dump(dict(result),stream,default=lambda value: value.tolist(),allow_nan=False)
                domain,repair=bounded_rational_core_mesh(result['vertices'],result['faces'],
                    length_scale=frozen_length_scale)
                with (output/f'boolean-{len(excluded)+1:03d}-representation.json').open('x') as stream:
                    json.dump(repair,stream,allow_nan=False)
                event('CORE_BOUNDED_REPRESENTATION',part=name,
                    **{k:v for k,v in repair.items() if k not in ('retained_original_face_indices','removed_face_witnesses')})
            else:
                vertices=np.asarray([[float(gmpy2.mpq(x)) for x in row] for row in result['vertices']])
                domain,repair=canonical_core_mesh(vertices,np.asarray(result['faces'],dtype=np.int32))
                event('CORE_ZERO_DISPLACEMENT_CANONICALIZATION',part=name,**repair)
            excluded.append(name)
            # A durable core checkpoint is numerical evidence, not a Nef certificate.
            packet=dict(schema='noodles-core-binary64-domain-v1',
                profile='CGAL_CORE_BOUNDED_ALIAS_CHAIN_V1' if bounded_alias_chain else 'CGAL_CORE_BINARY64_V1',
                vertices=domain.vertices.tolist(),faces=domain.faces.tolist(),excluded_parts=excluded,
                manifest_sha256=sha(case/'manifest.json'),geometry_admission=False,optimization_admission=False)
            target=output/f'domain-{len(excluded):03d}.json'
            with target.open('x') as stream:json.dump(packet,stream)
            event('CORE_CGAL_EXCLUSION_FINISHED',part=name,excluded=len(excluded),triangles=len(domain.faces))
        domain.export(output/'design-core.obj')
        with (output/'complete.json').open('x') as stream:
            json.dump(dict(core_domain_exclusions_complete=True,excluded_parts=excluded,
                seconds=time.monotonic()-started,geometry_admission=False,optimization_admission=False),stream)
        event('CORE_DOMAIN_COMPLETE_NOT_INDEPENDENTLY_AUDITED',excluded=len(excluded))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--previous-core',type=Path)
    parser.add_argument('--bounded-alias-chain',action='store_true')
    args=parser.parse_args();main(args.case,args.output,args.previous_core,args.bounded_alias_chain)
