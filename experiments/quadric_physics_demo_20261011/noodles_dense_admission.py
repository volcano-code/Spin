"""Explicit dense-v2 input adapter and honest broad-phase resource screen."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import trimesh
from noodles_admission import admit,validate_manifest


def resource_and_overlap_screen(case):
    case=Path(case);manifest=json.loads((case/'manifest.json').read_text())
    validate_manifest(manifest,'dense-v2')
    meshes={}
    for name,entry in manifest['parts'].items():
        path=case/entry['model_obj']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']['model']:
            raise ValueError('STALE_SCREEN_SOURCE:'+name)
        meshes[name]=trimesh.load(path,force='mesh',process=False)
    names=[n for n in meshes if n.startswith('tube-')]
    overlap=[]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            if np.all(np.minimum(meshes[a].bounds[1],meshes[b].bounds[1])-np.maximum(meshes[a].bounds[0],meshes[b].bounds[0])>=0):
                overlap.append([a,b])
    fixed_count=len(manifest['construction']['fixed_parts'])
    vertices=sum(len(m.vertices) for m in meshes.values());faces=sum(len(m.faces) for m in meshes.values())
    calls=4*fixed_count+1
    result=dict(schema='dense-noodles-broad-phase-resource-screen-v1',status='ESTIMATE_AND_BROAD_PHASE_ONLY',
        spec_sha256=manifest['spec_sha256'],tube_count=len(names),fixed_part_count=fixed_count,
        source_vertices=vertices,source_triangles=faces,source_binary64_array_bytes=vertices*3*8+faces*3*8,
        planned_exact_boolean_calls=calls,
        tube_pair_aabb_overlap_count=len(overlap),tube_pair_aabb_overlap_candidates=overlap,
        overlap_semantics='SOURCE_SOLIDS_MAY_OVERLAP; ONLY_EXACT_SET_UNION_CAN_VALIDATE_RESULT',
        individual_source_native_validity='PENDING_BOUNDED_ADMISSION',whole_union_validated=False,solver_admitted=False,
        time_reference=dict(v1_measured_seconds=229.,v1_boolean_calls=29,
            constant_cost_per_call_scenario_seconds=229.*calls/29,
            scope='SCENARIO_NOT_WALL_TIME_UPPER_BOUND; DENSE_INTERSECTION_COMPLEXITY_UNMEASURED'),
        peak_memory_upper_bound='UNKNOWN_EXACT_COORDINATE_AND_INTERSECTION_GROWTH',
        resource_protocol=dict(memory_bytes=24*1024**3,swap_bytes=0,maximum_seconds=3600),
        operation='INPUT_GEOMETRY_ONLY_NO_FROZEN_SOLVER_REQUEST_OR_SIMULATION')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path);parser.add_argument('--screen-only',action='store_true')
    args=parser.parse_args()
    if args.screen_only:print(json.dumps(resource_and_overlap_screen(args.input),indent=2,allow_nan=False))
    else:
        if args.output is None:parser.error('--output is required for admission')
        admit(args.input,args.output,input_version='dense-v2')
