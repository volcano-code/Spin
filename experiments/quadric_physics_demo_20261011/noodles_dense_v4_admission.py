"""Explicit v4 source admission only; no optimizer or simulation dispatch."""
import argparse
import json
from pathlib import Path
from noodles_admission import admit,validate_manifest


def resource_screen(case):
    record=json.loads((Path(case)/'manifest.json').read_text())
    validate_manifest(record,'dense-v4')
    parts=record['parts'];fixed=len(record['construction']['fixed_parts'])
    triangles=sum(p['audit']['triangles'] for p in parts.values())
    vertices=sum(p['audit']['vertices'] for p in parts.values())
    return dict(schema='dense-noodles-v4-resource-screen-v1',status='ESTIMATE_NOT_ADMISSION',
        spec_sha256=record['spec_sha256'],tube_count=len(parts)-4,fixed_parts=fixed,
        source_triangles=triangles,source_vertices=vertices,
        source_binary64_array_bytes=24*(triangles+vertices),planned_exact_boolean_calls=4*fixed+1,
        native_source_self_intersection='NOT_CHECKED',native_whole_union='NOT_CHECKED',
        peak_memory_upper_bound='UNKNOWN_INTERSECTION_AND_RATIONAL_GROWTH',
        wall_time_upper_bound='UNKNOWN_DENSE_INTERSECTION_COMPLEXITY',
        solver_admitted=False,optimization_invocations=0,simulation_invocations=0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path);parser.add_argument('--screen-only',action='store_true')
    args=parser.parse_args()
    if args.screen_only:print(json.dumps(resource_screen(args.input),indent=2,allow_nan=False))
    else:
        if args.output is None:parser.error('--output required')
        admit(args.input,args.output,input_version='dense-v4')
