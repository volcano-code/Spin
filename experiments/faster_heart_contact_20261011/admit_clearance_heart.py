"""Independent exact admission of the new distance-domain candidate."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
import trimesh
from admit_rounded_input import construct, exact_rows, save_region
from rounded_input import SPEC, capacity_bounds, exact_moments, publish
from clearance_source_binding import check_candidate_authority_inputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--old-authority', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-native-sha256', required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    began = time.monotonic()
    with (args.output/'events.jsonl').open('x') as events:
        def event(stage):
            row=dict(stage=stage, elapsed_seconds=time.monotonic()-began)
            events.write(json.dumps(row)+'\n'); events.flush(); os.fsync(events.fileno())
            print(json.dumps(row), flush=True)
        try:
            terminal=json.loads((args.candidate/'terminal.json').read_text())
            construction=json.loads((args.candidate/'construction.json').read_text())
            if terminal['status'] != 'CANDIDATE_PREPARED_NOT_ADMITTED':
                raise ValueError('INCOMPLETE_CANDIDATE')
            for name, sha in terminal['files'].items():
                if hashlib.sha256((args.candidate/name).read_bytes()).hexdigest() != sha:
                    raise ValueError('STALE_CANDIDATE: '+name)
            old=json.loads((args.old_authority/'admission.json').read_text())
            old_config=json.loads((args.old_authority.parent/'config.json').read_text())
            check_candidate_authority_inputs(args.candidate,construction,old_config['identity'],SPEC)
            outer_path=args.old_authority/'outer.rational.json'
            if hashlib.sha256(outer_path.read_bytes()).hexdigest() != old['files'][outer_path.name]:
                raise ValueError('STALE_OUTER_AUTHORITY')
            import _make_it_stand as native
            native_sha=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()
            if native_sha != args.expected_native_sha256:
                raise ValueError('NATIVE_IDENTITY_MISMATCH')
            from make_it_stand.verification.spin_problem import MeshSnapshot, SpinArtifactProblem
            event('LOAD_UNCHANGED_OUTER_AND_NEW_DOMAIN')
            saved=json.loads(outer_path.read_text()); outer=saved['vertices'],saved['faces']
            candidate=trimesh.load(args.candidate/'inner-source.obj', force='mesh', process=False)
            design=exact_rows(candidate.vertices),candidate.faces.tolist()
            # Check the candidate itself: intersection must not silently repair it.
            outer_snapshot=MeshSnapshot.from_rationals(*outer)
            design_snapshot=MeshSnapshot.from_rationals(*design)
            pin_mesh=trimesh.load(args.candidate/'pin.obj',force='mesh',process=False)
            pin_snapshot=MeshSnapshot.from_rationals(exact_rows(pin_mesh.vertices),pin_mesh.faces.tolist())
            pin_relation=native.relate_exact_material_regions(design_snapshot.exact_region(),
                                                             pin_snapshot.exact_region())
            publish(args.output/'pin-exclusion.json',dict(pin_relation))
            if not pin_relation['completed'] or pin_relation['positive_volume_overlap']:
                raise ValueError('FIXED_PIN_EXCLUSION_NOT_CERTIFIED')
            event('EXACT_FIXED_MATERIAL_CONSTRUCTION')
            fixed=construct(native,outer,design,'difference')
            save_region(args.output/'fixed.obj',fixed)
            save_region(args.output/'design.obj',design)
            save_region(args.output/'outer.obj',outer)
            event('EXACT_PARTITION_AND_DOMAIN_CHECK')
            fixed_snapshot=MeshSnapshot.from_rationals(*fixed)
            problem=SpinArtifactProblem(original=outer_snapshot,
                fixed_regions=[('fixed-heart-and-pin',fixed_snapshot,1.)],
                design_regions=[('distance-eroded-heart',design_snapshot,1.)],
                length_unit='25mm',mass_unit='18.75g',
                target=dict(mode='SPINNING_TOP',axis=[0.,0.,1.],origin=[0.,0.,0.],support_points=[]),
                minimum_outer_wall=construction['minimum_clearance_model'])
            event('EXACT_BOUNDARY_CLEARANCE')
            clearance=native.check_exact_boundary_clearance(outer_snapshot.exact_region(),
                design_snapshot.exact_region(),construction['minimum_clearance_model'])
            publish(args.output/'clearance.json',clearance)
            event('INDEPENDENT_EXACT_MOMENTS_AND_CAPACITY')
            moments={name:exact_moments(*region) for name,region in
                     [('outer',outer),('fixed',fixed),('design',design)]}
            delta=[a-b-c for a,b,c in zip(moments['outer'],moments['fixed'],moments['design'])]
            if any(delta):raise ValueError('PARTITION_MOMENTS_MISMATCH')
            capacity=capacity_bounds(*fixed,*design)
            publish(args.output/'capacity.json',capacity)
            passed=clearance['completed'] and clearance['at_least_minimum'] is True
            status=('INPUT_GEOMETRY_BLOCKED' if not passed else
                    'INPUT_CAPACITY_EXCLUDED' if capacity['status']=='EXCLUDED_BY_FIRST_MOMENT_BOUND'
                    else 'INPUT_ADMITTED_FEASIBILITY_UNKNOWN')
            publish(args.output/'admission.json',dict(status=status,
                problem_digest=problem.input_digest, native_sha256=native_sha,
                outer_authority_source_sha256=old['files'][outer_path.name],
                construction_sha256=hashlib.sha256((args.candidate/'construction.json').read_bytes()).hexdigest(),
                partition_moment_delta_exact=[str(x) for x in delta],
                independent_moments_exact={name:[str(x) for x in values] for name,values in moments.items()},
                optimization_invocations=0,simulation_invocations=0,
                files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.glob('*')
                       if p.suffix in ('.obj','.json')}))
            publish(args.output/'terminal.json',dict(status=status,completed=True,
                optimization_invocations=0,simulation_invocations=0,
                elapsed_seconds=time.monotonic()-began))
            event(status)
        except Exception as exc:
            publish(args.output/'failure.json',dict(type=type(exc).__name__,message=str(exc),
                completed=False,elapsed_seconds=time.monotonic()-began))
            raise


if __name__ == '__main__':main()
