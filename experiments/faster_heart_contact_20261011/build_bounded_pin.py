"""Explicit pin intersection companion; original pin and request stay unchanged."""
import argparse
import json
from pathlib import Path
import time
import gmpy2

from admit_bounded_request import save,text
from bounded_pin_correspondence import boundary_distance_certificate,numerical_pin_gate
from pin_role_witness import load_pin_boundary
from rational_request_view import sha


def run(failed_job,output,expected_native_sha256):
    import _make_it_stand as native
    import numpy as np
    from rounded_input import exact_moments
    from make_it_stand.verification.spin_problem import MeshSnapshot
    failed_job,output=Path(failed_job),Path(output)
    output.mkdir(parents=True,exist_ok=False)
    config=json.loads((failed_job/'config.json').read_text());identities=config['identity']
    keys=('authority-outer','authority-fixed','candidate-fixed.rational.json','candidate-design.rational.json',
          'pin','source_binding','original_admission','candidate','native')
    inputs={key:identities[key] for key in keys}
    began=time.monotonic();stage='IDENTITY'
    with (output/'events.jsonl').open('x') as events:
        def event(name):
            nonlocal stage
            import os
            stage=name;events.write(json.dumps(dict(stage=stage,elapsed_seconds=time.monotonic()-began))+'\n')
            events.flush();os.fsync(events.fileno());print(stage,flush=True)
        try:
            if sha(native.__file__)!=expected_native_sha256 or inputs['native']['sha256']!=expected_native_sha256:
                raise ValueError('NATIVE_IDENTITY_MISMATCH')
            for info in inputs.values():
                if sha(info['path'])!=info['sha256']:raise ValueError('STALE_PIN_CORRESPONDENCE_SOURCE')
            binding=json.loads(Path(inputs['source_binding']['path']).read_text())
            if binding.get('source_identity_verified')is not True or binding['frozen_source']['pin.obj']!=inputs['pin']['sha256'] or \
                    binding['admission_sha256']!=inputs['original_admission']['sha256']:
                raise ValueError('ORIGINAL_PIN_AUTHORITY_BINDING_MISMATCH')
            pin=load_pin_boundary(inputs['pin']['path'])
            fixed=json.loads(Path(inputs['candidate-fixed.rational.json']['path']).read_text())
            original_fixed=json.loads(Path(inputs['authority-fixed']['path']).read_text())
            design=json.loads(Path(inputs['candidate-design.rational.json']['path']).read_text())
            outer=json.loads(Path(inputs['authority-outer']['path']).read_text())
            ov=[[gmpy2.mpq(x) for x in row] for row in outer['vertices']]
            length=max(max(p[k] for p in ov)-min(p[k] for p in ov) for k in range(3))
            save(output/'identity.json',dict(schema='bounded-pin-companion-input-v1',inputs=inputs,
                failed_config_sha256=sha(failed_job/'config.json'),worker_sha256=sha(__file__),
                proof_algorithm_sha256=sha(Path(__file__).with_name('bounded_pin_correspondence.py')),
                original_authority_unchanged=True,request_geometry_unchanged=True,
                minimum_outer_wall_model='.015',normalization_length_authority_exact=text(length)))
            event('ORIGINAL_PIN_AND_FIXED_REGIONS')
            ps=MeshSnapshot.from_rationals(pin['vertices'],pin['faces'])
            fs=MeshSnapshot.from_rationals(fixed['vertices'],fixed['faces'])
            af=MeshSnapshot.from_rationals(original_fixed['vertices'],original_fixed['faces'])
            original_relation=native.relate_exact_material_regions(ps.exact_region(),af.exact_region())
            save(output/'original-pin-in-authority-fixed.json',original_relation)
            if not original_relation['completed'] or not original_relation['a_subset_b']:
                raise ValueError('ORIGINAL_PIN_NOT_IN_AUTHORITY_FIXED')
            event('SAME_FIXED_INTERSECTION_PIN')
            result=native.libigl_exact_boolean_construction_from_rationals(pin['vertices'],np.asarray(pin['faces'],dtype=np.int32),
                fixed['vertices'],np.asarray(fixed['faces'],dtype=np.int32),'intersection')
            if not result['succeeded']:raise ValueError('EXACT_PIN_INTERSECTION_FAILED')
            derived=dict(schema='bounded-fixed-pin-intersection-v1',vertices=[[text(gmpy2.mpq(str(x))) for x in row] for row in result['vertices']],
                faces=np.asarray(result['faces'],dtype=np.int32).tolist(),
                role='BOUNDED_DERIVED_PIN_NOT_ORIGINAL_PIN',construction='ORIGINAL_PIN_INTERSECT_SAME_DERIVED_FIXED',
                original_pin_sha256=inputs['pin']['sha256'],derived_fixed_sha256=inputs['candidate-fixed.rational.json']['sha256'])
            save(output/'derived-pin.rational.json',derived)
            ds=MeshSnapshot.from_rationals(derived['vertices'],derived['faces'])
            design_snapshot=MeshSnapshot.from_rationals(design['vertices'],design['faces'])
            relations=dict(derived_pin_in_fixed=native.relate_exact_material_regions(ds.exact_region(),fs.exact_region()),
                derived_pin_in_original_pin=native.relate_exact_material_regions(ds.exact_region(),ps.exact_region()),
                derived_pin_design=native.relate_exact_material_regions(ds.exact_region(),design_snapshot.exact_region()))
            save(output/'native-relations.json',relations)
            if any(not r['completed'] for r in relations.values()) or \
                    not relations['derived_pin_in_fixed']['a_subset_b'] or not relations['derived_pin_in_original_pin']['a_subset_b'] or \
                    relations['derived_pin_design']['positive_volume_overlap']:
                raise ValueError('ACTUAL_DERIVED_PIN_ROLE_FAILED')
            event('INDEPENDENT_MOMENTS_AND_ACTUAL_REMOVED_VOLUME')
            am=exact_moments(pin['vertices'],pin['faces']);dm=exact_moments(derived['vertices'],derived['faces'])
            removed=gmpy2.mpq(relations['derived_pin_in_original_pin']['symmetric_difference_volume_exact'])
            event('WHOLE_BOUNDARY_EXACT_CONVEX_COVER')
            boundary=boundary_distance_certificate(pin,derived,text(length/10**6))
            save(output/'whole-boundary-cover.json',boundary)
            numeric=numerical_pin_gate(boundary,[text(x) for x in am],[text(x) for x in dm],text(removed),text(length))
            save(output/'independent-numerical-gate.json',numeric)
            for info in inputs.values():
                if sha(info['path'])!=info['sha256']:raise ValueError('PIN_SOURCE_CHANGED_DURING_CERTIFICATION')
            status='PIN_CORRESPONDENCE_CERTIFIED' if numeric['status']=='CERTIFIED' else 'PIN_CORRESPONDENCE_'+numeric['status']
            save(output/'terminal.json',dict(schema='bounded-pin-companion-admission-v1',status=status,
                pin_correspondence_certified=numeric['status']=='CERTIFIED',
                authority_moments_exact=[text(x) for x in am],derived_moments_exact=[text(x) for x in dm],
                actual_removed_volume_exact=text(removed),
                geometry_admission_granted=False,solver_admission_granted=False,
                requires='FRESH_WHOLE_REQUEST_GATE_WITH_EXPLICIT_BOUNDED_PIN_COMPANION',
                inputs=inputs,files={p.name:sha(p) for p in output.iterdir() if p.is_file()},
                optimization_invocations=0,simulation_invocations=0,elapsed_seconds=time.monotonic()-began))
        except Exception as exc:
            save(output/'failure.json',dict(stage=stage,type=type(exc).__name__,message=str(exc),inputs=inputs,
                geometry_admission_granted=False,optimization_invocations=0,simulation_invocations=0,
                elapsed_seconds=time.monotonic()-began))
            raise


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--failed-job',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--expected-native-sha256',required=True)
    args=parser.parse_args();run(args.failed_job,args.output,args.expected_native_sha256)
