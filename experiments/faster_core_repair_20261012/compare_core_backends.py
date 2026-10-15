"""Frozen same-operand component comparison; no optimization success claims."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from phase_timing import PhaseTiming

ROOT = Path(__file__).resolve().parents[2]


def safe(value):
    if isinstance(value, np.ndarray): return safe(value.tolist())
    if isinstance(value, np.generic): return safe(value.item())
    if isinstance(value, dict): return {str(k):safe(v) for k,v in value.items()}
    if isinstance(value, (tuple,list)): return [safe(v) for v in value]
    return value


def main(output, probe):
    import trimesh
    import gmpy2
    from make_it_stand import spin
    sys.path.insert(0, str(ROOT/'experiments/spin_it_reproduction/strict'))
    from run_faster_reference_boolean import _write_obj, _read_obj, _summarize_vertex_construction_provenance
    from make_it_stand.verification.faster_material import read_exact_boolean_sidecar_coordinates
    from make_it_stand.verification.spin_problem import MeshSnapshot
    from make_it_stand.verification.cherchi_motion_differential import (
        normalize_closed_component_orientation, recover_construction_vertex_motion,
        differentiate_polyhedral_raw_moments)
    output.mkdir(parents=True, exist_ok=False)
    def save(name, value):
        with (output/name).open('x') as stream: json.dump(safe(value), stream, allow_nan=False)
    native = spin._native_module()
    domain = trimesh.creation.box(extents=(2.,2.,2.))
    outer = trimesh.creation.box(extents=(2.2,2.2,2.2))
    fixed = native.libigl_exact_boolean_construction(outer.vertices,outer.faces,
        domain.vertices,domain.faces,'difference')
    if not fixed['succeeded']: raise ValueError('FIXED_MATERIAL_CONSTRUCTION_FAILED')
    fixed_vertices = np.asarray([[float(gmpy2.mpq(x)) for x in row] for row in fixed['vertices']])
    fixed_mesh = trimesh.Trimesh(fixed_vertices, fixed['faces'], process=False)
    if not fixed_mesh.is_watertight or not fixed_mesh.is_winding_consistent or fixed_mesh.volume<=0:
        raise ValueError('FIXED_MATERIAL_INVALID')
    fixed_moments = spin.integrate_raw_mass_moments(fixed_vertices.tolist(), np.asarray(fixed['faces']).tolist(), .75)
    request = spin.SpinItFasterRequest(outer_vertices=outer.vertices.tolist(),outer_faces=outer.faces.tolist(),
        design_vertices=domain.vertices.tolist(),design_faces=domain.faces.tolist(), fixed_shell_moments=fixed_moments,
        material_density=2., mode='SPINNING_TOP',center_of_mass_weight=1.)
    settings = spin.SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
        icosphere_subdivisions=2,maximum_icosphere_subdivisions=2,retain_fixed_normal_provenance=True,
        cap_radius_policy='GROW_TO_CERTIFIED_V1',exact_binary64_self_intersection=True,
        exact_binary64_source_boundary=True)
    # These are explicit diagnostic states, never supplied as seeded solver starts.
    states = dict(ellipsoid=(-.55,.02,-.03,.01,1.,1.3,.8,0.,0.,0.),
        one_sheet=(-.15,.02,-.03,.01,1.,1.2,-.8,0.,0.,0.),
        two_sheet=(.15,.02,-.03,.01,1.,1.2,-.8,0.,0.,0.))
    identity = dict(schema='faster-core-same-operands-v1',request=asdict(request),settings=asdict(settings),
        states=states, fixed_material=dict(vertices=fixed_vertices,faces=fixed['faces'],density=.75),
        native_path=native.__file__,native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        probe_path=str(probe),probe_sha256=hashlib.sha256(probe.read_bytes()).hexdigest(),
        optimization_admission=False)
    save('frozen-inputs.json', identity)
    _, native_request = spin._native_request(request)
    rows = []
    for name, coefficients in states.items():
        q = np.asarray(coefficients);q /= np.linalg.norm(q)
        state_dir = output/name;state_dir.mkdir()
        for repetition in range(4):
            timer = PhaseTiming()
            configured = spin._native_evaluator_settings(native, settings)
            configured.diagnostic_stage_observer = timer
            started = time.perf_counter()
            try:
                evaluated = native.evaluate_spin_it_faster_quadric(native_request,q,configured)
                elapsed = time.perf_counter()-started
                record = dict(state=name,repetition=repetition,cgal_value_seconds=elapsed,phases=timer.finish(),
                    derivative_defined=evaluated['derivative_defined'],diagnostic=evaluated['derivative_diagnostic'])
                if repetition == 0:
                    sv=np.asarray(evaluated['interface_vertices']);sf=np.asarray(evaluated['interface_faces'])
                    _write_obj(state_dir/'design.obj',domain.vertices,domain.faces)
                    _write_obj(state_dir/'source.obj',sv,sf)
                    _write_obj(state_dir/'cgal-material.obj',evaluated['material_vertices'],evaluated['material_faces'])
                    save(name+'-cgal-evaluation.json',evaluated)
                    # Separate diagnostic call, explicitly counted: obtain the production Jacobian.
                    full_timer=PhaseTiming();begin=time.perf_counter()
                    full=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings,
                        diagnostic_stage_observer=full_timer)
                    record['separate_fullrank_diagnostic_seconds']=time.perf_counter()-begin
                    record['fullrank_phases']=full_timer.finish()
                    save(name+'-fullrank-evaluation.json',asdict(full))
                    operation='intersection' if evaluated['quadric_solid_is_positive_region'] else 'subtraction'
                    tri=native.SpinTriangulationSettings()
                    for key in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                                'retain_fixed_normal_provenance'):
                        setattr(tri,key,getattr(settings,key))
                    source=native.evaluate_spin_it_faster_source_motion(q,domain.vertices,tri)
                    if not source['defined'] or not np.array_equal(source['faces'],sf) or not np.array_equal(source['vertices'],sv):
                        raise ValueError('SOURCE_MOTION_DOES_NOT_MATCH_ACTUAL_BOOLEAN_OPERAND')
                    for cherchi_repetition in range(4):
                        material=state_dir/f'cherchi-{cherchi_repetition}.obj'
                        begin=time.perf_counter()
                        completed=subprocess.run([str(probe),operation,str(state_dir/'design.obj'),
                            str(state_dir/'source.obj'),str(material),'--dump-exact-rational'],
                            capture_output=True,text=True,timeout=300)
                        cherchi_seconds=time.perf_counter()-begin
                        if completed.returncode: raise ValueError('CHERCHI_BOOLEAN_FAILED:'+completed.stderr)
                        cv,cf=_read_obj(material)
                        rational_coordinates,_=read_exact_boolean_sidecar_coordinates(cv,Path(str(material)+'.exact.tsv'))
                        snapshot=MeshSnapshot.from_rationals(rational_coordinates,cf)
                        exact,faces=snapshot.arrays()
                        faces,orientation=normalize_closed_component_orientation(exact,faces,rational_vertices=snapshot.rational_vertices)
                        moments=np.asarray(spin.integrate_raw_mass_moments(exact.tolist(),faces.tolist(),request.material_density))
                        target=np.asarray(evaluated['variable_raw_moments'])
                        item=dict(state=name,repetition=cherchi_repetition,operation=operation,
                            cherchi_process_boolean_and_export_seconds=cherchi_seconds,
                            raw_moment_max_error=float(np.max(np.abs(moments-target))),orientation=orientation,
                            entity_equivalence='NOT_YET_CERTIFIED',derivative_comparison='NOT_YET_CHECKED')
                        if cherchi_repetition == 0:
                            begin=time.perf_counter()
                            construction=_summarize_vertex_construction_provenance(
                                (state_dir/'design.obj',state_dir/'source.obj'),material,Path(str(material)+'.construction.tsv'))
                            velocity=np.asarray(source['vertex_velocity']).reshape(len(sv),3,10)
                            motion=recover_construction_vertex_motion(construction['records'],domain.vertices,sv,velocity,exact)
                            jac=differentiate_polyhedral_raw_moments(exact,faces,motion,density=request.material_density,
                                rational_vertices=snapshot.rational_vertices)
                            item['cherchi_provenance_derivative_seconds']=time.perf_counter()-begin
                            item['jacobian_max_error']=None if full.coefficient_jacobian is None else float(np.max(np.abs(jac-np.asarray(full.coefficient_jacobian))))
                            item['derivative_comparison']='CHECKED_SAME_OPERANDS'
                            save(name+'-cherchi-jacobian.json',jac)
                        rows.append(item);save(f'{name}-cherchi-{cherchi_repetition}-record.json',item)
                rows.append(record);save(f'{name}-cgal-{repetition}-record.json',record)
            except Exception as exc:
                failure=dict(state=name,repetition=repetition,exception_type=type(exc).__name__,exception=str(exc),phases=timer.finish())
                rows.append(failure);save(f'{name}-failure-{repetition}.json',failure)
                break
    save('report.json',dict(rows=rows,optimization_admission=False,performance_admission=False,
        note='CGAL value and Cherchi process/export costs are separate; no paper timing equivalence claimed'))


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    args=p.parse_args();main(args.output,args.probe)
