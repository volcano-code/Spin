"""Run the existing Faster core on the frozen handled Heart, without deep audit.

This is a separately identified binary64 geometry profile. It does not upgrade
the deferred rational derivative batch or grant independent material acceptance.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time


def core_settings(original):
    if original.get('derivative_method')!='DISCRETE_FULL_RANK_QUADRIC':
        raise ValueError('core runner requires actual DISCRETE_FULL_RANK_QUADRIC')
    configured=dict(original)
    configured.update(rational_source_coordinates=False,compact_rational_source_coordinates=False,
        exact_binary64_material_boundary=False,exact_binary64_self_intersection=True)
    return configured


def paper_core_settings(original):
    """Explicit paper source plus exact material geometry, with motion gates intact."""
    configured=core_settings(original)
    configured.update(paper_projective_hyperboloid=True,exact_binary64_material_boundary=True)
    return configured


def main(source,output,seed,geometry_profile='legacy'):
    import numpy as np
    import trimesh
    from make_it_stand.spin import (SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig,
        evaluate_spin_it_faster_fullrank_quadric)
    from rounded_solver_journal import run_durable_fullrank,save,json_value
    source,output=Path(source),Path(output)
    output.mkdir(parents=True,exist_ok=False)
    frozen=json.loads(source.read_text())
    for path,expected in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:
            raise ValueError('FROZEN_INPUT_CHANGED:'+path)
    selected=next(row for row in frozen['requests'] if row['settings']['deterministic_seed']==seed)
    payload=dict(selected);payload['settings']=SpinSolverConfig(**selected['settings'])
    request=SpinItFasterRequest(**payload)
    if geometry_profile not in ('legacy','paper-exact'):
        raise ValueError('UNKNOWN_CORE_GEOMETRY_PROFILE')
    profile='FASTER_PAPER_CORE_EXACT_BOUNDARY_V2' if geometry_profile=='paper-exact' else 'FASTER_CORE_BINARY64_DISCRETE_V1'
    configured=paper_core_settings(frozen['evaluator_settings']) if geometry_profile=='paper-exact' else core_settings(frozen['evaluator_settings'])
    settings=SpinEvaluatorConfig(**configured)
    save(output/'core-policy.json',dict(profile=profile,
        request=selected,evaluator_settings=asdict(settings),source_request_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        units=frozen['units'],derivative_batch='DEFERRED_BY_USER',
        independent_material_audit='DEFERRED_BY_USER',material_admission=False))
    started=time.monotonic()
    result=run_durable_fullrank(request,settings,output/'run')
    accepted=sum(row.get('accepted_geodesic_step') is not None and row['accepted_geodesic_step']>0
                 for row in result['trace'])
    summary=dict(profile=profile,status=result['status'],
        seed=seed,solver_seconds=time.monotonic()-started,function_evaluations=result['function_evaluations'],
        accepted_updates=accepted,raw_equation_residual=result.get('raw_equation_residual'),
        physical_relative_residual=result.get('physical_relative_residual'),
        independent_material_audit='DEFERRED_BY_USER',material_admission=False,
        physics_admission=False,units=frozen['units'])
    if result['status']=='CONVERGED' and accepted>0:
        final=evaluate_spin_it_faster_fullrank_quadric(request,result['quadric_coefficients'],evaluator_settings=settings)
        save(output/'final-evaluation.json',asdict(final))
        material=trimesh.Trimesh(final.material_vertices,final.material_faces,process=False)
        quick=dict(watertight=bool(material.is_watertight),winding_consistent=bool(material.is_winding_consistent),
            positive_volume=bool(material.volume>0),positive_triangle_areas=bool(np.all(material.area_faces>0)),
            full_self_intersection_check=False)
        summary['basic_material_checks']=quick
        summary['total_raw_moments']=list(final.total_raw_moments)
        summary['quadric_coefficients']=json_value(result['quadric_coefficients'])
        scale=frozen['units']['model_length_m']
        material.apply_scale(scale)
        material.export(output/'optimized-inner-preview.obj')
        before=trimesh.Trimesh(request.outer_vertices,request.outer_faces,process=False)
        before.apply_scale(scale);before.export(output/'before.obj')
        fixed_path=next(Path(path) for path in frozen['input_sources'] if Path(path).name=='fixed.binary64.json')
        fixed_packet=json.loads(fixed_path.read_text())
        fixed=trimesh.Trimesh(fixed_packet['vertices'],fixed_packet['faces'],process=False)
        fixed.apply_scale(scale);fixed.export(output/'fixed.obj')
        trimesh.util.concatenate([fixed,material]).export(output/'assembly-preview.obj')
        summary['preview_role']='ASSEMBLY_COMPONENTS_NOT_FUSED_OR_INDEPENDENTLY_AUDITED'
    save(output/'summary.json',summary)
    print(json.dumps(json_value(summary)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--request',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--seed',type=int,choices=(0,1),required=True)
    parser.add_argument('--geometry-profile',choices=('legacy','paper-exact'),default='legacy')
    args=parser.parse_args();main(args.request,args.output,args.seed,args.geometry_profile)
