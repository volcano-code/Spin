"""Freeze the APPROVED layout/recipe for split jobs; no production map calls.

The geometry-owned native context is captured at attempt_start, then explicitly
aborted before initial material evaluation. This is context/identity preparation,
not optimizer success, a known-solution initializer, or derivative evidence.
"""
import argparse
from pathlib import Path

from production_point_job import require_point_resource_envelope,save
from make_it_stand.verification.quadric_point_jobs import (
    _read,file_sha256,freeze_point_manifest,verify_point_sources)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--registration',type=Path,required=True)
    parser.add_argument('--registration-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--job-envelope',type=Path)
    args=parser.parse_args()
    if file_sha256(args.registration)!=args.registration_sha256:
        raise ValueError('original registration identity mismatch')
    require_point_resource_envelope(args.job_envelope)
    prior=_read(args.registration)
    if prior.get('schema')!='quadric-production-five-point-preregistration-v1' or \
       prior.get('production_queries')!=81 or prior.get('scientific_tolerance')!=1e-6 or \
       prior.get('source_interval_width')!='1/'+str(2**384) or prior.get('optimizer_invocations')!=0:
        raise ValueError('not the approved original 81-point registration')
    # Preserve scientific input/registration assets; old code/native hashes are
    # historical and intentionally replaced by explicit new dependency binding.
    assets={}
    for path,digest in prior['sources'].items():
        candidate=Path(path)
        if candidate.suffix=='.json':
            if file_sha256(candidate)!=digest:
                raise ValueError('frozen scientific source changed')
            assets['historical_asset_'+str(len(assets))]=str(candidate.resolve())
    if not assets: raise ValueError('no authenticated original scientific assets')
    from make_it_stand import spin
    from make_it_stand.verification import _spin2014_authority as checker
    from make_it_stand.verification import mass_moments,rational_text,quadric_point_jobs,surface_derivatives,quadric_branch
    native=spin._native_module()
    request=dict(prior['request']);request['settings']=spin.SpinSolverConfig(**request['settings'])
    public_request=spin.SpinItFasterRequest(**request)
    _,actual_request=spin._native_request(public_request)
    config=spin.SpinEvaluatorConfig(**prior['evaluator'])
    if config.icosphere_subdivisions!=2 or config.maximum_icosphere_subdivisions!=2 or \
       not config.rational_source_coordinates:
        raise ValueError('original recipe does not provide frozen level2 authority sources')
    evaluator=spin._native_evaluator_settings(native,config)
    fixed_density=public_request.fixed_shell_density or public_request.material_density
    script=Path(__file__).resolve()
    paths=dict(assets,input=str(args.registration.resolve()),native=native.__file__,
        runner=str(script.with_name('production_point_job.py')),
        aggregator=str(script.with_name('aggregate_production_points.py')),
        freezer=str(script),checker=checker.__file__,mapping_interface=spin.__file__,
        mass_moments=mass_moments.__file__,rational_text=rational_text.__file__,
        manifest_checker=quadric_point_jobs.__file__,surface_derivatives=surface_derivatives.__file__,
        quadric_branch=quadric_branch.__file__)
    fingerprints={name:file_sha256(path) for name,path in paths.items()}
    seen=[]
    class ContextOnlyStop(RuntimeError): pass
    def stop(event):
        if event['type']!='attempt_start' or event['evaluation_index']!=0:
            raise ValueError('context preparation reached an unexpected evaluation event')
        seen.append(event)
        raise ContextOnlyStop('context only: production mapping not evaluated')
    objective=native.Spin2014QuadricObjectiveSettings()
    if public_request.mode=='YO_YO': objective.objective=native.Spin2014Objective.YOYO
    try:
        native.solve_spin2014_authority_quadric_hollowing(actual_request,objective,fixed_density,
            evaluator,evaluation_budget_per_attempt=1,event_sink=stop)
    except ContextOnlyStop:
        pass
    if len(seen)!=1: raise ValueError('missing geometry-owned pre-evaluation context identity')
    mapping=dict(coefficient_policy='COEFFICIENT_MPFR1024_V1',precision_bits=1024,normal_level=2,
        request=prior['request'],evaluator=prior['evaluator'],fixed_material_density=fixed_density,
        expected_request_identity=seen[0]['request_identity'],expected_context_identity=seen[0]['context_identity'],
        source_paths=paths,context_capture_scope='ATTEMPT_START_ABORTED_BEFORE_INITIAL_MATERIAL',
        context_capture_actual_material_evaluations=0,original_registration_sha256=args.registration_sha256)
    manifest=freeze_point_manifest(prior['all_raw_trials'][0],prior['original_base_raw_steps'],
                                   sources=fingerprints,mapping_settings=mapping)
    if [p['raw_coefficients'] for p in manifest['points']]!=prior['all_raw_trials'] or manifest['rows']!=prior['rows']:
        raise ValueError('split layout altered approved point positions/rows')
    verify_point_sources(manifest,paths)
    save(args.output,manifest)
    print(manifest['identity'],flush=True)
    return 0


if __name__=='__main__': raise SystemExit(main())
