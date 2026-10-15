"""Bind an admitted handled Heart to Faster's existing production request."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import gmpy2

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze(admission,candidate,out):
    from make_it_stand.spin import SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig
    admission,candidate,out=map(Path,(admission,candidate,out))
    proof=json.loads((admission/'admission.json').read_text())
    if proof.get('geometry_admission_granted') is not True:
        raise ValueError('HEART_WHOLE_INPUT_NOT_ADMITTED')
    for path,expected in proof['source_files'].items():
        if sha(path)!=expected:raise ValueError('GEOMETRY_ADMISSION_SOURCE_CHANGED:'+path)
    geometry={name:json.loads((candidate/(name+'.binary64.json')).read_text()) for name in ('outer','fixed','design')}
    bounds=json.loads((admission/'fixed-independent-bounds.json').read_text())
    fixed=tuple(float(gmpy2.mpq(q)) for q in bounds['independent_derived_moments_exact'])
    if len(fixed)!=10 or fixed[0]<=0:raise ValueError('INVALID_INDEPENDENT_FIXED_MOMENTS')
    requests=[]
    for seed in (0,1):
        request=SpinItFasterRequest(outer_vertices=geometry['outer']['vertices'],outer_faces=geometry['outer']['faces'],
            design_vertices=geometry['design']['vertices'],design_faces=geometry['design']['faces'],
            fixed_shell_moments=fixed,mode='SPINNING_TOP',target_origin=(0.,0.,0.),
            material_density=1.,center_of_mass_weight=1.,spinning_top_objective='INERTIA_AND_LOWERING',
            settings=SpinSolverConfig(deterministic_seed=seed,equation_profile='DUAL_MONOMIAL_RAW_V1',
                maximum_restarts=8,maximum_function_evaluations_per_restart=60,capture_evaluation_payload=True))
        requests.append(asdict(request))
    settings=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
        icosphere_subdivisions=2,maximum_icosphere_subdivisions=2,
        retain_fixed_normal_provenance=True,rational_source_coordinates=True,
        exact_binary64_source_boundary=True,
        cap_radius_policy='GROW_TO_CERTIFIED_V1')
    sources={str(path.resolve()):sha(path) for path in [admission/'admission.json',
        admission/'fixed-independent-bounds.json',candidate/'outer.binary64.json',
        candidate/'design.binary64.json',candidate/'fixed.binary64.json',Path(__file__)]}
    payload=dict(schema='rounded-heart-faster-frozen-request-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',
        requests=requests,evaluator_settings=asdict(settings),input_sources=sources,
        geometry_problem_digest=proof['derived_problem_digest'],
        units=dict(model_length_m=.025,model_mass_kg=.01875,density_kg_m3=1200.,
                   objective_coordinate_contract='INPUT_MODEL_UNITS_TARGET_ORIGIN_ZERO_V1'),
        initialization='REQUEST_SEED_ONLY_NO_EXPLICIT_QUADRICS',
        derivative_admission_granted=False,solver_admission_granted=False,material_admission_granted=False)
    out.mkdir(parents=True,exist_ok=False)
    (out/'request.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    (out/'identity.json').write_text(json.dumps(dict(request_sha256=sha(out/'request.json'),
        schema=payload['schema'],derivative_admission_granted=False),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('admission','candidate','output'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();freeze(args.admission,args.candidate,args.output)
