"""Bind inspected actual SI input to the existing normalized Faster interface."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
from bound_faster_input import normalize_bound_geometry
from make_it_stand.spin import SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
source=HERE/'noodles-balanced-shared-export-job-02/output'
manifest_path=source/'manifest.json';manifest=json.loads(manifest_path.read_text())
original_path=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5/manifest.json'
original=json.loads(original_path.read_text());spec=original['spec']
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
if manifest['model_manifest_sha256']!=sha(original_path) or not manifest['derived_input_partition_passed']:
    raise ValueError('WHOLE_INPUT_NOT_BOUND_OR_PARTITIONED')
if spec['target_axis']!=[0.,0.,1.] or spec['objective']!='FASTER_INERTIA_RATIO_PLUS_CENTER_OF_MASS_HEIGHT':
    raise ValueError('FROZEN_OBJECTIVE_OR_AXIS_CHANGED')
meshes={};sources={str(manifest_path):sha(manifest_path),str(original_path):sha(original_path)}
for name in ('outer','design','fixed'):
    entry=manifest['files'][name];path=source/entry['name']
    if sha(path)!=entry['sha256'] or entry['length_unit']!='m' or entry['role']!='INPUT_NOT_OPTIMIZED':
        raise ValueError('SAVED_DERIVED_INPUT_CHANGED:'+name)
    sources[str(path)]=sha(path);meshes[name]=trimesh.load(path,force='mesh',process=False)
target=np.asarray(spec['target_origin_model'])*spec['model_unit_m']
bound=normalize_bound_geometry(meshes['outer'],meshes['design'],meshes['fixed'],target,spec['density_kg_m3'])
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
requests=[]
for seed in spec['seeds']:
    request=SpinItFasterRequest(outer_vertices=bound['meshes']['outer'].vertices.tolist(),
        outer_faces=bound['meshes']['outer'].faces.tolist(),
        design_vertices=bound['meshes']['design'].vertices.tolist(),design_faces=bound['meshes']['design'].faces.tolist(),
        fixed_shell_moments=bound['fixed_shell_moments'],target_origin=bound['target_origin'],
        mode='SPINNING_TOP',material_density=1.,original_units_per_normalized_unit=bound['normalization']['L'],
        center_of_mass_weight=spec['center_of_mass_weight'],spinning_top_objective='INERTIA_AND_LOWERING',
        settings=SpinSolverConfig(deterministic_seed=seed,equation_profile='DUAL_MONOMIAL_RAW_V1',
            maximum_iterations=60,maximum_restarts=spec['maximum_attempts_per_seed'],
            maximum_function_evaluations_per_restart=spec['maximum_evaluations_per_attempt'],
            equation_tolerance=spec['raw_residual_tolerance'],capture_evaluation_payload=True))
    requests.append(asdict(request))
settings=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',
    icosphere_subdivisions=2,maximum_icosphere_subdivisions=2,retain_fixed_normal_provenance=True,
    rational_source_coordinates=False,exact_binary64_self_intersection=True,
    exact_binary64_source_boundary=True,exact_binary64_material_boundary=True,
    paper_projective_hyperboloid=True,cap_radius_policy='GROW_TO_CERTIFIED_V1')
fixed=bound['meshes']['fixed']
with (output/'fixed.binary64.json').open('x') as stream:
    json.dump(dict(vertices=fixed.vertices.tolist(),faces=fixed.faces.tolist()),stream,allow_nan=False)
sources[str(output/'fixed.binary64.json')]=sha(output/'fixed.binary64.json')
sources[str(Path(__file__))]=sha(Path(__file__));sources[str(HERE/'bound_faster_input.py')]=sha(HERE/'bound_faster_input.py')
with (output/'request.json').open('x') as stream:
    json.dump(dict(schema='balanced-noodles-faster-frozen-request-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',
        requests=requests,evaluator_settings=asdict(settings),input_sources=sources,
        normalization=bound['normalization'],backend='CGAL',coefficient_policy='BINARY64_CORE_V1',
        units=dict(model_length_m=bound['normalization']['L'],model_mass_kg=bound['normalization']['mass_unit_kg'],
            density_kg_m3=spec['density_kg_m3'],objective_coordinate_contract=spec['weight_coordinate_contract']),
        alpha=.5,t_max=.5,initialization='REQUEST_AND_SEED_ONLY_NO_EXPLICIT_QUADRICS',
        derivative_admission=False,solver_admission=False,optimization_invocations=0),stream,allow_nan=False)
print('Frozen normalized request only; production map admission pending',flush=True)
