"""Read-only consistency audit of the existing flower delivery and replay.

Reuses recorded exact geometry checks; does not rerun optimization or confer
manufacturing qualification. Writes a new audit report, never edits delivery.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
FLOWER = HERE.parent/'quadric_physics_demo_20261011/paper-flower-method-v2'


def require(condition, code):
    if not condition:
        raise ValueError(code)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def check_physics_matrix(summary, properties):
    require(summary.get('actual_integration') is True, 'MISSING_ACTUAL_INTEGRATION')
    scenes = summary['scenes']
    expected = set(itertools.product(('original', 'optimized'), ('ideal', 'flexible'),
                                     range(8), (.0005, .00025)))
    observed = [(row['variant'], row['kind'], row['phase_index'], row['dt']) for row in scenes]
    require(len(observed) == 64 and set(observed) == expected, 'INCOMPLETE_SCENE_MATRIX')
    support = scenes[0]['support']
    require(support['frequency_hz'] == 10. and support['damping_ratio'] == .1, 'SUPPORT_PROTOCOL_CHANGED')
    require(properties['original']['inertia_kg_m2'] != properties['optimized']['inertia_kg_m2'],
            'VARIANT_INERTIAS_IDENTICAL')
    maximum_inertia_error = 0.
    for row in scenes:
        require(row['support'] == support, 'SUPPORT_PARAMETERS_NOT_SHARED')
        require(row['actuator_count'] == 0, 'CONTINUOUS_DRIVE_PRESENT')
        require(abs(row['phase']-row['phase_index']*math.pi/4) <= 1e-14, 'PHASE_PROTOCOL_CHANGED')
        actual = np.asarray(row['loaded_inertia']); expected_inertia = np.asarray(properties[row['variant']]['inertia_kg_m2'])
        require(actual.shape == (3, 3) and np.isfinite(actual).all(), 'INVALID_LOADED_INERTIA')
        error = float(np.max(np.abs(actual-expected_inertia))/np.max(np.abs(expected_inertia)))
        require(error <= 1e-12, 'FULL_MATERIAL_INERTIA_MISMATCH')
        maximum_inertia_error = max(maximum_inertia_error, error)
        for field in ('peak_displacement', 'force_mean', 'energy_drift', 'energy_balance_error'):
            require(math.isfinite(row[field]), 'NONFINITE_PHYSICS_STATISTIC')
        if row['kind'] == 'ideal':
            mass = properties[row['variant']]['mass_kg']
            center = properties[row['variant']]['center_m']
            force = mass*math.hypot(center[0], center[1])*(2*math.pi*35)**2
            require(abs(row['force_mean']-force) <= 1e-12*max(force, 1e-12), 'IDEAL_BEARING_FORCE_MISMATCH')
    return dict(scene_count=len(scenes), full_inertia_bound=True, shared_support_parameters=True,
                maximum_relative_loaded_inertia_error=maximum_inertia_error)


def audit(root):
    root = Path(root).resolve()
    load = lambda name: json.loads((root/name).read_text())
    checked_files = {}
    def check(path, expected):
        path = Path(path).resolve()
        observed = sha(path)
        require(observed == expected, 'FILE_IDENTITY_MISMATCH:'+str(path))
        checked_files[str(path)] = observed
    def manifest_files(base, files):
        base = (root/base).resolve()
        for name, expected in files.items():
            path = (base/name).resolve()
            require(path.is_relative_to(base), 'MANIFEST_PATH_OUTSIDE_DELIVERY')
            check(path, expected)
    first, repeat = load('solve-job-01/output/application.json'), load('solve-job-02/output/application.json')
    require(first['scientific_hash'] == repeat['scientific_hash'], 'OPTIMIZATION_REPEAT_MISMATCH')
    for job in ('solve-job-01', 'solve-job-02', 'saved-independent-job-03', 'delivery-job-03',
                'mujoco-blender-v1/physics-job-01', 'mujoco-blender-v1/render-matched-job-02',
                'mujoco-blender-v1/encode-matched-job-01'):
        resource, config = load(job+'/result.json'), load(job+'/config.json')
        require(resource['status'] == 'PROCESS_EXITED' and resource['returncode'] == 0 and
                resource['unit'] == config['unit'], 'DELIVERY_RESOURCE_JOB_NOT_COMPLETE:'+job)
    require(first['cases'][0]['summary']['strict_terminal_numerics']['all_runs_numerically_qualified'] and
            first['cases'][0]['seed_sequence'] == [0, 1], 'REGISTERED_OPTIMIZATION_EVIDENCE_INCOMPLETE')
    saved = load('saved-independent-job-03/output/saved-independent-check.json')
    require(saved['full_material_gate'] == 'PASSED' and saved['simulation_admission'] is True and
            saved['accepted_steps'] > 0, 'INDEPENDENT_MATERIAL_GATE_NOT_PASSED')
    for path, expected in saved['sources'].items():
        check(path, expected)
    require(all(row['completed'] and row['valid'] for row in saved['exact_region_facts'].values()),
            'EXACT_REGION_EVIDENCE_INVALID')
    require(all(row['completed'] and row['complete_partition'] and row['no_positive_volume_overlap']
                for row in saved['partitions'].values()), 'PARTITION_EVIDENCE_INVALID')
    require(max(max(values) for values in saved['moment_errors'].values()) <= 1e-6,
            'SAVED_MATERIAL_MOMENTS_OUTSIDE_GATE')
    preparation = load('prepare-job-02/output/output-preparation.json')
    require(preparation['derived_geometry_policy'] == 'BOUNDED_SHARED_RATIONAL_ARRANGEMENT_V1',
            'UNREGISTERED_OBJ_REPRESENTATION')
    require(preparation['circle_center'] == saved['circle_center'] and
            preparation['circle_radius'] == saved['circle_radius'], 'OPTIMIZATION_CIRCLE_CHANGED')
    for row in preparation['files'].values():
        check(root/'prepare-job-02/output'/row['file'], row['sha256'])
        require(row['watertight'] and row['winding_consistent'] and row['positive_volume'] and
                row['exact_zero_area_face_count'] == 0, 'OBJ_TOPOLOGY_EVIDENCE_INVALID')
    delivery = load('delivery-job-03/output/delivery-manifest.json')
    check(root/'prepare-job-02/output/output-preparation.json', delivery['material_preparation_sha256'])
    check(root/'saved-independent-job-03/output/saved-independent-check.json', delivery['saved_check_sha256'])
    manifest_files('delivery-job-03/output', delivery['files'])
    manifest_files('delivery-job-03/output', delivery['archives'])
    table = load('table-scene-v1/delivery-manifest.json')
    manifest_files('table-scene-v1', table['files'])
    assembly = load('table-scene-v1/assembly-check.json')
    for name, row in assembly.items():
        check(root/'table-scene-v1'/name, row['sha256'])
        require(row['coordinates_and_faces_unchanged'] and row['watertight_components'], 'ASSEMBLY_CHANGED')
    require({'retained', 'top', 'bottom'}.issubset(assembly['flower-complete.obj']['parts']), 'ASSEMBLY_LAYER_MISSING')
    archive_path = root/'table-scene-v1/flower-complete-obj.zip'
    with zipfile.ZipFile(archive_path) as archive:
        for name in ('flower-complete.obj', 'flower-rotor.obj', 'physical-properties.json'):
            matches = [key for key in archive.namelist() if key == name or key.endswith('/'+name)]
            require(len(matches) == 1, 'ASSEMBLY_ARCHIVE_MEMBER_MISSING_OR_AMBIGUOUS:'+name)
            require(hashlib.sha256(archive.read(matches[0])).hexdigest() == sha(root/'table-scene-v1'/name),
                    'ASSEMBLY_ARCHIVE_DIFFERS_FROM_DELIVERY:'+name)
    # The page provides separate assembly and separated-parts downloads. Each
    # archive is checked against its own delivered files, without requiring
    # redundant copies of all parts in the complete-assembly archive.
    with zipfile.ZipFile(root/'delivery-job-03/output/final-assembly-obj.zip') as archive:
        for name in delivery['files']:
            if not name.startswith('final-assembly/'):
                continue
            require(name in archive.namelist(), 'SEPARATED_PART_ARCHIVE_MEMBER_MISSING:'+name)
            require(hashlib.sha256(archive.read(name)).hexdigest() == delivery['files'][name],
                    'SEPARATED_PART_ARCHIVE_DIFFERS_FROM_DELIVERY:'+name)
    physics_root = root/'mujoco-blender-v1/physics-job-01/output'
    summary = load('mujoco-blender-v1/physics-job-01/output/summary.json')
    properties = load('table-scene-v1/physical-properties.json')
    physics = check_physics_matrix(summary, properties)
    for row in summary['scenes']:
        check(physics_root/row['file'], row['sha256'])
        check(physics_root/Path(row['file']).with_suffix('.xml'), row['xml_sha256'])
    with zipfile.ZipFile(root/'mujoco-blender-v1/mujoco-results.zip') as archive:
        require(hashlib.sha256(archive.read('physics/summary.json')).hexdigest() == sha(physics_root/'summary.json'),
                'PHYSICS_ARCHIVE_SUMMARY_MISMATCH')
        for row in summary['scenes']:
            for name, expected in ((row['file'], row['sha256']),
                                   (str(Path(row['file']).with_suffix('.xml')), row['xml_sha256'])):
                require(hashlib.sha256(archive.read('physics/'+name)).hexdigest() == expected,
                        'PHYSICS_ARCHIVE_SCENE_MISMATCH:'+name)
    render_path = root/'mujoco-blender-v1/render-matched-job-02/output/render.json'
    render = json.loads(render_path.read_text())
    video = load('mujoco-blender-v1/encode-matched-job-01/output/video.json')
    check(render_path, video['render_sha256'])
    check(physics_root/'summary.json', video['physics_sha256'])
    check(root/'table-scene-v1/table-scene.json', video['material_sha256'])
    require(render['physics_sha256'] == video['physics_sha256'] and render['material_sha256'] == video['material_sha256'] and
            render['renders_actual_states'] and render['vibration_amplification'] == 1 and
            render['material_profile'] == video['material_profile'] == 'MATCHED_COMPARISON_V1' and
            render['selected_phase'] == video['selected_phase'] == 0 and render['playback_rates'] == [1, .05],
            'RENDER_TRAJECTORY_OR_PRESENTATION_MISMATCH')
    manifest_files('mujoco-blender-v1/render-matched-job-02/output', render['files'])
    check(root/'mujoco-blender-v1/encode-matched-job-01/output/flower-matched-comparison.mp4', video['video_sha256'])
    return dict(schema='flower-current-delivery-identity-audit-v1', passed=True,
        optimization_repeat_identity_verified=True, stored_exact_material_gate_verified=True,
        complete_assembly_has_three_layers=True, physics=physics,
        physics_download_contains_all_64_trajectories_and_scenes=True,
        matched_video_uses_bound_trajectory=True, checked_file_count=len(checked_files), checked_files=checked_files,
        author_original_cad=False, manufacturing_approved=False,
        scope='Current files and stored exact material/dynamics evidence consistency; no new geometry or calibration claim')


if __name__ == '__main__':
    report = audit(FLOWER)
    with Path(sys.argv[1]).open('x') as stream:
        json.dump(report, stream, sort_keys=True, allow_nan=False)
    print(json.dumps({key: value for key, value in report.items() if key != 'checked_files'}), flush=True)
