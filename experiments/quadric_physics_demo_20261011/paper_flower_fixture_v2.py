"""Preregistered full editable middle between two unchanged fixed covers.

The v1 metric-eroded web remains historical evidence. No saved optimized
circle, desired cut location, or post-solve rotation enters this request.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys

from paper_flower_fixture import (build_fixture as build_legacy_fixture,
    boundary_loops, exact_extrusion_moments, polygons)


def build_fixture(source):
    from shapely.geometry import Polygon
    problem, record, geometry = build_legacy_fixture(source)
    geometry['design'] = geometry['original']
    geometry['fixed_middle'] = Polygon()
    layers = [layer for layer in record['layers']
              if layer['id'] != 'fixed-middle-boundary-rim']
    problem.update(design_region_loops=boundary_loops(geometry['design']),
        fixed_extrusion_layers=[{k: v for k, v in layer.items() if k != 'id'}
                               for layer in layers])
    fixed = [F(0)] * 10
    for layer in layers:
        fixed = [a + b for a, b in zip(fixed, exact_extrusion_moments(layer))]
    record.update(schema='paper-flower-method-input-v2', layers=layers,
        case_id='paper-image-derived-flower-open-middle-v2',
        rim_width=0., rim_rule='no fixed middle hole-boundary web; full original middle editable',
        parent_input_policy='v1 retained unchanged; v2 uses original unoptimized profile',
        fixed_exact_rational_moments=[str(x) for x in fixed],
        fixed_model_mass_moments=[float(x) for x in fixed],
        mechanical_strength_status='NOT_VERIFIED',
        input_revision_reason='remove unintended protected web that crossed the optimized cavity')
    record['admission'].update(design_component_count=len(polygons(geometry['design'])),
        design_area=geometry['design'].area, fixed_middle_area=0.,
        design_rim_overlap_area=0., middle_partition_symmetric_difference_area=0.,
        minimum_design_boundary_clearance=0.,
        middle_boundary_policy='DESIGN_MEETS_FIXED_COVERS_AND_ORIGINAL_SIDE_BOUNDARY')
    return problem, record, geometry


def freeze(source, out):
    import trimesh
    from paper_flower_output import constrained_boundary_prism, write_roundtrip_obj
    problem, record, geometry = build_fixture(source)
    layers = record['layers']
    meshes = {}
    layer_geometry = {'bottom': geometry['original'], 'top': geometry['original'],
        'bearing-seat': geometry['seat'], 'outer-race': geometry['outer_race']}
    for layer in layers:
        meshes[layer['id']] = trimesh.util.concatenate([
            constrained_boundary_prism(polygon, layer['z_interval'], record['meters_per_model_unit'])
            for polygon in polygons(layer_geometry[layer['id']])])
    for name, interval in [('original-flower', [-.08, .08]), ('design', [-.06, .06])]:
        meshes[name] = constrained_boundary_prism(geometry['original'], interval,
                                                 record['meters_per_model_unit'])
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    checks = {}
    for name, mesh in meshes.items():
        mesh.apply_scale(record['meters_per_model_unit'])
        path = out / (name + '.obj')
        write_roundtrip_obj(mesh, path)
        loaded = trimesh.load(path, force='mesh', process=False)
        if not loaded.is_watertight or not loaded.is_winding_consistent or loaded.volume <= 0:
            raise ValueError('OBJ_RELOAD_INVALID: ' + name)
        checks[name] = dict(watertight=True, winding_consistent=True, positive_volume=True,
            volume_m3=float(loaded.volume), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    record['obj_checks'] = checks
    manifest = dict(schema='faster-planar-method-fixture-v1', scope='RECONSTRUCTED_METHOD_VALIDATION',
        cases=[dict(id=record['case_id'], mode='FIDGET_SPINNER', asset_class='RECONSTRUCTED',
                    planar_problem=problem)])
    for name, value in [('flower-cases.json', manifest), ('input.json', record)]:
        with (out / name).open('x') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write('\n')
    with (out / 'input-sha256.json').open('x') as stream:
        json.dump({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(out.iterdir()) if p.is_file()}, stream, indent=2)
        stream.write('\n')
    print(json.dumps(record['admission'], indent=2))
    return record


if __name__ == '__main__':
    freeze(Path(sys.argv[1]), Path(sys.argv[2]))
