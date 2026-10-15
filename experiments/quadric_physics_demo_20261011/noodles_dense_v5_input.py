"""Gentle nonfolding sweeps on the unchanged dense v3 layout and envelope.

Short tubes receive proportionally smaller waves: wave curvature scales as
amplitude / length squared. The large oscillating lateral arch is replaced
by a single parabolic path through the internal connection. This is a
preregistered input repair, not a response to an optimized material shape.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import numpy as np
import noodles_dense_v3_input as base
from noodles_input import mesh_audit, write_obj

SPEC = deepcopy(base.SPEC)
SPEC.update(schema='dense-noodles-reconstruction-v5',
    parent_source_spec_sha256=base.identity(),
    tube_surface_arch_y_model=0.,
    sweep_rule='PARABOLIC_LATERAL_PATH_WITH_LENGTH_SQUARED_SCALED_WAVES_V1',
    curvature_radius_screen_limit=.45,
    visual_change_purpose='Repair witnessed tube body self-intersection while keeping dense horizontal layout, tube radius and full envelope')


def tube_layout():
    return base.tube_layout()


def dense_centerline(index):
    layout = tube_layout()
    row = layout[index]
    u = np.linspace(0., 1., SPEC['tube_samples'])
    s = 2*u - 1
    scale = (row['half_length_model'] / max(r['half_length_model'] for r in layout))**2
    phase = SPEC['tube_wave_phase_increment']*index
    x = s*row['half_length_model']
    y = row['y_model']*s**2 + SPEC['tube_wave_y_model']*scale*np.sin(4*np.pi*u)*np.sin(np.pi*u)
    z = row['z_model'] + SPEC['tube_wave_z_model']*scale*np.sin(3*np.pi*u+phase)*np.sin(np.pi*u)
    return np.column_stack((x,y,z))


def identity():
    return hashlib.sha256(json.dumps(dict(spec=SPEC,tube_layout=tube_layout()),
        sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def components():
    parts = base.components()
    for index in range(len(tube_layout())):
        parts[f'tube-{index:03d}-fixed'] = base.rounded_tube(
            dense_centerline(index), SPEC['tube_radius_model'], SPEC['tube_sections'])
    return parts


def prepare(directory):
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=False)
    try:
        parts=components()
        entries={}
        for name,mesh in parts.items():
            audit=mesh_audit(mesh)
            if not audit['closed'] or not audit['winding_consistent'] or audit['components']!=1 or mesh.volume<=0 or audit['minimum_triangle_area_model2']<=0:
                raise ValueError('INVALID_SOURCE_PART: '+name)
            write_obj(directory/(name+'-model.obj'),mesh)
            si=mesh.copy();si.vertices*=SPEC['model_unit_m']
            write_obj(directory/(name+'-si.obj'),si)
            entries[name]=dict(audit=audit,model_obj=name+'-model.obj',si_obj=name+'-si.obj',
                sha256={s:hashlib.sha256((directory/(name+'-'+s+'.obj')).read_bytes()).hexdigest() for s in ('model','si')})
        record=dict(schema='dense-noodles-source-parts-v5',spec=SPEC,spec_sha256=identity(),
            tube_layout=tube_layout(),parts=entries,
            reference_image_sha256=hashlib.sha256(base.REFERENCE.read_bytes()).hexdigest() if base.REFERENCE.exists() else None,
            construction=dict(fixed_parts=['hub-fixed','pin-fixed']+[n for n in parts if n.startswith('tube-')],
                outer_operation='UNION(outer-ellipsoid, ALL fixed_parts)',
                design_operation='DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',
                fixed_material_operation='DIFFERENCE(outer, design)'),
            connections='EXACT_WHOLE_UNION_PENDING',
            preview_role='ACTUAL_INITIAL_MATERIAL_NOT_OPTIMIZED',
            status='SOURCE_PARTS_PREPARED_PENDING_EXACT_ADMISSION',
            optimization_invocations=0,simulation_invocations=0,
            creator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            tube_mesh_helper_sha256=hashlib.sha256(Path(base.__file__).with_name('noodles_input.py').read_bytes()).hexdigest(),
            inherited_geometry_generator_sha256=hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest())
        with (directory/'manifest.json').open('x') as stream:
            json.dump(record,stream,indent=2,allow_nan=False)
        base.preview(directory,parts)
        return record
    except Exception as error:
        with (directory/'failure.json').open('x') as stream:
            json.dump(dict(status='SOURCE_PREPARATION_FAILED',exception=type(error).__name__,
                           message=str(error),optimization_invocations=0),stream,indent=2)
        raise


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--prepare',type=Path,required=True)
    prepare(parser.parse_args().prepare)
