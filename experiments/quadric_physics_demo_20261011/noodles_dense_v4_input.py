"""New v4 input: actual smaller core, same frozen v3 dense fixed tube bundle.

No rendering hide or optimizer endpoint is used. The core size is a separately
registered reconstruction choice; v3's oversized-core preview remains intact.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import numpy as np
from noodles_input import mesh_audit,write_obj
import noodles_dense_v3_input as base

SPEC=deepcopy(base.SPEC)
SPEC.update(schema='dense-noodles-reconstruction-v4',
    parent_source_spec_sha256=base.identity(),
    tube_envelope_radii_model=list(base.SPEC['ellipsoid_radii_model']),
    core_ellipsoid_radii_model=[.62,.48,.46],ellipsoid_radii_model=[.62,.48,.46],
    visual_change_purpose='Actual smaller core; visible initial contour formed by fixed horizontal tubes, not hidden by renderer',
    core_change_source='PREREGISTERED_INPUT_GEOMETRY_PURPOSE_NOT_OPTIMIZATION_RESULT',
    actual_initial_material='UNION(actual core ellipsoid, all actual fixed tube/hub/pin solids)')


def tube_layout():return base.tube_layout()


def identity():
    return hashlib.sha256(json.dumps(dict(spec=SPEC,tube_layout=tube_layout()),sort_keys=True,
        separators=(',',':'),allow_nan=False).encode()).hexdigest()


def components():
    parts=base.components();center=np.asarray(SPEC['ellipsoid_center_model'])
    scale=np.asarray(SPEC['core_ellipsoid_radii_model'])/np.asarray(base.SPEC['ellipsoid_radii_model'])
    for name in ('outer-ellipsoid','inner-ellipsoid-candidate'):
        parts[name].vertices=center+(parts[name].vertices-center)*scale
    return parts


def preview(directory,parts):
    # The same renderer consumes ALL the actual v4 part OBJ files, including
    # the smaller core. Only the header text is changed; no mesh is hidden.
    base.preview(directory,parts)
    from PIL import Image,ImageDraw,ImageFont
    canvas=Image.open(directory/'input-preview.png').convert('RGB');draw=ImageDraw.Draw(canvas)
    draw.rectangle((0,0,1599,48),fill='#f4f5f7')
    font=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',23)
    draw.text((20,12),'面条 v4｜实际缩小核心＋79根横向圆帽细管＋握柄与底尖；未优化',font=font,fill='#243249')
    canvas.save(directory/'input-preview.png')


def prepare(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    try:
        parts=components();entries={}
        for name,mesh in parts.items():
            audit=mesh_audit(mesh)
            if not audit['closed'] or not audit['winding_consistent'] or audit['components']!=1 or mesh.volume<=0 or audit['minimum_triangle_area_model2']<=0:
                raise ValueError('INVALID_SOURCE_PART: '+name)
            write_obj(directory/(name+'-model.obj'),mesh)
            si=mesh.copy();si.vertices*=SPEC['model_unit_m'];write_obj(directory/(name+'-si.obj'),si)
            entries[name]=dict(audit=audit,model_obj=name+'-model.obj',si_obj=name+'-si.obj',
                sha256={s:hashlib.sha256((directory/(name+'-'+s+'.obj')).read_bytes()).hexdigest() for s in ('model','si')})
        manifest=dict(schema='dense-noodles-source-parts-v4',spec=SPEC,spec_sha256=identity(),tube_layout=tube_layout(),
            reference_image_sha256=hashlib.sha256(base.REFERENCE.read_bytes()).hexdigest() if base.REFERENCE.exists() else None,
            parts=entries,construction=dict(fixed_parts=['hub-fixed','pin-fixed']+[n for n in parts if n.startswith('tube-')],
                outer_operation='UNION(outer-ellipsoid, ALL fixed_parts)',
                design_operation='DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',
                fixed_material_operation='DIFFERENCE(outer, design)'),
            connections='SAME_V3_INTERNAL_TENDRIL_HUB_CONNECTIONS; EXACT_WHOLE_UNION_PENDING',
            preview_role='ALL_ACTUAL_INITIAL_MATERIAL_COMPONENTS_NOT_OPTIMIZED_BOUNDARY',
            parent_v3_ellipsoid_radii_model=list(base.SPEC['ellipsoid_radii_model']),
            current_core_radii_model=SPEC['core_ellipsoid_radii_model'],
            hidden_geometry_assumptions=SPEC['hidden_geometry_assumptions'],
            status='SOURCE_PARTS_PREPARED_PENDING_EXACT_ADMISSION',optimization_invocations=0,simulation_invocations=0,
            creator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            tube_mesh_helper_sha256=hashlib.sha256(Path(__file__).with_name('noodles_input.py').read_bytes()).hexdigest(),
            inherited_geometry_generator_sha256=hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest())
        (directory/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
        preview(directory,parts)
        return manifest
    except Exception as exc:
        (directory/'failure.json').write_text(json.dumps(dict(status='SOURCE_PREPARATION_FAILED',
            type=type(exc).__name__,message=str(exc),optimization_invocations=0,simulation_invocations=0),indent=2)+'\n')
        raise


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',type=Path,required=True)
    prepare(parser.parse_args().prepare)
