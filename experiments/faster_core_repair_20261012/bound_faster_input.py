"""The existing bbox normalization convention, fed only saved actual geometry."""
from pathlib import Path
import sys
import numpy as np
import trimesh
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'faster_heart_contact_20261011'))
from rounded_input import exact_moments

def normalize_bound_geometry(outer,design,fixed,target_origin,density):
    density=float(density);target=np.asarray(target_origin,dtype=float)
    if not np.isfinite(density) or density<=0 or target.shape!=(3,) or not np.isfinite(target).all():
        raise ValueError('INVALID_FASTER_INPUT_UNITS_OR_TARGET')
    low,high=outer.bounds;length=float(max(high-low));center=low+(high-low)*.5
    if not np.isfinite(length) or length<=0:raise ValueError('INVALID_FASTER_NORMALIZATION')
    meshes={name:trimesh.Trimesh((mesh.vertices-center)/length,mesh.faces.copy(),process=False)
            for name,mesh in [('outer',outer),('design',design),('fixed',fixed)]}
    raw=exact_moments(meshes['fixed'].vertices,meshes['fixed'].faces)
    return dict(meshes=meshes,fixed_shell_moments=tuple(float(x) for x in raw),
        target_origin=((target-center)/length).tolist(),material_density=1.,
        normalization=dict(origin=center.tolist(),L=length,R=density,mass_unit_kg=density*length**3,
            frame='INPUT_NORMALIZED_V1',authority='IMPLEMENTATION_DEFAULT'),solver_admission=False)
