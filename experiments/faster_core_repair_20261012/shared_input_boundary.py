"""Share the identical derived internal interface; no automatic containment gate."""
import numpy as np
import trimesh

def assemble_fixed_boundary(outer,design):
    for mesh in (outer,design):
        if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume<=0:
            raise ValueError('INVALID_SHARED_INPUT_BOUNDARY')
    return trimesh.Trimesh(np.vstack((outer.vertices,design.vertices)),
        np.vstack((outer.faces,design.faces[:,::-1]+len(outer.vertices))),process=False)
