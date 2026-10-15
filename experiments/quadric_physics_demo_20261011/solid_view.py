"""3D display derivatives of saved planar curves; never changes solver material."""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'render_deps'))

def extrude(loops,low,high):
    import mapbox_earcut,trimesh
    from matplotlib.path import Path as P
    rings=[np.asarray(x,float) for x in loops]
    def area(q):return np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1))/2
    shells=[q for q in rings if area(q)>0];holes=[q for q in rings if area(q)<0]
    parts=[]
    for outer in shells:
        group=[outer]+[h for h in holes if P(outer).contains_point(h[0])]
        xy=np.vstack(group);n=len(xy);ends=np.cumsum([len(x) for x in group]).astype(np.uint32)
        top=mapbox_earcut.triangulate_float64(xy,ends).reshape(-1,3)
        vertices=np.vstack([np.column_stack([xy,np.full(n,low)]),np.column_stack([xy,np.full(n,high)])])
        faces=[*top[:,::-1],*(top+n)];offset=0
        for ring in group:
            for j in range(len(ring)):
                a=offset+j;b=offset+(j+1)%len(ring);faces.extend([[a,b,b+n],[a,b+n,a+n]])
            offset+=len(ring)
        mesh=trimesh.Trimesh(vertices,faces,process=False)
        if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume<=0:raise ValueError('INVALID_DISPLAY_EXTRUSION')
        parts.append(mesh)
    if not parts:raise ValueError('EMPTY_DISPLAY_EXTRUSION')
    return trimesh.util.concatenate(parts)

def entity_loops(entity):
    if 'loops' in entity:return entity['loops']
    loops=[]
    for cycle in entity['boundary_cycles']:
        points=[]
        for s in cycle:
            if s['type']=='LINE_SEGMENT':points.append(s['start'])
            else:
                t=s['start_angle']+np.linspace(0,s['sweep'],max(4,int(abs(s['sweep'])*160)),endpoint=False)
                points.extend(np.array(s['center'])+s['radius']*np.column_stack([np.cos(t),np.sin(t)]))
        loops.append(points)
    return loops
