"""Reconstructed Utah Teapot experiment helpers (not author-original input)."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'experiments/heart_stem_contact'))
from stem_geometry import boolean, write_mesh

STRATEGIES = {'inertia': 'TOP_INERTIA_ONLY', 'lowering': 'TOP_LOWERING_ONLY',
              'mass_lowering': 'TOP_MASS_AND_LOWERING'}
SOURCE_SHA256 = '1b5396fedd74b577e32cef41146582c2f2e1a050d5b4915193c0ac1ad4187ed4'


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')


def file_hashes(paths):
    return {str(Path(p).resolve()): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def check_files(files):
    for name, expected in files.items():
        if not Path(name).is_file() or hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ValueError('STALE_RESULT: '+name)


def canonical_indices(mesh):
    """Weld identical binary64 positions; discard only repeated-index faces.

    No epsilon, averaging, area threshold, or displacement is permitted.
    CGAL rational intersections can coincide when converted to binary64.
    This does not repair non-collapsed collinear faces or invalid topology.
    """
    v, ix = np.unique(np.asarray(mesh.vertices), axis=0, return_inverse=True)
    f = ix[np.asarray(mesh.faces)]
    keep = (f[:,0]!=f[:,1]) & (f[:,0]!=f[:,2]) & (f[:,1]!=f[:,2])
    return trimesh.Trimesh(v, f[keep], process=False), dict(
        identical_vertices_merged=len(mesh.vertices)-len(v),collapsed_faces_removed=int((~keep).sum()),
        max_vertex_displacement=0.,tolerance=0.)


def cap_planar_boundaries(mesh):
    """Explicit convex planar caps; preserve every source vertex and triangle.

    This is reconstruction, not a silent repair. Reject other loop geometries.
    Final exact self-intersection/region checks remain independently required.
    """
    v = np.asarray(mesh.vertices).copy(); f = np.asarray(mesh.faces).copy()
    if not np.isfinite(v).all() or not mesh.is_winding_consistent:
        raise ValueError('INVALID_SOURCE_WINDING')
    counts = np.bincount(mesh.edges_unique_inverse)
    if (counts > 2).any(): raise ValueError('NONMANIFOLD_SOURCE')
    boundary = mesh.edges[counts[mesh.edges_unique_inverse] == 1]
    outgoing = {}; incoming = {}
    for a, b in boundary:
        if int(a) in outgoing or int(b) in incoming: raise ValueError('BRANCHED_BOUNDARY')
        outgoing[int(a)] = int(b); incoming[int(b)] = int(a)
    if set(outgoing) != set(incoming): raise ValueError('OPEN_BOUNDARY_CHAIN')
    verts = v.tolist(); faces = f.tolist(); records = []
    scale = max(float(np.linalg.norm(np.ptp(v, axis=0))), 1e-30)
    while outgoing:
        start = min(outgoing); loop = []; cur = start
        while True:
            loop.append(cur); cur = outgoing.pop(cur)
            if cur == start: break
            if cur not in outgoing: raise ValueError('BROKEN_BOUNDARY_LOOP')
        pts = v[loop]; center = pts.mean(axis=0)
        _, _, vh = np.linalg.svd(pts-center, full_matrices=False)
        residual = float(np.max(np.abs((pts-center)@vh[-1])))
        if residual > 1e-10*scale: raise ValueError('NONPLANAR_BOUNDARY')
        turns = np.cross(np.roll(pts,-1,axis=0)-pts, np.roll(pts,-2,axis=0)-np.roll(pts,-1,axis=0))@vh[-1]
        if (turns > 1e-12*scale**2).any() and (turns < -1e-12*scale**2).any():
            raise ValueError('NONCONVEX_BOUNDARY')
        k = len(verts); verts.append(center.tolist())
        faces.extend([[loop[(i+1)%len(loop)], a, k] for i, a in enumerate(loop)])
        records.append(dict(vertices=loop, center=center.tolist(), planarity_error=residual))
    result = trimesh.Trimesh(verts, faces, process=False)
    if not result.is_watertight or not result.is_winding_consistent or result.volume <= 0:
        raise ValueError('CLOSURE_FAILED')
    return result, dict(capped_loops=len(records), added_triangles=len(faces)-len(f), loops=records)


def shaft_mesh():
    angles = np.linspace(-np.pi/2, 0, 7)
    profile = [(0,0)]+[(.025*np.cos(t), .025+.025*np.sin(t)) for t in angles[1:]]
    profile += [(.06,.13), (.06,1.71), (0,1.71)]
    mesh = trimesh.creation.revolve(np.array(profile), sections=32)
    if mesh.volume < 0: mesh.invert()
    return mesh


def strategy_settings(name):
    from make_it_stand.verification.source_spin_material import SourceSolverSettings, AuditedSlqpSettings
    if name not in STRATEGIES: raise ValueError('unknown strategy')
    return SourceSolverSettings(objective=STRATEGIES[name], initial_level=3,
        maximum_optimization_level=5, maximum_adaptive_passes=12, maximum_leaf_count=12000,
        maximum_output_cuboids=12000, initial_fill=.5, initial_phi=.2,
        inertia_weight=1., center_weight=1., laplacian_weight=.001,
        slqp=AuditedSlqpSettings(maximum_iterations=240))
