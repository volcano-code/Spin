"""Conservative distance-domain candidate; independent exact admission required.

For an exact 1-Lipschitz distance field, interpolation differs from the true
field by at most the cell diagonal. A second cell diagonal protects the
interior of each marching-cubes triangle. Numerical mesh queries are NOT exact
distance oracles, so this construction deliberately never grants admission.
"""
import numpy as np


def extract_distance_domain(distance, low, high, *, spacing, clearance,
                            chunk_size=8192, event=None):
    import trimesh
    from skimage.measure import marching_cubes
    low, high = np.asarray(low, float), np.asarray(high, float)
    if (low.shape != (3,) or high.shape != (3,) or
        not np.isfinite([*low, *high, spacing, clearance]).all() or
        spacing <= 0 or clearance <= 0 or np.any(high <= low) or chunk_size <= 0):
        raise ValueError('INVALID_DISTANCE_DOMAIN_PARAMETERS')
    axes = [np.linspace(a, b, int(np.ceil((b-a)/spacing))+1)
            for a, b in zip(low, high)]
    steps = np.array([axis[1]-axis[0] for axis in axes])
    shape = tuple(len(axis) for axis in axes)
    count = int(np.prod(shape))
    field = np.empty(count, float)
    for start in range(0, count, chunk_size):
        indices = np.unravel_index(np.arange(start, min(start+chunk_size, count)), shape)
        points = np.column_stack([axis[index] for axis, index in zip(axes, indices)])
        values = np.asarray(distance(points), float)
        if values.shape != (len(points),) or not np.isfinite(values).all():
            raise ValueError('NONFINITE_DISTANCE_OR_WRONG_SHAPE')
        field[start:start+len(points)] = values
        if event is not None:
            event(dict(stage='DISTANCE_GRID', completed=start+len(points), total=count))
    field = field.reshape(shape)
    guard = 2*float(np.linalg.norm(steps))
    level = clearance+guard
    if field.max() <= level:
        raise ValueError('EMPTY_DISTANCE_DOMAIN_CANDIDATE')
    if any(np.any(np.take(field, index, axis=axis) >= level)
           for axis in range(3) for index in (0, -1)):
        raise ValueError('DOMAIN_INTERSECTS_GRID_BOUNDARY')
    # skimage internally represents field positions in float32; include this
    # only as a recorded candidate error, never as an exact certificate.
    vertices, faces, _, _ = marching_cubes(field, level=level, spacing=steps,
                                          allow_degenerate=False)
    mesh = trimesh.Trimesh(vertices=vertices+low, faces=faces, process=False)
    if mesh.volume < 0:
        mesh.invert()
    if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume <= 0:
        raise ValueError('INVALID_DISTANCE_DOMAIN_BOUNDARY')
    return mesh, dict(schema='conservative-distance-domain-candidate-v1',
        requested_clearance=clearance, actual_grid_spacing=steps.tolist(),
        grid_shape=list(shape), interpolation_and_triangle_guard=guard,
        extraction_level=level, distance_assumption='EXACT_1_LIPSCHITZ_ONLY_FOR_ANALYTIC_BOUND',
        numerical_mesh_distance_is_certificate=False, independent_exact_admission=False,
        solver_admission=False, simulation_admission=False)
