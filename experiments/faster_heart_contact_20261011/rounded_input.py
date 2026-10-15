"""Frozen rounded reconstruction and independent first-moment capacity bounds.

No optimized quadric, cavity template or previous solver endpoint is read here.
The displayed mesh is distinct from the authoritative rational construction.
"""
import hashlib
import json
from pathlib import Path
import sys

import gmpy2
import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments/reconstructed_heart_top'))
from heart_top import heart_mesh, integrate as signed_moments_array

SPEC = dict(schema='rounded-heart-reconstruction-input-v1',
    classification='RECONSTRUCTED_METHOD_VALIDATION', author_original=False,
    samples=64, layers=24, depth_factor=1.6, angle_degrees=-35.,
    body_centroid_model=[.1, 0., 1.12], inner_scale=.90,
    grip_radius_model=.11, grip_center_z_model=2.30,
    model_unit_m=.025, density_kg_m3=1200.,
    target_origin_model=[0., 0., 0.], target_axis=[0., 0., 1.],
    minimum_outer_clearance_model=.015,
    objective='FASTER_INERTIA_RATIO_PLUS_CENTER_OF_MASS_HEIGHT',
    center_of_mass_weight=1., weight_coordinate_contract='FASTER_NORMALIZED_REQUEST',
    backend='CGAL', normal_subdivision=2, seeds=[0, 1],
    maximum_attempts_per_seed=8, maximum_evaluations_per_attempt=60,
    raw_residual_tolerance=1e-12, alpha=.5, t_max=.5,
    parameters_source='PREREGISTERED_RECONSTRUCTION_NOT_AUTHOR_PARAMETERS')


def signed_moments(mesh):
    return signed_moments_array(mesh.vertices, mesh.faces)


def rounded_components():
    """Actual connected union is a later exact Boolean, not a mesh concatenation."""
    vertices, faces = heart_mesh(SPEC['samples'], SPEC['layers'])
    vertices[:, 1] *= SPEC['depth_factor']
    vertices -= signed_moments_array(vertices, faces)[1:4] / signed_moments_array(vertices, faces)[0]
    angle = np.deg2rad(SPEC['angle_degrees'])
    rotation = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0],
                         [-np.sin(angle), 0, np.cos(angle)]])
    center = np.asarray(SPEC['body_centroid_model'])
    vertices = vertices @ rotation.T + center
    body = trimesh.Trimesh(vertices, faces, process=False)
    inner = trimesh.Trimesh(center + SPEC['inner_scale'] * (vertices - center),
                            faces.copy(), process=False)
    profile = [(0., -.16)]
    for angle in np.linspace(-np.pi/2, 0, 7)[1:]:
        profile.append((.025*np.cos(angle), -.135+.025*np.sin(angle)))
    radius, height = SPEC['grip_radius_model'], SPEC['grip_center_z_model']
    start = -np.arccos(.065/radius)
    profile.extend([(.065, -.03), (.065, height+radius*np.sin(start))])
    for angle in np.linspace(start, np.pi/2, 17)[1:-1]:
        profile.append((radius*np.cos(angle), height+radius*np.sin(angle)))
    profile.append((0., height+radius))
    pin = trimesh.creation.revolve(np.array(profile), sections=32)
    if pin.volume < 0:
        pin.invert()
    return body, inner, pin


def _q(value):
    if isinstance(value, (float, np.floating)):
        return gmpy2.mpq(float(value))
    return gmpy2.mpq(value)


def _det(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) +
            a[1]*(b[2]*c[0]-b[0]*c[2]) + a[2]*(b[0]*c[1]-b[1]*c[0]))


def _half_surface_first(vertices, faces, axis, sign):
    """Exact first moment on sign*x>=0; the omitted plane cap contributes zero.

    The plane passes through the target origin. This computes a necessary bound
    for arbitrary subsets of a valid represented solid, not quadric feasibility.
    """
    total = gmpy2.mpq(0)
    for face in faces:
        poly = [vertices[int(i)] for i in face]
        clipped = []
        for a, b in zip(poly, poly[1:]+poly[:1]):
            inside, other = sign*a[axis] >= 0, sign*b[axis] >= 0
            if inside:
                clipped.append(a)
            if inside != other:
                t = -a[axis]/(b[axis]-a[axis])
                clipped.append(tuple(a[k]+t*(b[k]-a[k]) for k in range(3)))
        for j in range(1, len(clipped)-1):
            a, b, c = clipped[0], clipped[j], clipped[j+1]
            total += _det(a, b, c)*(a[axis]+b[axis]+c[axis])/24
    return total


def exact_moments(vertices, faces):
    """Independent oriented-tetrahedron formula, with ten unrounded rationals."""
    v = [tuple(_q(x) for x in row) for row in vertices]
    s = [gmpy2.mpq(0) for _ in range(10)]
    for face in faces:
        a, b, c = [v[int(i)] for i in face]
        volume = _det(a, b, c)/6
        sums = [a[k]+b[k]+c[k] for k in range(3)]
        s[0] += volume
        for k in range(3):
            s[1+k] += volume*sums[k]/4
            s[4+k] += volume*(sums[k]**2+a[k]**2+b[k]**2+c[k]**2)/20
        for i, (j, k) in enumerate(((1, 2), (0, 2), (0, 1))):
            s[7+i] += volume*(sums[j]*sums[k]+a[j]*a[k]+b[j]*b[k]+c[j]*c[k])/20
    return s


def capacity_bounds(fixed_vertices, fixed_faces, design_vertices, design_faces):
    fixed = exact_moments(fixed_vertices, fixed_faces)
    design = [tuple(_q(x) for x in row) for row in design_vertices]
    result = dict(schema='rounded-heart-first-moment-capacity-v1',
        scope='NECESSARY_CONTINUOUS_SUBSET_BOUND_GIVEN_VALID_INPUT_PARTITION',
        target_origin=[0., 0., 0.], target_axis=[0., 0., 1.], density_model=1.,
        feasible_witness=False, sealed_certificate='NOT_ISSUED')
    excluded = False
    for axis, label in enumerate(('x', 'y')):
        low = _half_surface_first(design, design_faces, axis, -1)
        high = _half_surface_first(design, design_faces, axis, 1)
        total_low, total_high = fixed[1+axis]+low, fixed[1+axis]+high
        rejected = total_low > 0 or total_high < 0
        excluded |= rejected
        result[label] = dict(fixed_first_exact=str(fixed[1+axis]),
            design_min_exact=str(low), design_max_exact=str(high),
            total_min_exact=str(total_low), total_max_exact=str(total_high),
            total_min=float(total_low), total_max=float(total_high), excludes_zero=rejected)
    result['status'] = ('EXCLUDED_BY_FIRST_MOMENT_BOUND' if excluded else
                        'UNKNOWN_NOT_EXCLUDED_BY_FIRST_MOMENT_BOUNDS')
    return result


def write_mesh(path, mesh):
    """Binary64 display/source serialization with full roundtrip precision."""
    with Path(path).open('x') as out:
        for row in mesh.vertices:
            out.write('v ' + ' '.join(format(float(x), '.17g') for x in row) + '\n')
        for row in mesh.faces:
            out.write('f ' + ' '.join(str(int(x)+1) for x in row) + '\n')


def publish(path, value):
    with Path(path).open('x') as out:
        json.dump(value, out, indent=2, allow_nan=False)


def preview(folder, body, inner, pin, *, filename='input-preview.png'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(12, 8), facecolor='#f5f7fb')
    for i, (title, elev, azim) in enumerate((('Front / fixed tilted pose', 0, -90),
        ('Side / rounded depth', 0, 0), ('Perspective / connected-input preview', 18, -65),
        ('Design-domain preview / not an optimized cavity', 8, -90))):
        ax = fig.add_subplot(2, 2, i+1, projection='3d')
        meshes = [(body, '#c55694', .22 if i == 3 else 1.), (pin, '#d9c69c', 1.)]
        if i == 3:
            meshes.append((inner, '#e9bc57', .8))
        for mesh, color, alpha in meshes:
            triangles = mesh.vertices[mesh.faces]*25
            from matplotlib.colors import to_rgb
            normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
            normals /= np.maximum(np.linalg.norm(normals, axis=1)[:, None], 1e-30)
            light = np.asarray([-.3, -.8, .6]); light /= np.linalg.norm(light)
            shade = .48 + .52*np.maximum(normals @ light, 0)
            colors = np.asarray(to_rgb(color))[None, :] * shade[:, None]
            collection = Poly3DCollection(triangles, facecolors=colors, edgecolors='none', alpha=alpha)
            ax.add_collection3d(collection)
        ax.set(xlim=(-30, 30), ylim=(-22, 22), zlim=(-5, 62),
               xlabel='x [mm]', ylabel='y [mm]', zlabel='z [mm]', title=title)
        ax.set_box_aspect((60, 44, 67)); ax.view_init(elev, azim)
        ax.plot([0, 0], [0, 0], [-4, 60], color='#247c9a', linewidth=1)
    fig.suptitle('Rounded asymmetric Heart v2 — RECONSTRUCTED INPUT ONLY', fontsize=15)
    destination = folder/filename
    if destination.exists():
        raise FileExistsError(destination)
    fig.tight_layout(); fig.savefig(destination, dpi=160); plt.close(fig)


def prepare(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    body, inner, pin = rounded_components()
    publish(folder/'spec.json', SPEC)
    for name, mesh in (('heart-body.obj', body), ('inner-source.obj', inner), ('pin.obj', pin)):
        write_mesh(folder/name, mesh)
    preview(folder, body, inner, pin)
    sources = {name: hashlib.sha256((folder/name).read_bytes()).hexdigest()
               for name in ('spec.json', 'heart-body.obj', 'inner-source.obj', 'pin.obj')}
    publish(folder/'preparation.json', dict(status='SOURCE_GEOMETRY_PREPARED_NOT_ADMITTED',
        sources=sources, implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        body_geometric_moments=signed_moments(body).tolist(),
        body_vertices=len(body.vertices), body_faces=len(body.faces),
        optimization='NOT_RUN', simulation='NOT_RUN'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.prepare)
