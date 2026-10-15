"""Frozen reconstructed Noodles source parts; exact union is a separate gate.

No optimization output or yellow paper material is used to generate this input.
Source OBJ vertices represent a finite triangulation, not exact smooth tubes.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import trimesh

SPEC = {
    'schema': 'reconstructed-noodles-input-v1',
    'classification': 'RECONSTRUCTED_METHOD_VALIDATION',
    'author_original': False,
    'model_unit_m': .025, 'density_kg_m3': 1200.,
    'target_origin_model': [0., 0., -.16], 'target_axis': [0., 0., 1.],
    'ellipsoid_center_model': [0., 0., 1.12],
    'ellipsoid_radii_model': [.62, .53, .64], 'inner_scale': .90,
    'ellipsoid_subdivisions': 2,
    'tube_angles_degrees': [0., 74., 147., 219., 292.],
    'tube_reaches_model': [.90, .84, .93, .88, .86],
    'tube_radius_model': .075, 'tube_samples': 41, 'tube_sections': 16,
    'tube_start_radius_model': .22, 'tube_start_z_model': .96,
    'tube_end_z_model': 1.31, 'tube_twist_radians': .75,
    'tube_angle_bend_radians': .12,
    'hub_radius_model': .13, 'hub_interval_model': [.48, 1.42],
    'shaft_radius_model': .065, 'tip_radius_model': .025,
    'tip_bottom_model': -.16, 'grip_radius_model': .115,
    'grip_center_z_model': 2.12, 'shaft_sections': 32,
    'objective': 'FASTER_INERTIA_RATIO_PLUS_CENTER_OF_MASS_HEIGHT',
    'center_of_mass_weight': 1.,
    'weight_coordinate_contract': 'FASTER_NORMALIZED_REQUEST',
    'backend': 'CGAL', 'normal_subdivision': 2, 'seeds': [0, 1],
    'maximum_attempts_per_seed': 8, 'maximum_evaluations_per_attempt': 60,
    'raw_residual_tolerance': 1e-12, 'alpha': .5, 't_max': .5,
    'parameters_source': 'PREREGISTERED_CONSTRUCTION_NOT_AUTHOR_PARAMETERS',
}


def spec_identity():
    return hashlib.sha256(json.dumps(SPEC, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def centerline(index):
    u = np.linspace(0., 1., SPEC['tube_samples'])
    angle = np.deg2rad(SPEC['tube_angles_degrees'][index])
    theta = angle + SPEC['tube_twist_radians']*u + SPEC['tube_angle_bend_radians']*np.sin(np.pi*u)
    radius = SPEC['tube_start_radius_model'] + SPEC['tube_reaches_model'][index]*np.sin(np.pi*u)
    z = (SPEC['tube_start_z_model'] +
         (SPEC['tube_end_z_model']-SPEC['tube_start_z_model'])*u)
    return np.column_stack((radius*np.cos(theta), radius*np.sin(theta), z))


def tube_mesh(points, radius, sections):
    """Parallel-transport circular rings and oriented end caps."""
    points = np.asarray(points, dtype=float)
    tangent = np.gradient(points, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    seed = np.eye(3)[np.argmin(np.abs(tangent[0]))]
    normal = seed - seed.dot(tangent[0])*tangent[0]
    normal /= np.linalg.norm(normal)
    normals = [normal]
    for a, b in zip(tangent[:-1], tangent[1:]):
        axis = np.cross(a, b); sine = np.linalg.norm(axis); cosine = np.clip(a.dot(b), -1., 1.)
        if sine > 1e-14:
            axis /= sine
            normal = normal*cosine + np.cross(axis, normal)*sine + axis*axis.dot(normal)*(1-cosine)
        normal -= normal.dot(b)*b; normal /= np.linalg.norm(normal)
        normals.append(normal.copy())
    normals = np.asarray(normals); binormals = np.cross(tangent, normals)
    angles = np.arange(sections)*2*np.pi/sections
    rings = points[:, None, :] + radius*(np.cos(angles)[None, :, None]*normals[:, None, :] +
                                         np.sin(angles)[None, :, None]*binormals[:, None, :])
    vertices = np.vstack((rings.reshape(-1, 3), points[[0, -1]]))
    faces = []
    for i in range(len(points)-1):
        for j in range(sections):
            a = i*sections+j; b = i*sections+(j+1)%sections
            c = (i+1)*sections+j; d = (i+1)*sections+(j+1)%sections
            faces.extend(((a, b, d), (a, d, c)))
    cap0, cap1 = len(vertices)-2, len(vertices)-1
    for j in range(sections):
        k = (j+1)%sections
        faces.extend(((cap0, k, j), (cap1, (len(points)-1)*sections+j, (len(points)-1)*sections+k)))
    mesh = trimesh.Trimesh(vertices, np.asarray(faces), process=False)
    if mesh.volume < 0: mesh.invert()
    return mesh


def components():
    sphere = trimesh.creation.icosphere(subdivisions=SPEC['ellipsoid_subdivisions'])
    sphere.vertices = sphere.vertices*np.asarray(SPEC['ellipsoid_radii_model']) + SPEC['ellipsoid_center_model']
    inner = sphere.copy()
    center = np.asarray(SPEC['ellipsoid_center_model'])
    inner.vertices = center + SPEC['inner_scale']*(inner.vertices-center)
    hub = trimesh.creation.cylinder(radius=SPEC['hub_radius_model'],
        height=np.diff(SPEC['hub_interval_model'])[0], sections=32)
    hub.apply_translation([0, 0, np.mean(SPEC['hub_interval_model'])])
    profile = [(0., SPEC['tip_bottom_model'])]
    for a in np.linspace(-np.pi/2, 0, 7)[1:]:
        profile.append((SPEC['tip_radius_model']*np.cos(a), -.135+SPEC['tip_radius_model']*np.sin(a)))
    r = SPEC['grip_radius_model']; z = SPEC['grip_center_z_model']
    start = -np.arccos(SPEC['shaft_radius_model']/r)
    profile += [(SPEC['shaft_radius_model'], -.03), (SPEC['shaft_radius_model'], z+r*np.sin(start))]
    profile += [(r*np.cos(a), z+r*np.sin(a)) for a in np.linspace(start, np.pi/2, 17)[1:-1]]
    profile.append((0., z+r))
    pin = trimesh.creation.revolve(np.asarray(profile), sections=SPEC['shaft_sections'])
    if pin.volume < 0: pin.invert()
    parts = {'outer-ellipsoid': sphere, 'inner-ellipsoid-candidate': inner, 'hub-fixed': hub, 'pin-fixed': pin}
    for i in range(len(SPEC['tube_angles_degrees'])):
        parts[f'tube-{i}-fixed'] = tube_mesh(centerline(i), SPEC['tube_radius_model'], SPEC['tube_sections'])
    return parts


def mesh_audit(mesh):
    triangles = mesh.vertices[mesh.faces]
    areas = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]), axis=1)/2
    return dict(vertices=len(mesh.vertices), triangles=len(mesh.faces),
        closed=bool(mesh.is_watertight), winding_consistent=bool(mesh.is_winding_consistent),
        components=len(mesh.split()), signed_volume_model3=float(mesh.volume),
        minimum_triangle_area_model2=float(areas.min()),
        self_intersection='NOT_CERTIFIED_PENDING_EXACT_SOURCE_ADMISSION')


def write_obj(path, mesh):
    with Path(path).open('x') as stream:
        for v in mesh.vertices: stream.write('v '+' '.join(format(float(x), '.17g') for x in v)+'\n')
        for face in mesh.faces: stream.write('f '+' '.join(str(int(i)+1) for i in face)+'\n')


def prepare(directory):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=False)
    parts = components(); entries = {}
    for name, mesh in parts.items():
        audit = mesh_audit(mesh)
        if not audit['closed'] or not audit['winding_consistent'] or audit['components'] != 1 or mesh.volume <= 0:
            raise ValueError(f'INVALID_CLOSED_SOURCE_PART:{name}')
        write_obj(directory/f'{name}-model.obj', mesh)
        si = mesh.copy(); si.vertices *= SPEC['model_unit_m']; write_obj(directory/f'{name}-si.obj', si)
        entries[name] = dict(audit=audit, model_obj=f'{name}-model.obj', si_obj=f'{name}-si.obj',
            sha256={suffix: hashlib.sha256((directory/f'{name}-{suffix}.obj').read_bytes()).hexdigest() for suffix in ('model', 'si')})
    fixed = ['hub-fixed', 'pin-fixed']+[f'tube-{i}-fixed' for i in range(5)]
    manifest = dict(schema='noodles-source-parts-manifest-v1', spec=SPEC, spec_sha256=spec_identity(), parts=entries,
        construction=dict(fixed_parts=fixed, outer_operation='UNION(outer-ellipsoid, ALL fixed_parts)',
            design_operation='DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))',
            fixed_material_operation='DIFFERENCE(outer, design)'),
        source_connectedness='ALL_PARTS_HAVE_POSITIVE_OVERLAP_WITH_CENTRAL_ELLIPSOID; EXACT_UNION_PENDING',
        status='SOURCE_PARTS_PREPARED_PENDING_EXACT_UNION_PARTITION_AND_CAPACITY',
        optimization_invocations=0, simulation_invocations=0, manufacturing_approved=False,
        creator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (directory/'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False))
    preview(directory, parts)
    return manifest


def preview(directory, parts):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(12, 7))
    for index, (title, elev, azim) in enumerate((('Front source parts', 0, -90), ('Side source parts', 0, 0),
            ('Perspective / exact union pending', 25, -55), ('Editable candidate / fixed parts excluded later', 15, -65))):
        ax = fig.add_subplot(2, 2, index+1, projection='3d')
        for name, mesh in parts.items():
            if name == 'inner-ellipsoid-candidate' and index != 3: continue
            color = '#e6b44c' if name == 'inner-ellipsoid-candidate' else '#488a9b' if 'tube' in name else '#b95886' if 'ellipsoid' in name else '#c9cbd0'
            alpha = .16 if index == 3 and name == 'outer-ellipsoid' else .8
            ax.add_collection3d(Poly3DCollection(mesh.vertices[mesh.faces], facecolors=color, edgecolors='none', alpha=alpha))
        ax.set(xlim=(-1.3, 1.3), ylim=(-1.3, 1.3), zlim=(-.2, 2.3), title=title)
        ax.set_box_aspect((2.6, 2.6, 2.5)); ax.view_init(elev, azim)
        ax.set_xlabel('x [model]'); ax.set_ylabel('y [model]'); ax.set_zlabel('z [model]')
    fig.suptitle('Reconstructed Noodles v1 — source preview; no optimization or physics yet')
    fig.tight_layout(); fig.savefig(directory/'input-preview.png', dpi=150); plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(); parser.add_argument('--prepare', type=Path, required=True)
    prepare(parser.parse_args().prepare)
