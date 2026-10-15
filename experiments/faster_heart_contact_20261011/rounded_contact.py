"""Actual PyBullet contact for one preregistered scene.

The caller must admit saved material and collision proxies before invoking this
backend. This module measures motion and contacts; it does not grant an optimized
material seal or claim a calibrated real stopping time.
"""
import hashlib
from pathlib import Path
import sys
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'experiments/heart_stem_contact'))
from bullet_contact_common_frame import create_body
from contact_matrix import scenarios


def simulate_scene(outer_file, collision_files, props, outer_vertices_si, scene, protocol):
    import pybullet as p
    if scene not in scenarios(protocol):
        raise ValueError('UNREGISTERED_CONTACT_SCENE')
    mass, center, inertia, diagonal, axes = props
    if (not np.isfinite(np.r_[mass, center, inertia.ravel(), diagonal, axes.ravel()]).all() or
            mass <= 0 or np.min(diagonal) <= 0 or np.linalg.det(axes) < 0 or
            not np.allclose(axes.T@axes, np.eye(3), atol=1e-12) or
            not np.allclose(axes@np.diag(diagonal)@axes.T, inertia, rtol=1e-12, atol=1e-18)):
        raise ValueError('INVALID_FULL_INERTIA_REPRESENTATION')
    files = [Path(outer_file), *(Path(x) for x in collision_files)]
    hashes = {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in files}
    displayed = trimesh.load(outer_file, force='mesh', process=False)
    if not np.array_equal(displayed.vertices, np.asarray(outer_vertices_si)):
        raise ValueError('OUTER_COLLISION_DISPLAY_COORDINATE_MISMATCH')
    dt = scene['dt_seconds']
    steps = int(round(protocol['duration_seconds']/dt))
    client = p.connect(p.DIRECT)
    if client < 0:
        raise RuntimeError('PYBULLET_CONNECTION_FAILED')
    try:
        p.setRealTimeSimulation(0, physicsClientId=client)
        p.setPhysicsEngineParameter(fixedTimeStep=dt, numSolverIterations=150,
            deterministicOverlappingPairs=1, physicsClientId=client)
        p.setGravity(0, 0, -protocol['gravity_m_s2'], physicsClientId=client)
        floor_shape = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.4, .4, .02], physicsClientId=client)
        floor = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=floor_shape,
            basePosition=[0, 0, -.02], physicsClientId=client)
        p.changeDynamics(floor, -1, lateralFriction=protocol['lateral_friction'],
            restitution=protocol['restitution'], physicsClientId=client)
        azimuth = np.deg2rad(scene['azimuth_degrees'])
        rotation = Rotation.from_rotvec(np.deg2rad(protocol['initial_tilt_degrees'])*
            np.array([np.cos(azimuth), np.sin(azimuth), 0.])).as_matrix()
        clearance = .0002
        origin = np.array([0., 0., clearance-(np.asarray(outer_vertices_si)@rotation.T)[:, 2].min()])
        body = create_body(p, client, outer_file, collision_files, props, rotation, origin,
                           protocol['nominal_spin_hz'])
        if p.getNumConstraints(physicsClientId=client) != 0:
            raise ValueError('UNEXPECTED_WORLD_CONSTRAINT')
        info = p.getDynamicsInfo(body, -1, physicsClientId=client)
        position, quat = p.getBasePositionAndOrientation(body, physicsClientId=client)
        principal_world = np.array(p.getMatrixFromQuaternion(quat)).reshape(3, 3)
        if (info[0] != mass or not np.allclose(principal_world@np.diag(info[2])@principal_world.T,
                rotation@inertia@rotation.T, rtol=1e-12, atol=1e-18) or
                not np.allclose(position, origin+rotation@center, atol=1e-12, rtol=0)):
            raise ValueError('PYBULLET_DID_NOT_LOAD_ACTUAL_MATERIAL_INERTIA')
        closest = p.getClosestPoints(body, floor, .01, physicsClientId=client)
        gap = min((x[8] for x in closest), default=None)
        if gap is None or abs(gap-clearance) > .0001:
            raise ValueError(f'COLLISION_MODEL_BLOCKER: measured initial gap {gap}, expected {clearance}')
        margins = [dict(link=link, margin_m=p.getDynamicsInfo(body, link, physicsClientId=client)[-1])
                   for link in range(-1, p.getNumJoints(body, physicsClientId=client))]
        records = []; maximum_tilt = 0.; maximum_penetration = 0.
        stride = max(1, int(round(1/(200*dt))))
        for step in range(steps+1):
            if step:
                p.stepSimulation(physicsClientId=client)
            position, quat = p.getBasePositionAndOrientation(body, physicsClientId=client)
            linear, angular = p.getBaseVelocity(body, physicsClientId=client)
            if not np.isfinite(np.r_[position, quat, linear, angular]).all():
                raise ValueError('NONFINITE_CONTACT_STATE')
            world_principal = np.array(p.getMatrixFromQuaternion(quat)).reshape(3, 3)
            world_model = world_principal@axes.T
            axis = world_model[:, 2]
            tilt = float(np.rad2deg(np.arccos(np.clip(axis[2], -1, 1))))
            maximum_tilt = max(maximum_tilt, tilt)
            contacts = p.getContactPoints(body, floor, physicsClientId=client)
            maximum_penetration = max(maximum_penetration, max((-x[8] for x in contacts), default=0.))
            if step % stride == 0 or step == steps:
                local_angular = world_principal.T@np.asarray(angular)
                energy = (.5*mass*np.dot(linear, linear)+.5*np.dot(diagonal*local_angular, local_angular)+
                          mass*protocol['gravity_m_s2']*position[2])
                records.append(dict(time_s=step*dt, center_position_m=list(position),
                    principal_quaternion=list(quat), model_quaternion=Rotation.from_matrix(world_model).as_quat().tolist(),
                    model_origin_m=(np.asarray(position)-world_model@center).tolist(),
                    linear_velocity_m_s=list(linear), angular_velocity_rad_s=list(angular),
                    tilt_degrees=tilt, axial_spin_hz=float(np.dot(angular, axis)/(2*np.pi)),
                    energy_J=float(energy), contacts=[dict(link=x[3], on_body_m=list(x[5]),
                        on_floor_m=list(x[6]), normal=list(x[7]), distance_m=x[8], normal_force_N=x[9],
                        lateral_force_1_N=x[10], lateral_direction_1=list(x[11]),
                        lateral_force_2_N=x[12], lateral_direction_2=list(x[13])) for x in contacts]))
        if any(hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest for path, digest in hashes.items()):
            raise ValueError('COLLISION_SOURCE_CHANGED_DURING_SIMULATION')
        return dict(schema='rounded-heart-actual-contact-scene-v1', scene=scene,
            status='PHYSICAL_TRAJECTORY_RECORDED_NOT_CALIBRATED', duration_seconds=protocol['duration_seconds'],
            actual_steps=steps, engine_api_version=p.getAPIVersion(), input_sha256=hashes,
            loaded_mass_kg=float(mass), loaded_center_m=center.tolist(), loaded_inertia_kg_m2=inertia.tolist(),
            runtime_collision_margins=margins, initial_contact_gap_m=gap,
            maximum_tilt_degrees=maximum_tilt, maximum_penetration_m=maximum_penetration,
            final_tilt_degrees=records[-1]['tilt_degrees'], final_spin_hz=records[-1]['axial_spin_hz'],
            records=records, world_constraint_count=0, persistent_drive_torque=False)
    finally:
        p.disconnect(client)
