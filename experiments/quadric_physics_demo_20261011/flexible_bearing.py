"""Unpowered rotor on two translating spring/damper support stages.

This engineering abstraction measures support vibration; it is not a calibrated
human-hand or rolling-contact bearing model. The same frozen support constants
are used for original and optimized rotors.
"""
import math
import numpy as np


def support_protocol(original_rotor_mass, moving_clamp_mass):
    if not math.isfinite(original_rotor_mass + moving_clamp_mass) or min(original_rotor_mass, moving_clamp_mass) <= 0:
        raise ValueError('INVALID_SUPPORT_MASS')
    frequency, damping = 10., .1
    # Each carriage carries half the actual moving clamp. The x stage moves
    # both carriages; the y stage moves only the second carriage and the rotor.
    stage_mass = moving_clamp_mass / 2
    effective = np.array([moving_clamp_mass + original_rotor_mass,
                          stage_mass + original_rotor_mass])
    omega = 2*math.pi*frequency
    return dict(schema='flexible-bearing-support-v1', frequency_hz=frequency,
                damping_ratio=damping, original_rotor_mass_kg=original_rotor_mass,
                moving_clamp_mass_kg=moving_clamp_mass, carriage_mass_kg=stage_mass,
                effective_reference_mass_kg=effective.tolist(),
                stiffness_N_m=(effective*omega**2).tolist(),
                damping_N_s_m=(2*damping*effective*omega).tolist(),
                scope='ENGINEERING_LINEAR_SUPPORT_NOT_CALIBRATED_HAND')


def simulate(mass, center, inertia, *, protocol, dt, duration, hz, phase):
    import pybullet as p
    from scipy.spatial.transform import Rotation
    c, I = np.asarray(center, float), np.asarray(inertia, float)
    if c.shape != (3,) or I.shape != (3, 3) or not np.isfinite(np.r_[mass, c, I.ravel(), dt, duration, hz, phase]).all():
        raise ValueError('INVALID_MASS_PROPERTIES_OR_PROTOCOL')
    if mass <= 0 or dt <= 0 or duration <= 0 or not np.allclose(I, I.T, atol=1e-18, rtol=1e-12):
        raise ValueError('INVALID_MASS_PROPERTIES_OR_PROTOCOL')
    values, axes = np.linalg.eigh(I)
    if min(values) <= 0:
        raise ValueError('INVALID_INERTIA')
    if np.linalg.det(axes) < 0:
        axes[:, 0] *= -1
    expected = support_protocol(protocol['original_rotor_mass_kg'], protocol['moving_clamp_mass_kg'])
    if expected != protocol:
        raise ValueError('SUPPORT_PROTOCOL_IDENTITY_MISMATCH')
    k, d = np.array(protocol['stiffness_N_m']), np.array(protocol['damping_N_s_m'])
    stage_mass = protocol['carriage_mass_kg']
    client = p.connect(p.DIRECT)
    try:
        p.setGravity(0, 0, -9.81, physicsClientId=client)
        p.setPhysicsEngineParameter(fixedTimeStep=dt, numSolverIterations=150,
                                   deterministicOverlappingPairs=1, physicsClientId=client)
        body = p.createMultiBody(baseMass=0, baseCollisionShapeIndex=-1,
            linkMasses=[stage_mass, stage_mass, mass], linkCollisionShapeIndices=[-1]*3,
            linkVisualShapeIndices=[-1]*3, linkPositions=[[0, 0, 0]]*3,
            linkOrientations=[[0, 0, 0, 1]]*3,
            linkInertialFramePositions=[[0, 0, 0], [0, 0, 0], c.tolist()],
            linkInertialFrameOrientations=[[0, 0, 0, 1], [0, 0, 0, 1], Rotation.from_matrix(axes).as_quat().tolist()],
            linkParentIndices=[0, 1, 2],
            linkJointTypes=[p.JOINT_PRISMATIC, p.JOINT_PRISMATIC, p.JOINT_REVOLUTE],
            linkJointAxis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], physicsClientId=client)
        for link in range(3):
            p.changeDynamics(body, link, maxJointVelocity=max(1000., 4*abs(2*math.pi*hz)), physicsClientId=client)
            p.changeDynamics(body, link, linearDamping=0, angularDamping=0, jointDamping=0,
                             localInertiaDiagonal=(values.tolist() if link == 2 else [stage_mass*1e-6]*3),
                             physicsClientId=client)
            p.setJointMotorControl2(body, link, p.VELOCITY_CONTROL, force=0, physicsClientId=client)
        p.enableJointForceTorqueSensor(body, 2, True, physicsClientId=client)
        p.resetJointState(body, 2, phase, 2*math.pi*hz, physicsClientId=client)
        loaded = p.getDynamicsInfo(body, 2, physicsClientId=client)
        loaded_I = axes @ np.diag(loaded[2]) @ axes.T
        np.testing.assert_allclose(loaded_I, I, rtol=1e-12, atol=1e-18)
        records, damping_loss, peak, peak_error = [], 0., 0., 0.
        initial_energy = None
        for step in range(round(duration/dt)+1):
            if step:
                q = np.array([p.getJointState(body, j, physicsClientId=client)[0] for j in (0, 1)])
                v = np.array([p.getJointState(body, j, physicsClientId=client)[1] for j in (0, 1)])
                damping_loss += float(np.dot(d, v*v))*dt
                for j, force in enumerate(-k*q-d*v):
                    p.setJointMotorControl2(body, j, p.TORQUE_CONTROL, force=float(force), physicsClientId=client)
                p.stepSimulation(physicsClientId=client)
            states = [p.getJointState(body, j, physicsClientId=client) for j in range(3)]
            q, v = np.array([s[0] for s in states[:2]]), np.array([s[1] for s in states[:2]])
            angle, omega, reaction, torque = states[2]
            energy = .5*float(np.dot(k, q*q))
            for link, link_mass in enumerate((stage_mass, stage_mass, mass)):
                ls = p.getLinkState(body, link, computeLinkVelocity=1, physicsClientId=client)
                linear, angular = np.array(ls[6]), np.array(ls[7])
                orientation = Rotation.from_quat(ls[1]).as_matrix()
                diagonal = values if link == 2 else np.full(3, stage_mass*1e-6)
                local_w = orientation.T @ angular
                energy += .5*link_mass*float(linear@linear) + .5*float(np.dot(diagonal, local_w*local_w))
            if initial_energy is None:
                initial_energy = energy
            error = abs(energy+damping_loss-initial_energy)/max(initial_energy, np.finfo(float).tiny)
            peak_error = max(peak_error, error)
            peak = max(peak, float(np.linalg.norm(q)))
            if not np.isfinite(np.r_[q, v, angle, omega, energy, reaction]).all():
                raise ValueError('NONFINITE_SUPPORT_STATE')
            records.append(dict(t=step*dt, support_xy_m=q.tolist(), support_velocity_m_s=v.tolist(),
                angle=angle, omega=omega, reaction=list(reaction), motor_torque=torque,
                energy_J=energy, damping_loss_J=damping_loss, relative_energy_balance_error=error))
        return dict(schema='flexible-bearing-actual-trajectory-v1', support_protocol=protocol,
            scope='UNPOWERED_REVOLUTE_ROTOR_ON_TWO_LINEAR_SUPPORT_STAGES',
            mass_kg=mass, center_m=c.tolist(), loaded_inertia_kg_m2=loaded_I.tolist(),
            dt_seconds=dt, duration_seconds=duration, initial_hz=hz, initial_phase_rad=phase,
            rotor_motor_enabled=False, maximum_support_displacement_m=peak,
            maximum_relative_energy_balance_error=peak_error, records=records)
    finally:
        p.disconnect(client)
