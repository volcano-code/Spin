"""Unpowered ideal revolute bearing; actual PyBullet states, not animation."""
import numpy as np

def simulate(mass,center,inertia,*,dt,duration,hz,phase=0.):
    import pybullet as p
    from scipy.spatial.transform import Rotation
    c=np.asarray(center,float);I=np.asarray(inertia,float)
    if mass<=0 or dt<=0 or duration<=0 or not np.isfinite(np.r_[mass,c,I.ravel(),dt,duration,hz,phase]).all():
        raise ValueError('INVALID_MASS_PROPERTIES')
    values,U=np.linalg.eigh(I)
    if min(values)<=0 or not np.allclose(I,I.T):raise ValueError('INVALID_INERTIA')
    if np.linalg.det(U)<0:U[:,0]*=-1
    client=p.connect(p.DIRECT)
    try:
        p.setGravity(0,0,-9.81,physicsClientId=client)
        p.setPhysicsEngineParameter(fixedTimeStep=dt,numSolverIterations=150,physicsClientId=client)
        body=p.createMultiBody(baseMass=0,baseCollisionShapeIndex=-1,
            linkMasses=[mass],linkCollisionShapeIndices=[-1],linkVisualShapeIndices=[-1],
            linkPositions=[[0,0,0]],linkOrientations=[[0,0,0,1]],
            linkInertialFramePositions=[c.tolist()],
            linkInertialFrameOrientations=[Rotation.from_matrix(U).as_quat().tolist()],
            linkParentIndices=[0],linkJointTypes=[p.JOINT_REVOLUTE],linkJointAxis=[[0,0,1]],physicsClientId=client)
        p.changeDynamics(body,0,maxJointVelocity=max(1000.,4*abs(2*np.pi*hz)),physicsClientId=client)
        p.changeDynamics(body,0,localInertiaDiagonal=values.tolist(),linearDamping=0,angularDamping=0,
                         jointDamping=0,physicsClientId=client)
        p.setJointMotorControl2(body,0,p.VELOCITY_CONTROL,force=0,physicsClientId=client)
        p.enableJointForceTorqueSensor(body,0,True,physicsClientId=client)
        p.resetJointState(body,0,phase,2*np.pi*hz,physicsClientId=client)
        loaded=p.getDynamicsInfo(body,0,physicsClientId=client)
        loadedI=U@np.diag(loaded[2])@U.T
        np.testing.assert_allclose(loadedI,I,rtol=1e-12,atol=1e-18)
        J=float(I[2,2]+mass*(c[0]**2+c[1]**2));E0=.5*J*(2*np.pi*hz)**2
        rows=[];maxerror=0
        for step in range(round(duration/dt)+1):
            if step:p.stepSimulation(physicsClientId=client)
            angle,w,reaction,torque=p.getJointState(body,0,physicsClientId=client)
            E=.5*J*w*w;maxerror=max(maxerror,abs(E/E0-1))
            rows.append(dict(t=step*dt,angle=angle,omega=w,reaction=list(reaction),motor_torque=torque,energy=E))
        return dict(records=rows,mass=mass,center=c.tolist(),loaded_inertia=loadedI.tolist(),axis_inertia=J,
                    maximum_relative_energy_drift=maxerror,dt=dt,duration=duration,initial_phase_rad=phase,
                    scope='IDEAL_FIXED_REVOLUTE_BEARING_NO_CONTACT_NO_MOTOR_NO_STOPPING_TIME')
    finally:p.disconnect(client)
