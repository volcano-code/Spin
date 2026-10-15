"""Free-contact PyBullet experiment; never uses a fixed joint or maintained spin."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/reconstructed_heart_top'))
from heart_top import integrate,properties
from scipy.spatial.transform import Rotation


def material_properties(vertices,faces,density):
    m,c,Ip=properties(integrate(vertices,faces)*density)
    Ic=Ip-m*((c@c)*np.eye(3)-np.outer(c,c))
    d,U=np.linalg.eigh(Ic)
    if not np.isfinite(d).all() or d.min()<=0: raise ValueError('invalid central inertia')
    if np.linalg.det(U)<0: U[:,0]*=-1
    return m,c,Ic,d,U


def create_body(p,client,outer_file,collision_files,props,R,origin,hz,color=(.85,.57,.12,1)):
    import trimesh
    m,c,I,d,U=props
    local_rotation=Rotation.from_matrix(U.T).as_quat().tolist()
    local_position=(-U.T@c).tolist()
    n=len(collision_files)
    if not n: raise ValueError('missing convex collision parts')
    # The first convex part is the tip-bearing shaft. This Bullet build ignores
    # link collisionMargin updates; only the base is reliably zero-margin.
    # Remaining heart links retain measured padding (reported, not hidden).
    # Fixed links assemble ONE free rigid object; none attach to the ground.
    collisions=[]
    for file in collision_files:
        mesh=trimesh.load(file,force='mesh',process=False)
        vertices=(np.asarray(mesh.vertices)-c)@U
        collisions.append(p.createCollisionShape(p.GEOM_MESH,vertices=vertices.tolist(),physicsClientId=client))
    visual=p.createVisualShape(p.GEOM_MESH,fileName=str(outer_file),meshScale=[1,1,1],
        visualFramePosition=local_position,visualFrameOrientation=local_rotation,
        rgbaColor=color,physicsClientId=client)
    k=n-1
    body=p.createMultiBody(baseMass=m,baseCollisionShapeIndex=collisions[0],baseVisualShapeIndex=visual,
        basePosition=(np.asarray(origin)+R@c).tolist(),baseOrientation=Rotation.from_matrix(R@U).as_quat().tolist(),
        linkMasses=[0.]*k,linkCollisionShapeIndices=collisions[1:],linkVisualShapeIndices=[-1]*k,
        linkPositions=[[0,0,0]]*k,linkOrientations=[[0,0,0,1]]*k,
        linkInertialFramePositions=[[0,0,0]]*k,linkInertialFrameOrientations=[[0,0,0,1]]*k,
        linkParentIndices=[0]*k,linkJointTypes=[p.JOINT_FIXED]*k,linkJointAxis=[[0,0,0]]*k,
        useMaximalCoordinates=False,physicsClientId=client)
    p.changeDynamics(body,-1,mass=float(m),localInertiaDiagonal=d.tolist(),lateralFriction=.8,
        rollingFriction=0,spinningFriction=0,restitution=0,linearDamping=0,angularDamping=0,
        collisionMargin=0,maxJointVelocity=1000,physicsClientId=client)
    for link in range(k):
        p.changeDynamics(body,link,lateralFriction=.8,rollingFriction=0,spinningFriction=0,
            restitution=0,collisionMargin=0,physicsClientId=client)
        # Bullet auto-generates collision visuals even for visual index -1.
        # Hide those approximations; the base visual is the actual full outer.
        p.changeVisualShape(body,link,rgbaColor=[0,0,0,0],physicsClientId=client)
    omega=R@np.array([0,0,2*np.pi*hz])
    p.resetBaseVelocity(body,linearVelocity=np.cross(omega,R@c).tolist(),
        angularVelocity=omega.tolist(),physicsClientId=client)
    return body


def simulate_contact(outer_file,collision_files,props,outer_vertices,*,hz,tilt,azimuth,
                     duration=8.,dt=1/2000,clearance=.0002,gui=False,frame_callback=None,
                     dissipative_defaults=False):
    import pybullet as p
    import time
    if duration<=0 or dt<=0 or not np.isfinite([hz,tilt,azimuth,duration,dt,clearance]).all():
        raise ValueError('invalid simulation parameters')
    client=p.connect(p.GUI if gui else p.DIRECT)
    if client<0: raise RuntimeError('PyBullet connection failed')
    try:
        p.setRealTimeSimulation(0,physicsClientId=client)
        p.setPhysicsEngineParameter(fixedTimeStep=dt,numSolverIterations=150,
            deterministicOverlappingPairs=1,physicsClientId=client)
        p.setGravity(0,0,-9.81,physicsClientId=client)
        floor_shape=p.createCollisionShape(p.GEOM_BOX,halfExtents=[.4,.4,.02],physicsClientId=client)
        floor_visual=p.createVisualShape(p.GEOM_BOX,halfExtents=[.4,.4,.02],
            rgbaColor=[.45,.29,.15,1],physicsClientId=client)
        floor=p.createMultiBody(baseMass=0,baseCollisionShapeIndex=floor_shape,baseVisualShapeIndex=floor_visual,
            basePosition=[0,0,-.02],physicsClientId=client)
        p.changeDynamics(floor,-1,lateralFriction=.8,restitution=0,physicsClientId=client)
        R=Rotation.from_rotvec(np.deg2rad(tilt)*np.array([np.cos(azimuth),np.sin(azimuth),0])).as_matrix()
        origin=np.array([0,0,clearance-(outer_vertices@R.T)[:,2].min()])
        body=create_body(p,client,outer_file,collision_files,props,R,origin,hz)
        if dissipative_defaults:
            for link in range(-1,p.getNumJoints(body,physicsClientId=client)):
                p.changeDynamics(body,link,rollingFriction=.001,spinningFriction=.001,
                    linearDamping=.04,angularDamping=.04,physicsClientId=client)
        if p.getNumConstraints(physicsClientId=client)!=0: raise ValueError('unexpected fixed constraint')
        m,c,I,d,U=props
        close=p.getClosestPoints(body,floor,.01,physicsClientId=client)
        initial_gap=min((q[8] for q in close),default=None)
        support_ok=initial_gap is not None and abs(initial_gap-clearance)<=.0001
        margins=[dict(link=i,margin_m=p.getDynamicsInfo(body,i,physicsClientId=client)[-1])
                 for i in range(-1,p.getNumJoints(body,physicsClientId=client))]
        if not support_ok:
            raise ValueError(f'COLLISION_MODEL_BLOCKER: actual initial gap {initial_gap}, expected {clearance}; margins={margins}')
        p.resetDebugVisualizerCamera(.11,0,-12,[0,0,.030],physicsClientId=client)
        if gui:
            p.addUserDebugText('Reconstructed heart / FREE CONTACT / no maintained torque',[-.05,0,.085],
                textSize=1.1,textColorRGB=[0,0,0],physicsClientId=client)
        samples=[]; max_tilt=0.; max_penetration=0.; max_count=0; max_energy=0.; heart_contact=False
        steps=int(round(duration/dt)); stride=max(1,round(1/(200*dt))); start=time.monotonic()
        for step in range(steps+1):
            if step: p.stepSimulation(physicsClientId=client)
            pos,q=p.getBasePositionAndOrientation(body,physicsClientId=client)
            linear,angular=p.getBaseVelocity(body,physicsClientId=client)
            if not np.isfinite(np.r_[pos,q,linear,angular]).all():
                raise ValueError('NONFINITE_PHYSICS_STATE: rejected, not a stability result')
            Q=np.array(p.getMatrixFromQuaternion(q)).reshape(3,3)
            axis=Q@U.T@np.array([0,0,1]); angle=float(np.rad2deg(np.arccos(np.clip(axis[2],-1,1))))
            wlocal=Q.T@angular
            energy=float(.5*m*np.dot(linear,linear)+.5*np.dot(d*wlocal,wlocal)+m*9.81*pos[2])
            if step==0: e0=energy
            max_energy=max(max_energy,energy); max_tilt=max(max_tilt,angle)
            contacts=p.getContactPoints(body,floor,physicsClientId=client)
            heart_contact=heart_contact or any(cp[3]>=0 for cp in contacts)
            max_count=max(max_count,len(contacts))
            max_penetration=max(max_penetration,max((-cp[8] for cp in contacts),default=0))
            if step%stride==0 or step==steps:
                samples.append(dict(time=step*dt,position=list(pos),principal_quaternion=list(q),
                    model_quaternion=Rotation.from_matrix(Q@U.T).as_quat().tolist(),
                    angular_velocity=list(angular),linear_velocity=list(linear),tilt_degrees=angle,
                    spin_hz=float(np.dot(angular,axis)/(2*np.pi)),energy_J=energy,
                    contacts=[dict(link=cp[3],on_body=list(cp[5]),on_table=list(cp[6]),normal=list(cp[7]),
                        distance=cp[8],normal_force=cp[9],friction_1=cp[10],friction_direction_1=list(cp[11]),
                        friction_2=cp[12],friction_direction_2=list(cp[13])) for cp in contacts]))
                if frame_callback is not None: frame_callback(p,client,step*dt)
                if gui:
                    delay=step*dt-(time.monotonic()-start)
                    if delay>0: time.sleep(min(delay,.05))
        return dict(protocol='PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT',duration=duration,dt=dt,hz=hz,
            initial_tilt=tilt,azimuth=azimuth,dissipative_defaults=dissipative_defaults,
            maximum_tilt_degrees=max_tilt,final_tilt_degrees=samples[-1]['tilt_degrees'],
            maximum_penetration_m=max_penetration,maximum_contact_count=max_count,
            initial_gap_m=initial_gap,collision_initial_support_passed=support_ok,
            runtime_collision_margins=margins,heart_proxy_contact_occurred=heart_contact,
            validation_status='OBSERVATIONAL_NOT_VALIDATED',
            collision_scope='zero-margin shaft base; padded V-HACD heart links; fallen-contact approximate',
            relative_energy_increase=max(0.,(max_energy-e0)/max(abs(e0),1e-30)),
            constraint_count=p.getNumConstraints(physicsClientId=client),samples=samples,
            final_spin_hz=samples[-1]['spin_hz'],api_version=p.getAPIVersion(),elapsed=time.monotonic()-start)
    finally: p.disconnect(client)
