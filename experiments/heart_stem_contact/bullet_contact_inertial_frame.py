"""Candidate version: place the material frame on the body, not inside hulls.

Opt-in only until geometry and dynamic contracts are verified. Historical
constructors and evidence remain unchanged.
"""
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation

def create_body(p,client,outer_file,collision_files,props,R,origin,hz,color=(.85,.57,.12,1)):
    if not collision_files:raise ValueError('missing convex collision parts')
    m,c,I,d,U=props
    shapes=[p.createCollisionShape(p.GEOM_MESH,
        vertices=np.asarray(trimesh.load(file,force='mesh',process=False).vertices).tolist(),
        physicsClientId=client) for file in collision_files]
    # The base-collider update reaches immediate compound children whereas
    # the link-collider update may only change the compound's own margin.
    # Initialize the shared shape before attaching it to a fixed link.
    for shape in shapes:
        temporary=p.createMultiBody(baseMass=0,baseCollisionShapeIndex=shape,physicsClientId=client)
        try:p.changeDynamics(temporary,-1,collisionMargin=0,physicsClientId=client)
        finally:p.removeBody(temporary,physicsClientId=client)
    visual=p.createVisualShape(p.GEOM_MESH,fileName=str(outer_file),rgbaColor=color,physicsClientId=client)
    k=len(shapes)-1
    body=p.createMultiBody(baseMass=m,baseCollisionShapeIndex=shapes[0],baseVisualShapeIndex=visual,
        basePosition=np.asarray(origin).tolist(),baseOrientation=Rotation.from_matrix(R).as_quat().tolist(),
        baseInertialFramePosition=np.asarray(c).tolist(),baseInertialFrameOrientation=Rotation.from_matrix(U).as_quat().tolist(),
        linkMasses=[0.]*k,linkCollisionShapeIndices=shapes[1:],linkVisualShapeIndices=[-1]*k,
        linkPositions=[[0,0,0]]*k,linkOrientations=[[0,0,0,1]]*k,
        linkInertialFramePositions=[[0,0,0]]*k,linkInertialFrameOrientations=[[0,0,0,1]]*k,
        linkParentIndices=[0]*k,linkJointTypes=[p.JOINT_FIXED]*k,linkJointAxis=[[0,0,0]]*k,
        useMaximalCoordinates=False,physicsClientId=client)
    p.changeDynamics(body,-1,mass=float(m),localInertiaDiagonal=d.tolist(),maxJointVelocity=1000,physicsClientId=client)
    for link in range(-1,k):
        p.changeDynamics(body,link,lateralFriction=.8,rollingFriction=0,spinningFriction=0,
            restitution=0,linearDamping=0,angularDamping=0,collisionMargin=0,physicsClientId=client)
    omega=R@np.array([0,0,2*np.pi*hz])
    p.resetBaseVelocity(body,linearVelocity=np.cross(omega,R@c).tolist(),angularVelocity=omega.tolist(),physicsClientId=client)
    return body
