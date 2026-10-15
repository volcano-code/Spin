"""Versioned common-frame collision constructor; legacy evidence stays unchanged.

This corrects material-frame-dependent convex hull construction. It does not
certify hull fidelity, contact integration, wall thickness or physical behavior.
"""
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation

def create_body(p,client,outer_file,collision_files,props,R,origin,hz,color=(.85,.57,.12,1)):
    m,c,I,d,U=props
    if not collision_files:
        raise ValueError('missing convex collision parts')
    local_rotation=Rotation.from_matrix(U.T).as_quat().tolist()
    local_position=(-U.T@c).tolist()
    shapes=[]
    for file in collision_files:
        vertices=np.asarray(trimesh.load(file,force='mesh',process=False).vertices)
        # Hull input is invariant under changes to material COM/principal axes.
        shapes.append(p.createCollisionShape(p.GEOM_MESH,vertices=vertices.tolist(),
            collisionFramePosition=local_position,collisionFrameOrientation=local_rotation,physicsClientId=client))
    visual=p.createVisualShape(p.GEOM_MESH,fileName=str(outer_file),meshScale=[1,1,1],
        visualFramePosition=local_position,visualFrameOrientation=local_rotation,rgbaColor=color,physicsClientId=client)
    k=len(shapes)-1
    body=p.createMultiBody(baseMass=m,baseCollisionShapeIndex=shapes[0],baseVisualShapeIndex=visual,
        basePosition=(np.asarray(origin)+R@c).tolist(),baseOrientation=Rotation.from_matrix(R@U).as_quat().tolist(),
        linkMasses=[0.]*k,linkCollisionShapeIndices=shapes[1:],linkVisualShapeIndices=[-1]*k,
        linkPositions=[[0,0,0]]*k,linkOrientations=[[0,0,0,1]]*k,
        linkInertialFramePositions=[[0,0,0]]*k,linkInertialFrameOrientations=[[0,0,0,1]]*k,
        linkParentIndices=[0]*k,linkJointTypes=[p.JOINT_FIXED]*k,linkJointAxis=[[0,0,0]]*k,
        useMaximalCoordinates=False,physicsClientId=client)
    p.changeDynamics(body,-1,mass=float(m),localInertiaDiagonal=d.tolist(),lateralFriction=.8,
        rollingFriction=0,spinningFriction=0,restitution=0,linearDamping=0,angularDamping=0,
        collisionMargin=0,maxJointVelocity=1000,physicsClientId=client)
    for link in range(k):
        # Retain the legacy padding behavior; it remains a separate audit gate.
        p.changeDynamics(body,link,lateralFriction=.8,rollingFriction=0,spinningFriction=0,
            restitution=0,collisionMargin=0,physicsClientId=client)
        p.changeVisualShape(body,link,rgbaColor=[0,0,0,0],physicsClientId=client)
    omega=R@np.array([0,0,2*np.pi*hz])
    p.resetBaseVelocity(body,linearVelocity=np.cross(omega,R@c).tolist(),angularVelocity=omega.tolist(),physicsClientId=client)
    return body
