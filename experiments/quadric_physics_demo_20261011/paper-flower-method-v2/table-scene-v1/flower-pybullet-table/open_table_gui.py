"""Open the exported URDF in PyBullet GUI. Close the window to exit."""
import json
from pathlib import Path
import time
import math
import pybullet as p

root = Path(__file__).resolve().parent
packet = json.loads((root/'table-scene.json').read_text())
client = p.connect(p.GUI)
try:
    p.setGravity(0, 0, -9.81)
    p.setPhysicsEngineParameter(fixedTimeStep=.00025, numSolverIterations=150)
    body = p.loadURDF(str(root/'flower-table.urdf'), useFixedBase=True,
                      flags=p.URDF_USE_INERTIA_FROM_FILE)
    for index in range(p.getNumJoints(body)):
        name = p.getJointInfo(body, index)[1].decode().removesuffix('_bearing')
        # Set velocity limits first: some Bullet builds reset link inertia here.
        before = p.getDynamicsInfo(body, index)
        p.changeDynamics(body, index, maxJointVelocity=1000)
        p.changeDynamics(body, index, localInertiaDiagonal=before[2],
                         linearDamping=0, angularDamping=0, jointDamping=0)
        p.setJointMotorControl2(body, index, p.VELOCITY_CONTROL, force=0)
        p.resetJointState(body, index, 0, 2*math.pi*packet['initial_hz'])
    p.resetDebugVisualizerCamera(.19, 0, -48, [0, 0, .738])
    start = time.monotonic()
    for step in range(32000):
        if not p.isConnected():
            break
        p.stepSimulation()
        if step % 40 == 0:
            time.sleep(max(0, start + step*.00025 - time.monotonic()))
    while p.isConnected():
        time.sleep(.1)
finally:
    if p.isConnected():
        p.disconnect()
