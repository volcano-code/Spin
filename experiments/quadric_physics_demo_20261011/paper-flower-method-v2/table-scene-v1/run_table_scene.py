"""Portable, unpowered PyBullet table-mounted bearing experiment and replay."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_bundle(bundle):
    packet = json.loads((bundle / 'table-scene.json').read_text())
    for name, digest in packet['files'].items():
        if sha(bundle / name) != digest:
            raise ValueError('BUNDLE_FILE_CHANGED:' + name)
    return packet


def build_world(bundle, packet, dt, *, visuals):
    import pybullet as p
    from scipy.spatial.transform import Rotation
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, -9.81, physicsClientId=client)
    p.setPhysicsEngineParameter(fixedTimeStep=dt, numSolverIterations=150, physicsClientId=client)
    def box(half, position, color):
        col = p.createCollisionShape(p.GEOM_BOX, halfExtents=half, physicsClientId=client)
        vis = p.createVisualShape(p.GEOM_BOX, halfExtents=half, rgbaColor=color, physicsClientId=client) if visuals else -1
        return p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col, baseVisualShapeIndex=vis,
                                 basePosition=position, physicsClientId=client)
    box([1, 1, .01], [0, 0, -.01], [.81, .84, .86, 1])
    box([.16, .10, .015], [0, 0, .705], [.58, .36, .19, 1])
    for x in (-.13, .13):
        for y in (-.075, .075):
            box([.008, .008, .345], [x, y, .345], [.28, .30, .32, 1])
    bodies, diagnostics = {}, {}
    for variant, x in [('original', -.050), ('optimized', .050)]:
        axis = [x, 0, .738]
        box([.009, .009, .005], [x, 0, .725], [.24, .27, .30, 1])
        box([.003, .003, .004], [x, 0, .734], [.42, .45, .48, 1])
        value = packet['properties'][variant]
        mass, c, inertia = value['mass_kg'], np.array(value['center_m']), np.array(value['inertia_kg_m2'])
        eigenvalues, basis = np.linalg.eigh(inertia)
        if mass <= 0 or min(eigenvalues) <= 0 or not np.isfinite(inertia).all():
            raise ValueError('INVALID_MASS_PROPERTIES')
        if np.linalg.det(basis) < 0:
            basis[:, 0] *= -1
        names = packet['rotor_parts'][variant]
        def visual(parts):
            if not visuals:
                return -1
            colors = [([.90, .37, .62, 1] if variant == 'original' else [.96, .72, .17, 1])
                      if name in ('design', 'retained') else [.44, .48, .52, 1] for name in parts]
            return p.createVisualShapeArray(shapeTypes=[p.GEOM_MESH]*len(parts),
                fileNames=[str(bundle / 'parts' / (name + '.obj')) for name in parts],
                meshScales=[[1, 1, 1]]*len(parts), rgbaColors=colors,
                visualFramePositions=[[0, 0, 0]]*len(parts),
                visualFrameOrientations=[[0, 0, 0, 1]]*len(parts), physicsClientId=client)
        body = p.createMultiBody(baseMass=0, basePosition=axis,
            baseVisualShapeIndex=visual(packet['clamp_parts']), baseCollisionShapeIndex=-1,
            linkMasses=[mass], linkCollisionShapeIndices=[-1], linkVisualShapeIndices=[visual(names)],
            linkPositions=[[0, 0, 0]], linkOrientations=[[0, 0, 0, 1]],
            linkInertialFramePositions=[c.tolist()],
            linkInertialFrameOrientations=[Rotation.from_matrix(basis).as_quat().tolist()],
            linkParentIndices=[0], linkJointTypes=[p.JOINT_REVOLUTE], linkJointAxis=[[0, 0, 1]],
            physicsClientId=client)
        p.changeDynamics(body, 0, maxJointVelocity=1000, physicsClientId=client)
        p.changeDynamics(body, 0, localInertiaDiagonal=eigenvalues.tolist(), linearDamping=0,
            angularDamping=0, jointDamping=0, physicsClientId=client)
        p.setJointMotorControl2(body, 0, p.VELOCITY_CONTROL, force=0, physicsClientId=client)
        p.enableJointForceTorqueSensor(body, 0, True, physicsClientId=client)
        p.resetJointState(body, 0, 0, 2*math.pi*packet['initial_hz'], physicsClientId=client)
        actual = p.getDynamicsInfo(body, 0, physicsClientId=client)
        loaded = basis @ np.diag(actual[2]) @ basis.T
        np.testing.assert_allclose(loaded, inertia, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(actual[3], c, rtol=1e-12, atol=1e-18)
        bodies[variant] = body
        diagnostics[variant] = dict(mass_kg=actual[0], center_m=list(actual[3]),
            loaded_inertia_kg_m2=loaded.tolist(), axis_world_m=axis,
            axis_inertia=float(inertia[2, 2] + mass*np.dot(c[:2], c[:2])),
            analytic_imbalance_force_N=float(mass*np.linalg.norm(c[:2])*(2*math.pi*packet['initial_hz'])**2))
    return client, bodies, diagnostics


def simulate(bundle, output, packet):
    import pybullet as p
    output.mkdir(parents=True, exist_ok=False)
    scenes = []
    for dt in packet['time_steps_seconds']:
        client, bodies, diagnostics = build_world(bundle, packet, dt, visuals=False)
        try:
            steps = round(packet['duration_seconds']/dt)
            arrays = {name: np.empty((steps + 1, 11)) for name in bodies}
            for step in range(steps + 1):
                if step:
                    p.stepSimulation(physicsClientId=client)
                for name, body in bodies.items():
                    angle, omega, reaction, torque = p.getJointState(body, 0, physicsClientId=client)
                    energy = .5*diagnostics[name]['axis_inertia']*omega**2
                    arrays[name][step] = [step*dt, angle, omega, *reaction, torque, energy]
            for name, data in arrays.items():
                file = f'{name}-dt{dt}.npz'
                np.savez_compressed(output/file, records=data,
                    columns=['t', 'angle', 'omega', 'Fx', 'Fy', 'Fz', 'Tx', 'Ty', 'Tz', 'motor_torque', 'energy'])
                force = np.linalg.norm(data[1:, 3:5], axis=1)
                row = dict(variant=name, dt=dt, steps=steps, file=file, sha256=sha(output/file),
                    maximum_motor_torque_Nm=float(abs(data[:, 9]).max()),
                    relative_energy_drift=float(abs(data[:, 10]/data[0, 10]-1).max()),
                    mean_horizontal_reaction_N=float(force.mean()),
                    final_hz=float(data[-1, 2]/(2*math.pi)), **diagnostics[name])
                row['analytic_force_error_N'] = abs(row['mean_horizontal_reaction_N']-row['analytic_imbalance_force_N'])
                if row['maximum_motor_torque_Nm'] != 0 or row['relative_energy_drift'] > 1e-6 or row['analytic_force_error_N'] > 1e-6:
                    raise ValueError('BEARING_BENCHMARK_FAILED:' + name)
                scenes.append(row)
            p.saveBullet(str(output/f'table-scene-dt{dt}.bullet'), physicsClientId=client)
        finally:
            p.disconnect(client)
    comparisons = {}
    for name in ('original', 'optimized'):
        coarse = np.load(output/f'{name}-dt0.0005.npz')['records']
        fine = np.load(output/f'{name}-dt0.00025.npz')['records'][::2]
        comparisons[name] = dict(max_angle_difference_rad=float(abs(coarse[:, 1]-fine[:, 1]).max()),
            max_hz_difference=float(abs(coarse[:, 2]-fine[:, 2]).max()/(2*math.pi)),
            mean_reaction_difference_N=float(abs(np.linalg.norm(coarse[1:, 3:5], axis=1).mean()-np.linalg.norm(fine[1:, 3:5], axis=1).mean())))
    summary = dict(schema='flower-pybullet-table-bearing-v1', scenes=scenes, comparisons=comparisons,
        source_config_sha256=sha(bundle/'table-scene.json'), actual_pybullet_integration=True,
        physics_scope=packet['physics_scope'], motor_drive=False, initial_phase=0,
        pybullet_api_version=p.getAPIVersion(), real_bearing_friction_calibrated=False)
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def render(bundle, output, packet):
    import pybullet as p
    from PIL import Image, ImageDraw, ImageFont
    summary = json.loads((output/'summary.json').read_text())
    if summary['source_config_sha256'] != sha(bundle/'table-scene.json'):
        raise ValueError('SIMULATION_CONFIG_MISMATCH')
    trajectories = {}
    for name in ('original', 'optimized'):
        entry = next(s for s in summary['scenes'] if s['variant'] == name and s['dt'] == .00025)
        path = output/entry['file']
        if sha(path) != entry['sha256']:
            raise ValueError('TRAJECTORY_CHANGED')
        trajectories[name] = np.load(path)['records']
    client, bodies, diagnostic = build_world(bundle, packet, .00025, visuals=True)
    font = ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc', 22)
    small = ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc', 17)
    width, height = 1080, 700
    projection = p.computeProjectionMatrixFOV(40, 540/540, .001, 3)
    def camera(eye, target, w, h, projection):
        return Image.fromarray(np.asarray(p.getCameraImage(w, h,
            p.computeViewMatrix(eye, target, [0, 0, 1]), projection,
            renderer=p.ER_TINY_RENDERER, shadow=1, lightDirection=[-1, -3, 5],
            physicsClientId=client)[2], dtype=np.uint8).reshape(h, w, 4)[:, :, :3])
    try:
        camera([.43, -.56, 1.05], [0, 0, .40], 1080, 700,
            p.computeProjectionMatrixFOV(44, 1080/700, .01, 4)).save(output/'table-overview.png')
        with (output/'video.log').open('xb') as log:
            process = subprocess.Popen(['ffmpeg', '-nostdin', '-n', '-f', 'rawvideo', '-pixel_format', 'rgb24',
                '-video_size', f'{width}x{height}', '-framerate', '24', '-i', 'pipe:0', '-an',
                '-c:v', 'libx264', '-threads', '4', '-preset', 'fast', '-crf', '20', '-pix_fmt', 'yuv420p',
                str(output/'flower-table-comparison.mp4')], stdin=subprocess.PIPE, stdout=log, stderr=log)
            try:
                for frame in range(432):
                    slow = frame >= 192
                    t = (frame-192)/480 if slow else frame/24
                    idx = min(round(t/.00025), len(trajectories['original'])-1)
                    for name, body in bodies.items():
                        data = trajectories[name][idx]
                        p.resetJointState(body, 0, data[1], data[2], physicsClientId=client)
                    image = Image.new('RGB', (width, height), '#f3f2ed')
                    draw = ImageDraw.Draw(image)
                    draw.text((20, 12), '桌面轴承支座 · PyBullet 实际轨迹', fill='#24313c', font=font)
                    for col, name in enumerate(('original', 'optimized')):
                        x = diagnostic[name]['axis_world_m'][0]
                        image.paste(camera([x+.042, -.070, .807], [x, 0, .735], 540, 540, projection), (540*col, 68))
                        draw.text((20+540*col, 44), '优化前' if name == 'original' else '优化后', fill='#24313c', font=font)
                        force = np.linalg.norm(trajectories[name][idx, 3:5])
                        draw.text((20+540*col, 612), f'{trajectories[name][idx,2]/(2*math.pi):.2f} Hz | 水平反力 {force:.3g} N', fill='#24313c', font=small)
                    draw.text((20, 648), f'物理时间 {t:.4f} s · '+('0.05倍慢放' if slow else '正常速度；高速转动存在帧采样效应'), fill='#24313c', font=small)
                    draw.text((20, 675), '中心轴承固定于桌上支座；无持续驱动。比较偏心载荷，不预测真实停转时间。', fill='#53606b', font=small)
                    if frame == 0:
                        image.save(output/'preview.png')
                    process.stdin.write(np.asarray(image).tobytes())
            finally:
                process.stdin.close()
            if process.wait(timeout=60):
                raise RuntimeError('VIDEO_ENCODING_FAILED')
        (output/'render.json').write_text(json.dumps(dict(video_sha256=sha(output/'flower-table-comparison.mp4'),
            simulation_summary_sha256=sha(output/'summary.json'), replay_from_saved_states=True,
            physical_intervals=[[0, 8], [0, .5]], playback_rates=[1, .05]), indent=2)+'\n')
    finally:
        p.disconnect(client)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundle', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output', type=Path, default=Path('table-result'))
    parser.add_argument('--render-only', action='store_true')
    args = parser.parse_args()
    packet = load_bundle(args.bundle)
    if not args.render_only:
        simulate(args.bundle, args.output, packet)
    render(args.bundle, args.output, packet)
