"""Create a single-file assembly and a portable authenticated table experiment."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import trimesh
from flower_run_physics import admit


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def combined_obj(path, names, meshes):
    offset = 0
    with path.open('x') as stream:
        stream.write('# Units: meters. Complete assembly; independently moving parts retain groups.\nmtllib flower.mtl\n')
        for name in names:
            mesh = meshes[name]
            stream.write(f'o {name}\ng {name}\nusemtl '+('yellow' if name == 'retained' else 'pink' if name == 'design' else 'grey')+'\n')
            for point in mesh.vertices:
                stream.write('v '+' '.join(format(float(x), '.17g') for x in point)+'\n')
            for face in mesh.faces:
                stream.write('f '+' '.join(str(int(x)+1+offset) for x in face)+'\n')
            offset += len(mesh.vertices)
    # Independently parse the resulting text. No welding or topology repair.
    vertices, faces = [], []
    for line in path.read_text().splitlines():
        if line.startswith('v '): vertices.append(list(map(float, line.split()[1:])))
        if line.startswith('f '): faces.append([int(x)-1 for x in line.split()[1:]])
    expected = trimesh.util.concatenate([meshes[name] for name in names])
    np.testing.assert_array_equal(vertices, expected.vertices)
    np.testing.assert_array_equal(faces, expected.faces)
    loaded = trimesh.load(path, force='mesh', process=False, maintain_order=True)
    if len(loaded.faces) != len(faces) or not loaded.is_watertight or not loaded.is_winding_consistent:
        raise ValueError('WHOLE_OBJ_RELOAD_FAILED')
    return dict(parts=names, vertices=len(vertices), faces=len(faces), sha256=sha(path),
        coordinates_and_faces_unchanged=True, watertight_components=True,
        representation='SINGLE_FILE_MULTI_PART_ASSEMBLY_NOT_BOOLEAN_FUSED_BEARING')


def prepare(case, output):
    material = case/'prepare-job-02/output'
    proof = case/'saved-independent-job-03/output/saved-independent-check.json'
    record, audit = admit(material, proof.parent, case, saved_material_check=proof)
    output.mkdir(parents=True, exist_ok=False)
    parts = output/'parts'; parts.mkdir()
    names = ['design', 'retained', *record['rotor_layers'], *record['clamp_layers']]
    meshes = {}
    for name in names:
        source = material/record['files'][name]['file']
        if sha(source) != record['files'][name]['sha256']:
            raise ValueError('MATERIAL_SOURCE_CHANGED')
        shutil.copyfile(source, parts/(name+'.obj'))
        meshes[name] = trimesh.load(source, force='mesh', process=False)
    rotor = ['retained', *record['rotor_layers']]
    clamps = list(record['clamp_layers'])
    (output/'flower.mtl').write_text('newmtl grey\nKd 0.44 0.48 0.52\nnewmtl yellow\nKd 0.96 0.72 0.17\nnewmtl pink\nKd 0.90 0.37 0.62\n')
    checks = {}
    for filename, selected in [('flower-complete.obj', rotor+clamps),
                               ('flower-rotor.obj', rotor),
                               ('flower-cutaway.obj', [x for x in rotor+clamps if x != 'top'])]:
        checks[filename] = combined_obj(output/filename, selected, meshes)
    shutil.copyfile(Path(__file__).with_name('flower_table_scene.py'), output/'run_table_scene.py')
    properties = {key: record[key] for key in ('original', 'optimized', 'clamp')}
    (output/'physical-properties.json').write_text(json.dumps(properties, indent=2)+'\n')
    config = dict(schema='flower-table-scene-input-v1', classification='RECONSTRUCTED_METHOD_VALIDATION',
        physics_scope='TABLE_MOUNTED_IDEAL_BEARING_NO_FREE_TABLE_CONTACT_OF_ROTOR',
        properties=properties, rotor_parts=dict(original=['design', *record['rotor_layers']], optimized=rotor),
        clamp_parts=clamps, initial_hz=35, duration_seconds=8,
        time_steps_seconds=[.0005, .00025], units='meters', gravity_m_s2=[0, 0, -9.81],
        rotor_collisions='DISABLED_BY_IDEAL_BEARING_MODEL_WITH_CLEARANCE_ABOVE_TABLE',
        support='STATIC_TABLE_AND_BOLTED_PEDESTAL', bearing_motor=False,
        input_proof_sha256=sha(proof), input_material_sha256=sha(material/'output-preparation.json'),
        files={str(p.relative_to(output)): sha(p) for p in sorted(parts.iterdir())})
    (output/'table-scene.json').write_text(json.dumps(config, indent=2)+'\n')
    (output/'assembly-check.json').write_text(json.dumps(checks, indent=2)+'\n')
    (output/'README.txt').write_text('''完整装配：flower-complete.obj（与 flower.mtl 放在同一目录）。
转动件：flower-rotor.obj。查看中间切口：flower-cutaway.obj（展示性移开上盖）。
单位：米。完整装配保留轴承、固定夹持件与转子的独立网格，未焊死轴承。
全部原部件顶点和三角形逐项保持；输入材料已独立检查，见 assembly-check.json。

桌面仿真：固定桌子及桌上支座，中心理想轴承允许花盘自由旋转。
使用各自真实质量、质心和完整惯量；35Hz 初速，无持续驱动或阻尼。
转子有桌面净空，理想轴承模型不模拟转子落桌碰撞或真实轴承摩擦。
因此比较轴承不平衡载荷，不以此预测停转时间；桌面自由接触是不同问题。
原始与优化各两个时间步，首个登记相位；完整64场景工程支撑结果仍在原案例页。

环境：Python 3、numpy、scipy、trimesh、pybullet、Pillow、ffmpeg。
重跑：python run_table_scene.py --bundle . --output my-table-result
输出：summary.json、两个时间步的轨迹NPZ、Bullet快照、桌面预览及MP4。
脚本重建场景是可移植复现入口；Bullet快照仅作本版本诊断附件。
视频从保存的PyBullet关节状态回放，含正常速度及0.05倍慢放。

与论文的关系：二维圆形切除使质心回到轴承中心，与论文6.2节方法目标一致。
本模型是图片重建；几何、固定区域和物理支撑并非作者原始模型或实测装置。
''', encoding='utf-8')
    return config


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.case.resolve(), args.output.resolve())
