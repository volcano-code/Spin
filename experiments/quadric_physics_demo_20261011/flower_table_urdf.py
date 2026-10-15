"""Export the same table scene as an independently loadable two-rotor URDF."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def vector(xs):
    return ' '.join(format(float(x), '.17g') for x in xs)


def export(bundle):
    data = json.loads((bundle/'table-scene.json').read_text())
    root = ET.Element('robot', name='reconstructed_flower_table_comparison')
    table = ET.SubElement(root, 'link', name='table_and_fixed_supports')
    fixed_inertial = ET.SubElement(table, 'inertial')
    ET.SubElement(fixed_inertial, 'mass', value='0')
    ET.SubElement(fixed_inertial, 'inertia', ixx='0', iyy='0', izz='0', ixy='0', ixz='0', iyz='0')
    def box(size, xyz, color):
        for role in ('visual', 'collision'):
            node = ET.SubElement(table, role)
            ET.SubElement(node, 'origin', xyz=vector(xyz))
            ET.SubElement(ET.SubElement(node, 'geometry'), 'box', size=vector(size))
            if role == 'visual':
                ET.SubElement(ET.SubElement(node, 'material', name='surface'), 'color', rgba=vector(color))
    def mesh(link, name, xyz, color):
        node = ET.SubElement(link, 'visual')
        ET.SubElement(node, 'origin', xyz=vector(xyz))
        ET.SubElement(ET.SubElement(node, 'geometry'), 'mesh', filename='parts/'+name+'.obj')
        ET.SubElement(ET.SubElement(node, 'material', name=name), 'color', rgba=vector(color))
    box([2, 2, .02], [0, 0, -.01], [.81, .84, .86, 1])
    box([.32, .20, .03], [0, 0, .705], [.58, .36, .19, 1])
    for x in (-.13, .13):
        for y in (-.075, .075):
            box([.016, .016, .69], [x, y, .345], [.28, .30, .32, 1])
    for variant, x in [('original', -.05), ('optimized', .05)]:
        box([.018, .018, .010], [x, 0, .725], [.24, .27, .30, 1])
        box([.006, .006, .008], [x, 0, .734], [.42, .45, .48, 1])
        for name in data['clamp_parts']:
            mesh(table, name, [x, 0, .738], [.44, .48, .52, 1])
        rotor = ET.SubElement(root, 'link', name=variant+'_rotor')
        value = data['properties'][variant]
        inertial = ET.SubElement(rotor, 'inertial')
        ET.SubElement(inertial, 'origin', xyz=vector(value['center_m']))
        ET.SubElement(inertial, 'mass', value=str(value['mass_kg']))
        tensor = value['inertia_kg_m2']
        ET.SubElement(inertial, 'inertia', **{name:format(tensor[i][j], '.17g')
            for name, i, j in [('ixx',0,0),('iyy',1,1),('izz',2,2),('ixy',0,1),('ixz',0,2),('iyz',1,2)]})
        for name in data['rotor_parts'][variant]:
            color = ([.90,.37,.62,1] if variant=='original' else [.96,.72,.17,1]) if name in ('design','retained') else [.44,.48,.52,1]
            mesh(rotor, name, [0,0,0], color)
        joint = ET.SubElement(root, 'joint', name=variant+'_bearing', type='continuous')
        ET.SubElement(joint, 'parent', link='table_and_fixed_supports')
        ET.SubElement(joint, 'child', link=variant+'_rotor')
        ET.SubElement(joint, 'origin', xyz=vector([x,0,.738]))
        ET.SubElement(joint, 'axis', xyz='0 0 1')
        ET.SubElement(joint, 'dynamics', damping='0', friction='0')
    ET.indent(root)
    path = bundle/'flower-table.urdf'
    if path.exists():
        raise ValueError('URDF_ALREADY_EXISTS')
    ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)
    # Reload with Bullet's real URDF importer, preserving the off-diagonal tensor.
    import pybullet as p
    import numpy as np
    client = p.connect(p.DIRECT)
    try:
        body = p.loadURDF(str(path.resolve()), useFixedBase=True,
            flags=p.URDF_USE_INERTIA_FROM_FILE, physicsClientId=client)
        checked = {}
        for index in range(p.getNumJoints(body, physicsClientId=client)):
            name = p.getJointInfo(body, index, physicsClientId=client)[1].decode().removesuffix('_bearing')
            actual = p.getDynamicsInfo(body, index, physicsClientId=client)
            basis = np.array(p.getMatrixFromQuaternion(actual[4])).reshape(3,3)
            tensor = basis@np.diag(actual[2])@basis.T
            expected = data['properties'][name]
            np.testing.assert_allclose(tensor, expected['inertia_kg_m2'], rtol=1e-10, atol=1e-18)
            np.testing.assert_allclose(actual[3], expected['center_m'], rtol=1e-12, atol=1e-18)
            np.testing.assert_allclose(actual[0], expected['mass_kg'], rtol=1e-12)
            checked[name] = dict(mass_kg=actual[0], loaded_inertia_kg_m2=tensor.tolist())
        (bundle/'urdf-reload-check.json').write_text(json.dumps(checked, indent=2)+'\n')
    finally:
        p.disconnect(client)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundle', type=Path, required=True)
    args = parser.parse_args()
    export(args.bundle)
