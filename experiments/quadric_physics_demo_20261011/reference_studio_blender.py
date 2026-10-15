"""Render unchanged supplied triangle meshes with studio lighting, in Blender."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import bpy
from mathutils import Vector

def material(role,opacity):
    colors={'fixed':(.045,.05,.055,1),'outer-envelope':(.25,.25,.25,1),
        'design':(.52,.012,.15,1),'design-candidate':(.52,.012,.15,1),'retained':(.92,.43,.008,1)}
    m=bpy.data.materials.new(role);m.use_nodes=True
    nodes=m.node_tree.nodes;links=m.node_tree.links;p=nodes.get('Principled BSDF')
    p.inputs['Base Color'].default_value=colors[role];p.inputs['Roughness'].default_value=.32
    p.inputs['Metallic'].default_value=.08 if role=='retained' else 0
    if opacity<1:
        t=nodes.new('ShaderNodeBsdfTransparent');mix=nodes.new('ShaderNodeMixShader')
        mix.inputs[0].default_value=opacity
        links.new(t.outputs[0],mix.inputs[1]);links.new(p.outputs[0],mix.inputs[2])
        links.new(mix.outputs[0],nodes.get('Material Output').inputs['Surface'])
    return m

def point_at(obj,target):obj.rotation_euler=(Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()

def render(entry,out):
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.device='CPU'
    scene.cycles.samples=48;scene.cycles.use_denoising=True
    scene.cycles.transparent_max_bounces=32;scene.render.threads_mode='FIXED';scene.render.threads=4
    scene.render.resolution_x=1440;scene.render.resolution_y=1200;scene.render.resolution_percentage=100
    scene.world=bpy.data.worlds.new('Studio');scene.world.use_nodes=True
    scene.world.node_tree.nodes['Background'].inputs[0].default_value=(.8,.8,.8,1)
    scene.world.node_tree.nodes['Background'].inputs[1].default_value=.65
    scene.view_settings.view_transform='AgX'
    points=[v for p in entry['parts'] for v in p['vertices']]
    lo=[min(v[k] for v in points) for k in range(3)];hi=[max(v[k] for v in points) for k in range(3)]
    center=Vector([(a+b)/2 for a,b in zip(lo,hi)]);scale=4/max(b-a for a,b in zip(lo,hi))
    for part in entry['parts']:
        if hashlib.sha256(Path(part['source']).read_bytes()).hexdigest()!=part['source_sha256']:
            raise ValueError('RENDER_SOURCE_IDENTITY_MISMATCH')
        vertices=[(Vector(v)-center)*scale for v in part['vertices']]
        mesh=bpy.data.meshes.new(part['name']);mesh.from_pydata(vertices,[],part['faces']);mesh.update()
        obj=bpy.data.objects.new(part['name'],mesh);scene.collection.objects.link(obj)
        mesh.materials.append(material(part['role'],part['opacity']))
        if part['smooth']:
            for p in mesh.polygons:p.use_smooth=True
    ground=(lo[2]-center.z)*scale-.16
    bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,ground));plane=bpy.context.object
    floor=bpy.data.materials.new('Ground');floor.diffuse_color=(.76,.76,.76,1);floor.use_nodes=True
    floor.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value=(.76,.76,.76,1)
    floor.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value=.8;plane.data.materials.append(floor)
    for pos,power,size in [((-3,-4,7),700,5),((4,-1,5),450,4),((1,5,6),600,3)]:
        light=bpy.data.lights.new('Softbox','AREA');light.energy=power;light.shape='DISK';light.size=size
        obj=bpy.data.objects.new('Softbox',light);scene.collection.objects.link(obj);obj.location=pos;point_at(obj,(0,0,0))
    cam=bpy.data.cameras.new('Camera');obj=bpy.data.objects.new('Camera',cam);scene.collection.objects.link(obj)
    obj.location=entry['camera'];point_at(obj,(0,0,0));cam.type='ORTHO'
    inverse=obj.rotation_euler.to_matrix().transposed()
    projected=[inverse@((Vector(v)-center)*scale) for v in points]
    width=max(v.x for v in projected)-min(v.x for v in projected)
    height=max(v.y for v in projected)-min(v.y for v in projected)
    cam.ortho_scale=max(width,height*scene.render.resolution_x/scene.render.resolution_y)*1.18
    scene.camera=obj
    scene.render.image_settings.file_format='PNG';scene.render.filepath=str(out/(entry['name']+'.png'))
    bpy.ops.render.render(write_still=True)

parser=argparse.ArgumentParser();parser.add_argument('--scene',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:]);args.output.mkdir(parents=True,exist_ok=False)
bundle=json.loads(args.scene.read_text())
for entry in bundle['scenes']:render(entry,args.output)
(args.output/'render-manifest.json').write_text(json.dumps(dict(schema='reference-studio-render-v1',
    scene_sha256=hashlib.sha256(args.scene.read_bytes()).hexdigest(),engine='BLENDER_CYCLES_CPU',
    display_only=True,optimization_invocations=0,simulation_invocations=0,
    images={e['name']:hashlib.sha256((args.output/(e['name']+'.png')).read_bytes()).hexdigest() for e in bundle['scenes']}),indent=2))
