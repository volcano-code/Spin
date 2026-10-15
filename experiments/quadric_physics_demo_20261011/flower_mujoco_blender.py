"""Cycles replay of actual MuJoCo states, preserving all three material layers."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import bpy
import numpy as np
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parent))
from flower_render_material_policy import require_material_profile, comparison_material_role


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def aim(obj,target): obj.rotation_euler=(Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()


def mat(name,color,metal=0,rough=.35):
    m=bpy.data.materials.new(name);m.use_nodes=True
    p=m.node_tree.nodes.get('Principled BSDF');p.inputs['Base Color'].default_value=(*color,1)
    p.inputs['Metallic'].default_value=metal;p.inputs['Roughness'].default_value=rough
    return m


def box(name,position,scale,material):
    bpy.ops.mesh.primitive_cube_add(size=1,location=position)
    o=bpy.context.object;o.name=name;o.scale=scale
    bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
    o.data.materials.append(material)
    bevel=o.modifiers.new('Soft manufactured edges','BEVEL');bevel.width=.05;bevel.segments=3
    return o


def main(bundle,physics,out,material_profile='PAPER_ROLES_V1'):
    require_material_profile(material_profile)
    out.mkdir(parents=True,exist_ok=False)
    packet=json.loads((bundle/'table-scene.json').read_text());summary=json.loads((physics/'summary.json').read_text())
    if summary['input_sha256']!=sha(bundle/'table-scene.json'):raise ValueError('PHYSICS_INPUT_CHANGED')
    for n,h in packet['files'].items():
        if sha(bundle/n)!=h:raise ValueError('MESH_CHANGED')
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=16
    scene.cycles.use_denoising=True;scene.render.threads_mode='FIXED';scene.render.threads=4
    scene.render.resolution_x=960;scene.render.resolution_y=600;scene.render.resolution_percentage=100
    scene.render.fps=24;scene.frame_start=1;scene.frame_end=288
    scene.world=bpy.data.worlds.new('Soft studio');scene.world.use_nodes=True
    scene.world.node_tree.nodes['Background'].inputs[0].default_value=(.6,.67,.75,1)
    scene.world.node_tree.nodes['Background'].inputs[1].default_value=.35
    scene.view_settings.view_transform='AgX'
    grey=mat('Fixed shell graphite',(.11,.13,.15));yellow=mat('Optimized polymer',(.93,.48,.045));pink=mat('Original design polymer',(.63,.055,.23))
    polymer=mat('Shared comparison polymer',(.93,.48,.045))
    metal=mat('Bearing steel',(.42,.46,.5),.8,.22);wood=mat('Oak table',(.32,.16,.065),0,.4)
    nodes=wood.node_tree.nodes;links=wood.node_tree.links
    tex=nodes.new('ShaderNodeTexNoise');tex.inputs['Scale'].default_value=5
    coord=nodes.new('ShaderNodeTexCoord');mapping=nodes.new('ShaderNodeVectorMath');mapping.operation='MULTIPLY';mapping.inputs[1].default_value=(1,28,3)
    ramp=nodes.new('ShaderNodeValToRGB');ramp.color_ramp.elements[0].color=(.08,.032,.012,1);ramp.color_ramp.elements[1].color=(.5,.28,.11,1)
    links.new(coord.outputs['Generated'],mapping.inputs[0]);links.new(mapping.outputs[0],tex.inputs['Vector']);links.new(tex.outputs['Fac'],ramp.inputs[0]);links.new(ramp.outputs[0],nodes.get('Principled BSDF').inputs['Base Color'])
    box('Table top',(0,0,-1.5),(32,20,3),wood)
    floor=mat('Room',(.5,.55,.6),0,.8)
    box('Floor',(0,0,-73),(250,250,1),floor)
    for x in (-13,13):
        for y in (-7.5,7.5):box('Table leg',(x,y,-37),(1.6,1.6,68),grey)
    moving={};rotors={};top_objects=[];bottom_objects=[]
    rows={}
    for variant,x in [('original',-5),('optimized',5)]:
        row=next(r for r in summary['scenes'] if r['kind']=='flexible' and r['variant']==variant and r['phase_index']==0 and r['dt']==.00025)
        if sha(physics/row['file'])!=row['sha256']:raise ValueError('TRAJECTORY_CHANGED')
        rows[variant]=np.load(physics/row['file'])['records']
        box('Fixed pedestal',(x,0,.5),(1.8,1.8,1),grey)
        box('Linear support housing',(x,0,1.25),(1.4,1.4,.5),metal)
        support=bpy.data.objects.new(variant+' support',None);scene.collection.objects.link(support);support.location=(x,0,1.8);moving[variant]=support
        rotor=bpy.data.objects.new(variant+' rotor',None);scene.collection.objects.link(rotor);rotor.parent=support;rotors[variant]=rotor
        for name in packet['rotor_parts'][variant]+packet['clamp_parts']:
            vertices=[];faces=[]
            for line in (bundle/'parts'/(name+'.obj')).read_text().splitlines():
                if line.startswith('v '):vertices.append(tuple(float(v)*100 for v in line.split()[1:4]))
                if line.startswith('f '):faces.append(tuple(int(v.split('/')[0])-1 for v in line.split()[1:]))
            mesh=bpy.data.meshes.new(name);mesh.from_pydata(vertices,[],faces);mesh.update()
            obj=bpy.data.objects.new(variant+' '+name,mesh);scene.collection.objects.link(obj)
            obj.parent=rotor if name in packet['rotor_parts'][variant] else support
            if material_profile=='MATCHED_COMPARISON_V1':
                material={'polymer':polymer,'metal':metal,'fixed':grey}[
                    comparison_material_role(name,packet['clamp_parts'])]
            else:
                material=yellow if name=='retained' else pink if name=='design' else metal if name in packet['clamp_parts']+['outer-race'] else grey
            mesh.materials.append(material)
            if name=='top':top_objects.append(obj)
            if name=='bottom':bottom_objects.append(obj)
        # All transforms below come from actual samples. Display length scale
        # is uniformly 100; support vibration is NOT artificially amplified.
        for frame in range(1,289):
            t=(frame-1)/24 if frame<=192 else (frame-193)*.05/24
            data=rows[variant][min(round(t/.00025),len(rows[variant])-1)]
            support.location=(x+100*data[1],100*data[2],1.8)
            rotor.rotation_euler=(0,0,float(data[3]))
            support.keyframe_insert(data_path='location',frame=frame)
            rotor.keyframe_insert(data_path='rotation_euler',frame=frame)
    for location,power,size in [((-8,-12,22),1600,12),((12,-4,14),1100,10),((0,12,18),1600,9)]:
        light=bpy.data.lights.new('Softbox','AREA');light.energy=power;light.shape='DISK';light.size=size
        obj=bpy.data.objects.new('Softbox',light);scene.collection.objects.link(obj);obj.location=location;aim(obj,(0,0,1))
    cam=bpy.data.cameras.new('Camera');obj=bpy.data.objects.new('Camera',cam);scene.collection.objects.link(obj)
    obj.location=(8,-20,20);aim(obj,(0,0,1.5));cam.type='ORTHO';cam.ortho_scale=19;scene.camera=obj
    scene.frame_set(1);scene.render.image_settings.file_format='PNG';scene.render.filepath=str(out/'preview.png')
    bpy.ops.render.render(write_still=True)
    # Separate labelled diagnostic image: explicitly move covers only for this
    # still, restore them before saving the actual simulation replay scene.
    for o in top_objects:o.location.z=1.8
    for o in bottom_objects:o.location.z=-1.0
    scene.render.filepath=str(out/'three-layer-exploded-display-only.png');bpy.ops.render.render(write_still=True)
    for o in top_objects+bottom_objects:o.location.z=0
    bpy.ops.wm.save_as_mainfile(filepath=str(out/'flower-mujoco-replay.blend'))
    frames=out/'frames';frames.mkdir();scene.render.filepath=str(frames/'frame-')
    bpy.ops.render.render(animation=True)
    (out/'render.json').write_text(json.dumps(dict(schema='mujoco-blender-replay-v1',
        physics_sha256=sha(physics/'summary.json'),material_sha256=sha(bundle/'table-scene.json'),
        blender_version=bpy.app.version_string,frames=288,fps=24,display_length_scale=100,
        vibration_amplification=1,intervals=[[0,8],[0,.2]],playback_rates=[1,.05],
        selected_phase=0,model='flexible',renders_actual_states=True,
        material_profile=material_profile,
        exploded_image_scope='DISPLAY_TRANSFORM_ONLY_NOT_SIMULATION',
        files={p.name:sha(p) for p in out.iterdir() if p.is_file()}),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,required=True);p.add_argument('--physics',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--material-profile',default='PAPER_ROLES_V1',choices=('PAPER_ROLES_V1','MATCHED_COMPARISON_V1'))
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);main(a.bundle,a.physics,a.output,a.material_profile)
