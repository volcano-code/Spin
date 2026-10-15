"""Illustrative bearing/cap assembly. No solver or dynamics acceptance."""
from pathlib import Path
import json
import hashlib
import sys

def assembly_spec():
    return dict(
        upper_cap=dict(radius=.19,z=[.12,.20],motion='FIXED_CLAMP'),
        lower_cap=dict(radius=.19,z=[-.20,-.12],motion='FIXED_CLAMP'),
        shaft=dict(radius=.052,z=[-.12,.12],motion='FIXED_CLAMP'),
        inner_race=dict(radius=.076,inner_radius=.055,z=[-.065,.065],motion='FIXED_CLAMP'),
        outer_race=dict(radius=.128,inner_radius=.104,z=[-.065,.065],motion='ROTATING_WITH_FLOWER'))

def main(source,out):
    import numpy as np
    import trimesh
    import pybullet as p
    from PIL import Image,ImageDraw,ImageFont
    from solid_view import extrude
    data=json.loads((source/'preview-input.json').read_text())
    inputs={name:hashlib.sha256((source/name).read_bytes()).hexdigest()
            for name in ['preview-input.json','bottom-preview.obj','middle-preview.obj','top-preview.obj']}
    out.mkdir(parents=True,exist_ok=False)
    parts=assembly_spec();files={}
    def circle(radius):
        a=np.arange(96)*2*np.pi/96
        return np.column_stack([radius*np.cos(a),radius*np.sin(a)]).tolist()
    for name,s in parts.items():
        loops=[circle(s['radius'])]
        if 'inner_radius' in s:loops.append(circle(s['inner_radius'])[::-1])
        mesh=extrude(loops,*s['z']);path=out/(name+'.obj');mesh.export(path,digits=17);files[name]=path
    # A separately labeled seat connects the raster-derived bore to the
    # illustrative circular bearing without changing any flower vertices.
    bore=data['loops'][data['bearing_loop_index']]
    seat=extrude([list(reversed(bore)),list(reversed(circle(.128)))],-.065,.065)
    seat.export(out/'bearing-seat.obj',digits=17);files['bearing_seat']=out/'bearing-seat.obj'
    parts['bearing_seat']=dict(motion='ROTATING_WITH_FLOWER',geometry='SOURCE_BORE_TO_CIRCULAR_RACE_SEAT')
    ball=trimesh.creation.icosphere(subdivisions=2,radius=.013)
    ball.export(out/'rolling-element.obj',digits=17)
    parts['rolling_elements']=dict(motion='ROLLING_RELATIVE_TO_BOTH_RACES',count=10,
        radius=.013,pitch_radius=.09,scope='ILLUSTRATIVE_NOT_CONTACT_VALIDATED')
    client=p.connect(p.DIRECT);bodies={}
    try:
        for name in ('bottom','middle','top'):
            files[name]=source/(name+'-preview.obj')
        for name,path in files.items():
            color=([.13,.43,.76,1] if name in ('upper_cap','lower_cap','shaft','inner_race') else
                   [.85,.65,.22,1] if name in ('outer_race','bearing_seat') else
                   [.78,.31,.57,1] if name=='middle' else [.43,.49,.54,1])
            visual=p.createVisualShape(p.GEOM_MESH,fileName=str(path),rgbaColor=color)
            bodies[name]=p.createMultiBody(baseMass=0,baseVisualShapeIndex=visual,baseCollisionShapeIndex=-1)
        rolling=[]
        for k in range(10):
            t=k*2*np.pi/10
            vis=p.createVisualShape(p.GEOM_MESH,fileName=str(out/'rolling-element.obj'),rgbaColor=[.8,.82,.86,1])
            rolling.append(p.createMultiBody(baseMass=0,baseVisualShapeIndex=vis,baseCollisionShapeIndex=-1,basePosition=[.09*np.cos(t),.09*np.sin(t),0]))
        def render(exploded):
            offsets={'top':.6,'bottom':-.6,'upper_cap':1.,'lower_cap':-1.} if exploded else {}
            for name,body in bodies.items():p.resetBasePositionAndOrientation(body,[0,0,offsets.get(name,0)],[0,0,0,1])
            view=p.computeViewMatrix([.3,-3.8,2.9] if not exploded else [.2,-4.5,3.3],[0,0,0],[0,0,1])
            projection=p.computeProjectionMatrixFOV(38,1,.01,20)
            pixels=p.getCameraImage(850,850,view,projection,renderer=p.ER_TINY_RENDERER,lightDirection=[-2,-3,5])[2]
            return Image.fromarray(np.asarray(pixels,np.uint8).reshape(850,850,4)[:,:,:3])
        assembled,exploded=render(False),render(True)
    finally:p.disconnect(client)
    canvas=Image.new('RGB',(1700,1030),'#f3f5f7');draw=ImageDraw.Draw(canvas)
    font='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
    title=ImageFont.truetype(font,32);small=ImageFont.truetype(font,25)
    canvas.paste(assembled,(0,70));canvas.paste(exploded,(850,70))
    draw.text((30,20),'花盘装配预览：已补上下夹持帽',font=title,fill='#18324b')
    draw.text((880,20),'展开图：夹持件与旋转件分开',font=title,fill='#18324b')
    draw.text((30,930),'蓝色：固定夹持帽、轴与内圈；金色：随花盘转动的轴承外圈和配合座。',font=small,fill='#18324b')
    draw.text((30,970),'参考图重建／尚未优化；轴承及帽为示意尺寸，不是制造认证或物理仿真结果。',font=small,fill='#18324b')
    canvas.save(out/'assembly-preview.png')
    assert all(hashlib.sha256((source/name).read_bytes()).hexdigest()==sha for name,sha in inputs.items())
    (out/'assembly.json').write_text(json.dumps(dict(schema='flower-grip-preview-v1',
        source=str(source),source_sha256=inputs,source_unchanged=True,
        parts=parts,meters_per_model_unit=.02,optimization_runs=0,
        simulation_runs=0,manufacturing_approved=False,
        mechanical_fit='ILLUSTRATIVE_CLEARANCES_NO_PRESS_FIT_OR_CONTACT_CERTIFICATION'),indent=2)+'\n')
    print(out/'assembly-preview.png')

if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
