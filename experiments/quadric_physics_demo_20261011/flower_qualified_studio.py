"""Studio views bound to the qualified v2 flower material and actual physics."""
import argparse
import json
from pathlib import Path
import sys
from reference_studio import obj_part,sha,HERE,ROOT


def main():
    parser=argparse.ArgumentParser()
    for name in ('case','saved-check','physics','output'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    for name in ('case','saved_check','physics','output'):
        setattr(args,name,getattr(args,name).resolve())
    from flower_result_figure import validate_status_sources
    status=validate_status_sources(args.case,saved_material_check=args.saved_check,physics_summary=args.physics/'summary.json')
    material=Path(status['material_preparation_path']).parent
    record=json.loads((material/'output-preparation.json').read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    scenes=[]
    for mode in ('input','cutaway','exploded','assembly'):
        parts=[]
        primary='design' if mode=='input' else 'retained'
        info=record['files'][primary]
        parts.append(obj_part(material/info['file'],info['sha256'],'design' if mode=='input' else 'retained',1.))
        for name in record['rotor_layers']:
            if mode=='cutaway' and name in ('top','bottom'):continue
            info=record['files'][name]
            dz=.016 if mode=='exploded' and name=='top' else -.016 if mode=='exploded' and name=='bottom' else 0
            parts.append(obj_part(material/info['file'],info['sha256'],'fixed',1.,opacity=.12 if mode=='input' else 1.,offset=(0,0,dz)))
        if mode=='assembly':
            for name in record['clamp_layers']:
                info=record['files'][name]
                parts.append(obj_part(material/info['file'],info['sha256'],'fixed',1.))
        scenes.append(dict(name=mode,status='OPTIMIZED_QUALIFIED_MATERIAL',parts=parts,
            camera=[0,-9,8],provenance=status))
    scene=args.output/'scenes.json'
    scene.write_text(json.dumps(dict(schema='qualified-flower-studio-scenes-v1',scenes=scenes,
        classification='RECONSTRUCTED_METHOD_VALIDATION',display_only=True)))
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    worker=HERE/'reference_studio_blender.py';blender=Path('/home/lin/.local/bin/blender').resolve()
    files=dict(preparer=Path(__file__),helper=HERE/'reference_studio.py',renderer=worker,
        validator=HERE/'flower_result_figure.py',physics_guard=HERE/'flower_run_physics.py',
        scene=scene,proof=args.saved_check,physics=args.physics/'summary.json',
        material=material/'output-preparation.json',blender=blender)
    for name,info in record['files'].items():files['mesh.'+name]=material/info['file']
    result=run_source_job(args.output/'render-job-01',[str(blender),'-b','-t','4','--python',str(worker),
        '--','--scene',str(scene),'--output','{job}/output'],files)
    if result.status=='PROCESS_EXITED' and result.returncode==0:
        rendered=args.output/'render-job-01/output'
        manifest=json.loads((rendered/'render-manifest.json').read_text())
        if manifest['scene_sha256']!=sha(scene) or set(manifest['images'])!=set(x['name'] for x in scenes):
            raise ValueError('STUDIO_OUTPUT_IDENTITY_MISMATCH')
        for name,digest in manifest['images'].items():
            if sha(rendered/(name+'.png'))!=digest:raise ValueError('STUDIO_IMAGE_MISSING_OR_CHANGED')
    print(json.dumps(dict(status=result.status,returncode=result.returncode)),flush=True)
    return 0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1


if __name__=='__main__':raise SystemExit(main())
