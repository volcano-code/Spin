"""Source-bound studio rendering, never modifies a solver input or result."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def obj_part(path,expected,role,scale,*,opacity=1.,offset=(0,0,0),smooth=False):
    path=Path(path)
    if sha(path)!=expected:raise ValueError('SOURCE_MESH_IDENTITY_MISMATCH')
    mesh=trimesh.load(path,force='mesh',process=False)
    return dict(name=path.stem,source=str(path.resolve()),source_sha256=expected,role=role,
        vertices=(np.asarray(mesh.vertices)*scale+offset).tolist(),faces=np.asarray(mesh.faces).tolist(),
        opacity=opacity,smooth=smooth,display_only=True)

def rational_part(path,role,scale,opacity):
    import gmpy2
    d=json.loads(path.read_text())
    return dict(name=path.stem,source=str(path),source_sha256=sha(path),role=role,
        vertices=[[float(gmpy2.mpq(x))*scale for x in row] for row in d['vertices']],
        faces=d['faces'],opacity=opacity,smooth=True,display_only=True,
        representation='BINARY64_RENDER_VIEW_OF_AUTHORITY_RATIONAL_NOT_SCIENTIFIC_OUTPUT')

def prepare(out):
    out.mkdir(parents=True,exist_ok=False)
    flower=HERE/'paper-flower-method-v1'
    from flower_result_figure import validate_status_sources
    status=validate_status_sources(flower,saved_material_check=flower/'saved-independent-job-03/output/saved-independent-check.json',
        physics_summary=flower/'physics-job-02/output/summary.json')
    material=Path(status['material_preparation_path']).parent
    record=json.loads((material/'output-preparation.json').read_text())
    parts=[]
    for name in ['top','bottom','fixed-middle-boundary-rim','bearing-seat','outer-race','retained']:
        info=record['files'][name]
        dz=.018 if name=='top' else -.018 if name=='bottom' else 0
        parts.append(obj_part(material/info['file'],info['sha256'],'retained' if name=='retained' else 'fixed',
            1.,offset=(0,0,dz)))
    scenes=[dict(name='flower',status='OPTIMIZED_QUALIFIED_MATERIAL',parts=parts,
        title='Flower — actual optimized material',legend='Grey: fixed layers   Yellow: retained material',
        note='Display-only layer separation; reconstructed method validation',camera=[0,-9,6],
        provenance=status)]
    heart=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-admission-job-01'
    config=json.loads((heart/'config.json').read_text());pin=config['identity']['candidate-pin.obj']
    parts=[rational_part(heart/'output/outer.rational.json','fixed',.025,.22),
        rational_part(heart/'output/design.rational.json','design',.025,1.),
        obj_part(pin['path'],pin['sha256'],'fixed',.025,smooth=True)]
    scenes.append(dict(name='heart',status='INPUT_PREVIEW_NOT_OPTIMIZED',parts=parts,
        title='Handled heart — actual input',legend='Transparent grey: outer envelope   Pink: design domain',
        note='Optimization pending derived-input admission; no optimized cavity shown',camera=[-.8,9,3.0],
        provenance=dict(authority_admission_sha256=sha(heart/'output/admission.json'))))
    # Reuse the existing full-envelope candidate rather than secretly enlarging v4.
    noodles=HERE/'noodles-dense-v3';manifest=json.loads((noodles/'manifest.json').read_text());parts=[]
    for name,info in manifest['parts'].items():
        if name=='outer-ellipsoid':role,opacity='outer-envelope',.045
        elif name=='inner-ellipsoid-candidate':role,opacity='design-candidate',.32
        else:role,opacity='fixed',1.
        parts.append(obj_part(noodles/info['model_obj'],info['sha256']['model'],role,.025,
            opacity=opacity,smooth=True))
    scenes.append(dict(name='noodles',status='INVALID_INPUT_CANDIDATE_PREVIEW_ONLY',parts=parts,
        title='Noodles — full-envelope candidate',legend='Grey: fixed tubes   Transparent pink: candidate design envelope',
        note='Raw source parts, not an admitted Boolean partition; tube self-intersection unresolved',
        camera=[4,-9,4],provenance=dict(manifest=str(noodles/'manifest.json'),manifest_sha256=sha(noodles/'manifest.json'),
        geometry_version='dense-v3',version_choice='EXISTING_FULL_ENVELOPE_REFERENCE_PREVIEW_NOT_NEW_SOLVER_INPUT')))
    bundle=dict(schema='reference-studio-scenes-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',
        optimization_invocations=0,simulation_invocations=0,scenes=scenes)
    (out/'scenes.json').write_text(json.dumps(bundle))
    return out/'scenes.json'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    out=args.output.resolve();scene=prepare(out)
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    blender=Path('/home/lin/.local/bin/blender').resolve()
    worker=HERE/'reference_studio_blender.py'
    result=run_source_job(out/'render-job-01',[str(blender),'-b','-t','4','--python',str(worker),
        '--','--scene',str(scene),'--output','{job}/output'],
        dict(scene=scene,renderer=worker,preparer=Path(__file__),blender=blender),
        maximum_wall_ns=3600_000_000_000,maximum_memory_bytes=24*1024**3)
    (out/'terminal.json').write_text(json.dumps(dict(status=result.status,returncode=result.returncode)))
    print(json.dumps(dict(status=result.status,returncode=result.returncode,output=str(out))),flush=True)
    if result.status!='PROCESS_EXITED' or result.returncode: return 1
    bundle=json.loads(scene.read_text());cards=[]
    for s in bundle['scenes']:
        cards.append(f'<section><h2>{s["title"]}</h2><img src="render-job-01/output/{s["name"]}.png"><p>{s["legend"]}</p><p>{s["status"]}: {s["note"]}</p></section>')
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>三案例 · 参考式三维渲染</title><style>body{font-family:system-ui;background:#eee;color:#222;max-width:1200px;margin:32px auto}section{background:white;padding:24px;margin:24px 0;border-radius:12px}img{width:100%}p{line-height:1.6}</style><h1>参考式三维渲染 · 重建方法验证</h1><p>花盘为真实优化材料。心形、面条为输入预览，不是优化成功或物理仿真。</p>'+''.join(cards))
    return 0

if __name__=='__main__':raise SystemExit(main())
