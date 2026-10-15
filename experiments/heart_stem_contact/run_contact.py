"""Original-only observations remain possible when optimized geometry is blocked."""
import argparse,json,hashlib,itertools,sys,time
from pathlib import Path
import numpy as np
import trimesh
from stem_geometry import ROOT
from bullet_contact import material_properties,simulate_contact
from run import write


def main():
    p=argparse.ArgumentParser(); p.add_argument('run',type=Path)
    p.add_argument('--output-name',default='pybullet-original'); p.add_argument('--model',choices=['original','carved'],default='original')
    p.add_argument('--collision-dir',default='collision-tip-base')
    p.add_argument('--gui-hz',type=float,default=35)
    p.add_argument('--gui-tilt',type=float,default=2)
    p.add_argument('--gui',action='store_true'); a=p.parse_args()
    out=a.run/a.output_name; out.mkdir(exist_ok=False)
    config=json.loads((a.run/'config.json').read_text()); manifest=json.loads((a.run/a.collision_dir/'manifest.json').read_text())
    if hashlib.sha256((a.run/'outer.obj').read_bytes()).hexdigest()!=manifest['source_outer_sha256']:
        raise ValueError('STALE_COLLIDER_SOURCE')
    if config['meters_per_model_unit']!=manifest['meters_per_model_unit']:
        raise ValueError('COLLIDER_SCALE_MISMATCH')
    if hashlib.sha256(Path(manifest['outer_file']).read_bytes()).hexdigest()!=manifest['outer_sha256']:
        raise ValueError('STALE_OUTER_VISUAL')
    for profile in manifest['profiles'].values():
        if not profile['near_upright_support_passed']: raise ValueError('COLLISION_MODEL_BLOCKER')
        if set(profile['files'])!=set(profile['mesh_sha256']): raise ValueError('UNBOUND_COLLIDER')
        for f,h in profile['mesh_sha256'].items():
            if hashlib.sha256(Path(f).read_bytes()).hexdigest()!=h: raise ValueError('STALE_COLLIDER')
    from make_it_stand.verification.spin_common import native
    n=native(); mesh=trimesh.load(a.run/'outer.obj',force='mesh',process=False)
    original=(np.asarray(mesh.vertices),np.asarray(mesh.faces))
    displayed=trimesh.load(manifest['outer_file'],force='mesh',process=False)
    if not (np.array_equal(displayed.faces,original[1]) and
            np.array_equal(displayed.vertices,original[0]*config['meters_per_model_unit'])):
        raise ValueError('COLLIDER_VISUAL_SOURCE_GEOMETRY_MISMATCH')
    if a.model=='original':
        material=original
        facts=n.make_exact_material_region(*original).facts
        if not facts['completed'] or not facts['valid']: raise ValueError('invalid input material')
        write(out/'material-admission.json',dict(scope='ORIGINAL_INPUT_NOT_OPTIMIZED',facts=facts))
    else:
        from make_it_stand.verification.spin_bundle import load_spin_artifact
        summary=json.loads((a.run/'summary.json').read_text())
        if not summary['sealed']: raise ValueError('BLOCKED_UNSEALED_CARVING: no raw mesh fallback')
        artifact=load_spin_artifact(a.run/'sealed',expected_problem_digest=summary['problem_digest'])
        if not all(np.array_equal(x,y) for x,y in zip(artifact.outer.arrays(),original)):
            raise ValueError('ARTIFACT_OUTER_COLLIDER_MISMATCH')
        material=artifact.final_material.arrays()
        write(out/'material-admission.json',dict(artifact_digest=artifact.scientific_digest))
    props=material_properties(material[0]*config['meters_per_model_unit'],material[1],config['density_kg_m3'])
    m,c,I,d,U=props
    write(out/'mass.json',dict(mass_kg=m,center_m=c.tolist(),central_inertia=I.tolist(),principal_moments=d.tolist(),
        principal_axes=U.tolist(),input_sha256=hashlib.sha256((a.run/'outer.obj').read_bytes()).hexdigest()))
    outer_si=original[0]*config['meters_per_model_unit']
    sources=[Path(__file__),Path(__file__).with_name('bullet_contact.py'),Path(n.__file__),
        ROOT/'experiments/reconstructed_heart_top/heart_top.py',a.run/a.collision_dir/'manifest.json',
        a.run/'config.json',a.run/'outer.obj']
    import pybullet
    sources.append(Path(pybullet.__file__))
    if a.gui: sources.append(Path(__file__).with_name('capture.py'))
    before={str(f.resolve()):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
    write(out/'source.json',before)
    def execute(label,*args,**kwargs):
        try:
            return simulate_contact(*args,**kwargs)
        except Exception as error:
            write(out/(label+'-failed.json'),dict(status='FAILED_NOT_A_STABILITY_RESULT',error=str(error)))
            raise
    if a.gui:
        profile=manifest['profiles']['robust']
        from capture import GuiCapture
        with GuiCapture(out) as capture:
            r=execute('gui',manifest['outer_file'],profile['files'],props,outer_si,
                hz=a.gui_hz,tilt=a.gui_tilt,azimuth=0,gui=True,frame_callback=capture)
        write(out/'gui-run.json',r)
        return
    # Whole predefined matrix. No selecting only a favorable azimuth.
    plan=list(itertools.product(['default','robust'],[6.,35.,45.],[-1,*range(8)]))
    write(out/'plan.json',dict(scenarios=[dict(preset=p,hz=h,direction=d) for p,h,d in plan],
        model=a.model,duration_s=8,dt=1/2000,dt_refinement=1/4000,
        scope='RECONSTRUCTED_FREE_CONTACT_UNCALIBRATED_FRICTION',phenomenological_damping=0,
        sliding_friction=.8,rolling_friction=0,spinning_friction=0,
        additional_dissipative_defaults='nominal35/45Hz both proxies: rolling/spinning=.001, damping=.04',
        coordinate_velocity_limit_rad_s=1000))
    rows=[]; repeats=[]; refinements=[]
    for preset,hz,direction in plan:
        profile=manifest['profiles'][preset]
        if not profile['near_upright_support_passed']: raise ValueError('COLLISION_MODEL_BLOCKER')
        for f,h in profile['mesh_sha256'].items():
            if hashlib.sha256(Path(f).read_bytes()).hexdigest()!=h: raise ValueError('STALE_COLLIDER')
        label=f'{preset}-{hz:g}hz-dir{direction}'
        print(label,flush=True)
        kw=dict(hz=hz,tilt=0 if direction<0 else 2,azimuth=0 if direction<0 else direction*np.pi/4)
        args=(manifest['outer_file'],profile['files'],props,outer_si)
        r=execute(label,*args,**kw); write(out/(label+'.json'),r)
        rows.append(dict(preset=preset,direction=direction,**{k:v for k,v in r.items() if k not in ['samples','elapsed']}))
        if direction==-1 and hz in [35,45]:
            repeated=execute(label+'-repeat',*args,**kw); write(out/(label+'-repeat.json'),repeated)
            eq=all(r[k]==repeated[k] for k in r if k!='elapsed')
            repeats.append(dict(case=label,scientific_bytes_equal=eq))
            fine=execute(label+'-refined',*args,**kw,dt=1/4000); write(out/(label+'-refined.json'),fine)
            error=float(np.max(np.abs(np.array([s['tilt_degrees'] for s in r['samples']])-
                                      np.array([s['tilt_degrees'] for s in fine['samples']]))))
            refinements.append(dict(case=label,maximum_angle_difference_degrees=error,
                base_maximum_tilt=r['maximum_tilt_degrees'],refined_maximum_tilt=fine['maximum_tilt_degrees']))
            dissipated=execute(label+'-dissipative-defaults',*args,**kw,dissipative_defaults=True)
            write(out/(label+'-dissipative-defaults.json'),dissipated)
    unchanged=all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in before.items())
    write(out/'report.json',dict(model=a.model,scenarios=rows,repeats=repeats,refinements=refinements,
        source_unchanged=unchanged,paper_reproduction=False,physical_print_validation=False,
        carving_sealed=json.loads((a.run/'summary.json').read_text())['sealed']))
    print('Recorded',len(rows),'unique observational contact scenarios; NOT a validated stability claim',flush=True)


if __name__=='__main__': main()
