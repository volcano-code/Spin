"""Seal and simulate the original solid; never substitute it for an optimization."""
import argparse
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import trimesh
from teapot import ROOT, STRATEGIES, write, file_hashes, check_files
from admission import check_collider_context
sys.path.insert(0,str(ROOT/'experiments/heart_stem_contact'))
from bullet_contact import material_properties, simulate_contact
from heart_top import integrate


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--collisions',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);out=a.output;start=time.monotonic()
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_artifact import SealedSpinArtifact
    from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
    from make_it_stand.verification.spin_common import native,FRAME,ORDER
    from make_it_stand.verification.source_spin_material import _normalize_model_moments
    meta=json.loads((a.input/'manifest.json').read_text());collider=json.loads((a.collisions/'manifest.json').read_text())
    hashes=file_hashes([Path(__file__),Path(__file__).with_name('admission.py'),Path(native().__file__),
        ROOT/'python/make_it_stand/verification/mass_moments.py',
        ROOT/'python/make_it_stand/verification/spin_artifact.py',
        ROOT/'experiments/heart_stem_contact/bullet_contact.py',ROOT/'experiments/reconstructed_heart_top/heart_top.py',
        a.input/'manifest.json',a.collisions/'manifest.json'])
    hashes.update(meta['geometry_hashes']);hashes.update(meta['source_files'])
    hashes.update(collider['outer_hash']);hashes.update(collider['source_hashes'])
    for profile in collider['profiles'].values():hashes.update(profile['hashes'])
    write(out/'source.json',dict(files=hashes))
    try:
        check_files(hashes);check_collider_context(collider,meta['problem_digest'])
        def arrays(name):
            m=trimesh.load(a.input/(name+'.obj'),force='mesh',process=False)
            return np.asarray(m.vertices),np.asarray(m.faces)
        outer,design,fixed=map(arrays,['outer','design','fixed'])
        problem=SpinArtifactProblem(original=outer,fixed_regions=[('fixed-shell-and-shaft',fixed,1.)],
            design_regions=[('body-interior',design,1.)],length_unit='25mm',mass_unit='18.75g',
            target={'mode':'SPINNING_TOP','axis':[0,0,1],'origin':[0,0,0],'support_points':[]},minimum_outer_wall=.015)
        if problem.input_digest!=meta['problem_digest']:raise ValueError('INPUT_MISMATCH')
        s=integrate(*outer)
        artifact=SealedSpinArtifact.create(problem,outer=outer,design_results=[('body-interior',design,None)],
            final_material=outer,newly_carved=None,
            solver_claim={'frame':FRAME,'order':list(ORDER),'moments':_normalize_model_moments(s,problem.normalization)})
        print('BASELINE SEAL CREATED; NO OPTIMIZATION CLAIM',flush=True)
        save_spin_artifact(artifact,out/'sealed')
        artifact=load_spin_artifact(out/'sealed',expected_problem_digest=problem.input_digest)
        props=material_properties(outer[0]*.025,outer[1],1200.)
        mass,center,Ic,d,U=props
        model=dict(label='original',mass_kg=mass,center_m=center.tolist(),central_inertia=Ic.tolist(),
            artifact_digest=artifact.scientific_digest)
        rows=[];trajectory_hashes={}
        for hz in [35.,45.]:
            for profile,proxy in collider['profiles'].items():
                name=f'original-{hz:g}hz-{profile}-dir0';print('BASELINE SIM',name,flush=True)
                r=simulate_contact(Path(collider['outer']),proxy['files'],props,outer[0]*.025,
                    hz=hz,tilt=2.,azimuth=0.,duration=8.,dt=1/2000)
                r['teapot_proxy_contact_occurred']=r.pop('heart_proxy_contact_occurred')
                r['collision_scope']='zero-margin shaft base; padded V-HACD teapot links; fallen contact approximate'
                write(out/(name+'.json'),r);trajectory_hashes.update(file_hashes([out/(name+'.json')]))
                rows.append(dict(name=name,label='original',hz=hz,profile=profile,direction=0,status='RECORDED',
                    path=str((out/(name+'.json')).resolve()),max_tilt=r['maximum_tilt_degrees'],
                    final_tilt=r['final_tilt_degrees'],body_contact=r['teapot_proxy_contact_occurred']))
        hashes.update(file_hashes([p for p in (out/'sealed').iterdir() if p.is_file()]))
        check_files(hashes);write(out/'source.json',dict(files=hashes))
        write(out/'report.json',dict(scope='ORIGINAL_SOLID_BASELINE_ONLY_NOT_OPTIMIZED_RESULT',protocol='contact',
            problem_digest=problem.input_digest,models=[model],runs=rows,trajectory_hashes=trajectory_hashes,
            cases=[dict(strategy=s,status='NOT_ADMITTED_BASELINE_ONLY') for s in STRATEGIES],
            source_unchanged=True,all_scenarios_repeated=False,refinement_directions=[],
            validation_status='OBSERVATIONAL_NOT_VALIDATED',elapsed_seconds=time.monotonic()-start))
    except Exception as e:
        write(out/'failure.json',dict(error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        raise


if __name__=='__main__':main()
