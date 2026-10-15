"""Read checked sealed materials; retain all requested trajectories and failures."""
import argparse
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import trimesh
from teapot import ROOT, STRATEGIES, write, file_hashes, check_files
from admission import check_case_context, check_collider_context
sys.path.insert(0,str(ROOT/'experiments/reconstructed_heart_top'))
from heart_top import integrate,properties,simulate
sys.path.insert(0,str(ROOT/'experiments/heart_stem_contact'))
from bullet_contact import material_properties,simulate_contact


def main():
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--collisions',type=Path)
    p.add_argument('--protocol',choices=['ideal','contact'],required=True);a=p.parse_args()
    dest=a.output;dest.mkdir(exist_ok=False);start=time.monotonic()
    from make_it_stand.verification.spin_bundle import load_spin_artifact
    source_paths=[Path(__file__),Path(__file__).with_name('admission.py'),ROOT/'experiments/reconstructed_heart_top/heart_top.py',
        ROOT/'experiments/heart_stem_contact/bullet_contact.py']
    inherited={};common_digest=None;common_settings=None
    models=[];case_status=[];original_added=False
    for strategy in STRATEGIES:
        folder=a.campaign/strategy
        if not (folder/'summary.json').is_file():
            case_status.append(dict(strategy=strategy,status='MISSING_RESULT'));continue
        result=json.loads((folder/'summary.json').read_text())
        upstream=json.loads((folder/'source.json').read_text())['files'];check_files(upstream);inherited.update(upstream)
        native_run=json.loads((folder/'native-run.json').read_text());config=json.loads((folder/'config.json').read_text())
        if common_digest is None:common_digest=result['problem_digest']
        check_case_context(native_run,config,result,strategy,common_digest)
        settings=dict(solver={k:v for k,v in config['solver'].items() if k!='objective'},
            geometry=config['geometry'],boundary=config['boundary'])
        if common_settings is None:common_settings=settings
        if settings!=common_settings:raise ValueError('STRATEGY_NUMERICAL_SETTINGS_DIFFER')
        source_paths.extend([folder/'summary.json',folder/'source.json',folder/'native-run.json',folder/'config.json'])
        if not result['sealed']:
            case_status.append(dict(strategy=strategy,status='UNSEALED_NOT_ADMITTED',solver_status=result['status']));continue
        artifact=load_spin_artifact(folder/'sealed',expected_problem_digest=result['problem_digest'])
        if artifact.scientific_digest!=native_run['artifact_scientific_digest']:raise ValueError('ARTIFACT_SOLVER_MISMATCH')
        source_paths.extend(p for p in (folder/'sealed').iterdir() if p.is_file())
        if not original_added:
            models.append(('original',artifact.problem.original.arrays(),artifact.problem.original.arrays(),artifact.scientific_digest))
            original_added=True
        models.append((strategy,artifact.final_material.arrays(),artifact.outer.arrays(),artifact.scientific_digest))
        case_status.append(dict(strategy=strategy,status='GEOMETRY_SEALED',artifact=artifact.scientific_digest))
        cv,cf=artifact.newly_carved.arrays()
        np.savez(dest/(strategy+'-geometry.npz'),outer_v=artifact.outer.arrays()[0],outer_f=artifact.outer.arrays()[1],
            material_v=artifact.final_material.arrays()[0],material_f=artifact.final_material.arrays()[1],cavity_v=cv,cavity_f=cf)
    if not models:raise ValueError('NO_SEALED_MATERIAL: no optimized physics or video admitted')
    sources=file_hashes(source_paths);sources.update(inherited)
    collider=None
    if a.protocol=='contact':
        if a.collisions is None:raise ValueError('contact needs collider manifest')
        collider=json.loads((a.collisions/'manifest.json').read_text())
        check_collider_context(collider,common_digest)
        sources.update(file_hashes([a.collisions/'manifest.json']))
        sources.update(collider['outer_hash']);sources.update(collider['source_hashes'])
        for profile in collider['profiles'].values():sources.update(profile['hashes'])
    write(dest/'source.json',dict(files=sources))
    rows=[];model_records=[];trajectory_hashes={}
    for label,material,outer,artifact_digest in models:
        s=integrate(*material);scale=.025;density=1200.
        si=s*density*scale**np.array([3,4,4,4,5,5,5,5,5,5])
        mass,center,pivot=properties(si)
        model_records.append(dict(label=label,mass_kg=mass,center_m=center.tolist(),pivot_inertia=pivot.tolist(),
            geometric_moments=s.tolist(),artifact_digest=artifact_digest))
        props=material_properties(material[0]*scale,material[1],density)
        for hz in [10.,20.,35.,45.]:
            # Contact matrix focuses on the same two higher speeds without
            # claiming that omitted speeds were simulated in this protocol.
            if a.protocol=='contact' and hz not in [35.,45.]:continue
            profiles=['ideal'] if collider is None else list(collider['profiles'])
            for profile in profiles:
                for direction in range(8):
                    name=f'{label}-{hz:g}hz-{profile}-dir{direction}'
                    print('SIM',name,flush=True)
                    try:
                        azimuth=direction*np.pi/4
                        if collider is None:
                            r=simulate(mass,center,pivot,[0,0,2*np.pi*hz],duration=8.,tilt_degrees=2.,azimuth=azimuth,max_step=1/480)
                            r['protocol']='IDEAL_FIXED_PIVOT_GRAVITY_NOT_TABLE_CONTACT'
                        else:
                            r=simulate_contact(Path(collider['outer']),collider['profiles'][profile]['files'],props,
                                outer[0]*scale,hz=hz,tilt=2.,azimuth=azimuth,duration=8.,dt=1/2000)
                            r['teapot_proxy_contact_occurred']=r.pop('heart_proxy_contact_occurred')
                            r['collision_scope']='zero-margin shaft base; padded V-HACD teapot links; fallen contact approximate'
                        write(dest/(name+'.json'),r)
                        trajectory_hashes.update(file_hashes([dest/(name+'.json')]))
                        row=dict(name=name,label=label,hz=hz,profile=profile,direction=direction,status='RECORDED',
                            path=str((dest/(name+'.json')).resolve()),problem_digest=common_digest,artifact_digest=artifact_digest,
                            max_tilt=r.get('max_tilt_degrees',r.get('maximum_tilt_degrees')),
                            final_tilt=r['final_tilt_degrees'])
                        if direction==0:
                            if collider is None:
                                fine=simulate(mass,center,pivot,[0,0,2*np.pi*hz],duration=8.,tilt_degrees=2.,azimuth=azimuth,max_step=1/4000)
                                error=float(np.max(np.abs(np.array(r['tilt_degrees'])-fine['tilt_degrees'])))
                            else:
                                fine=simulate_contact(Path(collider['outer']),collider['profiles'][profile]['files'],props,
                                    outer[0]*scale,hz=hz,tilt=2.,azimuth=azimuth,duration=8.,dt=1/4000)
                                error=float(np.max(np.abs(np.array([x['tilt_degrees'] for x in r['samples']])-
                                    np.interp([x['time'] for x in r['samples']],[x['time'] for x in fine['samples']],
                                        [x['tilt_degrees'] for x in fine['samples']]))))
                            write(dest/(name+'-refined.json'),fine)
                            trajectory_hashes.update(file_hashes([dest/(name+'-refined.json')]))
                            row.update(refinement_max_tilt_difference=error,refinement_passed=error<=.001)
                        rows.append(row)
                    except Exception as e:
                        failure=dict(name=name,status='FAILED',error=str(e),traceback=traceback.format_exc())
                        write(dest/(name+'-failure.json'),failure);rows.append(failure)
                    write(dest/'progress.json',dict(cases=case_status,models=model_records,runs=rows))
    check_files(sources)
    write(dest/'report.json',dict(scope='RECONSTRUCTED_TEAPOT_NOT_ORIGINAL_PAPER_EXPERIMENT',
        protocol=a.protocol,cases=case_status,models=model_records,runs=rows,source_unchanged=True,
        problem_digest=common_digest,trajectory_hashes=trajectory_hashes,
        validation_status='OBSERVATIONAL_NOT_FULLY_VALIDATED',elapsed_seconds=time.monotonic()-start,
        refinement_directions=[0],all_scenarios_repeated=False,physical_print_validation=False))


if __name__=='__main__':main()
