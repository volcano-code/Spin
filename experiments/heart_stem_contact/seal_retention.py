"""Retain critical cavity cells, reintegrate and seal a separate extension.

An archived actual solve is the candidate source, not a new optimizer run.
Post-retention geometry is not claimed stationary for the paper objective.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
import traceback
import numpy as np
import trimesh
from stem_geometry import ROOT,write_mesh
from heart_top import integrate
from run import write
from topology_retention import expand_cuboids,retain_critical_cells,boundary_mesh,box_moments


def check_archived_configuration(config,prior):
    pairs=[('geometry','source_settings'),('solver','solver_settings'),
           ('boundary','boundary_settings'),('minimum_wall','minimum_outer_wall')]
    if (any(config[a]!=prior[b] for a,b in pairs) or prior['density']!=1. or
            config['meters_per_model_unit']!=.025 or config['density_kg_m3']!=1200.):
        raise ValueError('ARCHIVED_CONFIGURATION_MISMATCH: settings or fixed physical units differ')


def top_constraints(moments,*,phi,root_side):
    s=np.asarray(moments,float)
    if s.shape!=(10,) or not np.isfinite(np.r_[s,phi,root_side]).all() or s[0]<=0 or root_side<=0:
        raise ValueError('invalid mass moments/constraint context')
    c=np.cos(phi); t=np.sin(phi)
    raw=np.array([s[1],s[2],s[8],s[7],c*t*(s[4]-s[5])+(c*c-t*t)*s[9]])
    normalized=raw/np.array([root_side**4]*2+[root_side**5]*3)
    error=float(np.abs(normalized).max())
    if error>1e-6: raise ValueError(f'REPAIRED_CONSTRAINT_FAILED: normalized residual {error}')
    return dict(raw=raw.tolist(),normalized=normalized.tolist(),maximum_normalized_residual=error,
                tolerance=1e-6,phi=phi,scope='TOP_PAPER_EQUALITIES_ONLY_NOT_KKT_OR_DYNAMIC_STABILITY')


def _subtract_mesh(outer,cavity):
    v,f=outer; cv,cf=cavity
    return np.vstack([v,cv]),np.vstack([f,cf[:,::-1]+len(v)])


def seal(source_run,out):
    from make_it_stand.verification.spin_common import native,digest,FRAME,ORDER
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_artifact import SealedSpinArtifact
    from make_it_stand.verification.spin_bundle import save_spin_artifact
    from make_it_stand.verification.source_spin_material import _normalize_model_moments
    source_run=source_run.resolve(); out=out.resolve(); out.mkdir(exist_ok=False)
    start=time.monotonic()
    inputs=[source_run/n for n in ['native-run.json','summary.json','config.json','outer.obj','design.obj',
                                  'fixed.obj','heart.obj','rod.obj']]
    sources=[Path(__file__),Path(__file__).with_name('topology_retention.py'),Path(__file__).with_name('stem_geometry.py'),
             ROOT/'experiments/reconstructed_heart_top/heart_top.py',Path(native().__file__)]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs+sources}
    write(out/'source.json',dict(files=hashes,scope='DECLARED_FILES_AND_LOADED_BINARY_NOT_FULL_BUILD_ATTESTATION'))
    try:
        prior=json.loads((source_run/'native-run.json').read_text()); prior_summary=json.loads((source_run/'summary.json').read_text())
        if digest(prior)!=prior_summary['scientific_digest']: raise ValueError('ARCHIVED_RUN_DIGEST_MISMATCH')
        config=json.loads((source_run/'config.json').read_text()); solve=prior['native_run']['solve']
        check_archived_configuration(config,prior)
        if solve['status']!='BINARY_MOMENT_VALIDATED': raise ValueError('candidate has no binary moment validation')
        if config['solver']['objective']!='TOP_MASS_AND_LOWERING': raise ValueError('unsupported candidate objective')
        for n in ['outer.obj','design.obj','fixed.obj','heart.obj','rod.obj','config.json']:
            shutil.copyfile(source_run/n,out/n)
        def mesh(name):
            m=trimesh.load(source_run/(name+'.obj'),force='mesh',process=False)
            return np.asarray(m.vertices),np.asarray(m.faces)
        outer,design,fixed=mesh('outer'),mesh('design'),mesh('fixed')
        problem=SpinArtifactProblem(original=outer,fixed_regions=[('shaft-and-shell',fixed,1.)],
            design_regions=[('design-excluding-shaft',design,1.)],length_unit='25mm',mass_unit='18.75g',
            target={'mode':'SPINNING_TOP','axis':[0,0,1],'origin':[0,0,0],'support_points':[]},
            minimum_outer_wall=config['minimum_wall'])
        if problem.input_digest!=prior['problem_digest']: raise ValueError('ARCHIVED_INPUT_GEOMETRY_MISMATCH')
        g=config['geometry']; grid=dict(depth=g['maximum_depth'],origin=g['root_origin'],side=g['root_side'])
        original_cells=expand_cuboids(solve['binary_cuboids']['cuboids'],**grid)
        cells,trace=retain_critical_cells(original_cells)
        write(out/'retention-trace.json',dict(profile='CONSERVATIVE_TOPOLOGY_RETENTION_V1',
            authority='IMPLEMENTATION_EXTENSION_NOT_PAPER_ROUNDING',grid=grid,
            source_run=str(source_run),source_scientific_digest=prior_summary['scientific_digest'],
            original_cell_count=len(original_cells),final_cell_count=len(cells),passes=trace,
            original_cells=[list(c) for c in sorted(original_cells)],final_cells=[list(c) for c in sorted(cells)],
            no_additional_carving=cells<=original_cells,optimizer_reexecuted=False,
            stationary_after_retention='NOT_ESTABLISHED',resolution_claim='EXACT_EXPANSION_NOT_NEW_OPTIMIZATION'))
        print('Retention',len(original_cells),'->',len(cells),'cavity cells; exact final validation',flush=True)
        cavity=boundary_mesh(cells,**grid); material=_subtract_mesh(outer,cavity); retained=_subtract_mesh(design,cavity)
        for name,arrays in [('cavity',cavity),('material',material)]:
            write_mesh(out/(name+'.obj'),trimesh.Trimesh(*arrays,process=False))
        original_moments=integrate(*outer)
        analytic=original_moments-box_moments(cells,**grid)
        actual=integrate(*material)
        scales=np.array([g['root_side']**3]+[g['root_side']**4]*3+[g['root_side']**5]*6)
        moment_error=float(np.max(np.abs(actual-analytic)/scales))
        if not np.isfinite(moment_error) or moment_error>1e-10: raise ValueError('REPAIRED_MOMENT_MISMATCH')
        constraints=top_constraints(actual,phi=solve['final_phi'],root_side=g['root_side'])
        claim=_normalize_model_moments(analytic,problem.normalization)
        artifact=SealedSpinArtifact.create(problem,outer=outer,
            design_results=[('design-excluding-shaft',retained,cavity)],final_material=material,newly_carved=cavity,
            solver_claim={'frame':FRAME,'order':list(ORDER),'moments':claim})
        # This API slot is a moment claim, not proof of post-retention optimization.
        # Its provenance is explicit in retention-trace.json and this summary.
        save_spin_artifact(artifact,out/'sealed')
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in hashes.items()):
            raise ValueError('STALE_RESULT')
        write(out/'summary.json',dict(status='SEALED_GEOMETRY_RETENTION_EXTENSION',sealed=True,
            problem_digest=problem.input_digest,scientific_digest=artifact.scientific_digest,
            author_original=False,paper_reproduction=False,physical_validation='NOT_VALIDATED',
            solver_stationarity_after_retention='NOT_ESTABLISHED',claim_provenance='RECOMPUTED_ANALYTIC_CUBOID_SUBTRACTION',
            source_unchanged=True,normalized_moment_error=moment_error,constraints=constraints,
            retained_cell_count=len(original_cells)-len(cells),
            removed_fraction=float(1-actual[0]/original_moments[0]),
            mass_kg=float(actual[0]*config['density_kg_m3']*config['meters_per_model_unit']**3),
            center_m=(actual[1:4]/actual[0]*config['meters_per_model_unit']).tolist(),
            artifact_facts=artifact.facts,elapsed_seconds=time.monotonic()-start))
        print('SEALED_GEOMETRY_RETENTION_EXTENSION',artifact.scientific_digest,flush=True)
        return artifact
    except Exception as error:
        write(out/'failure.json',dict(error=str(error),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        write(out/'summary.json',dict(status='RETENTION_FAILED',sealed=False,error=str(error)))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source_run',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    seal(a.source_run,a.output)
