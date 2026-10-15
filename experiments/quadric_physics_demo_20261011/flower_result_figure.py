"""Report the saved solve without constructing a new optimization result."""
import argparse
import hashlib
import json
from pathlib import Path

def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_status_sources(folder, *, saved_material_check=None, physics_summary=None):
    """Consume actual independent gates, never upgrade a plotting result itself."""
    folder=Path(folder)
    status=dict(physics_result=False,material_result=False)
    if physics_summary is not None and saved_material_check is None:
        raise ValueError('SAVED_MATERIAL_CHECK_REQUIRED_FOR_PHYSICS_STATUS')
    if saved_material_check is None:return status
    from flower_run_physics import admit
    saved=Path(saved_material_check)
    proof=json.loads(saved.read_text())
    preparations=[Path(path) for path in proof.get('sources',{}) if Path(path).name=='output-preparation.json']
    if len(preparations)!=1:raise ValueError('SAVED_CHECK_PREPARATION_IDENTITY_REQUIRED')
    preparation=preparations[0]
    record,audit=admit(preparation.parent,saved.parent,folder,saved_material_check=saved)
    status.update(material_result=True,saved_material_check_path=str(saved.resolve()),
        saved_material_check_sha256=_sha(saved),material_preparation_path=str(preparation.resolve()),
        material_preparation_sha256=_sha(preparation),
        independent_saved_balance=audit['independent_saved_balance'])
    if physics_summary is None:return status
    physics_path=Path(physics_summary)
    physics=json.loads(physics_path.read_text())
    sources=physics.get('sources',{})
    if sources.get('material_preparation_sha256')!=_sha(preparation) or \
            sources.get('saved_material_check_sha256')!=_sha(saved) or \
            Path(sources.get('saved_material_check_path','')).resolve()!=saved.resolve() or \
            physics.get('original')!=record['original'] or physics.get('optimized')!=record['optimized']:
        raise ValueError('PHYSICS_MATERIAL_IDENTITY_MISMATCH')
    if physics.get('schema')!='admitted-flower-actual-physics-matrix-v1' or \
            physics.get('models')!=['ideal','flexible'] or physics.get('phase_count')!=8 or \
            physics.get('time_steps_seconds')!=[.0005,.00025] or physics.get('duration_seconds')!=8. or \
            physics.get('expected_scenes')!=64 or physics.get('recorded_scenes')!=64 or \
            physics.get('motor_drive') is not False:
        raise ValueError('PHYSICS_SCENE_PROTOCOL_INCOMPLETE')
    expected={(model,variant,phase,dt) for model in ('ideal','flexible')
              for variant in ('original','optimized') for phase in range(8) for dt in (.0005,.00025)}
    rows=physics.get('scenes',[])
    identities=[(r.get('model'),r.get('variant'),r.get('phase_index'),r.get('dt_seconds')) for r in rows]
    if len(rows)!=64 or len(set(identities))!=64 or set(identities)!=expected:
        raise ValueError('PHYSICS_SCENE_MATRIX_INCOMPLETE')
    for row in rows:
        identity=f"{row['model']}-{row['variant']}-phase{row['phase_index']}-dt{row['dt_seconds']}"
        if row.get('id')!=identity or row.get('actual_steps')!=round(8/row['dt_seconds']):
            raise ValueError('PHYSICS_SCENE_IDENTITY_MISMATCH')
        trajectory=physics_path.parent/(identity+'.npz')
        if trajectory.is_symlink() or not trajectory.is_file() or _sha(trajectory)!=row.get('trajectory_sha256'):
            raise ValueError('PHYSICS_TRAJECTORY_IDENTITY_MISMATCH')
    status.update(physics_result=True,physics_summary_path=str(physics_path.resolve()),
        physics_summary_sha256=_sha(physics_path),recorded_scenes=len(rows),
        support_scope='IDEAL_BEARING_AND_FLEXIBLE_SUPPORT_ENGINEERING_MODELS_NOT_CALIBRATED_HAND')
    return status


def render(folder, output, *, saved_material_check=None, physics_summary=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Polygon as Patch
    from shapely.geometry import Polygon, Point
    from paper_flower_output import profile_geometry
    folder=Path(folder); output=Path(output)
    status=validate_status_sources(folder,saved_material_check=saved_material_check,physics_summary=physics_summary)
    source=folder/'solve-job-01/output/application.json'
    case=json.loads(source.read_text())['cases'][0]
    frozen=json.loads((folder/'input.json').read_text())
    successful=[a for a in case['attempts'] if a['seed']==0 and a['status']=='CONVERGED']
    if len(successful)!=1: raise ValueError('PREDECLARED_SEED_SUCCESS_REQUIRED')
    solved=successful[0]['solver_result']
    if solved['accepted_steps']<=0: raise ValueError('NO_ACCEPTED_OPTIMIZATION')
    original=Polygon(frozen['retained_original_loops'][0],frozen['retained_original_loops'][1:])
    design=profile_geometry(case['planar_problem']['design_region_loops'])
    cut=Point(solved['center']).buffer(solved['radius'],resolution=1024)
    retained=design.difference(cut); removed=design.intersection(cut)
    def draw(ax,geometry,color):
        pieces=[geometry] if geometry.geom_type=='Polygon' else geometry.geoms
        for p in pieces:
            ax.add_patch(Patch(p.exterior.coords,facecolor=color,edgecolor='#324152',linewidth=.45))
            for hole in p.interiors: ax.add_patch(Patch(hole.coords,facecolor='white',edgecolor='#324152',linewidth=.35))
    fig,axes=plt.subplots(1,3,figsize=(14,6))
    fig.subplots_adjust(top=.75,bottom=.12,left=.05,right=.99,wspace=.28)
    draw(axes[0],original,'#e2b7bd')
    draw(axes[1],original,'#dddfe5');draw(axes[1],retained,'#ddb342');draw(axes[1],removed,'#b9dcef')
    axes[1].add_patch(Circle(solved['center'],solved['radius'],fill=False,edgecolor='#1774a4',lw=1.8))
    draw(axes[2],removed,'#b9dcef')
    for ax,title in zip(axes,['Frozen three-layer flower outline','Actual middle-layer solve: retained / removed','Removed material (not preset)']):
        ax.set_aspect('equal');ax.set_xlim(-1.3,1.3);ax.set_ylim(-1.3,1.3);ax.set_title(title,fontsize=11)
        ax.axhline(0,color='#aeb7c1',lw=.4);ax.axvline(0,color='#aeb7c1',lw=.4)
        ax.plot(0,0,'k+',ms=9);ax.set_xlabel('Model x');ax.set_ylabel('Model y')
    title=('Spin-It Faster 2D disk / dogleg — RECONSTRUCTED_METHOD_VALIDATION\n'
        f"48 original holes; accepted steps={solved['accepted_steps']}; residual={successful[0]['normalized_equation_residual']:.3g}\n"
        'Displayed circle tessellation is illustrative. OBJ exact partition and case physics are not yet approved.')
    if status['material_result']:
        title=('Spin-It Faster 2D disk / dogleg — RECONSTRUCTED_METHOD_VALIDATION\n'
            f"Saved OBJ material gate: PASSED; accepted steps={solved['accepted_steps']}; residual={successful[0]['normalized_equation_residual']:.3g}\n"
            'Target: planar first-moment balance at bearing origin; circle display tessellation only.\n'
            +('64 actual scenes: ideal bearing + flexible support (engineering models, not calibrated hand).'
              if status['physics_result'] else 'Physics not supplied or approved in this figure.'))
    fig.suptitle(title,fontsize=11,y=.97)
    output.mkdir(parents=True,exist_ok=False)
    fig.savefig(output/'actual-optimization.png',dpi=160);plt.close(fig)
    (output/'source.json').write_text(json.dumps(dict(application_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        seed=0,run=successful[0]['run'],restart=successful[0]['restart'],center=solved['center'],radius=solved['radius'],
        display_only_circle_tessellation=True,**status),indent=2)+'\n')
    ledger=solved.get('evaluation_ledger',{}).get('solver_evaluations',[])
    values=[row for row in ledger if row.get('has_residual') is True]
    if values:
        fig,ax=plt.subplots(figsize=(9,4))
        ax.semilogy([row['evaluation_index'] for row in values],
                    [row['normalized_residual'] for row in values],'.',color='#607d8b',label='All residual evaluations / probes')
        accepted=[row for row in values if row.get('accepted_trial') is True]
        ax.semilogy([row['evaluation_index'] for row in accepted],
                    [row['normalized_residual'] for row in accepted],'o',color='#bd7832',label='Accepted trials')
        ax.set(xlabel='Actual solver evaluation index',ylabel='Normalized residual',
               title='Stored dogleg evaluation ledger — probes are not accepted steps')
        ax.legend();fig.tight_layout();fig.savefig(output/'actual-evaluation-ledger.png',dpi=160);plt.close(fig)
    (output/'evaluation-ledger.json').write_text(json.dumps(dict(
        application_sha256=_sha(source),solver_evaluations=ledger,
        residual_trace_available=bool(values)),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--saved-material-check',type=Path);p.add_argument('--physics-summary',type=Path)
    a=p.parse_args();render(a.folder,a.output,saved_material_check=a.saved_material_check,physics_summary=a.physics_summary)
