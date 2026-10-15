"""Independent algebra on actual saved fill/leaf moments, not geometry proof.

Top constraints are linear in moments at the *fixed final* Givens angle.
The decomposition does not reoptimize phi, change rounding, or admit a result.
"""
import math
import numpy as np


def rounding_diagnostics(record):
    solve=record['solve']
    if solve['status'] not in ('BINARY_MOMENT_VALIDATED','BINARY_VALIDATION_FAILED'):
        return dict(available=False,status=solve['status'],reason='NO_COMPLETED_BINARY_COMPARISON')
    settings=record['settings']
    if settings['solver_settings']['objective'] not in (
            'TOP_INERTIA_ONLY','TOP_LOWERING_ONLY','TOP_MASS_AND_LOWERING'):
        raise ValueError('TOP_ONLY_DIAGNOSTIC')
    beta=np.asarray(solve['relaxed_fill'],dtype=float)
    binary=np.asarray(solve['binary_fill'],dtype=float)
    leaves=solve['final_partition']
    moments=np.asarray([c['moments'] for c in leaves],dtype=float)
    outer=np.asarray(record['source']['outer_surface_moments'],dtype=float)
    claimed=np.asarray(solve['binary_moments'],dtype=float)
    if (beta.ndim!=1 or len(beta)==0 or binary.shape!=beta.shape or
            moments.shape!=(len(beta),10) or outer.shape!=(10,) or claimed.shape!=(10,)):
        raise ValueError('INCOMPLETE_ROUNDING_STATE')
    if not all(np.isfinite(v).all() for v in (beta,binary,moments,outer,claimed)):
        raise ValueError('NONFINITE_ROUNDING_STATE')
    if np.any((beta<0)|(beta>1)):
        raise ValueError('INVALID_FILL')
    if not np.array_equal(binary,(beta>=.5).astype(float)):
        raise ValueError('ROUNDING_MISMATCH')
    side=float(settings['source_settings']['root_side']);density=float(settings['density'])
    phi=float(solve['final_phi'])
    if not all(math.isfinite(v) for v in (side,density,phi)) or min(side,density)<=0:
        raise ValueError('INVALID_NORMALIZATION')
    scale=density*np.array([side**3]+[side**4]*3+[side**5]*6)
    if not np.isfinite(scale).all() or np.any(scale<=0):raise ValueError('INVALID_NORMALIZATION')
    weighted=moments*density
    relaxed=outer*density-beta@weighted
    actual=outer*density-binary@weighted
    if np.max(np.abs(actual-claimed)/scale)>1e-10:
        raise ValueError('MOMENT_MISMATCH')
    # Order is M, X,Y,Z, XX,YY,ZZ,YZ,XZ,XY; normalization is input-fixed.
    transform=np.zeros((5,10))
    transform[0,1]=transform[1,2]=transform[2,8]=transform[3,7]=1.
    c,s=math.cos(phi),math.sin(phi)
    transform[4,4]=c*s;transform[4,5]=-c*s;transform[4,9]=c*c-s*s
    constraint_scale=density*np.array([side**4]*2+[side**5]*3)
    per_cell=-(binary-beta)[:,None]*weighted
    per_constraint=(per_cell@transform.T)/constraint_scale
    rc=(transform@relaxed)/constraint_scale;bc=(transform@actual)/constraint_scale
    delta=per_constraint.sum(axis=0)
    if not np.allclose(delta,bc-rc,rtol=0,atol=1e-12):
        raise ValueError('DECOMPOSITION_MISMATCH')
    by_level={}
    for level in sorted({cell['level'] for cell in leaves}):
        mask=np.array([cell['level']==level for cell in leaves])
        by_level[str(level)]=dict(cells=int(mask.sum()),
            changed_cells=int(np.count_nonzero((binary-beta)[mask])),
            constraint_delta_normalized=per_constraint[mask].sum(axis=0).tolist())
    return dict(available=True,status=solve['status'],
        scope='TRACE_ROUNDING_DECOMPOSITION_NOT_GEOMETRY_VERIFICATION',
        relaxed_moments=relaxed.tolist(),binary_moments=actual.tolist(),
        raw_moment_delta=per_cell.sum(axis=0).tolist(),
        relaxed_constraints_normalized=rc.tolist(),binary_constraints_normalized=bc.tolist(),
        constraint_delta_normalized=delta.tolist(),
        normalized_relaxed_residual=float(np.max(np.abs(rc))),
        normalized_binary_residual=float(np.max(np.abs(bc))),
        rounding_changed_cells=int(np.count_nonzero(binary-beta)),by_level=by_level,
        total_cells=len(beta),
        max_per_cell_constraint_change=float(np.max(np.abs(per_constraint))),
        sum_absolute_constraint_changes=np.abs(per_constraint).sum(axis=0).tolist())
