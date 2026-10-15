"""Diagnostic replay of recorded SLQP metrics; never a replacement solver.

NumPy SVD replaces Eigen COD for the tiny Schur solve. Agreement with saved
accepted steps is measured, not presumed; disagreement invalidates attribution.
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import numpy as np
from make_it_stand.verification.source_trace_archive import read_archive

def inverse_action(vector,pairs):
    q=np.array(vector,dtype=float,copy=True)
    if not np.isfinite(q).all():raise ValueError('NONFINITE_VECTOR')
    alphas=[]
    for s,y in reversed(pairs):
        curvature=float(s@y)
        if not np.isfinite(curvature) or curvature<=0:raise ValueError('INVALID_CURVATURE')
        alpha=(s@q)/curvature
        q-=np.outer(y,alpha) if q.ndim==2 else y*alpha
        alphas.append(alpha)
    if pairs:
        s,y=pairs[-1]
        q*=np.clip((s@y)/(y@y),1e-6,1e6)
    for (s,y),alpha in zip(pairs,reversed(alphas)):
        beta=(y@q)/(s@y)
        q+=np.outer(s,alpha-beta) if q.ndim==2 else s*(alpha-beta)
    return q

def quasi_newton_direction(g,a,c,active,pairs):
    active=np.asarray(active)!=0
    v=np.column_stack((g,a.T));v[active]=0
    inv=inverse_action(v,pairs);inv[active]=0
    hg,ha=inv[:,0],inv[:,1:]
    lam=np.linalg.lstsq(a@ha,c-a@hg,rcond=1e-12)[0] if len(c) else np.empty(0)
    return -hg-ha@lam

def bound_step(x,d,bounded):
    x=np.asarray(x)[:bounded];d=np.asarray(d)[:bounded]
    if not np.isfinite(x).all() or not np.isfinite(d).all():raise ValueError('NONFINITE_BOUND_STEP')
    pos=d>0;neg=d<0
    return float(np.clip(min(1.,np.min((1-x[pos])/d[pos],initial=1.),
                             np.min(-x[neg]/d[neg],initial=1.)),0.,1.))

def hessian_action(vector,pairs):
    # Direct BFGS Hessian recursion, independent of the inverse two-loop pass.
    scale=1. if not pairs else float(np.clip((pairs[-1][0]@pairs[-1][1])/(pairs[-1][1]@pairs[-1][1]),1e-6,1e6))
    terms=[]
    def apply(v):
        q=np.array(v,dtype=float,copy=True)/scale
        for bs,y,sbs,sy in terms:
            q+=(np.outer(y,y@v)/sy-np.outer(bs,bs@v)/sbs) if q.ndim==2 else y*(y@v)/sy-bs*(bs@v)/sbs
        return q
    for s,y in pairs:
        bs=apply(s);sbs=float(s@bs);sy=float(s@y)
        if min(sbs,sy)<=0 or not np.isfinite([sbs,sy]).all():raise ValueError('INVALID_HESSIAN_CURVATURE')
        terms.append((bs,y,sbs,sy))
    return apply(vector)

def tangent_qp_direction(x,g,a,c,bounded,pairs):
    """Experimental tangent-cone QP, not wired into native optimization.

    Add/release faces using the same SPD metric. Cycles, inaccurate solves and
    large boundary sets remain UNKNOWN, never an infeasibility certificate.
    Interior box limits still require bound_step and the native merit search.
    """
    x=np.asarray(x,dtype=float);g=np.asarray(g,dtype=float)
    a=np.asarray(a,dtype=float);c=np.asarray(c,dtype=float)
    if (x.ndim!=1 or len(x)==0 or g.shape!=x.shape or c.ndim!=1 or
            a.shape!=(len(c),len(x)) or not 0<=bounded<=len(x) or
            not all(np.isfinite(v).all() for v in (x,g,a,c)) or
            np.any(x[:bounded]<0) or np.any(x[:bounded]>1)):
        raise ValueError('INVALID_TANGENT_QP_INPUT')
    lower=np.flatnonzero(x[:bounded]<=1e-10)
    upper=np.flatnonzero(x[:bounded]>=1-1e-10)
    faces={int(i):-1 for i in lower}|{int(i):1 for i in upper}
    if len(faces)>256:return dict(status='QP_UNKNOWN',reason='DIAGNOSTIC_FACE_LIMIT')
    working=sorted(faces);seen=set();hg=inverse_action(g,pairs)
    for iteration in range(4*len(faces)+4):
        key=tuple(working)
        if key in seen:return dict(status='QP_UNKNOWN',reason='WORKING_SET_CYCLE')
        seen.add(key)
        rows=np.zeros((len(working),len(x)))
        for j,i in enumerate(working):rows[j,i]=1
        constraints=np.vstack((a,rows));rhs=np.concatenate((c,np.zeros(len(working))))
        hc=inverse_action(constraints.T,pairs)
        lam=np.linalg.lstsq(constraints@hc,rhs-constraints@hg,rcond=1e-12)[0] if len(rhs) else np.empty(0)
        d=-hg-hc@lam
        # These entries are equality-constrained to exactly zero. Enforce that
        # representation, then check BOTH equality and stationarity again.
        # This does not clamp unconstrained outward components or accepted beta.
        d[working]=0.
        primal=float(np.max(np.abs(constraints@d+rhs),initial=0.))
        if not np.isfinite(d).all() or primal>1e-9:return dict(status='QP_UNKNOWN',reason='LINEARIZED_CONSTRAINT_ERROR')
        violations=[(faces[i]*d[i],i) for i in faces if i not in working and faces[i]*d[i]>1e-10]
        if violations:
            working=sorted(working+[max(violations)[1]]);continue
        dual=[(-faces[i]*lam[len(c)+j],i) for j,i in enumerate(working)
              if -faces[i]*lam[len(c)+j]>1e-10]
        if dual:
            working.remove(max(dual)[1]);continue
        stationarity=float(np.max(np.abs(g+hessian_action(d,pairs)+constraints.T@lam),initial=0.))
        if not np.isfinite(stationarity) or stationarity>1e-9:return dict(status='QP_UNKNOWN',reason='STATIONARITY_ERROR')
        return dict(status='QP_FIRST_ORDER_CHECKED',direction=d,working_set=working,
            iterations=iteration+1,equality_error=primal,stationarity=stationarity)
    return dict(status='QP_UNKNOWN',reason='WORKING_SET_LIMIT')

def secant_pair(previous,current):
    s=np.asarray(current['state'])-previous['state']
    old_lambda=np.asarray(previous['kkt']['lambda'])
    delta_a=np.asarray(current['constraint_jacobian'])-previous['constraint_jacobian']
    y=np.asarray(current['gradient'])-previous['gradient']+delta_a.T@old_lambda
    return s,y

def _checked_qp_direction(p,x,g,a,c,bounded,pairs,active_tolerance,residual_tolerance,finite_box=False,scaled_schur=False):
    """Check the saved tangent QP against independently reconstructed BFGS B.

    This verifies the local direction certificate, not global optimality or
    the nonlinear objective/geometry. Never substitute the unselected QN step.
    """
    if (not np.isfinite(active_tolerance) or not 0<=active_tolerance<.5 or
            not np.isfinite(residual_tolerance) or residual_tolerance<=0):
        raise ValueError('METRIC_INVALID_TOLERANCES')
    q=p['finite_box' if finite_box else 'metric_tangent']
    if (not all(q.get(k) is True for k in ('attempted','checked','selected','accepted')) or
            q.get('status')!=('BOX_QP_FIRST_ORDER_CHECKED' if finite_box else 'QP_FIRST_ORDER_CHECKED') or p['final_state'] or
            q.get('fallback_attempted') is not False or
            p['used_first_order_fallback'] or p['restoration']['accepted'] or
            not np.isfinite(p['accepted_step']) or p['accepted_step']<=0):
        raise ValueError('METRIC_ACCEPTANCE_CONFLICT')
    trials=p.get('objective_trials');count=q.get('trial_count')
    if (type(count) is not int or count<=0 or not isinstance(trials,list) or len(trials)!=count or
            any(t.get('accepted') is not False for t in trials[:-1]) or
            trials[-1].get('accepted') is not True or trials[-1].get('valid') is not True or
            trials[-1].get('reason')!='ACCEPTED' or trials[-1].get('step')!=p['accepted_step']):
        raise ValueError('METRIC_ACCEPTED_TRIAL_MISMATCH')
    d=np.asarray(q.get('direction'),dtype=float);lam=np.asarray(q.get('multipliers'),dtype=float)
    working=q.get('working_set')
    if (not isinstance(working,list) or any(type(i) is not int for i in working) or
            len(set(working))!=len(working) or any(i<0 or i>=(2*bounded if finite_box else bounded) for i in working) or
            d.shape!=x.shape or lam.shape!=(len(c)+len(working),) or
            not np.isfinite(d).all() or not np.isfinite(lam).all()):
        raise ValueError('METRIC_CERTIFICATE_DIMENSION_OR_VALUE')
    faces={i:(-1 if x[i]<=active_tolerance else 1) for i in range(bounded)
        if x[i]<=active_tolerance or x[i]>=1-active_tolerance}
    indices=[i//2 if finite_box else i for i in working]
    if finite_box:
        faces={i:(1 if i%2 else -1) for i in range(2*bounded)}
        if (len(set(indices))!=len(indices) or any(abs(d[i//2]-(i%2-x[i//2]))>residual_tolerance for i in working)):
            raise ValueError('BOX_WORKING_FACE_MISMATCH')
    elif any(i not in faces or d[i]!=0. for i in working):
        raise ValueError('METRIC_WORKING_SET_NOT_FIXED_BOUNDARY')
    # Native checks a full-rank Schur system before issuing CHECKED. Use an
    # independent SVD rank check; this is not bitwise Eigen COD replication.
    for fixed in (([],working) if finite_box else (list(faces),working)):
        constraints=np.zeros((len(a)+len(fixed),len(x)))
        constraints[:len(a)]=a
        for j,i in enumerate(fixed):constraints[len(a)+j,i//2 if finite_box else i]=1.
        if len(constraints):
            schur=constraints@inverse_action(constraints.T,pairs)
            if scaled_schur:
                diagonal=np.diag(schur)
                if not np.isfinite(diagonal).all() or np.any(diagonal<=0):
                    raise ValueError('METRIC_INVALID_SCALED_SCHUR')
                scale=1/np.sqrt(diagonal)
                schur=scale[:,None]*schur*scale[None,:]
            singular=np.linalg.svd(schur,compute_uv=False)
            if (len(constraints)>len(x) or not np.isfinite(singular).all() or
                    np.count_nonzero(singular>1e-12*singular[0])!=len(constraints)):
                raise ValueError('METRIC_RANK_DEFICIENT')
    stationarity=g+hessian_action(d,pairs)+a.T@lam[:len(c)]
    stationarity[indices]+=lam[len(c):]
    equality_error=float(np.max(np.abs(a@d+c),initial=0.))
    stationarity_error=float(np.max(np.abs(stationarity),initial=0.))
    residual=max(equality_error,stationarity_error)
    if not np.isfinite(residual) or residual>residual_tolerance:
        raise ValueError('METRIC_KKT_CERTIFICATE_DISAGREES')
    face_residual=max(0.,max((faces[i]*(d[i//2]-(i%2-x[i//2])) if finite_box else faces[i]*d[i] for i in faces),default=0.),
        max((-faces[i]*lam[len(c)+j] for j,i in enumerate(working)),default=0.))
    if face_residual>(min(active_tolerance,residual_tolerance) if finite_box else active_tolerance):
        raise ValueError('METRIC_FACE_KKT_DISAGREES')
    if finite_box:
        complementarity=max((abs(lam[len(c)+j]*(d[i//2]-(i%2-x[i//2])))
                             for j,i in enumerate(working)),default=0.)
        if not np.isfinite(complementarity) or complementarity>residual_tolerance:
            raise ValueError('BOX_COMPLEMENTARITY_DISAGREES')
        residual=max(residual,complementarity)
    trace=q.get('trace')
    if (not isinstance(trace,list) or not trace or trace[-1].get('action')!='CHECKED' or
            trace[-1].get('working_set')!=working or trace[-1].get('changed_index')!=-1):
        raise ValueError('METRIC_TERMINAL_TRACE_MISMATCH')
    for key,value in (('equality_residual',equality_error),('stationarity',stationarity_error)):
        recorded=q.get(key)
        if (not isinstance(recorded,(int,float)) or not np.isfinite(recorded) or recorded<0 or
                recorded>residual_tolerance or abs(recorded-value)>residual_tolerance or
                trace[-1].get(key)!=recorded):
            raise ValueError('METRIC_RECORDED_RESIDUAL_MISMATCH')
    if p['accepted_step']>bound_step(x,d,bounded)+1e-12:
        raise ValueError('METRIC_ACCEPTED_STEP_OUTSIDE_BOX')
    return d,max(residual,face_residual)

def checked_metric_direction(p,x,g,a,c,bounded,pairs,active_tolerance,residual_tolerance):
    return _checked_qp_direction(p,x,g,a,c,bounded,pairs,active_tolerance,residual_tolerance)

def checked_finite_box_direction(p,x,g,a,c,bounded,pairs,active_tolerance,residual_tolerance,scaled_schur=False):
    """Independent direct BFGS/KKT check of all finite displacement bounds."""
    return _checked_qp_direction(p,x,g,a,c,bounded,pairs,active_tolerance,residual_tolerance,True,scaled_schur)

def replay(path,probe_tangent=False):
    rows=[];previous=None;pairs=[];pending=None;pending_metric=False;pending_box=False;count=0
    for event in read_archive(path):
        count+=1;p=event['payload']
        if event['kind']=='header':
            settings=p['settings']['solver_settings'];capacity=settings['slqp']['lbfgs_memory']
        elif event['kind']=='pass':
            previous=None;pairs=[];pending=None;pending_metric=False;pending_box=False;bounded=p['bounded_variable_count']
            expected_n=len(p['final_state']) if 'final_state' in p else None
            expected_m=None
            row=dict(index=p['index'],status=p['status'],records=0,fallbacks=0,
                fallback_zero_bound_step=0,fallback_positive_bound_step=0,
                history_resets=0,accepted_step_comparisons=0,max_step_replay_error=0.,
                metric_accepted_step_comparisons=0,max_metric_kkt_replay_error=0.,
                finite_box_accepted_step_comparisons=0,max_finite_box_kkt_replay_error=0.,
                min_raw_scale=None,max_raw_scale=None,scale_clamp_count=0,examples=[],
                tangent_probe=[])
            rows.append(row)
        else:
            row['records']+=1
            x=np.asarray(p['state']);g=np.asarray(p['gradient']);c=np.asarray(p['constraints'])
            a=np.asarray(p['constraint_jacobian']);active=np.asarray(p['active_membership'])
            if x.ndim!=1 or c.ndim!=1:raise ValueError('TRACE_DIMENSION_MISMATCH')
            if expected_n is None:expected_n=len(x)
            if expected_m is None:expected_m=len(c)
            if (len(x)!=expected_n or len(x)==0 or not 0<=bounded<=len(x) or g.shape!=x.shape or
                    len(c)!=expected_m or a.shape!=(len(c),len(x)) or active.shape!=x.shape or
                    np.asarray(p['kkt']['lambda']).shape!=(len(c),)):
                raise ValueError('TRACE_DIMENSION_MISMATCH')
            if (not all(np.isfinite(v).all() for v in (x,g,a,c,active,p['kkt']['lambda'])) or
                    not np.isin(active,[0,1,2]).all() or np.any(active[bounded:]!=0)):
                raise ValueError('INVALID_TRACE_STATE')
            if previous is not None:
                if pending is not None:
                    if x.shape!=pending.shape:raise ValueError('PENDING_DIMENSION_MISMATCH')
                    error=float(np.max(np.abs(x-pending)))
                    if not np.isfinite(error):raise ValueError('NONFINITE_REPLAY_ERROR')
                    row['accepted_step_comparisons']+=1
                    row['metric_accepted_step_comparisons']+=int(pending_metric)
                    row['finite_box_accepted_step_comparisons']+=int(pending_box)
                    row['max_step_replay_error']=max(row['max_step_replay_error'],error)
                if previous['used_first_order_fallback'] or previous['restoration']['accepted']:pairs=[]
                if previous['curvature_pair_accepted']:
                    pairs.append(secant_pair(previous,p));pairs=pairs[-capacity:]
            if len(pairs)!=p['history_before']:raise ValueError('HISTORY_BEFORE_MISMATCH')
            reason=p['history_reset_reason']
            if reason.startswith(('NEW_SOLVE_OR_PARTITION','ACTIVE_MEMBERSHIP_CHANGED')):
                pairs=[];row['history_resets']+=1
            if len(pairs)!=p['lbfgs_history_size']:raise ValueError('HISTORY_APPLIED_MISMATCH')
            pending=None;pending_metric=False;pending_box=False
            metric_direction=None
            metric=p.get('metric_tangent',{})
            box=p.get('finite_box',{})
            box_policy=settings['slqp'].get('direction_policy') in ('FINITE_BOX_QP_V1','FINITE_BOX_PRIMAL_V2')
            if box and (metric or not box_policy):
                raise ValueError('BOX_POLICY_OR_EVENT_CONFLICT')
            if (box_policy and
                    p['accepted_step']>0 and not box):
                raise ValueError('BOX_MISSING_DIRECTION_EVIDENCE')
            if box.get('accepted') or (box_policy and p['accepted_step']>0 and
                    not p['used_first_order_fallback'] and not p['restoration']['accepted']):
                metric_direction,error=checked_finite_box_direction(p,x,g,a,c,bounded,pairs,
                    settings['slqp'].get('active_set_tolerance',1e-10),
                    settings['slqp'].get('finite_box_residual_tolerance',1e-9),
                    scaled_schur=settings['slqp'].get('direction_policy')=='FINITE_BOX_PRIMAL_V2')
                row['max_finite_box_kkt_replay_error']=max(row['max_finite_box_kkt_replay_error'],error)
            if metric.get('accepted') or (metric.get('selected') and p['accepted_step']>0 and
                    not p['used_first_order_fallback'] and not p['restoration']['accepted']):
                metric_direction,error=checked_metric_direction(p,x,g,a,c,bounded,pairs,
                    settings['slqp'].get('active_set_tolerance',1e-10),
                    settings['slqp'].get('metric_tangent_residual_tolerance',1e-9))
                row['max_metric_kkt_replay_error']=max(row['max_metric_kkt_replay_error'],error)
            if not p['final_state']:
                d=quasi_newton_direction(g,a,c,active,pairs)
                step=bound_step(x,d,bounded)
                if pairs:
                    s,y=pairs[-1];scale=float((s@y)/(y@y))
                    row['min_raw_scale']=scale if row['min_raw_scale'] is None else min(row['min_raw_scale'],scale)
                    row['max_raw_scale']=scale if row['max_raw_scale'] is None else max(row['max_raw_scale'],scale)
                    row['scale_clamp_count']+=int(not 1e-6<=scale<=1e6)
                if p['used_first_order_fallback']:
                    row['fallbacks']+=1
                    row['fallback_zero_bound_step' if step<=1e-14 else 'fallback_positive_bound_step']+=1
                    if len(row['examples'])<8:
                        row['examples'].append(dict(iteration=p['iteration'],qn_bound_step=step,
                            directional_derivative=float(g@d),linearized_residual=float(np.max(np.abs(a@d+c),initial=0.)),
                            history_size=len(pairs)))
                    if probe_tangent:
                        trial=tangent_qp_direction(x,g,a,c,bounded,pairs)
                        if 'direction' in trial:
                            td=trial.pop('direction')
                            trial.update(bound_step=bound_step(x,td,bounded),directional_derivative=float(g@td))
                        trial['iteration']=p['iteration'];row['tangent_probe'].append(trial)
                elif p['accepted_step']>0 and not p['restoration']['accepted']:
                    if metric_direction is not None:
                        d=metric_direction
                        pending_metric=not bool(box)
                        pending_box=bool(box)
                    pending=x+p['accepted_step']*d
                    pending[:bounded]=np.clip(pending[:bounded],0,1)
            previous=p
    with Path(path).open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
    if not rows or sum(r['accepted_step_comparisons'] for r in rows)==0:raise ValueError('NO_REPLAY_COMPARISONS')
    if max(r['max_step_replay_error'] for r in rows)>1e-9:raise ValueError('REPLAY_DISAGREES_WITH_NATIVE_STEPS')
    return dict(scope='DIAGNOSTIC_NUMPY_REPLAY_NOT_OPTIMIZATION_OR_PAPER_ACCEPTANCE',
        archive_sha256=digest,verified_records=count,passes=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--probe-tangent',action='store_true');args=p.parse_args()
    print(json.dumps(replay(args.archive,args.probe_tangent),indent=2,allow_nan=False))
