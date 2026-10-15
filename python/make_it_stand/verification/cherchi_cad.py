"""Explicit full-rank Cherchi terminal-to-analytic-CAD parameters.

No Boolean, optimizer, source relabelling or manufacturing approval occurs
here. The existing CAD exporter still validates the actual OCCT output against
every sealed material region. Rank-deficient carriers remain unsupported on
this path. Rehydrating a historical result additionally requires saved-run and
saved-material checks; a dictionary or a `passed` field is not this interface.
"""
from dataclasses import asdict, replace
from fractions import Fraction

import numpy as np

from .spin_common import digest, require, scalar


def full_rank_cad_carrier(coefficients):
    """Classify OCCT's bounded carrier and its sign, not the retained side.

Retained material is always p>=0. The carrier has a definite quadratic block
for ellipsoids, or signature (+,+,-) after orienting hyperboloids. Whether it
represents p>=0 determines Common versus Cut. Exact rational elimination
classifies the represented center level even arbitrarily close to a cone.
The eigensystem conditioning check matches the existing OCCT input envelope.
"""
    q=np.asarray(coefficients,dtype=float)
    require(q.shape==(10,) and np.isfinite(q).all(), 'CAD requires ten finite coefficients')
    require(abs(float(np.linalg.norm(q))-1.)<=1e-10, 'CAD coefficients must have unit norm')
    matrix=np.array(((q[4],q[9],q[8]),(q[9],q[5],q[7]),(q[8],q[7],q[6])))
    values=np.linalg.eigvalsh(matrix)
    threshold=1e-10*float(np.max(np.abs(matrix)))
    require(threshold>0. and np.all(np.abs(values)>threshold),
            'rank-deficient or ill-conditioned CAD carrier is unsupported')
    positive=int(np.count_nonzero(values>0.))
    rows=[[Fraction(float(x)) for x in row] + [-Fraction(float(q[i+1]))/2]
          for i,row in enumerate(matrix)]
    for i in range(3):
        pivot=next((j for j in range(i,3) if rows[j][i]),None)
        require(pivot is not None,'exactly singular quadratic block')
        rows[i],rows[pivot]=rows[pivot],rows[i]
        divisor=rows[i][i]; rows[i]=[x/divisor for x in rows[i]]
        for j in range(3):
            if j!=i:
                factor=rows[j][i]
                rows[j]=[a-factor*b for a,b in zip(rows[j],rows[i])]
    center_level=Fraction(float(q[0]))+sum(Fraction(float(q[i+1]))*rows[i][3]/2 for i in range(3))
    require(center_level!=0,'exact cone has no supported analytic CAD carrier')
    sign=1 if center_level>0 else -1
    if positive in (0,3):
        require((positive==3 and sign<0) or (positive==0 and sign>0),
                'empty ellipsoid has no supported analytic CAD carrier')
        kind='ELLIPSOID'; positive_side=positive==0
    else:
        positive_side=positive==1
        oriented_sign=-sign if positive_side else sign
        kind='ONE_SHEET_HYPERBOLOID' if oriented_sign<0 else 'TWO_SHEET_HYPERBOLOID'
    return {'quadric_type':kind,'quadric_solid_is_positive_region':positive_side,
            'center_level_sign':sign,'classification':'FULL_RANK_ANALYTIC_CARRIER_ONLY'}


def cherchi_cad_parameters(result):
    """Validate the genuine terminal contract without substituting CGAL geometry."""
    from .faster_material import (FasterMaterialResult,_prepare,_normalized_material_moments,
        _transform_moments,rotational_objective_definition,material_run_scientific_identity)
    from ..spin import SpinEvaluatorConfig,SpinSolverConfig
    require(isinstance(result,FasterMaterialResult),'expected strict FasterMaterialResult')
    identity=result.identity
    require(digest(material_run_scientific_identity(identity))==result.scientific_digest and result.status=='CONVERGED' and
            identity.get('status')=='CONVERGED' and result.seal_error is None and
            identity.get('seal_error') is None and result.artifact is not None,
            'Cherchi CAD requires an intact converged result and sealed material')
    artifact=result.artifact; problem=artifact.problem
    require(identity.get('operation')=='SOLVE_CHERCHI_EXTERNAL_GEOMETRY' and
            identity.get('problem_digest')==problem.input_digest and
            identity.get('artifact_scientific_digest')==artifact.scientific_digest,
            'Cherchi CAD operation/problem/material identity mismatch')
    require(identity.get('coefficient_frame')=='INPUT_NORMALIZED_TARGET_V1' and
            identity.get('coefficient_convention')=='NATIVE_MATRIX_ENTRY_P',
            'unsupported Cherchi CAD coefficient coordinates')
    _,request,rebuilt=_prepare(problem,SpinEvaluatorConfig(**identity['evaluator_settings']),
                              SpinSolverConfig(**identity['solver_settings']))
    for key in ('normalization','target_axes','target_origin','fixed_moments_input_normalized'):
        require(identity.get(key)==rebuilt[key],'Cherchi CAD locked '+key+' mismatch')
    objective=identity.get('objective_definition',{})
    require(request.mode in ('SPINNING_TOP','YO_YO'),
            'Cherchi CAD supports rotational terminal requests; use dedicated application adapters')
    expected=rotational_objective_definition(request.mode,
        objective.get('center_of_mass_weight'),objective.get('spinning_top_objective'))
    require(objective==expected,'Cherchi CAD objective differs from the locked request')
    request=replace(request,center_of_mass_weight=expected['center_of_mass_weight'],
                    spinning_top_objective=expected['spinning_top_objective'])
    require(identity.get('native_request_digest')==digest(asdict(request)),
            'Cherchi CAD request digest mismatch')
    run=identity.get('native_run',{})
    raw=scalar(run.get('raw_equation_residual'),'raw Eq.14 residual')
    mu=scalar(run.get('multiplier_mu'),'Eq.14 multiplier')
    require(run.get('status')=='CONVERGED' and run.get('equation_profile')=='DUAL_MONOMIAL_RAW_V1'
            and 0.<=raw<=1e-12 and mu>0.,'Cherchi CAD terminal Eq.14 acceptance failed')
    q=run.get('quadric_coefficients',())
    carrier=full_rank_cad_carrier(q)
    trace=run.get('trace',[]); binding=identity.get('external_terminal_binding',{})
    complement=identity.get('terminal_complement',{})
    require(isinstance(trace,list) and len(trace)>0 and
            trace[-1].get('quadric_coefficients')==q and
            trace[-1].get('material_geometry_hash')==binding.get('retained_geometry') and
            trace[-1].get('material_topology_hash')==binding.get('retained_topology') and
            all(binding.get(key) for key in ('retained_geometry','retained_topology',
                                             'complement_geometry','complement_topology')),
            'Cherchi CAD terminal coefficients/geometry/trace mismatch')
    require(complement.get('status')=='PASSED' and
            complement.get('geometry_hash')==binding.get('complement_geometry') and
            complement.get('topology_hash')==binding.get('complement_topology') and
            str(complement.get('backend','')).startswith(
                'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:') and
            complement.get('boolean_operation') in ('intersection','subtraction'),
            'Cherchi CAD complement binding failed')
    # Reintegrate the authenticated rational material, not its display mesh.
    norm=problem.normalization
    axes=np.asarray(request.target_axes); origin=np.asarray(request.target_origin)
    def target_moments(mesh,rho):
        if mesh.explicit_empty:return np.zeros(10)
        # Subtract the world origin on coordinates BEFORE integration. Shifting
        # rounded world second moments loses the small object covariance.
        model=_normalized_material_moments(mesh,rho,norm)
        return np.asarray(_transform_moments(model,axes.T,-axes.T@origin))
    fixed=sum((target_moments(r.mesh,r.density) for r in problem.fixed_regions),np.zeros(10))
    require(len(artifact.design_results)==1,'Cherchi CAD requires one design region')
    region,retained,_=artifact.design_results[0]
    require(region==problem.design_regions[0].region_id,'Cherchi CAD design region mismatch')
    variable=target_moments(retained,problem.design_regions[0].density)
    for key,want in (('fixed_moments_target',fixed),('variable_moments_target',variable),
                     ('total_moments_target',fixed+variable)):
        got=np.asarray(binding.get(key,()),dtype=float)
        require(got.shape==(10,) and np.isfinite(got).all() and
                np.max(np.abs(got-want))<=1e-10*max(abs(float((fixed+variable)[0])),1e-14),
                'Cherchi CAD independently integrated '+key+' mismatch')
    return {**carrier,'coefficients':tuple(float(x) for x in q),
            'source_operation':identity['operation'],'source_result_digest':result.scientific_digest,
            'source_scope':identity.get('scope'),'terminal_binding_digest':digest(binding)}
