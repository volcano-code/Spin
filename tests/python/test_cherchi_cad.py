"""CAD adapter contracts, NOT new Cherchi optimization successes.

The terminal fixture has a real, independently sealed two-density solid and a
synthetic solver ledger. This deliberately tests dispatch/identity/CAD in
isolation; original saved solver runs require their own evidence checks.
"""
from dataclasses import asdict, replace
import copy
from fractions import Fraction
import importlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from test_faster_material_bridge import problem
from make_it_stand import spin
from make_it_stand.verification import faster_material as bridge
from make_it_stand.verification.spin_common import digest

ROOT=Path(__file__).resolve().parents[2]
PROBE=Path(os.environ.get('MAKE_IT_STAND_OCCT_STEP_PROBE',
    ROOT/'build/occt-nearcone-20261002/occt_step_roundtrip'))

def unit_q(values):
    values=np.asarray(values,dtype=float)
    return tuple(float(x) for x in values/np.linalg.norm(values))

def terminal_contract_fixture(offset=(3.,-4.,5.)):
    p=problem(density=2.,fixed_density=.75,scale=2.,offset=offset,axis=(1.,0.,0.))
    _,request,identity=bridge._prepare(p,
        spin.SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC'),
        spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1'))
    request=replace(request,center_of_mass_weight=0.)
    identity['native_request_digest']=digest(asdict(request))
    identity['objective_definition']=bridge.rotational_objective_definition('YO_YO',0.,'INERTIA_AND_LOWERING')
    target_axes=np.asarray(request.target_axes); target_origin=np.asarray(request.target_origin)
    fixed=bridge._transform_moments(request.fixed_shell_moments,target_axes.T,-target_axes.T@target_origin)
    variable=spin.integrate_raw_mass_moments(request.design_vertices,request.design_faces,request.material_density)
    variable=bridge._transform_moments(variable,target_axes.T,-target_axes.T@target_origin)
    total=(np.asarray(fixed)+variable).tolist()
    empty=(np.empty((0,3)),np.empty((0,3),dtype=int))
    artifact=bridge._seal(p,request,(request.design_vertices,request.design_faces),empty,
        (request.outer_vertices,request.outer_faces),total)
    # (x-3)^2+y^2+z^2-1 > 0 throughout [-.25,.25]^3.
    q=unit_q((8.,-6.,0.,0.,1.,1.,1.,0.,0.,0.))
    native={'status':'CONVERGED','equation_profile':'DUAL_MONOMIAL_RAW_V1',
        'quadric_coefficients':list(q),'raw_equation_residual':0.,'multiplier_mu':1.,
        'trace':[{'quadric_coefficients':list(q),'material_geometry_hash':'contract-geometry',
                  'material_topology_hash':'contract-topology'}]}
    identity.update(operation='SOLVE_CHERCHI_EXTERNAL_GEOMETRY',native_run=native,
        scope='SYNTHETIC_ADAPTER_CONTRACT_NOT_SOLVER_RUN',
        external_terminal_binding={'retained_geometry':'contract-geometry',
            'retained_topology':'contract-topology','complement_geometry':'contract-empty',
            'complement_topology':'contract-empty-topology','fixed_moments_target':fixed,
            'variable_moments_target':variable,'total_moments_target':total},
        terminal_complement={'status':'PASSED','geometry_hash':'contract-empty',
            'topology_hash':'contract-empty-topology',
            'backend':'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:contract',
            'boolean_operation':'intersection'})
    return bridge._result('CONVERGED',identity,artifact)

class CherchiCadTests(unittest.TestCase):
    def api(self):
        spec=importlib.util.find_spec('make_it_stand.verification.cherchi_cad')
        self.assertIsNotNone(spec,'explicit Cherchi terminal CAD adapter is missing')
        return importlib.import_module(spec.name)

    def test_full_rank_carrier_side_and_exact_near_cone(self):
        classify=self.api().full_rank_cad_carrier
        fixtures=[((-1,0,0,0,1,1,1,0,0,0),'ELLIPSOID',False),
                  ((-1,0,0,0,1,1,-1,0,0,0),'ONE_SHEET_HYPERBOLOID',False),
                  ((1,0,0,0,1,1,-1,0,0,0),'TWO_SHEET_HYPERBOLOID',False),
                  ((-1e-36,0,0,0,1,1,-1,0,0,0),'ONE_SHEET_HYPERBOLOID',False),
                  ((1e-36,0,0,0,1,1,-1,0,0,0),'TWO_SHEET_HYPERBOLOID',False),
                  # q=(x-1)^2+y^2-z^2+2^-40, near-center cancellation.
                  ((1+2**-40,-2,0,0,1,1,-1,0,0,0),'TWO_SHEET_HYPERBOLOID',False)]
        for q,kind,side in fixtures:
            for sign in (1,-1):
                with self.subTest(q=q,sign=sign):
                    result=classify(unit_q(np.array(q)*sign))
                    self.assertEqual(result['quadric_type'],kind)
                    self.assertIs(result['quadric_solid_is_positive_region'],side if sign==1 else not side)

    def test_large_translation_does_not_destroy_normalized_terminal_moments(self):
        result=terminal_contract_fixture(offset=(1e8,-1e8,1e8))
        report=self.api().cherchi_cad_parameters(result)
        self.assertEqual(report['quadric_type'],'ELLIPSOID')
        self.assertEqual(report['source_result_digest'],result.scientific_digest)

    def test_unsupported_polynomials_are_not_coerced_into_cad(self):
        classify=self.api().full_rank_cad_carrier
        for q in ((0,0,0,0,1,1,-1,0,0,0),(1,0,0,0,1,1,1,0,0,0),
                  (0,0,0,1,1,1,0,0,0,0),(0,0,0,0,1,1e-15,-1,0,0,0)):
            with self.subTest(q=q),self.assertRaises(ValueError):classify(unit_q(q))
        with self.assertRaises(ValueError):classify([float('nan')]*10)

    def test_actual_cad_keeps_source_operation_frame_and_two_densities(self):
        result=terminal_contract_fixture(); before=result.identity
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            # A CGAL re-evaluation would silently substitute a different map.
            with patch.object(bridge,'evaluate_spin_it_faster_quadric',
                              side_effect=AssertionError('no CGAL fallback')):
                cad=bridge.export_spin_it_faster_quadric_material_step(result,
                    occt_executable=PROBE,brep_path=root/'m.brep',step_path=root/'m.step')
            self.assertEqual(result.identity,before)
            self.assertEqual(cad['source_operation'],'SOLVE_CHERCHI_EXTERNAL_GEOMETRY')
            self.assertEqual(cad['source_scope'],'SYNTHETIC_ADAPTER_CONTRACT_NOT_SOLVER_RUN')
            self.assertEqual(cad['material_region_ids'],['design','shell'])
            np.testing.assert_allclose(cad['entity_densities'],[1.,.375],rtol=0,atol=0)
            self.assertAlmostEqual(cad['raw_mass_moments']['step'][0],58/128,places=10)
            self.assertLessEqual(cad['maximum_normalized_material_moment_error'],1e-6)
            self.assertTrue((root/'m.step').is_file())

    def test_changed_terminal_request_and_complement_cannot_reach_cad(self):
        api=self.api();result=terminal_contract_fixture()
        changes=[('target_origin',[1.,0.,0.]),('target_axes',np.eye(3).tolist()),
                 ('fixed_moments_input_normalized',[0.]*10),
                 ('native_request_digest','0'*64),('normalization',{}),
                 ('problem_digest','0'*64)]
        for key,value in changes:
            identity=copy.deepcopy(result.identity);identity[key]=value
            altered=bridge._result('CONVERGED',identity,result.artifact)
            with self.subTest(key=key),self.assertRaises(ValueError):api.cherchi_cad_parameters(altered)
        for section,key,value in [('native_run','raw_equation_residual',1e-10),
            ('native_run','multiplier_mu',0.),('native_run','raw_equation_residual',float('inf')),
            ('terminal_complement','status','FAILED'),('terminal_complement','geometry_hash','wrong'),
            ('external_terminal_binding','total_moments_target',[0.]*10)]:
            identity=copy.deepcopy(result.identity);identity[section][key]=value
            with self.subTest(section=section,key=key),self.assertRaises(ValueError):
                altered=bridge._result('CONVERGED',identity,result.artifact)
                api.cherchi_cad_parameters(altered)

    def test_wrong_analytic_material_is_rejected_before_publication(self):
        result=terminal_contract_fixture(); identity=result.identity
        # This sphere actually removes material, unlike the sealed full domain.
        q=list(unit_q((-.01,0,0,0,1,1,1,0,0,0)))
        identity['native_run']['quadric_coefficients']=q
        identity['native_run']['trace'][-1]['quadric_coefficients']=q
        altered=bridge._result('CONVERGED',identity,result.artifact)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with self.assertRaisesRegex(ValueError,'normalized moment error'):
                bridge.export_spin_it_faster_quadric_material_step(altered,occt_executable=PROBE,
                    brep_path=root/'m.brep',step_path=root/'m.step')
            for name in ('m.brep','m.step','m.step.identity.json'):
                self.assertFalse((root/name).exists())

    def test_exact_retained_moments_do_not_revalidate_collapsed_display(self):
        from make_it_stand.verification.spin_problem import MeshSnapshot
        from make_it_stand.verification.rational_text import fraction_to_rational_text
        helper=getattr(bridge,'_normalized_material_moments',None)
        self.assertTrue(callable(helper),'CAD normalization must preserve authoritative rational coordinates')
        rows=((Fraction(0),Fraction(0),Fraction(1)),
              (Fraction(1,2),Fraction(0),Fraction(1)+Fraction(1,2**60)),
              (Fraction(1),Fraction(0),Fraction(1)),(Fraction(0),Fraction(1),Fraction(0)))
        faces=((0,2,1),(0,1,3),(0,3,2),(1,2,3))
        mesh=MeshSnapshot.from_rationals([[fraction_to_rational_text(x) for x in row] for row in rows],faces)
        with self.assertRaisesRegex(ValueError,'zero-area'):MeshSnapshot(mesh.arrays())
        normalized=[[x/2 for x in row] for row in rows]
        volume=Fraction(1,2**60)/48; rho=Fraction(5,3)
        sums=[sum(row[i] for row in normalized) for i in range(3)]
        second=[[rho*volume*(sums[i]*sums[j]+sum(row[i]*row[j] for row in normalized))/20
                 for j in range(3)] for i in range(3)]
        want=[float(x) for x in (rho*volume,*(rho*volume*x/4 for x in sums),
            second[0][0],second[1][1],second[2][2],second[1][2],second[0][2],second[0][1])]
        actual=helper(mesh,5.,{'L':2.,'R':3.,'origin':[0.,0.,0.]})
        np.testing.assert_allclose(actual,want,rtol=1e-14,atol=0.)

if __name__=='__main__':unittest.main()
