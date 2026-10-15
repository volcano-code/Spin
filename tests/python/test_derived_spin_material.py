"""Public derived exporter re-enters the real material checker, not trust moments."""
import importlib.util
from fractions import Fraction as Q
from pathlib import Path
import unittest
import numpy as np
import test_sealed_spin_artifact as seal_fixtures
from test_exact_material_region import box, combine
from make_it_stand.verification.spin_problem import MeshSnapshot


class DerivedMaterialTests(unittest.TestCase):
    def adapter(self):
        import make_it_stand.verification as api
        self.assertTrue(hasattr(api,'derive_float_material'),
                        'public derived material interface is missing')
        return api.derive_float_material

    def fixture(self):
        helper=seal_fixtures.SpinSealTests();helper.setUp()
        problem,outer,retained,carved=helper.hollow()
        return helper.seal(problem,outer,retained,carved,retained,carved)

    def test_real_material_recheck_keeps_cavity_density_and_no_manufacturing_upgrade(self):
        derive=self.adapter();authority=self.fixture()
        candidate,attachment=derive(authority)
        self.assertEqual(candidate.problem.input_digest,authority.problem.input_digest)
        self.assertEqual(candidate.outer,authority.outer)
        self.assertEqual(candidate.facts['physical']['mass'],7/8)
        self.assertEqual(candidate.facts['bounded_void_components'],1)
        self.assertEqual(candidate.facts['geometry_status'],'MESH_MATERIAL_VALIDATED')
        self.assertEqual(candidate.facts['manufacturing_admission'],'BLOCKED_UNKNOWN_WALL_THICKNESS')
        self.assertEqual(attachment['authority_scientific_digest'],authority.scientific_digest)
        self.assertEqual(attachment['derived_scientific_digest'],candidate.scientific_digest)
        self.assertEqual(attachment['scope'],'DERIVED_MESH_MATERIAL_NOT_MANUFACTURING_OR_SOLVER_APPROVAL')
        self.assertEqual(authority.scientific_digest,self.fixture().scientific_digest)

    def test_wrong_solver_claim_is_rejected_even_with_valid_boundaries(self):
        derive=self.adapter();authority=self.fixture()
        # Deliberately forged stale object: revalidation must use real geometry.
        import json
        claim=authority.solver_claim;claim['moments'][0]+=0.2
        object.__setattr__(authority,'_claim_json',json.dumps(claim))
        with self.assertRaisesRegex(ValueError,'moment|recomputed|claim'):
            derive(authority)

    def test_collapsed_exact_face_reaches_full_partition_without_moving_vertices(self):
        derive=self.adapter();helper=seal_fixtures.SpinSealTests();helper.setUp()
        outer=box();vertices,faces=outer
        exact=[[Q(float(x)) for x in row] for row in vertices]
        exact.append([Q(1,2),Q(1,2)+Q(1,2**54),Q(1)])
        triangles=np.vstack((faces[[i for i in range(len(faces)) if i not in (2,3)]],
                             [[4,6,8],[6,4,5],[4,8,7],[8,6,7]]))
        text=[[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in exact]
        retained=MeshSnapshot.from_rationals(text,triangles)
        self.assertFalse(helper.n.audit_exact_triangle_boundary(*retained.arrays(), 'AnyNonzero')['boundary_valid'])
        problem=helper.problem(outer)
        authority=helper.a.SealedSpinArtifact.create(problem,outer=outer,
            design_results=[('d',retained,None)],final_material=retained,newly_carved=None,
            solver_claim=helper.claim(problem,[(outer,1.)]))
        candidate,attachment=derive(authority)
        np.testing.assert_array_equal(candidate.final_material.arrays()[0],retained.arrays()[0])
        self.assertNotEqual(candidate.final_material.identity,retained.identity)
        self.assertEqual(candidate.facts['physical']['mass'],1.)
        self.assertEqual(candidate.facts['geometry_status'],'MESH_MATERIAL_VALIDATED')
        self.assertGreater(sum(len(row.get('proposal',{}).get('patches',[])) for row in attachment['records']),0)

    def test_real_nonzero_fixed_shell_keeps_two_densities(self):
        derive=self.adapter();helper=seal_fixtures.SpinSealTests();helper.setUp()
        outer=box();domain=box((.25,)*3,(.75,)*3)
        fixed=combine(outer,(domain[0],domain[1][:,::-1]))
        retained=box((.25,)*3,(.5,.75,.75))
        carved=box((.5,.25,.25),(.75,)*3)
        material=combine(outer,(carved[0],carved[1][:,::-1]))
        p=helper.problem(outer=outer,fixed=[('shell',fixed,.75)],design=[('d',domain,2.)],wall=.1)
        authority=helper.a.SealedSpinArtifact.create(p,outer=outer,
            design_results=[('d',retained,carved)],final_material=material,newly_carved=carved,
            solver_claim=helper.claim(p,[(fixed,.75),(retained,2.)]))
        candidate,_=derive(authority)
        self.assertEqual(candidate.facts['physical']['mass'],.78125)
        self.assertEqual([r['density'] for r in candidate.facts['regions']],[.75,2.])

    def test_unknown_or_unsealed_input_cannot_grant_a_candidate(self):
        derive=self.adapter()
        for value in (None,{'passed':True},box()):
            with self.subTest(value=type(value)),self.assertRaises((ValueError,TypeError)):
                derive(value)

    def test_derived_result_persists_and_reloads_under_original_problem_identity(self):
        import tempfile
        from make_it_stand.verification.spin_bundle import save_spin_artifact,load_spin_artifact
        candidate,attachment=self.adapter()(self.fixture())
        with tempfile.TemporaryDirectory(prefix='spin-staged-derived-') as temp:
            path=Path(temp)/'derived-material'
            save_spin_artifact(candidate,path)
            loaded=load_spin_artifact(path,expected_problem_digest=candidate.problem.input_digest)
            self.assertEqual(loaded.scientific_digest,attachment['derived_scientific_digest'])
            self.assertEqual(loaded.facts['physical']['mass'],7/8)
            with self.assertRaises(ValueError):load_spin_artifact(path,expected_problem_digest='0'*64)

    def test_rational_rounding_gap_is_not_hidden_by_nearly_equal_mass(self):
        derive=self.adapter();helper=seal_fixtures.SpinSealTests();helper.setUp()
        outer=box()
        vertices,faces=box((.25,)*3,(.75,)*3)
        # Use exact 1/3 for one design/cavity boundary; it cannot be binary64.
        rational=[[Q(float(x)) for x in row] for row in vertices]
        for row in rational:
            if row[0]==Q(1,4):row[0]=Q(1,3)
        text=lambda rows:[[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rows]
        cavity=MeshSnapshot.from_rationals(text(rational),faces)
        ov,of=outer
        shell_exact=[[Q(float(x)) for x in row] for row in ov]+rational
        shell_faces=np.vstack((of,faces[:,::-1]+len(ov)))
        shell=MeshSnapshot.from_rationals(text(shell_exact),shell_faces)
        problem=helper.problem(outer=outer,fixed=[('shell',shell,1.)],
            design=[('d',cavity,2.)],wall=.1)
        # No carving: analytical box shell + two-density inner cube.
        frame=problem.normalization
        rows=[]
        for lo,hi,rho in [([Q(0)]*3,[Q(1)]*3,Q(1)),
                          ([Q(1,3),Q(1,4),Q(1,4)],[Q(3,4)]*3,Q(1))]:
            lo=[(x-Q(o))/Q(frame['L']) for x,o in zip(lo,frame['origin'])]
            hi=[(x-Q(o))/Q(frame['L']) for x,o in zip(hi,frame['origin'])]
            volume=rho/Q(frame['R'])
            for a,b in zip(lo,hi):volume*=b-a
            mean=[(a+b)/2 for a,b in zip(lo,hi)]
            diag=[(a*a+a*b+b*b)/3 for a,b in zip(lo,hi)]
            rows.append([volume,*[volume*x for x in mean],*[volume*x for x in diag],
                         volume*mean[1]*mean[2],volume*mean[0]*mean[2],volume*mean[0]*mean[1]])
        claim={'frame':'INPUT_NORMALIZED_V1','order':helper.claim(helper.problem(),[(outer,1.)])['order'],
               'moments':[float(a+b) for a,b in zip(*rows)]}
        authority=helper.a.SealedSpinArtifact.create(problem,outer=outer,
            design_results=[('d',cavity,None)],final_material=outer,newly_carved=None,solver_claim=claim)
        with self.assertRaisesRegex(ValueError,'partition|outside|overlap|cover|subset'):
            derive(authority)

if __name__=='__main__':unittest.main()
