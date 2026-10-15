"""Output-only exact patch contracts; not manufacturing or paper acceptance."""
from fractions import Fraction as Q
import importlib.util
import unittest



class CoplanarCandidateTests(unittest.TestCase):
    def engine(self):
        name='make_it_stand.verification._coplanar_boundary'
        spec=importlib.util.find_spec(name)
        self.assertIsNotNone(spec,'production exact-patch proposal is missing')
        return importlib.import_module(name).propose_coplanar_diagonal_repairs

    def fixture(self):
        # Round-to-even collapses the small exact positive area at c.
        exact=[(Q(0),Q(0),Q(0)),(Q(1),Q(1),Q(0)),
               (Q(1,2),Q(1,2)+Q(1,2**54),Q(0)),(Q(1),Q(0),Q(0))]
        return exact,[[float(x) for x in row] for row in exact],[(0,1,2),(1,0,3)]

    def call(self, exact=None, display=None, faces=None, pairs=None, regions=None):
        e,v,f=self.fixture()
        return self.engine()(e if exact is None else exact,v if display is None else display,
            f if faces is None else faces,[(0,1)] if pairs is None else pairs,
            face_regions=['material','material'] if regions is None else regions)

    def test_repairs_rounding_collapsed_diagonal_without_changing_vertices_or_patch_boundary(self):
        e,v,f=self.fixture();before=(list(e),[list(row) for row in v],list(f))
        result=self.call(e,v,f)
        self.assertEqual(result['faces'],[(2,3,1),(2,0,3)])
        self.assertEqual((e,v,f),before)
        self.assertEqual(len(result['patches']),1)
        self.assertEqual(result['patches'][0]['face_indices'],[0,1])
        self.assertEqual(result['status'],'PROPOSED_REQUIRES_FULL_MATERIAL_AUDIT')

    def test_valid_display_is_not_retriangulated_for_cosmetic_reasons(self):
        e,_,f=self.fixture();e[2]=(Q(1,2),Q(3,4),Q(0))
        result=self.call(e,[[float(x) for x in row] for row in e],f)
        self.assertEqual(result['faces'],f);self.assertEqual(result['patches'],[])

    def test_rejects_non_coplanar_patch_even_below_binary64_resolution(self):
        e,_,f=self.fixture();e[3]=(e[3][0],e[3][1],Q(1,2**1080))
        result=self.call(e,[[float(x) for x in row] for row in e],f)
        self.assertEqual(result['faces'],f);self.assertEqual(result['patches'],[])

    def test_does_not_flip_concave_exact_quad_or_inconsistent_winding(self):
        e,_,f=self.fixture();e[2]=(Q(2),Q(2)+Q(1,2**54),Q(0))
        result=self.call(e,[[float(x) for x in row] for row in e],f)
        self.assertEqual(result['faces'],f);self.assertEqual(result['patches'],[])
        e,v,f=self.fixture();f[1]=(0,1,3)
        self.assertEqual(self.call(e,v,f)['patches'],[])

    def test_does_not_cross_material_interface_or_nonmanifold_existing_diagonal(self):
        self.assertEqual(self.call(regions=['shell','inner'])['patches'],[])
        e,v,f=self.fixture()
        for extra in ((0,1,3),(2,3,0)):
            faces=f+[extra]
            result=self.call(e,v,faces,regions=['material']*3)
            self.assertEqual(result['faces'],faces);self.assertEqual(result['patches'],[])

    def test_deduplicates_pair_order_and_preserves_distinct_patch_identity(self):
        result=self.call(pairs=[(1,0),(0,1),(1,0)])
        self.assertEqual(result,self.call())

    def test_rejects_lossy_coordinates_invalid_indices_and_labels(self):
        e,v,f=self.fixture()
        invalid=[dict(display=[[False,0.,0.],*v[1:]]),
                 dict(display=[[float('nan'),0.,0.],*v[1:]]),
                 dict(display=[[0.25,0.,0.],*v[1:]]),
                 dict(exact=[(0.,Q(0),Q(0)),*e[1:]]),
                 dict(faces=[(False,1,2),f[1]]),dict(pairs=[(-1,1)]),
                 dict(pairs=[(0,2)]),dict(regions=['material'])]
        for kwargs in invalid:
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):self.call(**kwargs)

    def test_invalid_exact_patch_is_not_blessed_by_a_float_repair(self):
        e,v,f=self.fixture();e[2]=(Q(1,2),Q(1,2),Q(0))
        self.assertEqual(self.call(e,v,f)['patches'],[])


    def test_multiple_independent_patches_preserve_all_vertices_and_hint_order(self):
        e,v,f=self.fixture()
        extra=[tuple(x+Q(2) if j==2 else x for j,x in enumerate(row)) for row in e]
        exact=e+extra;display=v+[[float(x) for x in row] for row in extra]
        faces=f+[tuple(i+4 for i in face) for face in f]
        result=self.call(exact,display,faces,pairs=[(3,2),(1,0)],regions=['m']*4)
        self.assertEqual(result['faces'],[(2,3,1),(2,0,3),(6,7,5),(6,4,7)])
        self.assertEqual([x['face_indices'] for x in result['patches']],[[0,1],[2,3]])
        self.assertEqual(result,self.call(exact,display,faces,pairs=[(0,1),(2,3)],regions=['m']*4))
        self.assertEqual(display,v+[[float(x) for x in row] for row in extra])

    def test_shared_face_proposal_requires_new_full_audit(self):
        e,v,f=self.fixture()
        e=e+[(Q(0),Q(0),Q(1))];v=v+[[0.,0.,1.]]
        faces=f+[(0,3,4)]
        result=self.call(e,v,faces,pairs=[(0,1),(1,2)],regions=['m']*3)
        self.assertEqual(len(result['patches']),1)
        self.assertEqual(result['faces'][2],faces[2])
        self.assertEqual(result['skipped_pairs'],[{'face_indices':[1,2],
            'reason':'OVERLAPPING_PROPOSAL_REQUIRES_NEW_WHOLE_AUDIT'}])

    def test_new_edge_from_first_patch_blocks_second_replacement(self):
        e,v,f=self.fixture()
        e=e+[e[0],e[1]];v=v+[v[0],v[1]]
        faces=f+[(4,5,2),(5,4,3)]
        result=self.call(e,v,faces,pairs=[(0,1),(2,3)],regions=['m']*4)
        self.assertEqual(len(result['patches']),1)
        self.assertEqual(result['faces'][2:],faces[2:])
        self.assertEqual(result['skipped_pairs'],[{'face_indices':[2,3],
            'reason':'REPLACEMENT_DIAGONAL_ALREADY_EXISTS'}])

if __name__=='__main__':unittest.main()
