"""Moving interpolation scalars are not discrete construction features."""
import copy
import unittest
import cherchi_core_map as core


class FeatureBranchTests(unittest.TestCase):
    def signature(self, records, faces):
        helper=getattr(core,'feature_branch_identity',None)
        self.assertTrue(callable(helper),'source-feature branch helper missing')
        return helper('intersection',faces,[[0,1,2]],'ellipsoid',0,[[0,1]],records)

    def records(self):
        return [dict(point_type='LPI',parameter_t=.2+i*.1,operands=[dict(status='UNIQUE',
            source_vertices=[dict(input_object=1,input_vertex=i)])]) for i in range(3)]

    def test_motion_scalar_not_a_branch(self):
        a=self.records();b=copy.deepcopy(a);b[0]['parameter_t']=.9
        self.assertEqual(self.signature(a,[[0,1,2]]),self.signature(b,[[0,1,2]]))

    def test_vertex_reindex_and_cyclic_face_rotation_not_a_branch(self):
        a=self.records();b=[a[2],a[0],a[1]]
        self.assertEqual(self.signature(a,[[0,1,2]]),self.signature(b,[[2,0,1]]))

    def test_changed_source_feature_is_a_branch(self):
        a=self.records();b=copy.deepcopy(a);b[0]['operands'][0]['source_vertices'][0]['input_vertex']=9
        self.assertNotEqual(self.signature(a,[[0,1,2]]),self.signature(b,[[0,1,2]]))

    def test_reversed_face_is_a_branch(self):
        a=self.records()
        self.assertNotEqual(self.signature(a,[[0,1,2]]),self.signature(a,[[0,2,1]]))

    def test_ambiguous_features_cannot_claim_correspondence(self):
        a=self.records();a[0]['operands'][0]['status']='AMBIGUOUS'
        with self.assertRaisesRegex(ValueError,'CORRESPONDENCE_UNAVAILABLE'):
            self.signature(a,[[0,1,2]])

    def test_duplicate_features_cannot_claim_correspondence(self):
        a=self.records();a[1]=copy.deepcopy(a[0])
        with self.assertRaisesRegex(ValueError,'DUPLICATE_BRANCH_VERTEX_FEATURE'):
            self.signature(a,[[0,1,2]])
