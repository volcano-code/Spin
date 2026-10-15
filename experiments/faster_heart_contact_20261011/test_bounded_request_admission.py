import unittest
from fractions import Fraction as Q


def tetra():
    return dict(vertices=[['0/1']*3,['1/1','0/1','0/1'],['0/1','1/1','0/1'],['0/1','0/1','1/1']],
                faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]])


class BoundedAdmissionTests(unittest.TestCase):
    def module(self):
        import admit_bounded_request
        return admit_bounded_request

    def test_explicit_bounded_pin_policy_requires_companion(self):
        self.assertTrue(callable(getattr(self.module(),'resolve_pin_policy',None)))
        with self.assertRaisesRegex(ValueError,'PIN_COMPANION_REQUIRED'):
            self.module().resolve_pin_policy('BOUNDED_PIN_INTERSECTION_V1',None,{})

    def test_original_strict_policy_never_implicitly_consumes_companion(self):
        self.assertTrue(callable(getattr(self.module(),'resolve_pin_policy',None)))
        self.assertIsNone(self.module().resolve_pin_policy('ORIGINAL_PIN_STRICT_V1',None,{}))
        with self.assertRaisesRegex(ValueError,'EXPLICIT_PIN_POLICY_REQUIRED'):
            self.module().resolve_pin_policy('ORIGINAL_PIN_STRICT_V1','terminal.json',{})

    def test_bounded_pin_policy_rejects_stale_or_hand_labelled_attachment(self):
        import tempfile,json
        from pathlib import Path
        self.assertTrue(callable(getattr(self.module(),'resolve_pin_policy',None)))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'terminal.json'
            path.write_text(json.dumps(dict(schema='bounded-pin-companion-admission-v1',
                status='PIN_CORRESPONDENCE_CERTIFIED',pin_correspondence_certified=True,
                inputs={},files={})))
            with self.assertRaisesRegex(ValueError,'SOURCE_OR_FILES_REQUIRED'):
                self.module().resolve_pin_policy('BOUNDED_PIN_INTERSECTION_V1',path,{})
            with self.assertRaises(FileNotFoundError):
                self.module().resolve_pin_policy('BOUNDED_PIN_INTERSECTION_V1',path.parent/'missing.json',{})

    def test_face_reordering_cannot_claim_same_correspondence(self):
        a=tetra(); b=tetra(); b['faces'][0],b['faces'][1]=b['faces'][1],b['faces'][0]
        with self.assertRaisesRegex(ValueError,'CONNECTIVITY'):
            self.module().check_pair(a,b,Q(1))

    def test_nonbinary64_derived_coordinates_are_rejected(self):
        a=tetra(); b=tetra(); b['vertices'][1][0]='1/3'
        with self.assertRaisesRegex(ValueError,'BINARY64'):
            self.module().check_pair(a,b,Q(1))

    def test_wrong_rounding_even_if_small_is_rejected(self):
        a=tetra(); b=tetra(); b['vertices'][1][0]='4503599627370497/4503599627370496'
        with self.assertRaisesRegex(ValueError,'MAPPING'):
            self.module().check_pair(a,b,Q(1))

    def test_rounding_above_binary64_midpoint_does_not_double_round(self):
        a=tetra(); b=tetra()
        q=Q(1)+Q(1,2**53)+Q(1,2**110)
        a['vertices'][1][0]=str(q.numerator)+'/'+str(q.denominator)
        rounded=Q.from_float(float(q))
        b['vertices'][1][0]=str(rounded.numerator)+'/'+str(rounded.denominator)
        report=self.module().check_pair(a,b,Q(1))
        self.assertTrue(report['numerical_bounds_within_tolerance'])

    def test_exact_pair_bound_is_not_native_admission(self):
        report=self.module().check_pair(tetra(),tetra(),Q(1))
        self.assertEqual(report['maximum_linf_shift_exact'],'0/1')
        self.assertEqual(report['conditional_symmetric_difference_volume_upper_bound_exact'],'0/1')
        self.assertFalse(report['entity_bound_applied'])
        self.assertFalse(report['solver_admission_granted'])
        self.assertTrue(report['numerical_bounds_within_tolerance'])

    def test_a_large_rounding_volume_bound_fails_numerical_gate(self):
        a=tetra(); b=tetra()
        # Giant geometry gets huge error relative to an independently fixed unit.
        a['vertices'][1][0]='100000000000000000001/1'
        b['vertices'][1][0]='100000000000000000000/1'
        report=self.module().check_pair(a,b,Q(1))
        self.assertFalse(report['numerical_bounds_within_tolerance'])

    def test_candidate_integer_text_compares_by_exact_value_not_spelling(self):
        actual=dict(maximum_linf_shift_exact='0/1',
            conditional_symmetric_difference_volume_upper_bound_exact='1/3',
            independent_authority_moments_exact=['0/1'],
            independent_derived_moments_exact=['1/2'],independent_moment_delta_exact=['1/2'])
        claimed=dict(maximum_linf_vertex_shift_exact='0',
            conditional_symmetric_difference_volume_upper_bound_exact='1/3',
            independent_authority_moments_exact=['0'],
            independent_derived_moments_exact=['1/2'],independent_moment_delta_exact=['1/2'])
        self.module().require_report_matches(actual,claimed)
        claimed['independent_moment_delta_exact']=['1/3']
        with self.assertRaisesRegex(ValueError,'CLAIM_MISMATCH'):
            self.module().require_report_matches(actual,claimed)


if __name__=='__main__':unittest.main()
