"""Independent rational differentiation, not a native solver success test."""
import unittest
from fractions import Fraction

class IndependentHessianTests(unittest.TestCase):
    def test_integer_divisor_cannot_introduce_float_reciprocal(self):
        from make_it_stand.verification.authority_hessians import _Dual
        x=_Dual(Fraction(3,7),[Fraction(1)]+[Fraction(0)]*9)
        y=x/3
        self.assertEqual(y.value,Fraction(1,7))
        self.assertEqual(y.derivative[0],Fraction(1,3))
        self.assertIsInstance(y.value,Fraction)

    def test_independent_hessians_match_three_step_gradient_differences(self):
        from make_it_stand.verification.authority_hessians import independent_authority_hessians
        from make_it_stand.verification._spin2014_authority import _independent_authority_physics_values
        raw=[Fraction(x) for x in ['2','1/5','-1/7','1','3','4','5','-1/9','1/11','1/13']]
        arithmetic=(Fraction,lambda x:x,None)
        texts=[str(x.numerator)+'/'+str(x.denominator) for x in raw]
        for name in ('YOYO','TOP_INERTIA_ONLY','TOP_LOWERING_ONLY','TOP_MASS_AND_LOWERING','STATIC_LOWERING'):
            settings=dict(objective=name,inertia_weight=2,center_weight=3,
                solver_mass_unit_in_declared_units=3,solver_length_unit_in_declared_units=2)
            expected=independent_authority_hessians(texts,settings)['objective_hessian']
            for column in range(10):
                for factor in (1,2,4):
                    h=Fraction(1,10000*factor)*max(1,abs(raw[column]))
                    plus=raw.copy();minus=raw.copy();plus[column]+=h;minus[column]-=h
                    a=_independent_authority_physics_values(plus,settings,arithmetic)['objective_gradient']
                    b=_independent_authority_physics_values(minus,settings,arithmetic)['objective_gradient']
                    for row,(av,bv) in enumerate(zip(a,b)):
                        wanted=expected[row][column]
                        self.assertLessEqual(abs((av-bv)/(2*h)-wanted)/max(1,abs(wanted)),Fraction(1,10**6))

    def test_native_hessians_match_independent_formulas_for_five_objectives(self):
        import _make_it_stand as native
        from make_it_stand.verification.authority_hessians import independent_authority_hessians
        moments=['2/1','1/5','-1/7','1/1','3/1','4/1','5/1','-1/9','1/11','1/13']
        for name in ('YOYO','TOP_INERTIA_ONLY','TOP_LOWERING_ONLY','TOP_MASS_AND_LOWERING','STATIC_LOWERING'):
            settings=dict(objective=name,inertia_weight=2,center_weight=3,
                solver_mass_unit_in_declared_units=3,solver_length_unit_in_declared_units=2)
            ns=native.Spin2014QuadricObjectiveSettings()
            for key,value in settings.items():
                setattr(ns,key,getattr(native.Spin2014Objective,value) if key=='objective' else value)
            actual=native.evaluate_spin2014_authority_moment_problem(ns,moments)
            expected=independent_authority_hessians(moments,settings)
            for key in ('objective_hessian','constraint_hessians'):
                matrices=[(actual[key],expected[key])] if key=='objective_hessian' else zip(actual[key],expected[key])
                for a,b in matrices:
                    for ar,br in zip(a,b):
                        for av,bv in zip(ar,br):
                            self.assertLessEqual(abs(Fraction(av)-bv),Fraction(1,10**100)*max(1,abs(bv)))

    def test_static_hessian_matches_hand_derivation(self):
        from make_it_stand.verification.authority_hessians import independent_authority_hessians
        result=independent_authority_hessians(['2/1','0/1','0/1','3/1','5/1','6/1','7/1','0/1','0/1','0/1'],
            dict(objective='STATIC_LOWERING',center_weight=2))
        expected=[[Fraction(0)]*10 for _ in range(10)]
        expected[0][0]=Fraction(27,4)
        expected[0][3]=expected[3][0]=Fraction(-3)
        expected[3][3]=Fraction(1)
        self.assertEqual(result['objective_hessian'],expected)

    def test_all_five_hessians_are_exactly_symmetric_and_backend_equal(self):
        from make_it_stand.verification.authority_hessians import independent_authority_hessians
        moments=['2/1','1/5','-1/7','1/1','3/1','4/1','5/1','-1/9','1/11','1/13']
        for name in ('YOYO','TOP_INERTIA_ONLY','TOP_LOWERING_ONLY','TOP_MASS_AND_LOWERING','STATIC_LOWERING'):
            settings=dict(objective=name,inertia_weight=2,center_weight=3,
                solver_mass_unit_in_declared_units=3,solver_length_unit_in_declared_units=2)
            a=independent_authority_hessians(moments,settings)
            b=independent_authority_hessians(moments,settings,arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
            self.assertEqual(a,b)
            for matrix in [a['objective_hessian']]+a['constraint_hessians']:
                self.assertEqual(matrix,[list(row) for row in zip(*matrix)])

if __name__=='__main__': unittest.main()
