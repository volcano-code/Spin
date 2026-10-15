import unittest
import gmpy2
import core_boolean_mesh as core


class RationalCoordinateTokens(unittest.TestCase):
    def token(self,x):
        fn=getattr(core,'binary64_rational_token',None)
        self.assertTrue(callable(fn),'canonical binary64 rational encoder missing')
        return fn(x)

    def test_zero_and_integers_keep_explicit_denominator(self):
        self.assertEqual(self.token(0.0),'0/1')
        self.assertEqual(self.token(-2.0),'-2/1')

    def test_binary64_value_not_decimal_reinterpretation(self):
        self.assertEqual(self.token(.5),'1/2')
        self.assertEqual(gmpy2.mpq(self.token(.1)),gmpy2.mpq(.1))

    def test_nonfinite_is_rejected(self):
        for x in (float('nan'),float('inf')):
            with self.assertRaises(ValueError):self.token(x)
