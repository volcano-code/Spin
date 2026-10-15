import unittest
import phase_timing


class TimingTests(unittest.TestCase):
    def test_nested_stages_measure_each_interval_and_count_calls(self):
        clock = iter((0., 1., 3., 5., 6., 7.))
        timer = phase_timing.PhaseTiming(lambda: next(clock))
        for name in ('outer-begin','inner-begin','inner-end','outer-end','outer-begin','outer-end'):
            timer(name)
        report = timer.finish()
        self.assertEqual(report['stages']['inner'], {'calls':1, 'seconds':2.})
        self.assertEqual(report['stages']['outer'], {'calls':2, 'seconds':6.})
        self.assertFalse(report['incomplete_stages'])

    def test_missing_end_remains_incomplete_instead_of_zero_cost(self):
        timer = phase_timing.PhaseTiming(lambda: 1.)
        timer('boolean-begin')
        self.assertEqual(timer.finish()['incomplete_stages'], ['boolean'])


if __name__ == '__main__': unittest.main()
