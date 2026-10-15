"""CTest adapter: empty, skipped and expected-failure suites cannot pass."""
from pathlib import Path
import sys
import unittest


def main():
    if len(sys.argv)!=3:
        raise SystemExit('usage: runner SOURCE_ROOT TEST_MODULE')
    root=Path(sys.argv[1]).resolve();module=sys.argv[2]
    if not module.startswith('test_') or not module.isidentifier():
        raise SystemExit('invalid test module')
    if not (root/'tests/python'/f'{module}.py').is_file():
        raise SystemExit('required test module missing')
    sys.path[:0]=[str(root),str(root/'tests/python')]
    loader=unittest.TestLoader();suite=loader.loadTestsFromName(module)
    planned=suite.countTestCases()
    if loader.errors or planned==0:
        raise SystemExit('required suite failed to load or contains no tests')
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if (not result.wasSuccessful() or result.testsRun!=planned or result.skipped
        or result.expectedFailures or result.unexpectedSuccesses):
        raise SystemExit('required suite did not execute every test successfully')


if __name__=='__main__':main()
