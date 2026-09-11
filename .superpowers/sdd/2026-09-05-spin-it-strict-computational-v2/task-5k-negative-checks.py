"""Mutation checks for the independent BBW checker, not native solver tests."""
import argparse
import copy
import importlib.util
from pathlib import Path

path = Path(__file__).with_name("task-5k-check-native.py")
spec = importlib.util.spec_from_file_location("checker", path)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)

parser = argparse.ArgumentParser()
parser.add_argument("native")
parser.add_argument("reference")
args = parser.parse_args()
native, reference = checker.load(args.native), checker.load(args.reference)
assert native["cells_per_axis"] == 4 and native["protocol"] == "INDEPENDENT_NORMALIZED"
checker.check(native, reference)


def rejects(label, mutation):
    changed = copy.deepcopy(native)
    mutation(changed)
    try:
        checker.check(changed, reference)
    except (AssertionError, ValueError, IndexError):
        print("REJECTED", label)
        return
    raise AssertionError("checker accepted " + label)


rejects("hand-set success with wrong raw weights", lambda x: x["raw_weights"][31].__setitem__(0, .5))
rejects("substituted node", lambda x: x["nodes"][0].__setitem__(0, .125))
rejects("duplicate handle", lambda x: x["handle_nodes"].__setitem__(1, x["handle_nodes"][0]))
rejects("normalization mislabelled coupled optimum", lambda x: x["final_assessment"].__setitem__("optimality_passed", True))
rejects("dropped actual backend call", lambda x: x["attempts"].pop())
rejects("nonfinite raw weight", lambda x: x["raw_weights"][0].__setitem__(0, float("nan")))
rejects("geometry-energy mismatch", lambda x: x["raw_assessment"].__setitem__("energy", 0))
print("7/7 negative checks rejected; SYNTHETIC_CHECKER_TEST_ONLY")
