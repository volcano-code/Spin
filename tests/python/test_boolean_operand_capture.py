import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest

path = Path(__file__).resolve().parents[2] / 'experiments/faster_core_repair_20261012/boolean_operand_capture.py'
spec = importlib.util.spec_from_file_location('boolean_operand_capture', path)
capture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(capture)


def test_persist_actual_operands_before_call_and_stop_before_entry():
    events = []
    def boolean(*args):
        events.append(('native', args))
        return 'result'
    proxy = capture.CapturingBridge(SimpleNamespace(boolean=boolean),
        lambda name, body: events.append(('saved', body)), lambda: {'coefficients': list(range(10))}, stop_at=2)
    operands = ([[1., 2., 3.]], [[0, 0, 0]], [[4., 5., 6.]], [[0, 0, 0]], 'difference')
    assert proxy.boolean(*operands) == 'result'
    assert [x[0] for x in events] == ['saved', 'native']
    assert events[0][1]['source_vertices'] == operands[2]
    with pytest.raises(capture.StopBeforeBoolean):
        proxy.boolean(*operands)
    assert [x[0] for x in events] == ['saved', 'native', 'saved']


def test_failed_durable_publish_prevents_native_entry():
    def fail(*args):
        raise OSError('disk full')
    proxy = capture.CapturingBridge(SimpleNamespace(boolean=lambda *args: pytest.fail('entered native')),
                                    fail, lambda: {})
    with pytest.raises(OSError):
        proxy.boolean([], [], [], [], 'intersection')
