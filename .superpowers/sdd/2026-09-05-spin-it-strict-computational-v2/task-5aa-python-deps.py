"""Observe modules actually imported while checking the actual export."""
import importlib.util
from pathlib import Path
import sys

path = Path(__file__).with_name('task-5aa-controller-check.py')
spec = importlib.util.spec_from_file_location('binary_checker', path)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
checker.check_file(sys.argv[1])
paths = {path.with_name('task-5y-controller-check.py').resolve(), Path(__file__).resolve(), path.resolve(),
         path.with_name('task-5w-controller-check.py').resolve(),
         path.with_name('task-5u-controller-check.py').resolve()}
for module in tuple(sys.modules.values()):
    value = getattr(module, '__file__', None)
    if value and Path(value).is_file(): paths.add(Path(value).resolve())
for item in sorted(paths): print(item)
