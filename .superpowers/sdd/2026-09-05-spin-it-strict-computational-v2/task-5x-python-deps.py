"""List code modules actually imported while checking an actual native export."""
import importlib.util
import json
from pathlib import Path
import runpy
import sys

directory = Path(__file__).resolve().parent
checker = runpy.run_path(str(directory/'task-5x-controller-check.py'))
rows = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
checker['verify_rows'](rows)
paths = {Path(__file__).resolve(), directory/'task-5x-controller-check.py'}
for module in tuple(sys.modules.values()):
    path = getattr(module, '__file__', None)
    if not path:
        continue
    if path.endswith('.pyc'):
        path = importlib.util.source_from_cache(path)
    path = Path(path).resolve()
    if path.is_file():
        paths.add(path)
for path in sorted(paths):
    print(path)
