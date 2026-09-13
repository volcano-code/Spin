"""Read-only catalogue inventory. Counts requirements, never completed work."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[3]
path = root/'experiments/spin_it_reproduction/strict/requirements_catalogue.py'
spec = importlib.util.spec_from_file_location('strict_catalogue_inventory',path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
catalogue = module.catalogue()
rows = catalogue['requirements']
ids = [r['requirement_id'] for r in rows]
assert len(ids) == len(set(ids))
physical = [r for r in rows if r['requirement_kind']=='PHYSICAL_MEASUREMENT']
computational = [r for r in rows if r not in physical]
print(json.dumps({
    'scope':'CATALOGUE_INVENTORY_NOT_IMPLEMENTATION_OR_COMPLETION_PERCENTAGE',
    'catalogue_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
    'total_requirements':len(rows),'computational_requirements':len(computational),
    'physical_out_of_scope':len(physical),
    'by_namespace':dict(sorted(Counter(r['requirement_id'].split('.')[0] for r in rows).items())),
    'by_kind':dict(sorted(Counter(r['requirement_kind'] for r in rows).items())),
    'requirements':rows,
},ensure_ascii=False,indent=2,allow_nan=False))
