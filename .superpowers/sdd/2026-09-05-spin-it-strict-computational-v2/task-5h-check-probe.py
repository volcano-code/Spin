"""Independent checks of actual native grid exports; not a BBW solve."""
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)

paths = [Path(value) for value in sys.argv[1:]]
require(len(paths) == 2, "two actual exports required")
for path in paths:
    raw = path.read_bytes()
    data = json.loads(raw)
    require(data["scope"] == "SYNTHETIC_GRID_INTERPOLATION_NOT_BBW", "wrong scope")
    nodes = data["nodes"]
    cells = data["cells"]
    edges = data["edges"]
    handles = data["handle_nodes"]
    require((len(cells), len(nodes), len(edges), len(handles)) == (1728,2197,6084,218), "counts")
    expected_nodes = list(itertools.product(range(2,15), repeat=3))
    require(nodes == [[c/16 for c in key] for key in expected_nodes], "actual sparse nodes")
    expected_cells = list(itertools.product(range(2,14), repeat=3))
    node_lookup = {key: index for index,key in enumerate(expected_nodes)}
    expected_edges = set()
    for key, actual in zip(expected_cells,cells):
        corners = [node_lookup[tuple(key[a]+((bit>>a)&1) for a in range(3))] for bit in range(8)]
        require(actual == corners, "actual cell connectivity/corner convention")
        for bit in range(8):
            for axis in range(3):
                if not (bit & (1<<axis)):
                    expected_edges.add(tuple(sorted((corners[bit],corners[bit|(1<<axis)]))))
    require(edges == [list(e) for e in sorted(expected_edges)], "edges exactly from cells")
    expected_handles = [node_lookup[key] for key in itertools.product(range(2,15,2),repeat=3)
                        if any(c in (2,14) for c in key)]
    require(handles == expected_handles, "all coarse handles retained in order")
    maximum_error = 0.0
    for query in data["queries"]:
        x,y,z = query["point"]
        # Deliberately brute-force independent oracle, not production sparse lookup.
        containing = [key for key in expected_cells if all(key[a]/16 <= query["point"][a] <= (key[a]+1)/16 for a in range(3))]
        require(containing and query["cell_key"] == list(min(containing)), "closed lex-first cell")
        expected = [1,x,y,z,x*y,x*z,y*z,x*y*z,2]
        require(len(query["values"]) == len(expected), "field columns")
        for actual, want in zip(query["values"], expected):
            require(math.isfinite(actual) and abs(actual-want) <= 1e-12, "field reconstruction")
            maximum_error = max(maximum_error, abs(actual-want))
    require(len(data["queries"]) >= 5, "missing query checks")
    require(data["outside_rejected"] is True, "nextafter outside accepted")
    print(json.dumps({"path":str(path), "sha256":hashlib.sha256(raw).hexdigest(),
        "maximum_field_error":maximum_error,"actual_counts":[len(cells),len(nodes),len(edges),len(handles)],
        "scope":data["scope"]},sort_keys=True))
require(paths[0].read_bytes() == paths[1].read_bytes(), "repeated actual exports differ")
print("PASS: actual complete native grid, handle map, field transfer and repeated export checked independently")
