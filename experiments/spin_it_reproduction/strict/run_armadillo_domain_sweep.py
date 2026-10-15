#!/usr/bin/env python3
"""Audit a pre-registered Armadillo tetrahedral core sweep.

This runner is deliberately an asset/topology diagnostic.  It never invokes
an optimizer and never turns a threshold that happens to look better into a
paper design domain.  The tested thresholds are fixed before reading the
result: tetrahedra at dual-graph distance >= 1..10 from the explicit source
boundary.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
import hashlib
import json
from pathlib import Path


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, indent=2, separators=(",", ": ")) + "\n"


def read_medit(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    sections: dict[str, list[list[str]]] = {}
    i = 0
    while i < len(lines):
        token = lines[i].strip()
        if token in {"Vertices", "Tetrahedra"}:
            count = int(lines[i + 1].strip())
            rows = [line.split() for line in lines[i + 2 : i + 2 + count]]
            if len(rows) != count:
                raise ValueError(f"truncated {token}")
            sections[token] = rows
            i += count + 2
        else:
            i += 1
    if set(sections) != {"Vertices", "Tetrahedra"}:
        raise ValueError("Medit mesh must contain Vertices and Tetrahedra")
    vertices = [tuple(map(float, row[:3])) for row in sections["Vertices"]]
    tetrahedra = [tuple(int(x) - 1 for x in row[:4]) for row in sections["Tetrahedra"]]
    return vertices, tetrahedra


def topology_for_threshold(
    tetrahedra: list[tuple[int, int, int, int]],
    faces: dict[tuple[int, int, int], list[int]],
    adjacency: list[list[int]],
    distance: list[int],
    threshold: int,
) -> dict[str, object]:
    selected = {i for i, layer in enumerate(distance) if layer >= threshold}
    component = [-1] * len(tetrahedra)
    components: list[list[int]] = []
    for seed in sorted(selected):
        if component[seed] >= 0:
            continue
        cid = len(components)
        component[seed] = cid
        q = deque([seed])
        cells: list[int] = []
        while q:
            current = q.popleft()
            cells.append(current)
            for nxt in adjacency[current]:
                if nxt in selected and component[nxt] < 0:
                    component[nxt] = cid
                    q.append(nxt)
        components.append(cells)

    edge_faces: dict[tuple[int, int], list[tuple[int, int, int]]] = defaultdict(list)
    vertex_links: dict[int, list[tuple[int, int]]] = defaultdict(list)
    boundary_faces = 0
    for face, incident in faces.items():
        if sum(tet in selected for tet in incident) != 1:
            continue
        boundary_faces += 1
        for k in range(3):
            edge = tuple(sorted((face[k], face[(k + 1) % 3])))
            edge_faces[edge].append(face)
        for k, center in enumerate(face):
            vertex_links[center].append(
                tuple(sorted((face[(k + 1) % 3], face[(k + 2) % 3])))
            )

    bad_edges = sum(len(incidents) > 2 for incidents in edge_faces.values())
    bad_vertices = 0
    for links in vertex_links.values():
        degree: dict[int, int] = defaultdict(int)
        graph: dict[int, set[int]] = defaultdict(set)
        for a, b in links:
            degree[a] += 1
            degree[b] += 1
            graph[a].add(b)
            graph[b].add(a)
        if not degree or any(value != 2 for value in degree.values()):
            bad_vertices += 1
            continue
        seen = {next(iter(degree))}
        todo = list(seen)
        while todo:
            current = todo.pop()
            for nxt in graph[current]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        if len(seen) != len(degree):
            bad_vertices += 1

    return {
        "threshold": threshold,
        "selected_tetrahedra": len(selected),
        "components": len(components),
        "largest_component_tetrahedra": max(map(len, components), default=0),
        "boundary_triangles": boundary_faces,
        "boundary_nonmanifold_edges": bad_edges,
        "boundary_bad_vertex_links": bad_vertices,
        "manifold_boundary": bad_edges == 0 and bad_vertices == 0,
    }


def audit(mesh: Path, minimum: int, maximum: int) -> dict[str, object]:
    vertices, tetrahedra = read_medit(mesh)
    faces: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for ti, (a, b, c, d) in enumerate(tetrahedra):
        for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            faces[tuple(sorted(face))].append(ti)
    if any(len(incident) > 2 for incident in faces.values()):
        raise ValueError("source tetrahedral complex has face incidence > 2")
    boundary = {incident[0] for incident in faces.values() if len(incident) == 1}
    adjacency = [[] for _ in tetrahedra]
    for incident in faces.values():
        if len(incident) == 2:
            a, b = incident
            adjacency[a].append(b)
            adjacency[b].append(a)
    distance = [-1] * len(tetrahedra)
    queue = deque(sorted(boundary))
    for tet in boundary:
        distance[tet] = 0
    while queue:
        current = queue.popleft()
        for nxt in adjacency[current]:
            if distance[nxt] < 0:
                distance[nxt] = distance[current] + 1
                queue.append(nxt)
    thresholds = [topology_for_threshold(tetrahedra, faces, adjacency, distance, threshold)
                  for threshold in range(minimum, maximum + 1)]
    return {
        "schema": "armadillo-design-domain-layer-sweep-v1",
        "classification": "SOURCE_MODEL_ONLY",
        "formal_paper_case_credit": False,
        "optimizer_invoked": False,
        "runner": {
            "path": str(Path(__file__).resolve()),
            "sha256": hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        },
        "asset": {
            "path": str(mesh),
            "sha256": hashlib.sha256(mesh.read_bytes()).hexdigest(),
            "vertices": len(vertices),
            "tetrahedra": len(tetrahedra),
            "source_boundary_tetrahedra": len(boundary),
        },
        "pre_registered_rule": {
            "distance": "dual face-adjacency distance from explicit source boundary",
            "thresholds": list(range(minimum, maximum + 1)),
            "no_component_deletion": True,
            "no_optimizer": True,
        },
        "thresholds": thresholds,
        "conclusion": (
            "BLOCKED_NO_MANIFOLD_CORE_IN_PRE_REGISTERED_RANGE"
            if not any(row["manifold_boundary"] for row in thresholds)
            else "REQUIRES_REVIEW_MANIFOLD_THRESHOLD_FOUND"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-threshold", type=int, default=1)
    parser.add_argument("--maximum-threshold", type=int, default=10)
    args = parser.parse_args()
    if args.minimum_threshold < 1 or args.maximum_threshold < args.minimum_threshold:
        parser.error("invalid pre-registered threshold range")
    mesh = args.mesh.resolve(strict=True)
    result = audit(mesh, args.minimum_threshold, args.maximum_threshold)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(canonical(result), encoding="utf-8")
    print(canonical(result), end="")
    return 0 if result["conclusion"] == "BLOCKED_NO_MANIFOLD_CORE_IN_PRE_REGISTERED_RANGE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
