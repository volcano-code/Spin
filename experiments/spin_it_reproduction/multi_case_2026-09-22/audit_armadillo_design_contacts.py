#!/usr/bin/env python3
"""Recompute the frozen three-face-layer Armadillo design mask and local contacts."""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
import hashlib
import json
from pathlib import Path


def read_sections(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    sections = {}
    i = 0
    while i < len(lines):
        token = lines[i].strip()
        if token in {"Vertices", "Tetrahedra"}:
            count = int(lines[i + 1].strip())
            rows = []
            for line in lines[i + 2 : i + 2 + count]:
                fields = line.split()
                rows.append(fields)
            if len(rows) != count:
                raise ValueError(f"truncated {token}: expected {count}, got {len(rows)}")
            sections[token] = rows
            i += count + 2
        else:
            i += 1
    if set(sections) != {"Vertices", "Tetrahedra"}:
        raise ValueError("Medit mesh lacks Vertices or Tetrahedra section")
    vertices = [tuple(map(float, row[:3])) for row in sections["Vertices"]]
    tetrahedra = [tuple(int(x) - 1 for x in row[:4]) for row in sections["Tetrahedra"]]
    return vertices, tetrahedra


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--output", type=Path, help="write complete deterministic JSON evidence")
    parser.add_argument("--summary", action="store_true", help="print counts instead of all witnesses")
    args = parser.parse_args()
    mesh_path = args.mesh.resolve()
    asset_hash = hashlib.sha256(mesh_path.read_bytes()).hexdigest()
    vertices, tetrahedra = read_sections(mesh_path)

    faces = defaultdict(list)
    for ti, tet in enumerate(tetrahedra):
        a, b, c, d = tet
        for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            faces[tuple(sorted(face))].append(ti)
    if any(len(adjacent) > 2 for adjacent in faces.values()):
        raise ValueError("source tetrahedral complex itself has a face with incidence >2")
    boundary_tets = {adjacent[0] for adjacent in faces.values() if len(adjacent) == 1}
    distance = [-1] * len(tetrahedra)
    queue = deque()
    for ti in sorted(boundary_tets):
        distance[ti] = 0
        queue.append(ti)
    adjacency = [[] for _ in tetrahedra]
    for attached in faces.values():
        if len(attached) == 2:
            a, b = attached
            adjacency[a].append(b)
            adjacency[b].append(a)
    while queue:
        current = queue.popleft()
        for nxt in adjacency[current]:
            if distance[nxt] == -1:
                distance[nxt] = distance[current] + 1
                queue.append(nxt)
    design = {i for i, layer in enumerate(distance) if layer >= 3}

    component = [-1] * len(tetrahedra)
    components = []
    for seed in sorted(design):
        if component[seed] != -1:
            continue
        cid = len(components)
        component[seed] = cid
        cells = []
        q = deque([seed])
        while q:
            current = q.popleft()
            cells.append(current)
            for nxt in adjacency[current]:
                if nxt in design and component[nxt] == -1:
                    component[nxt] = cid
                    q.append(nxt)
        components.append(cells)

    boundary_faces = []
    for face, incident in faces.items():
        selected = [ti for ti in incident if ti in design]
        if len(selected) == 1:
            boundary_faces.append((face, selected[0]))

    edge_faces = defaultdict(list)
    vertex_links = defaultdict(list)
    for face, ti in boundary_faces:
        for k in range(3):
            edge = tuple(sorted((face[k], face[(k + 1) % 3])))
            edge_faces[edge].append((face, ti))
        for k, center in enumerate(face):
            others = tuple(sorted((face[(k + 1) % 3], face[(k + 2) % 3])))
            vertex_links[center].append(others)

    nonmanifold_edges = []
    component_pair_contacts = defaultdict(int)
    for edge, incidents in sorted(edge_faces.items()):
        if len(incidents) <= 2:
            continue
        cids = sorted({component[ti] for _, ti in incidents})
        for i, ca in enumerate(cids):
            for cb in cids[i + 1 :]:
                component_pair_contacts[(ca, cb)] += 1
        p, q = (vertices[index] for index in edge)
        nonmanifold_edges.append(
            {
                "vertices_zero_based": list(edge),
                "endpoints": [list(p), list(q)],
                "length": sum((p[k] - q[k]) ** 2 for k in range(3)) ** 0.5,
                "boundary_triangle_incidence": len(incidents),
                "incident_design_components": cids,
                "incident_boundary_faces": [
                    {
                        "vertices_zero_based": list(face),
                        "design_tetrahedron_zero_based": ti,
                        "source_boundary_distance": distance[ti],
                        "face_component": component[ti],
                    }
                    for face, ti in sorted(incidents)
                ],
            }
        )

    bad_vertex_links = []
    for vertex, link_edges in sorted(vertex_links.items()):
        degree = defaultdict(int)
        graph = defaultdict(set)
        for a, b in link_edges:
            degree[a] += 1
            degree[b] += 1
            graph[a].add(b)
            graph[b].add(a)
        if not degree or any(d != 2 for d in degree.values()):
            bad_vertex_links.append({"vertex_zero_based": vertex, "link_node_degrees": sorted(degree.values())})
            continue
        seen = set()
        q = [next(iter(degree))]
        seen.add(q[0])
        while q:
            for nxt in graph[q.pop()]:
                if nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
        if len(seen) != len(degree):
            bad_vertex_links.append(
                {"vertex_zero_based": vertex, "link_node_degrees": sorted(degree.values()), "link_cycles": "multiple"}
            )

    result = {
        "schema": "armadillo-design-contact-topology-audit-v1",
        "asset_sha256": asset_hash,
        "rule": "face adjacency distance >= 3 from explicit source boundary",
        "source_vertices": len(vertices),
        "source_tetrahedra": len(tetrahedra),
        "source_boundary_tetrahedra": len(boundary_tets),
        "design_tetrahedra": len(design),
        "design_face_connected_components": len(components),
        "component_tetrahedron_counts": sorted(map(len, components), reverse=True),
        "largest_component_tetrahedra": max(map(len, components)),
        "design_boundary_triangles": len(boundary_faces),
        "design_boundary_nonmanifold_edges": len(nonmanifold_edges),
        "design_boundary_bad_vertex_links": len(bad_vertex_links),
        "nonmanifold_edges_with_multiple_face_components": sum(
            1 for edge in nonmanifold_edges if len(edge["incident_design_components"]) > 1
        ),
        "component_pair_edge_contacts": [
            {"components": list(pair), "shared_nonmanifold_edges": count}
            for pair, count in sorted(component_pair_contacts.items())
        ],
        "nonmanifold_edge_witnesses": nonmanifold_edges,
        "bad_vertex_link_witnesses": bad_vertex_links,
    }
    if result["design_tetrahedra"] != 11315 or result["design_boundary_nonmanifold_edges"] != 82:
        raise SystemExit("recomputed candidate no longer matches the frozen diagnostic counts")
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    if args.summary:
        summary = {key: value for key, value in result.items() if key not in {"nonmanifold_edge_witnesses", "bad_vertex_link_witnesses"}}
        summary["nonmanifold_edge_witness_count"] = len(result["nonmanifold_edge_witnesses"])
        summary["bad_vertex_link_count"] = len(result["bad_vertex_link_witnesses"])
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
