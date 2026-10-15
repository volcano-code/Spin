"""Exact diagnostic recovery of source mesh facets for Boolean output faces."""

from __future__ import annotations

from pathlib import Path


def recover_face_sources(input_meshes: tuple[Path, ...], output_mesh: Path,
                         exact_vertex_rows: list[list[str]],
                         face_label_rows: list[list[str]]) -> list[dict[str, object]]:
    """Return all exact input facets containing each output triangle.

    This is a diagnostic candidate recovery, not a unique derivative map: a
    face on coincident input sheets can legitimately have multiple sources.
    Coordinates from input OBJ files are converted through binary64 first to
    match the values consumed by the Cherchi reader.
    """
    from fractions import Fraction

    Vector = tuple[Fraction, Fraction, Fraction]

    def parse_binary_rational(token: str) -> Fraction:
        if token == "0/1":
            return Fraction(0)
        negative = token.startswith("-")
        if negative:
            token = token[1:]
        numerator, denominator = token.split("/")

        def integer(value: str) -> int:
            return int(value[2:], 2) if value.startswith("0b") else int(value)

        result = Fraction(integer(numerator), integer(denominator))
        return -result if negative else result

    def load_vertices_and_faces(path: Path) -> tuple[list[Vector], list[tuple[int, int, int]]]:
        vertices: list[Vector] = []
        faces: list[tuple[int, int, int]] = []
        for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "v":
                if len(fields) < 4:
                    raise ValueError(f"{path}:{line_number}: malformed vertex")
                vertices.append(tuple(Fraction.from_float(float(value)) for value in fields[1:4]))
            elif fields[0] == "f":
                if len(fields) != 4:
                    raise ValueError(f"{path}:{line_number}: expected a triangle")
                indices = tuple(int(value.split("/")[0]) - 1 for value in fields[1:4])
                if any(index < 0 or index >= len(vertices) for index in indices):
                    raise ValueError(f"{path}:{line_number}: unsupported or invalid face index")
                faces.append(indices)
        return vertices, faces

    def subtract(a: Vector, b: Vector) -> Vector:
        return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]

    def dot(a: Vector, b: Vector) -> Fraction:
        return sum((a[i] * b[i] for i in range(3)), Fraction(0))

    def cross(a: Vector, b: Vector) -> Vector:
        return (a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def contains(triangle: tuple[Vector, Vector, Vector], point: Vector) -> bool:
        a, b, c = triangle
        ab, ac, ap = subtract(b, a), subtract(c, a), subtract(point, a)
        if dot(cross(ab, ac), ap) != 0:
            return False
        d00, d01, d11 = dot(ab, ab), dot(ab, ac), dot(ac, ac)
        d20, d21 = dot(ap, ab), dot(ap, ac)
        denominator = d00 * d11 - d01 * d01
        if denominator == 0:
            return False
        u = (d11 * d20 - d01 * d21) / denominator
        v = (d00 * d21 - d01 * d20) / denominator
        return u >= 0 and v >= 0 and u + v <= 1

    if len(face_label_rows) != len(load_vertices_and_faces(output_mesh)[1]):
        raise ValueError("face-label rows do not align with output mesh faces")
    output_obj_vertices, output_faces = load_vertices_and_faces(output_mesh)
    if len(exact_vertex_rows) != len(output_obj_vertices):
        raise ValueError("exact coordinate rows do not align with output OBJ vertices")
    for index, row in enumerate(exact_vertex_rows):
        if len(row) < 11 or int(row[0]) != index or row[1] != "MATCHED":
            raise ValueError(f"exact coordinate row {index} is missing, reordered, or unsupported")
    if len(face_label_rows) != len(output_faces):
        raise ValueError("face-label rows do not align with output mesh faces")
    for index, row in enumerate(face_label_rows):
        if len(row) != 2 or int(row[0]) != index:
            raise ValueError(f"face-label row {index} is missing or reordered")
        if int(row[1], 16) == 0:
            raise ValueError(f"face-label row {index} has an empty input-object mask")
    exact_by_obj_coordinate: dict[Vector, Vector] = {}
    for row in exact_vertex_rows:
        if len(row) < 11:
            raise ValueError("exact vertex row lacks the OBJ coordinate key")
        obj_coordinate = tuple(Fraction.from_float(float(token)) for token in row[8:11])
        exact_coordinate = tuple(parse_binary_rational(token) for token in row[2:5])
        previous = exact_by_obj_coordinate.setdefault(obj_coordinate, exact_coordinate)
        if previous != exact_coordinate:
            raise ValueError("one OBJ coordinate maps to multiple exact coordinates")
    try:
        output_vertices = [exact_by_obj_coordinate[vertex] for vertex in output_obj_vertices]
    except KeyError as error:
        raise ValueError("output OBJ vertex is absent from exact-coordinate sidecar") from error
    input_data = [load_vertices_and_faces(path) for path in input_meshes]
    records: list[dict[str, object]] = []
    for face_index, (output_face, label_row) in enumerate(zip(output_faces, face_label_rows)):
        output_triangle = tuple(output_vertices[index] for index in output_face)
        label_mask = int(label_row[1], 16)
        sources: list[dict[str, int]] = []
        for object_index, (vertices, faces) in enumerate(input_data):
            for source_face_index, source_face in enumerate(faces):
                source_triangle = tuple(vertices[index] for index in source_face)
                if all(contains(source_triangle, point) for point in output_triangle):
                    sources.append({"input_object": object_index,
                                    "input_triangle": source_face_index})
        candidate_mask = 0
        for source in sources:
            candidate_mask |= 1 << source["input_object"]
        status = ("UNMATCHED" if not sources else
                  "LABEL_MISMATCH" if candidate_mask & label_mask == 0 else
                  "UNIQUE" if len(sources) == 1 else "AMBIGUOUS")
        records.append({"face_index": face_index,
                        "input_object_mask": label_mask,
                        "source_triangles": sources,
                        "status": status})
    return records
