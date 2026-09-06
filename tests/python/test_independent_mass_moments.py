from __future__ import annotations

import math
import os
from pathlib import Path
import subprocess
import sys
import unittest

from make_it_stand.verification import (
    DensityRegion,
    GeometryValidationError,
    TriangleMesh,
    aggregate_mass_properties,
    compute_mass_properties,
    integrate_geometric_moments,
)


TETRAHEDRON = TriangleMesh(
    vertices=(
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 3.0),
    ),
    faces=((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)),
)


def box_mesh(minimum, maximum, *, inward=False) -> TriangleMesh:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    vertices = (
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    if inward:
        faces = tuple((face[0], face[2], face[1]) for face in faces)
    return TriangleMesh(vertices, faces)


def combine_meshes(*meshes: TriangleMesh) -> TriangleMesh:
    vertices = []
    faces = []
    for mesh in meshes:
        offset = len(vertices)
        vertices.extend(mesh.vertices)
        faces.extend(tuple(index + offset for index in face) for face in mesh.faces)
    return TriangleMesh(tuple(vertices), tuple(faces))


def transformed_mesh(mesh: TriangleMesh, transform) -> TriangleMesh:
    return TriangleMesh(tuple(transform(vertex) for vertex in mesh.vertices), mesh.faces)


class IndependentMassMomentOracleTests(unittest.TestCase):
    def assertVectorAlmostEqual(self, actual, expected, places=12):
        self.assertEqual(len(actual), len(expected))
        for got, want in zip(actual, expected):
            self.assertAlmostEqual(got, want, places=places)

    def assertMatrixAlmostEqual(self, actual, expected, places=12):
        self.assertEqual(len(actual), len(expected))
        for got_row, want_row in zip(actual, expected):
            self.assertVectorAlmostEqual(got_row, want_row, places=places)

    def test_hand_derived_tetrahedron_moments_and_inertia(self) -> None:
        geometric = integrate_geometric_moments(TETRAHEDRON)
        self.assertAlmostEqual(geometric.volume, 1.0, places=12)
        self.assertVectorAlmostEqual(
            geometric.first_moment, (1.0 / 4.0, 1.0 / 2.0, 3.0 / 4.0)
        )
        self.assertMatrixAlmostEqual(
            geometric.second_moment,
            (
                (1.0 / 10.0, 1.0 / 10.0, 3.0 / 20.0),
                (1.0 / 10.0, 2.0 / 5.0, 3.0 / 10.0),
                (3.0 / 20.0, 3.0 / 10.0, 9.0 / 10.0),
            ),
        )

        mass = compute_mass_properties(TETRAHEDRON)
        self.assertAlmostEqual(mass.mass, 1.0, places=12)
        self.assertVectorAlmostEqual(mass.center_of_mass, (1.0 / 4.0, 1.0 / 2.0, 3.0 / 4.0))
        self.assertMatrixAlmostEqual(
            mass.central_inertia,
            (
                (39.0 / 80.0, 1.0 / 40.0, 3.0 / 80.0),
                (1.0 / 40.0, 3.0 / 8.0, 3.0 / 40.0),
                (3.0 / 80.0, 3.0 / 40.0, 3.0 / 16.0),
            ),
        )

    def test_rigid_transform_rotates_non_diagonal_central_inertia(self) -> None:
        # A 90-degree rotation about z followed by translation (4, -3, 2).
        transformed = transformed_mesh(
            TETRAHEDRON,
            lambda vertex: (-vertex[1] + 4.0, vertex[0] - 3.0, vertex[2] + 2.0),
        )
        mass = compute_mass_properties(transformed)
        self.assertVectorAlmostEqual(mass.center_of_mass, (3.5, -2.75, 2.75))
        self.assertMatrixAlmostEqual(
            mass.central_inertia,
            (
                (3.0 / 8.0, -1.0 / 40.0, -3.0 / 40.0),
                (-1.0 / 40.0, 39.0 / 80.0, 3.0 / 80.0),
                (-3.0 / 40.0, 3.0 / 80.0, 3.0 / 16.0),
            ),
        )

    def test_large_translation_does_not_cancel_central_inertia(self) -> None:
        translated = transformed_mesh(
            TETRAHEDRON,
            lambda vertex: (
                vertex[0] + 1.0e8,
                vertex[1] - 2.0e8,
                vertex[2] + 3.0e8,
            ),
        )
        mass = compute_mass_properties(translated)
        self.assertVectorAlmostEqual(
            (
                mass.center_of_mass[0] - 1.0e8,
                mass.center_of_mass[1] + 2.0e8,
                mass.center_of_mass[2] - 3.0e8,
            ),
            (1.0 / 4.0, 1.0 / 2.0, 3.0 / 4.0),
            places=7,
        )
        self.assertMatrixAlmostEqual(
            mass.central_inertia,
            (
                (39.0 / 80.0, 1.0 / 40.0, 3.0 / 80.0),
                (1.0 / 40.0, 3.0 / 8.0, 3.0 / 40.0),
                (3.0 / 80.0, 3.0 / 40.0, 3.0 / 16.0),
            ),
            places=10,
        )

    def test_box_and_concentric_reversed_cavity_have_analytic_moments(self) -> None:
        solid = compute_mass_properties(box_mesh((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)))
        self.assertAlmostEqual(solid.mass, 8.0, places=12)
        self.assertMatrixAlmostEqual(
            solid.second_moment,
            ((8.0 / 3.0, 0.0, 0.0), (0.0, 8.0 / 3.0, 0.0), (0.0, 0.0, 8.0 / 3.0)),
        )

        hollow_mesh = combine_meshes(
            box_mesh((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)),
            box_mesh((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5), inward=True),
        )
        hollow = compute_mass_properties(hollow_mesh)
        self.assertAlmostEqual(hollow.mass, 7.0, places=12)
        self.assertVectorAlmostEqual(hollow.center_of_mass, (0.0, 0.0, 0.0))
        self.assertMatrixAlmostEqual(
            hollow.second_moment,
            ((31.0 / 12.0, 0.0, 0.0), (0.0, 31.0 / 12.0, 0.0), (0.0, 0.0, 31.0 / 12.0)),
        )
        self.assertMatrixAlmostEqual(
            hollow.central_inertia,
            ((31.0 / 6.0, 0.0, 0.0), (0.0, 31.0 / 6.0, 0.0), (0.0, 0.0, 31.0 / 6.0)),
        )

    def test_offset_reversed_cavity_shifts_center_and_inertia(self) -> None:
        mesh = combine_meshes(
            box_mesh((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)),
            box_mesh((-0.25, -0.5, -0.5), (0.75, 0.5, 0.5), inward=True),
        )
        mass = compute_mass_properties(mesh)
        self.assertAlmostEqual(mass.mass, 7.0, places=12)
        self.assertVectorAlmostEqual(mass.center_of_mass, (-1.0 / 28.0, 0.0, 0.0))
        self.assertMatrixAlmostEqual(
            mass.central_inertia,
            (
                (31.0 / 6.0, 0.0, 0.0),
                (0.0, 107.0 / 21.0, 0.0),
                (0.0, 0.0, 107.0 / 21.0),
            ),
        )

    def test_multiple_densities_are_applied_once_during_aggregation(self) -> None:
        left = box_mesh((-1.5, -0.5, -0.5), (-0.5, 0.5, 0.5))
        right = box_mesh((0.5, -0.5, -0.5), (1.5, 0.5, 0.5))
        mass = aggregate_mass_properties(
            (DensityRegion(left, 2.0), DensityRegion(right, 4.0))
        )
        self.assertAlmostEqual(mass.mass, 6.0, places=12)
        self.assertVectorAlmostEqual(mass.first_moment, (2.0, 0.0, 0.0))
        self.assertVectorAlmostEqual(mass.center_of_mass, (1.0 / 3.0, 0.0, 0.0))
        self.assertMatrixAlmostEqual(
            mass.central_inertia,
            ((1.0, 0.0, 0.0), (0.0, 19.0 / 3.0, 0.0), (0.0, 0.0, 19.0 / 3.0)),
        )

    def test_scale_and_density_covariance_is_stable_across_units(self) -> None:
        density = 2.5
        for scale in (1.0e-6, 1.0e6):
            with self.subTest(scale=scale):
                mesh = transformed_mesh(
                    TETRAHEDRON,
                    lambda vertex, scale=scale: tuple(scale * value for value in vertex),
                )
                mass = compute_mass_properties(mesh, density=density)
                self.assertAlmostEqual(mass.mass / (density * scale**3), 1.0, places=10)
                self.assertVectorAlmostEqual(
                    tuple(value / scale for value in mass.center_of_mass),
                    (1.0 / 4.0, 1.0 / 2.0, 3.0 / 4.0),
                    places=10,
                )
                inertia_factor = density * scale**5
                self.assertMatrixAlmostEqual(
                    tuple(
                        tuple(value / inertia_factor for value in row)
                        for row in mass.central_inertia
                    ),
                    (
                        (39.0 / 80.0, 1.0 / 40.0, 3.0 / 80.0),
                        (1.0 / 40.0, 3.0 / 8.0, 3.0 / 40.0),
                        (3.0 / 80.0, 3.0 / 40.0, 3.0 / 16.0),
                    ),
                    places=10,
                )

    def test_public_checker_import_does_not_require_native_extension(self) -> None:
        root = Path(__file__).resolve().parents[2]
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(root / "python")
        completed = subprocess.run(
            (
                sys.executable,
                "-c",
                "import sys; "
                "sys.modules['_make_it_stand'] = None; "
                "sys.modules['make_it_stand._make_it_stand'] = None; "
                "from make_it_stand.verification import TriangleMesh; "
                "print(TriangleMesh.__name__)",
            ),
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, msg=completed.stderr)
        self.assertEqual(completed.stdout.strip(), "TriangleMesh")


class IndependentMassMomentValidationTests(unittest.TestCase):
    def assertRejected(self, mesh: TriangleMesh, message_fragment: str) -> None:
        with self.assertRaisesRegex(GeometryValidationError, message_fragment):
            integrate_geometric_moments(mesh)

    def test_invalid_face_indices_and_triangles_are_rejected(self) -> None:
        vertices = TETRAHEDRON.vertices
        cases = {
            "out-of-range": TriangleMesh(vertices, ((0, 2, 4), *TETRAHEDRON.faces[1:])),
            "non-integer": TriangleMesh(vertices, ((0, 2, 1.5), *TETRAHEDRON.faces[1:])),
            "repeats": TriangleMesh(vertices, ((0, 2, 2), *TETRAHEDRON.faces[1:])),
            "duplicates": TriangleMesh(vertices, (*TETRAHEDRON.faces, (1, 2, 0))),
            "degenerate": TriangleMesh(
                (*vertices, (0.5, 0.0, 0.0)),
                ((0, 4, 1), *TETRAHEDRON.faces[1:]),
            ),
        }
        for message, mesh in cases.items():
            with self.subTest(message=message):
                self.assertRejected(mesh, message)

    def test_nonfinite_vertices_are_rejected(self) -> None:
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value):
                vertices = list(TETRAHEDRON.vertices)
                vertices[1] = (value, 0.0, 0.0)
                self.assertRejected(TriangleMesh(tuple(vertices), TETRAHEDRON.faces), "not finite")

    def test_open_boundary_and_inconsistent_winding_are_rejected(self) -> None:
        self.assertRejected(
            TriangleMesh(TETRAHEDRON.vertices, TETRAHEDRON.faces[:-1]),
            "not manifold",
        )
        flipped = list(TETRAHEDRON.faces)
        flipped[0] = (flipped[0][0], flipped[0][2], flipped[0][1])
        self.assertRejected(
            TriangleMesh(TETRAHEDRON.vertices, tuple(flipped)),
            "inconsistent face orientation",
        )

    def test_two_closed_fans_pinched_at_one_vertex_are_rejected(self) -> None:
        vertices = TETRAHEDRON.vertices + (
            (-1.0, 0.0, 0.0),
            (0.0, -2.0, 0.0),
            (0.0, 0.0, -3.0),
        )
        second_faces = tuple(
            (face[0] if face[0] == 0 else face[0] + 3,
             face[2] if face[2] == 0 else face[2] + 3,
             face[1] if face[1] == 0 else face[1] + 3)
            for face in TETRAHEDRON.faces
        )
        self.assertRejected(
            TriangleMesh(vertices, TETRAHEDRON.faces + second_faces),
            "non-manifold.*multiple fans",
        )

    def test_empty_zero_negative_material_and_bad_density_are_rejected(self) -> None:
        with self.assertRaisesRegex(GeometryValidationError, "no vertices"):
            integrate_geometric_moments(TriangleMesh((), ()))
        reversed_faces = tuple((a, c, b) for a, b, c in TETRAHEDRON.faces)
        self.assertRejected(
            TriangleMesh(TETRAHEDRON.vertices, reversed_faces),
            "zero or negative total oriented volume",
        )
        cancelled = combine_meshes(
            TETRAHEDRON,
            TriangleMesh(TETRAHEDRON.vertices, reversed_faces),
        )
        self.assertRejected(cancelled, "zero or negative total oriented volume")
        for density in (0.0, -1.0, math.nan, math.inf):
            with self.subTest(density=density):
                with self.assertRaisesRegex(GeometryValidationError, "density must be finite and positive"):
                    compute_mass_properties(TETRAHEDRON, density=density)
        with self.assertRaisesRegex(GeometryValidationError, "at least one density region"):
            aggregate_mass_properties(())


if __name__ == "__main__":
    unittest.main()
