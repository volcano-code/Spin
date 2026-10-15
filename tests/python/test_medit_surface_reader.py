from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from make_it_stand.verification.medit_surface import read_medit_surface


class MeditSurfaceReaderTests(unittest.TestCase):
    def write_mesh(self, text: str) -> Path:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        path = Path(temporary.name) / "mesh.mesh"
        path.write_text(text, encoding="ascii")
        return path

    def test_extracts_triangle_surface_without_reordering_or_changing_coordinates(self):
        path = self.write_mesh("""MeshVersionFormatted 1
Dimension 3
Vertices
5
0 0 0 1
1 0 0 1
0 1 0 1
0 0 1 1
9 9 9 1
Triangles
4
1 3 2 1
1 2 4 1
2 3 4 1
3 1 4 1
Tetrahedra
1
1 2 3 4 1
End
""")

        surface = read_medit_surface(path)

        np.testing.assert_array_equal(surface.vertices, [
            [0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
        ])
        np.testing.assert_array_equal(surface.faces, [
            [0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3],
        ])
        self.assertEqual(surface.medit_vertex_count, 5)
        self.assertEqual(surface.surface_vertex_count, 4)
        self.assertEqual(surface.triangle_count, 4)
        self.assertEqual(surface.excluded_unreferenced_vertex_count, 1)

    def test_rejects_out_of_range_triangle_indices(self):
        path = self.write_mesh("""MeshVersionFormatted 1
Dimension 3
Vertices
3
0 0 0 1
1 0 0 1
0 1 0 1
Triangles
1
1 2 4 1
End
""")

        with self.assertRaisesRegex(ValueError, "triangle index outside vertex range"):
            read_medit_surface(path)

    def test_rejects_duplicate_required_sections(self):
        path = self.write_mesh("""MeshVersionFormatted 1
Dimension 3
Vertices
3
0 0 0 1
1 0 0 1
0 1 0 1
Vertices
3
0 0 0 1
1 0 0 1
0 1 0 1
Triangles
1
1 2 3 1
End
""")

        with self.assertRaisesRegex(ValueError, "exactly one Vertices and Triangles section"):
            read_medit_surface(path)


if __name__ == "__main__":
    unittest.main()
