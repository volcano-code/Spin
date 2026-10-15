"""Source-bound feature-map admission, before any complete thickness claim."""
from __future__ import annotations

import importlib
import unittest


class FeatureThicknessAttachmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = importlib.import_module(
            "make_it_stand.verification.feature_thickness_contract"
        )

    @staticmethod
    def identity():
        return {
            "explicit_empty": False,
            "vertices": {"dtype": "<f8", "shape": [8, 3], "sha256": "a" * 64},
            "faces": {"dtype": "<i4", "shape": [12, 3], "sha256": "b" * 64},
        }

    def source_map(self):
        return {
            "schema": "feature-thickness-source-map-v1",
            "surface": self.identity(),
            "source_authority": "PUBLISHED_EXACT_SPECIFICATION",
            "physical_units": {
                "model_length_unit": "mm",
                "model_to_mm": 1.0,
                "minimum_thickness_mm": 1.2,
                "process_authority": "DECLARED_PROCESS_SPECIFICATION",
            },
            "patches": [
                {"id": "outer", "source_ids": [0, 1], "role": "ORIGINAL_SURFACE"},
                {"id": "cavity", "source_ids": [2], "role": "CAVITY_SURFACE"},
            ],
            "germ_equivalence": [["outer:ridge-left", "outer:ridge-right"]],
            "coverage": {
                "boundary_germs": "DECLARED_COMPLETE",
                "opposing_pair_domains": "DECLARED_COMPLETE",
            },
        }

    def test_missing_map_cannot_admit_complete_feature_thickness(self):
        """Removing source semantics must leave an artifact unqualified."""
        result = self.module.bind_feature_thickness_source_map(self.identity(), None)
        self.assertEqual(result["state"], "UNKNOWN_MISSING_SOURCE_FEATURE_MAP")
        self.assertFalse(result["eligible_for_geometric_proof"])
        self.assertEqual(result["complete_wall_thickness"], "UNKNOWN")

    def test_declared_map_is_bound_but_not_self_promoted_to_a_pass(self):
        """A caller's map declaration cannot substitute all-surface proof/review."""
        result = self.module.bind_feature_thickness_source_map(
            self.identity(), self.source_map()
        )
        self.assertEqual(result["state"], "PENDING_INDEPENDENT_REVIEW")
        self.assertTrue(result["eligible_for_geometric_proof"])
        self.assertEqual(result["complete_wall_thickness"], "UNKNOWN")
        self.assertEqual(result["minimum_thickness_mm"], 1.2)

    def test_map_for_another_surface_is_rejected_before_any_claim(self):
        """Changing source bytes must invalidate the attached feature semantics."""
        other = self.identity()
        other["vertices"]["sha256"] = "c" * 64
        with self.assertRaisesRegex(ValueError, "source surface identity"):
            self.module.bind_feature_thickness_source_map(other, self.source_map())

    def test_public_verification_api_exports_the_fail_closed_binder(self):
        """Removing the public binder would leave callers no safe admission path."""
        verification = importlib.import_module("make_it_stand.verification")
        self.assertIs(
            verification.bind_feature_thickness_source_map,
            self.module.bind_feature_thickness_source_map,
        )


if __name__ == "__main__":
    unittest.main()
