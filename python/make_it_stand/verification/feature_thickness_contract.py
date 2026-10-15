"""Admission of source-level feature semantics for future wall-thickness proofs.

This module deliberately does not measure a wall.  It binds the semantics a
complete all-surface checker would need to the exact bytes of its source
surface, and refuses to turn a caller declaration into a manufacturing pass.
"""
from __future__ import annotations

import re

from .spin_common import canonical, digest, require, scalar, text_label


_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_AUTHORITIES = {"AUTHOR_ORIGINAL", "PUBLISHED_EXACT_SPECIFICATION"}
_ROLES = {"ORIGINAL_SURFACE", "CAVITY_SURFACE", "MATERIAL_ISLAND_SURFACE"}


def _identity(value, label):
    require(isinstance(value, dict) and set(value) == {"explicit_empty", "vertices", "faces"},
            label + ": expected mesh identity")
    require(value["explicit_empty"] is False, label + ": source surface cannot be empty")
    for name, dtype in (("vertices", "<f8"), ("faces", "<i4")):
        record = value[name]
        require(isinstance(record, dict) and set(record) == {"dtype", "shape", "sha256"},
                label + ": invalid " + name + " identity")
        require(record["dtype"] == dtype, label + ": unexpected " + name + " dtype")
        require(isinstance(record["shape"], list) and len(record["shape"]) == 2 and
                all(type(x) is int and x >= 0 for x in record["shape"]),
                label + ": invalid " + name + " shape")
        require(record["shape"][1] == 3, label + ": " + name + " must be triplets")
        require(isinstance(record["sha256"], str) and _DIGEST.fullmatch(record["sha256"]) is not None,
                label + ": invalid " + name + " digest")
    require(value["vertices"]["shape"][0] > 0 and value["faces"]["shape"][0] > 0,
            label + ": source surface cannot be empty")
    return canonical(value)


def _source_map(value, expected_identity):
    require(isinstance(value, dict), "source feature map: expected object")
    require(set(value) == {"schema", "surface", "source_authority", "physical_units", "patches",
                           "germ_equivalence", "coverage"},
            "source feature map: missing or unknown fields")
    require(value["schema"] == "feature-thickness-source-map-v1", "source feature map: unsupported schema")
    require(_identity(value["surface"], "source surface identity") == expected_identity,
            "source surface identity does not match the sealed surface")
    require(value["source_authority"] in _AUTHORITIES,
            "source feature map: source authority cannot support strict admission")

    units = value["physical_units"]
    require(isinstance(units, dict) and set(units) == {
        "model_length_unit", "model_to_mm", "minimum_thickness_mm", "process_authority"},
        "source feature map: invalid physical unit declaration")
    text_label(units["model_length_unit"], "model length unit")
    scalar(units["model_to_mm"], "model-to-mm", positive=True)
    scalar(units["minimum_thickness_mm"], "minimum thickness mm", positive=True)
    text_label(units["process_authority"], "process authority")

    require(isinstance(value["patches"], list) and bool(value["patches"]),
            "source feature map: nonempty patches required")
    patch_ids = set()
    for patch in value["patches"]:
        require(isinstance(patch, dict) and set(patch) == {"id", "source_ids", "role"},
                "source feature map: invalid patch")
        patch_id = text_label(patch["id"], "patch ID")
        require(patch_id not in patch_ids, "source feature map: duplicate patch ID")
        patch_ids.add(patch_id)
        require(patch["role"] in _ROLES, "source feature map: unsupported patch role")
        require(isinstance(patch["source_ids"], list) and bool(patch["source_ids"]) and
                all(type(x) is int and x >= 0 for x in patch["source_ids"]) and
                len(set(patch["source_ids"])) == len(patch["source_ids"]),
                "source feature map: invalid patch source IDs")

    require(isinstance(value["germ_equivalence"], list),
            "source feature map: germ equivalence must be a list")
    for group in value["germ_equivalence"]:
        require(isinstance(group, list) and len(group) >= 2 and
                all(isinstance(germ, str) and ":" in germ for germ in group),
                "source feature map: invalid germ equivalence group")
        require(len(set(group)) == len(group), "source feature map: duplicate germ")
        require(all(germ.split(":", 1)[0] in patch_ids for germ in group),
                "source feature map: germ names unknown patch")

    coverage = value["coverage"]
    require(isinstance(coverage, dict) and coverage == {
        "boundary_germs": "DECLARED_COMPLETE",
        "opposing_pair_domains": "DECLARED_COMPLETE",
    }, "source feature map: incomplete coverage declaration")
    return canonical(value)


def bind_feature_thickness_source_map(surface_identity, source_map):
    """Bind a declaration to a source surface without granting a thickness pass.

    ``None`` is the only representation for a missing map.  A structurally
    valid map is merely eligible for a future independent geometric proof;
    neither its coverage declaration nor a caller-supplied review flag can
    promote ``complete_wall_thickness`` from ``UNKNOWN``.
    """
    surface = _identity(surface_identity, "sealed surface identity")
    if source_map is None:
        result = {
            "schema": "feature-thickness-source-map-attachment-v1",
            "surface": surface,
            "state": "UNKNOWN_MISSING_SOURCE_FEATURE_MAP",
            "eligible_for_geometric_proof": False,
            "complete_wall_thickness": "UNKNOWN",
            "review_status": "NOT_SUBMITTED",
            "source_map": None,
            "minimum_thickness_mm": None,
        }
    else:
        source = _source_map(source_map, surface)
        result = {
            "schema": "feature-thickness-source-map-attachment-v1",
            "surface": surface,
            "state": "PENDING_INDEPENDENT_REVIEW",
            "eligible_for_geometric_proof": True,
            "complete_wall_thickness": "UNKNOWN",
            "review_status": "NOT_APPROVED",
            "source_map": source,
            "minimum_thickness_mm": float(source_map["physical_units"]["minimum_thickness_mm"]),
        }
    result["attachment_digest"] = digest(result)
    return result
