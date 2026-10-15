#!/usr/bin/env python3
"""Bounded task batches for the approved complete rank-two source plan.

Only explicitly implemented modes are accepted. Exit zero is process evidence,
not source/material/derivative or paper acceptance. Never promotes a shared build.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
PYTHON = Path('/home/lin/anaconda3/bin/python3')
CMAKE = Path('/usr/bin/cmake')
CTEST = Path('/usr/bin/ctest')
NINJA = Path('/home/lin/anaconda3/bin/ninja')
PLAN = ROOT / 'docs/superpowers/plans/2026-10-05-rank-two-complete-source-implementation.md'
SPEC = ROOT / 'docs/superpowers/specs/2026-10-05-rank-two-complete-source-design.md'
FACTOR_BUILD = ROOT / '.superpowers/sdd/2026-10-05-rank-two-complete-source-implementation/flint-feasibility-build-02'
FACTOR_ROOT = FACTOR_BUILD / 'flint-3.3.1'
CERTIFICATION_TIMING_SOURCE = FACTOR_BUILD.parent / 'complete-finite-coverage-factor-retry-04/stdout.log'
FIBER_CAPTURE_AUDIT = FACTOR_BUILD.parent / 'weighted-caller-independent-audit-03'
MODES = {
    'actual-parent-affine': (600,24,['rank-two-actual-parent-affine-regression'],
        '^spin-faster-rank-two-actual-parent-affine$'),
    'actual-parent-affine-integration': (600,24,['rank-two-actual-parent-affine-regression',
        'rank-two-carrier-regression','rank-two-algebraic-regression','rank-two-inventory-material-mesh-regression'],
        '^spin-faster-rank-two-(actual-parent-affine|carrier-kernels|algebraic|inventory-material-mesh)$'),
    'complete-owner-material': (3600,24,['rank-two-complete-owner-material-regression'],
        '^spin-faster-rank-two-complete-owner-material$'),
    'inventory-material-integration': (600,24,['rank-two-inventory-material-mesh-regression',
        'rank-two-automatic-material-mesh-regression','rank-two-whole-interior-inventory-regression'],
        '^spin-faster-rank-two-(inventory-material-mesh|automatic-material-mesh|whole-interior-inventory)$'),
    'inventory-material-mesh': (600,24,['rank-two-inventory-material-mesh-regression'],
        '^spin-faster-rank-two-inventory-material-mesh$'),
    'whole-interior-inventory': (600,24,['rank-two-whole-interior-inventory-regression'],
        '^spin-faster-rank-two-whole-interior-inventory$'),
    'whole-interior-integration': (600,24,['rank-two-whole-interior-inventory-regression',
        'rank-two-interior-inventory-regression','rank-two-shared-finite-seams-regression'],
        '^spin-faster-rank-two-(whole-interior-inventory|interior-inventory|shared-finite-seams)$'),
    'automatic-material-mesh': (600,24,['rank-two-automatic-material-mesh-regression'],
        '^spin-faster-rank-two-automatic-material-mesh$'),
    'automatic-material-integration': (600,24,['rank-two-automatic-material-mesh-regression',
        'rank-two-interior-inventory-regression','rank-two-strict-material-mesh-regression',
        'rank-two-interior-mesh-regression','rank-two-material-lifts-regression'],
        '^spin-faster-rank-two-(automatic-material-mesh|interior-inventory|strict-material-mesh|interior-mesh|material-lifts)$'),
    'interior-inventory': (300,24,['rank-two-interior-inventory-regression'],
        '^spin-faster-rank-two-interior-inventory$'),
    'interior-inventory-integration': (600,24,['rank-two-interior-inventory-regression',
        'rank-two-moving-interior-regression','rank-two-strict-material-mesh-regression',
        'rank-two-material-lifts-regression','rank-two-interior-mesh-regression'],
        '^spin-faster-rank-two-(interior-inventory|moving-interior|strict-material-mesh|material-lifts|interior-mesh)$'),
    'strict-material-mesh': (300,24,['rank-two-strict-material-mesh-regression'],
        '^spin-faster-rank-two-strict-material-mesh$'),
    'strict-material-annulus': (600,24,['rank-two-strict-material-mesh-regression'],
        '^spin-faster-rank-two-strict-material-annulus$'),
    'strict-material-integration': (600,24,['rank-two-strict-material-mesh-regression',
        'rank-two-material-lifts-regression','rank-two-moving-mesh-edge-regression',
        'rank-two-interior-mesh-regression','rank-two-finite-mesh-scope-regression'],
        '^spin-faster-rank-two-(strict-material-mesh|material-lifts|moving-mesh-edge|interior-mesh|finite-mesh-scope)$'),
    'strict-material-complete': (600,24,['rank-two-strict-material-mesh-regression',
        'rank-two-material-lifts-regression','rank-two-moving-mesh-edge-regression',
        'rank-two-interior-mesh-regression','rank-two-finite-mesh-scope-regression'],
        '^spin-faster-rank-two-(strict-material-mesh|strict-material-annulus|material-lifts|moving-mesh-edge|interior-mesh|finite-mesh-scope)$'),
    'finite-mesh-scope': (300,24,['rank-two-finite-mesh-scope-regression'],
        '^spin-faster-rank-two-finite-mesh-scope$'),
    'finite-mesh-scope-integration': (300,24,['rank-two-finite-mesh-scope-regression',
        'rank-two-chord-sweep-regression','rank-two-covered-finite-mesh-regression',
        'rank-two-finite-cell-topology-regression'],
        '^spin-faster-rank-two-(finite-mesh-scope|chord-sweep|covered-finite-mesh|finite-cell-topology)$'),
    'finite-mesh-scope-affected': (300,24,['rank-two-finite-mesh-scope-regression',
        'rank-two-chord-incidence-regression','rank-two-registry-rotation-regression',
        'rank-two-covered-finite-mesh-regression','rank-two-finite-cell-topology-regression'],
        '^spin-faster-rank-two-(finite-mesh-scope|chord-incidence|registry-rotation|covered-finite-mesh|finite-cell-topology)$'),
    'complete-finite-mesh': (3600,24,['rank-two-complete-finite-mesh-regression'],
        '^spin-faster-rank-two-complete-finite-mesh$'),
    'strip-memo-consumer': (300,24,['rank-two-regular-registry-regression'],
        '^spin-faster-rank-two-regular-registry$'),
    'strip-memo-kernel': (300,24,['rank-two-strip-memo-regression','rank-two-arc-strip-regression'],
        '^spin-faster-rank-two-(strip-memo|arc-strip)$'),
    'strip-memo-integration': (300,24,['rank-two-strip-memo-regression','rank-two-arc-strip-regression',
        'rank-two-registry-edge-regression','rank-two-regular-registry-regression',
        'rank-two-refinement-reuse-regression','rank-two-complete-endpoint-binding-regression'],
        '^spin-faster-rank-two-(strip-memo|arc-strip|registry-edge|regular-registry|refinement-reuse|complete-endpoint-binding)$'),
    'weighted-root-comparison': (300,24,['rank-two-weighted-profile'],
        '^spin-faster-rank-two-weighted-profile$'),
    'weighted-root-optimized': (300,24,['rank-two-weighted-optimized-profile'],
        '^spin-faster-rank-two-weighted-optimized-profile$'),
    'native-fiber-resultant-comparison': (300,24,['rank-two-fiber-factor-resultant-profile','rank-two-fiber-factor-resultant-reverse-profile'],
        '^spin-faster-rank-two-fiber-factor-resultant(-reverse)?-profile$'),
    'phase-profiler-utility': (300,24,['rank-two-phase-profile-regression'],
        '^spin-faster-rank-two-phase-profile$'),
    'native-fiber-factor-comparison': (300,24,['rank-two-fiber-native-profile','rank-two-fiber-decomposed-profile','rank-two-fiber-optimized-profile'],
        '^spin-faster-rank-two-fiber-(native|decomposed|optimized)-profile$'),
    'native-fiber-optimized-profile': (300,24,['rank-two-fiber-optimized-profile'],
        '^spin-faster-rank-two-fiber-optimized-profile$'),
    'native-fiber-profile': (300,24,['rank-two-fiber-native-profile','rank-two-fiber-decomposed-profile'],
        '^spin-faster-rank-two-fiber-(native|decomposed)-profile$'),
    'coverage-timing-validation': (300,24,['rank-two-phase-profile-regression',
        'rank-two-complete-endpoint-binding-regression','rank-two-endpoint-timing-profile'],
        '^spin-faster-rank-two-(phase-profile|complete-endpoint-binding|endpoint-timing-profile)$'),
    'complete-finite-coverage-timing': (300,24,['rank-two-complete-finite-coverage-timing-profile'],
        '^spin-faster-rank-two-complete-finite-coverage-timing-profile$'),
    'native-certification-profile': (300,24,['rank-two-certification-profile'],
        '^spin-faster-rank-two-certification-profile$'),
    'native-handle-reuse': (300,24,['rank-two-carrier-regression'],
        '^spin-faster-rank-two-carrier-kernels$'),
    'native-handle-reuse-integration': (300,24,['rank-two-carrier-regression','rank-two-algebraic-regression',
        'rank-two-complete-endpoint-binding-regression','rank-two-refinement-reuse-regression'],
        '^spin-faster-rank-two-(carrier-kernels|algebraic|complete-endpoint-binding|refinement-reuse)$'),
    'refinement-reuse': (300,24,['rank-two-refinement-reuse-regression'],
        '^spin-faster-rank-two-refinement-reuse$'),
    'refinement-reuse-integration': (300,24,['rank-two-refinement-reuse-regression','rank-two-registry-reuse-regression',
        'rank-two-complete-endpoint-binding-regression','rank-two-finite-coverage-continuation-regression'],
        '^spin-faster-rank-two-(refinement-reuse|registry-reuse|complete-endpoint-binding|finite-coverage-continuation)$'),
    'registry-reuse': (300,24,['rank-two-registry-reuse-regression'],
        '^spin-faster-rank-two-registry-reuse$'),
    'complete-finite-coverage': (300,24,['rank-two-complete-finite-coverage-regression'],
        '^spin-faster-rank-two-complete-finite-coverage$'),
    'complete-endpoint-binding': (300,24,['rank-two-complete-endpoint-binding-regression'],
        '^spin-faster-rank-two-complete-endpoint-binding$'),
    'complete-endpoint-binding-integration': (300,24,['rank-two-complete-endpoint-binding-regression',
        'rank-two-complete-endpoint-import-regression','rank-two-finite-coverage-continuation-regression',
        'rank-two-mixed-intervals-regression'],
        '^spin-faster-rank-two-(complete-endpoint-binding|complete-endpoint-import|finite-coverage-continuation|mixed-intervals)$'),
    'complete-endpoint-import': (300,24,['rank-two-complete-endpoint-import-regression'],
        '^spin-faster-rank-two-complete-endpoint-import$'),
    'complete-endpoint-import-integration': (300,24,['rank-two-complete-endpoint-import-regression',
        'rank-two-projection-endpoint-import-regression','rank-two-all-seam-endpoints-regression',
        'rank-two-finite-endpoint-inventory-regression','rank-two-endpoint-neighborhood-regression',
        'rank-two-finite-coverage-continuation-regression','rank-two-arc-registry-regression'],
        '^spin-faster-rank-two-(complete-endpoint-import|projection-endpoint-import|all-seam-endpoints|finite-endpoint-inventory|endpoint-neighborhood|finite-coverage-continuation|arc-registry)$'),
    'projection-endpoint-import': (300,24,['rank-two-projection-endpoint-import-regression'],
        '^spin-faster-rank-two-projection-endpoint-import$'),
    'projection-endpoint-import-integration': (300,24,['rank-two-projection-endpoint-import-regression',
        'rank-two-all-seam-endpoints-regression','rank-two-finite-coverage-continuation-regression',
        'rank-two-arc-registry-regression'],
        '^spin-faster-rank-two-(projection-endpoint-import|all-seam-endpoints|finite-coverage-continuation|arc-registry)$'),
    'finite-endpoint-inventory': (300,24,['rank-two-finite-endpoint-inventory-regression'],
        '^spin-faster-rank-two-finite-endpoint-inventory$'),
    'finite-endpoint-inventory-integration': (300,24,['rank-two-finite-endpoint-inventory-regression',
        'rank-two-all-seam-endpoints-regression','rank-two-endpoint-neighborhood-regression',
        'rank-two-finite-coverage-continuation-regression','rank-two-arc-strip-regression',
        'rank-two-carrier-regression'],
        '^spin-faster-rank-two-(finite-endpoint-inventory|all-seam-endpoints|endpoint-neighborhood|finite-coverage-continuation|arc-strip|carrier-kernels)$'),
    'all-seam-endpoints': (300,24,['rank-two-all-seam-endpoints-regression'],
        '^spin-faster-rank-two-all-seam-endpoints$'),
    'all-seam-endpoints-integration': (300,24,['rank-two-all-seam-endpoints-regression',
        'rank-two-seam-endpoint-regression','rank-two-finite-coverage-continuation-regression'],
        '^spin-faster-rank-two-(all-seam-endpoints|seam-endpoint|finite-coverage-continuation)$'),
    'finite-coverage-continuation': (300,24,['rank-two-finite-coverage-continuation-regression'],
        '^spin-faster-rank-two-finite-coverage-continuation$'),
    'finite-coverage-continuation-integration': (300,24,['rank-two-finite-coverage-continuation-regression',
        'rank-two-finite-arc-coverage-regression','rank-two-seam-endpoint-regression'],
        '^spin-faster-rank-two-(finite-coverage-continuation|finite-arc-coverage|seam-endpoint)$'),
    'seam-endpoint': (300,24,['rank-two-seam-endpoint-regression'],
        '^spin-faster-rank-two-seam-endpoint$'),
    'seam-endpoint-integration': (300,24,['rank-two-seam-endpoint-regression',
        'rank-two-seam-registry-regression','rank-two-neighborhood-registry-regression',
        'rank-two-seam-refinement-regression'],
        '^spin-faster-rank-two-(seam-endpoint|seam-registry|neighborhood-registry|seam-refinement)$'),
    'seam-refinement': (300,24,['rank-two-seam-refinement-regression'],
        '^spin-faster-rank-two-seam-refinement$'),
    'seam-refinement-integration': (300,24,['rank-two-seam-refinement-regression',
        'rank-two-seam-registry-regression','rank-two-registry-subdivision-regression',
        'rank-two-mixed-intervals-regression'],
        '^spin-faster-rank-two-(seam-refinement|seam-registry|registry-subdivision|mixed-intervals)$'),
    'seam-registry-integration': (300,24,['rank-two-seam-registry-regression',
        'rank-two-arc-registry-regression','rank-two-shared-finite-seams-regression'],
        '^spin-faster-rank-two-(seam-registry|arc-registry|shared-finite-seams)$'),
    'seam-registry': (300,24,['rank-two-seam-registry-regression'],
                           '^spin-faster-rank-two-seam-registry$'),
    'shared-finite-seams': (300,24,['rank-two-shared-finite-seams-regression'],
                           '^spin-faster-rank-two-shared-finite-seams$'),
    'finite-cell-topology': (300,24,['rank-two-finite-cell-topology-regression'],
                           '^spin-faster-rank-two-finite-cell-topology$'),
    'algebraic-proof-reuse': (300,24,['rank-two-algebraic-regression'],
                           '^spin-faster-rank-two-algebraic$'),
    'annulus-all-edges': (300,24,['rank-two-annulus-edge-repair-regression'],
                       '^spin-faster-rank-two-annulus-all-edges$'),
    'annulus-asymmetric-edge': (300,24,['rank-two-annulus-edge-repair-regression'],
                       '^spin-faster-rank-two-annulus-asymmetric-edge$'),
    'annulus-edge-repair': (300,24,['rank-two-annulus-edge-repair-regression'],
                       '^spin-faster-rank-two-annulus-edge-repair$'),
    'material-lifts': (300,24,['rank-two-material-lifts-regression'],
                       '^spin-faster-rank-two-material-lifts$'),
    'moving-mesh-edge': (300,24,['rank-two-moving-mesh-edge-regression'],
                       '^spin-faster-rank-two-moving-mesh-edge$'),
    'interior-mesh': (300,24,['rank-two-interior-mesh-regression'],
                    '^spin-faster-rank-two-interior-mesh$'),
    'moving-interior': (300,24,['rank-two-moving-interior-regression'],
                    '^spin-faster-rank-two-moving-interior$'),
    'covered-finite-mesh': (300,24,['rank-two-covered-finite-mesh-regression'],
                    '^spin-faster-rank-two-covered-finite-mesh$'),
    'affine-dependency': (120,24,['rank-two-affine-dependency-regression'],
                    '^spin-faster-rank-two-affine-dependency$'),
    'finite-arc-coverage': (300,24,['rank-two-finite-arc-coverage-regression'],
                    '^spin-faster-rank-two-finite-arc-coverage$'),
    'finite-triangulation-predicate-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-finite-partition-regression','rank-two-finite-predicate-regression',
                    'rank-two-finite-triangulation-regression'],
                    '^spin-faster-rank-two-(algebraic|finite-partition|finite-predicate|finite-triangulation)$'),
    'finite-predicate': (120,24,['rank-two-finite-predicate-regression'],
                    '^spin-faster-rank-two-finite-predicate$'),
    'finite-triangulation': (300,24,['rank-two-finite-triangulation-regression'],
                    '^spin-faster-rank-two-finite-triangulation$'),
    'finite-triangulation-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-finite-partition-regression','rank-two-finite-triangulation-regression'],
                    '^spin-faster-rank-two-(algebraic|finite-partition|finite-triangulation)$'),
    'finite-partition': (300,24,['rank-two-finite-partition-regression'],
                    '^spin-faster-rank-two-finite-partition$'),
    'finite-partition-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-arc-registry-regression','rank-two-finite-partition-regression'],
                    '^spin-faster-rank-two-(algebraic|arc-registry|finite-partition)$'),
    'mixed-intervals': (300,24,['rank-two-mixed-intervals-regression'],
                    '^spin-faster-rank-two-mixed-intervals$'),
    'mixed-intervals-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-registry-subdivision-regression','rank-two-registry-edge-regression',
                    'rank-two-mixed-arc-refinement-regression','rank-two-mixed-intervals-regression'],
                    '^spin-faster-rank-two-(algebraic|registry-subdivision|registry-edge|mixed-arc-refinement|mixed-intervals)$'),
    'mixed-arc-refinement': (300,24,['rank-two-mixed-arc-refinement-regression'],
                    '^spin-faster-rank-two-mixed-arc-refinement$'),
    'mixed-arc-refinement-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-registry-subdivision-regression','rank-two-registry-edge-regression',
                    'rank-two-mixed-arc-refinement-regression'],
                    '^spin-faster-rank-two-(algebraic|registry-subdivision|registry-edge|mixed-arc-refinement)$'),
    'registry-subdivision': (300,24,['rank-two-registry-subdivision-regression'],
                    '^spin-faster-rank-two-registry-subdivision$'),
    'registry-subdivision-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-arc-registry-regression','rank-two-neighborhood-registry-regression',
                    'rank-two-registry-subdivision-regression'],
                    '^spin-faster-rank-two-(algebraic|arc-registry|neighborhood-registry|registry-subdivision)$'),
    'neighborhood-registry': (300,24,['rank-two-neighborhood-registry-regression'],
                    '^spin-faster-rank-two-neighborhood-registry$'),
    'neighborhood-registry-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-arc-registry-regression','rank-two-endpoint-neighborhood-regression',
                    'rank-two-neighborhood-registry-regression'],
                    '^spin-faster-rank-two-(algebraic|arc-registry|endpoint-neighborhood|neighborhood-registry)$'),
    'endpoint-neighborhood-integration': (300,24,['rank-two-algebraic-regression',
                    'rank-two-endpoint-fan-regression','rank-two-endpoint-neighborhood-regression'],
                    '^spin-faster-rank-two-(algebraic|endpoint-fan|endpoint-neighborhood)$'),
    'endpoint-neighborhood': (300,24,['rank-two-endpoint-neighborhood-regression'],
                             '^spin-faster-rank-two-endpoint-neighborhood$'),
    'chord-sweep': (300, 24, ['rank-two-chord-sweep-regression'],
                         '^spin-faster-rank-two-chord-sweep$'),
    'chord-sweep-integration': (300, 24, ['rank-two-algebraic-regression',
                    'rank-two-chord-incidence-regression', 'rank-two-registry-rotation-regression',
                    'rank-two-chord-sweep-regression'],
                    '^spin-faster-rank-two-(algebraic|chord-incidence|registry-rotation|chord-sweep)$'),
    'regular-registry-integration': (300, 24, ['rank-two-algebraic-regression',
                    'rank-two-registry-edge-regression', 'rank-two-arc-refinement-regression',
                    'rank-two-registry-rotation-regression', 'rank-two-regular-registry-regression'],
                    '^spin-faster-rank-two-(algebraic|registry-edge|arc-refinement|registry-rotation|regular-registry)$'),
    'regular-registry': (300, 24, ['rank-two-regular-registry-regression'],
                         '^spin-faster-rank-two-regular-registry$'),
    'registry-rotation': (300, 24, ['rank-two-registry-rotation-regression'],
                          '^spin-faster-rank-two-registry-rotation$'),
    'registry-rotation-integration': (300, 24, ['rank-two-algebraic-regression',
                    'rank-two-arc-registry-regression', 'rank-two-chord-incidence-regression',
                    'rank-two-registry-rotation-regression'],
                    '^spin-faster-rank-two-(algebraic|arc-registry|chord-incidence|registry-rotation)$'),
    'chord-incidence-integration': (300, 24, ['rank-two-algebraic-regression',
                                         'rank-two-carrier-regression',
                                         'rank-two-arc-strip-regression',
                                         'rank-two-arc-sampling-regression',
                                         'rank-two-endpoint-arc-regression',
                                         'rank-two-endpoint-fan-regression',
                                         'rank-two-arc-registry-regression',
                                         'rank-two-registry-edge-regression',
                                         'rank-two-arc-refinement-regression',
                                         'rank-two-chord-incidence-regression'],
                    '^spin-faster-rank-two-(algebraic|carrier-kernels|arc-strip|arc-sampling|endpoint-arc|endpoint-fan|arc-registry|registry-edge|arc-refinement|chord-incidence)$'),
    'chord-incidence': (300, 24, ['rank-two-chord-incidence-regression'],
                       '^spin-faster-rank-two-chord-incidence$'),
    'arc-refinement-integration': (300, 24, ['rank-two-algebraic-regression',
                                         'rank-two-carrier-regression',
                                         'rank-two-arc-strip-regression',
                                         'rank-two-arc-sampling-regression',
                                         'rank-two-endpoint-arc-regression',
                                         'rank-two-endpoint-fan-regression',
                                         'rank-two-arc-registry-regression',
                                         'rank-two-registry-edge-regression',
                                         'rank-two-arc-refinement-regression'],
                    '^spin-faster-rank-two-(algebraic|carrier-kernels|arc-strip|arc-sampling|endpoint-arc|endpoint-fan|arc-registry|registry-edge|arc-refinement)$'),
    'arc-refinement': (300, 24, ['rank-two-arc-refinement-regression'],
                     '^spin-faster-rank-two-arc-refinement$'),
    'registry-edge-integration': (300, 24, ['rank-two-algebraic-regression',
                                         'rank-two-carrier-regression',
                                         'rank-two-arc-strip-regression',
                                         'rank-two-arc-sampling-regression',
                                         'rank-two-endpoint-arc-regression',
                                         'rank-two-endpoint-fan-regression',
                                         'rank-two-arc-registry-regression',
                                         'rank-two-registry-edge-regression'],
                    '^spin-faster-rank-two-(algebraic|carrier-kernels|arc-strip|arc-sampling|endpoint-arc|endpoint-fan|arc-registry|registry-edge)$'),
    'registry-edge': (300, 24, ['rank-two-registry-edge-regression'],
                     '^spin-faster-rank-two-registry-edge$'),
    'arc-registry-integration': (300, 24, ['rank-two-algebraic-regression',
                                         'rank-two-carrier-regression',
                                         'rank-two-arc-strip-regression',
                                         'rank-two-arc-sampling-regression',
                                         'rank-two-endpoint-arc-regression',
                                         'rank-two-endpoint-fan-regression',
                                         'rank-two-arc-registry-regression'],
                    '^spin-faster-rank-two-(algebraic|carrier-kernels|arc-strip|arc-sampling|endpoint-arc|endpoint-fan|arc-registry)$'),
    'algebraic-degree40': (120, 24, ['rank-two-algebraic-regression'],
                     '^spin-faster-rank-two-algebraic-degree40$'),
    'arc-registry-profile': (120, 24, ['rank-two-arc-registry-regression'],
                     '^spin-faster-rank-two-arc-registry-profile$'),
    'arc-registry': (300, 24, ['rank-two-arc-registry-regression'],
                     '^spin-faster-rank-two-arc-registry$'),
    'endpoint-fan-integration': (300, 24, ['rank-two-carrier-regression',
                                         'rank-two-arc-strip-regression',
                                         'rank-two-arc-sampling-regression',
                                         'rank-two-endpoint-arc-regression',
                                         'rank-two-endpoint-fan-regression'],
                    '^spin-faster-rank-two-(carrier-kernels|arc-strip|arc-sampling|endpoint-arc|endpoint-fan)$'),
    'endpoint-fan': (300, 24, ['rank-two-endpoint-fan-regression'],
                     '^spin-faster-rank-two-endpoint-fan$'),
    'endpoint-integration': (300, 24, ['rank-two-carrier-regression',
                                     'rank-two-arc-strip-regression',
                                     'rank-two-arc-sampling-regression',
                                     'rank-two-endpoint-arc-regression'],
                    '^spin-faster-rank-two-(carrier-kernels|arc-strip|arc-sampling|endpoint-arc)$'),
    'endpoint-arc': (300, 24, ['rank-two-endpoint-arc-regression'],
                     '^spin-faster-rank-two-endpoint-arc$'),
    'arc-strip-mapping-profile': (120, 24, ['rank-two-arc-sampling-regression'],
                                '^spin-faster-rank-two-arc-strip-mapping-profile$'),
    'carrier-local-integration': (300, 24, ['rank-two-carrier-regression',
                                           'rank-two-arc-strip-regression',
                                           'rank-two-arc-sampling-regression'],
                            '^spin-faster-rank-two-(carrier-kernels|arc-strip|arc-sampling)$'),
    'arc-strip': (300, 24, ['rank-two-arc-strip-regression'],
                    '^spin-faster-rank-two-arc-strip$'),
    'carrier-foundations': (300, 24, ['rank-two-carrier-regression', 'rank-two-arc-sampling-regression'],
                            '^spin-faster-rank-two-(carrier-kernels|arc-sampling)$'),
    'arc-sampling-profile': (120, 24, ['rank-two-arc-sampling-regression'],
                             '^spin-faster-rank-two-arc-sampling-profile$'),
    'arc-sampling': (300, 24, ['rank-two-arc-sampling-regression'],
                     '^spin-faster-rank-two-arc-sampling$'),
    'carrier-kernels': (300, 24, ['rank-two-carrier-regression'],
                        '^spin-faster-rank-two-carrier-kernels$'),
    'edge-scopes': (300, 24, ['rank-two-frozen-edges-regression'],
                    '^spin-faster-rank-two-frozen-edges-scopes$'),
    'observer-profile': (120, 24, ['rank-two-whole-arrangement-regression'],
                         '^spin-faster-rank-two-whole-arrangement-profile$'),
    'observer': (300, 24, ['rank-two-whole-arrangement-regression'],
                 '^spin-faster-rank-two-whole-arrangement$'),
    'edge-stages': (300, 24, ['rank-two-frozen-edges-regression'],
                    '^spin-faster-rank-two-frozen-edges-stages$'),
    'edge-first': (900, 24, ['rank-two-frozen-edges-regression'],
                   '^spin-faster-rank-two-frozen-edges-first$'),
    'edge-all': (3600, 24, ['rank-two-frozen-edges-regression'],
                 '^spin-faster-rank-two-frozen-edges-all$'),
    'triangulator': (3600, 24, ['quadric-triangulator-regression'],
                     '^spin-faster-quadric-triangulator$'),
    'backend': (600, 4, ['rank-two-algebraic-backend-probe'],
                '^spin-faster-rank-two-backend$'),
    'native': (3600, 24, ['rank-two-source-input-regression',
                        'rank-two-algebraic-regression',
                        'rank-two-arrangement-regression',
                        'rank-two-normal-mesh-regression',
                        'rank-two-whole-arrangement-regression',
                        'rank-two-algebraic-backend-probe'],
               '^spin-faster-rank-two-(source-input|backend|algebraic|arrangement|normal-mesh|whole-arrangement)$'),
}
MODES.update({
    f'edge-state-{state}': (900, 24, ['rank-two-frozen-edges-regression'],
                          f'^spin-faster-rank-two-frozen-edges-state-{state}$')
    for state in (0, 1, 2, 18, 56, 57, 58, 59)
})


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def publish(path, data):
    with Path(path).open('x', encoding='utf-8') as out:
        json.dump(data, out, sort_keys=True, indent=2, allow_nan=False)
        out.write('\n')


def configuration(job,mode=None):
    values = {
        'CMAKE_MAKE_PROGRAM': NINJA, 'CMAKE_CXX_COMPILER': '/usr/bin/c++',
        'CMAKE_BUILD_TYPE': 'Release', 'CMAKE_EXPORT_COMPILE_COMMANDS': 'ON',
        'LIBIGL_SKIP_DOWNLOAD': 'ON', 'MAKE_IT_STAND_LIBIGL_DIR': ROOT / 'libigl',
        'MAKE_IT_STAND_BUILD_VIEWER': 'OFF', 'MAKE_IT_STAND_BUILD_REACHABILITY': 'OFF',
        'MAKE_IT_STAND_BUILD_SPIN_IT_2014': 'OFF', 'MAKE_IT_STAND_BUILD_SPIN_IT_FASTER': 'ON',
        'MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS': 'ON',
        'MAKE_IT_STAND_BUILD_PYTHON': 'OFF', 'BUILD_TESTING': 'ON',
        'Python_EXECUTABLE': PYTHON, 'Python3_EXECUTABLE': PYTHON,
        'Boost_DIR': '/home/lin/anaconda3/lib/cmake/Boost-1.82.0',
        'GMP_INCLUDE_DIR': '/usr/include/x86_64-linux-gnu',
        'GMP_LIBRARIES': '/usr/lib/x86_64-linux-gnu/libgmp.so',
        'MPFR_INCLUDE_DIR': '/usr/include',
        'MPFR_LIBRARIES': '/usr/lib/x86_64-linux-gnu/libmpfr.so',
        'MAKE_IT_STAND_OBSERVE_COMPILERS': 'ON',
        'MAKE_IT_STAND_COMPILER_EVIDENCE_DIR': job / 'compiler-observation',
        'MAKE_IT_STAND_COMPILER_TIMEOUT': '600',
        'MAKE_IT_STAND_RANK_TWO_FLINT_ROOT': FACTOR_ROOT,
        'MAKE_IT_STAND_BUILD_RANK_TWO_DIAGNOSTIC_PROFILE':
            'ON' if mode in ('weighted-root-comparison','weighted-root-optimized','native-certification-profile','native-fiber-profile','native-fiber-optimized-profile','native-fiber-factor-comparison','native-fiber-resultant-comparison') else 'OFF',
        'MAKE_IT_STAND_BUILD_RANK_TWO_COVERAGE_PROFILE':
            'ON' if mode in ('phase-profiler-utility','coverage-timing-validation','complete-finite-coverage-timing') else 'OFF',
    }
    return [str(CMAKE), '-S', str(ROOT), '-B', str(job / 'build'), '-G', 'Ninja',
            *[f'-D{k}={v}' for k, v in values.items()]]


def identities(mode=None):
    # Effective compiler dependencies are additionally captured by the observer.
    files = {str(p): p for p in (PLAN, SPEC, Path(__file__).resolve(), ROOT/'CMakeLists.txt')}
    if mode in ('actual-parent-affine','actual-parent-affine-integration'):
        for name in ('diagnosis.json','result.json','complete.json'):
            path=FACTOR_BUILD.parent/'r188-owner-failure-audit-07'/name
            files[str(path)]=path
    if mode=='native-certification-profile':
        files[str(CERTIFICATION_TIMING_SOURCE)]=CERTIFICATION_TIMING_SOURCE
    if mode in ('weighted-root-comparison','weighted-root-optimized','native-fiber-profile','native-fiber-optimized-profile','native-fiber-factor-comparison','native-fiber-resultant-comparison'):
        for name in ('config.json','result.json','complete.json','stdout.log','stderr.log','mathematical-capture.json'):
            path=FIBER_CAPTURE_AUDIT/name;files[str(path)]=path
    # The approved optional backend is part of these new strict builds. Bind
    # its actual archive, provenance receipt, generated headers and sources.
    for path in [FACTOR_BUILD/'build-receipt.json',FACTOR_ROOT/'libflint.a',
                 FACTOR_ROOT/'COPYING',FACTOR_ROOT/'COPYING.LESSER',
                 FACTOR_BUILD.parent/'flint-3.3.1.tar.xz']:
        files[str(path)]=path
    for path in (FACTOR_ROOT/'src').rglob('*'):
        if path.is_file() and path.suffix in {'.h','.c','.in'}:
            files[str(path)]=path
    normals = ROOT/'.superpowers/sdd/2026-09-29-seven-day-faster-sprint/rank-two-normal-boundary/diagnostic-01/bundle/normals.tsv'
    files[str(normals)] = normals
    edges = ROOT/'.superpowers/sdd/2026-09-29-seven-day-faster-sprint/rank-two-edge-certificate/green-03/bundle/native.tsv'
    files[str(edges)] = edges
    for folder, suffixes in [('cmake', {'.cmake'}), ('include', {'.h', '.hpp'}),
                             ('src', {'.cpp', '.h', '.hpp'}),
                             ('python/make_it_stand/verification', {'.py'})]:
        for path in (ROOT / folder).rglob('*'):
            if path.is_file() and path.suffix in suffixes:
                files[str(path)] = path
    for path in (ROOT / 'tests').glob('rank_two_*'):
        if path.is_file():
            files[str(path)] = path
    for fixture in (ROOT/'tests/fixtures').glob('rank_two_*'):
        if fixture.is_file():
            files[str(fixture)] = fixture
    legacy_test = ROOT / 'tests/quadric_triangulator_regression.cpp'
    files[str(legacy_test)] = legacy_test
    for path in [PYTHON, CMAKE, CTEST, NINJA, Path('/usr/bin/c++'),
                 Path('/usr/bin/ld'), Path('/usr/bin/ar'),
                 Path('/usr/lib/x86_64-linux-gnu/libgmp.so'),
                 Path('/usr/lib/x86_64-linux-gnu/libmpfr.so'),
                 ROOT/'libigl/CMakeLists.txt',
                 ROOT/'libigl/cmake/LibiglDownloadExternal.cmake']:
        files[str(path)] = path
    return files


def worker(mode, job):
    _, _, targets, pattern = MODES[mode]
    cfg = json.loads((job / 'config.json').read_text())
    snapshot = job / 'input-snapshot'
    snapshot.mkdir()
    saved = {}
    for entry in cfg['identity'].values():
        source = Path(entry['path'])
        if source.is_relative_to(ROOT):
            if source in saved:
                if saved[source] != entry['sha256']:
                    raise ValueError('CONFLICTING_INPUT_IDENTITY')
                continue
            destination = snapshot / source.relative_to(ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with source.open('rb') as src, destination.open('xb') as dst:
                for block in iter(lambda: src.read(1024 * 1024), b''):
                    dst.write(block)
            if digest(destination) != entry['sha256']:
                raise ValueError('STALE_SNAPSHOT_INPUT')
            saved[source] = entry['sha256']
    publish(job / 'git-context.json', {
        'head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'branch': subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
        'status': subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=ROOT,text=True),
        'plan_sha256': digest(PLAN), 'spec_sha256': digest(SPEC),
    })
    if mode in ('actual-parent-affine','actual-parent-affine-integration'):
        capture_dir=FACTOR_BUILD.parent/'r188-owner-failure-audit-07'
        result=json.loads((capture_dir/'result.json').read_text())
        complete=json.loads((capture_dir/'complete.json').read_text())
        if (result['status'],result['returncode'])!=('PROCESS_EXITED',0) or complete['result_sha256']!=digest(capture_dir/'result.json'):
            raise ValueError('ACTUAL_PARENT_CAPTURE_AUDIT_INVALID')
        capture=json.loads((capture_dir/'diagnosis.json').read_text())
        if capture['edge']!=[5,8] or capture['degree_product']!=1152:
            raise ValueError('ACTUAL_PARENT_CAPTURE_SCOPE_MISMATCH')
        # Only ORIGINAL immutable root certificates, never selected factors
        # or a solved midpoint, cross the native test boundary.
        with (job/'rank-two-r188-failed-affine.txt').open('x') as out:
            out.write('RANK_TWO_R188_FAILED_AFFINE_V1\n')
            for parent in capture['roots']:
                out.write(str(len(parent['coefficients']))+'\n'+'\n'.join(parent['coefficients'])+'\n')
                out.write(f"{parent['ordinal']}\n{parent['lower']}\n{parent['upper']}\n{parent['identity']}\n")
        if (job/'rank-two-r188-failed-affine.txt').read_bytes()!=(ROOT/'tests/fixtures/rank_two_r188_failed_affine.txt').read_bytes():
            raise ValueError('REGISTERED_ACTUAL_PARENT_FIXTURE_DIFFERS_FROM_CAPTURE')
    if mode in ('weighted-root-comparison','weighted-root-optimized','native-fiber-profile','native-fiber-optimized-profile','native-fiber-factor-comparison','native-fiber-resultant-comparison'):
        audit=json.loads((FIBER_CAPTURE_AUDIT/'result.json').read_text())
        if (audit['status'],audit['returncode'])!=('PROCESS_EXITED',0):
            raise ValueError('SELECTED_ACTUAL_FIBER_AUDIT_NOT_PASSED')
        if json.loads((FIBER_CAPTURE_AUDIT/'complete.json').read_text())['result_sha256']!=digest(FIBER_CAPTURE_AUDIT/'result.json'):
            raise ValueError('SELECTED_ACTUAL_FIBER_AUDIT_RECEIPT_MISMATCH')
        selected=json.loads((FIBER_CAPTURE_AUDIT/'mathematical-capture.json').read_text())
        if selected['schema']!='r143-caller-capture-audit-v1' or selected['selected']['product']<64:
            raise ValueError('SELECTED_ACTUAL_FIBER_SCOPE_MISMATCH')
        publish(job/'fiber-selection.json',{'schema':'r143-native-fiber-selection-v1',
            'audit_sha256':digest(FIBER_CAPTURE_AUDIT/'mathematical-capture.json'),
            'input_events_sha256':selected['events_sha256'],'selection':selected['selected'],
            'projection_y':selected['projection_y'],'scope':'SAVED_ACTUAL_NATIVE_BRANCH_NOT_WHOLE_COVERAGE'})
        with (job/'rank-two-fiber-input.txt').open('x') as out:
            out.write(f"RANK_TWO_ACTUAL_FIBER_V1\n{selected['saved_arc']}\n{selected['saved_branch']}\n{selected['saved_curve']}\n")
            for case in (selected['selected']['output'],selected['projection_y']):
                out.write(str(len(case['coefficients']))+'\n'+'\n'.join(case['coefficients'])+'\n')
                out.write(f"{case['ordinal']}\n{case['lower']}\n{case['upper']}\n{case['identity']}\n")
        if mode in ('weighted-root-comparison','weighted-root-optimized'):
            capture=selected['selected']
            with (job/'rank-two-weighted-input.txt').open('x') as out:
                out.write(f"RANK_TWO_ACTUAL_AFFINE_V1\n{capture['fraction']}\n{capture['bits']}\n")
                for case in (capture['left'],capture['right'],capture['output']):
                    out.write(str(len(case['coefficients']))+'\n'+'\n'.join(case['coefficients'])+'\n')
                    out.write(f"{case['ordinal']}\n{case['lower']}\n{case['upper']}\n{case['identity']}\n")
    if mode=='native-certification-profile':
        selected={};total=0
        with CERTIFICATION_TIMING_SOURCE.open() as source:
            for number,line in enumerate(source,1):
                if 'algebraic-root:validated-certificate ' not in line:
                    continue
                identity,ordinal,lower,upper,serialized=line.split('algebraic-root:validated-certificate ',1)[1].split()
                coefficients=serialized.rstrip(',').split(',');degree=len(coefficients)-1
                raw='rank-two-algebraic-root-v1;'+ordinal+';'+''.join(f'{len(c)}:{c};' for c in coefficients)
                if hashlib.sha256(raw.encode()).hexdigest()!=identity:
                    raise ValueError('CORRUPT_CAPTURED_ROOT_IDENTITY')
                selected[degree]={'line':number,'identity':identity,'ordinal':int(ordinal),
                    'lower':lower,'upper':upper,'coefficients':coefficients}
                total+=1
        if total!=8016 or sorted(selected)!=[1,2,4,8,16,32,64,72,144]:
            raise ValueError('MISSING_ACTUAL_DEGREE_CLASS')
        publish(job/'certification-selection.json',{
            'schema':'native-certification-timing-selection-v1',
            'policy':'LAST_ACTUAL_CERTIFICATE_OF_EVERY_OBSERVED_DEGREE',
            'source_sha256':digest(CERTIFICATION_TIMING_SOURCE),'total_certificates':total,
            'cases':[selected[k] for k in sorted(selected)]})
        with (job/'rank-two-profile-input.txt').open('x') as out:
            out.write(str(len(selected))+'\n')
            for degree in sorted(selected):
                case=selected[degree];out.write(str(len(case['coefficients']))+'\n')
                out.write('\n'.join(case['coefficients'])+'\n')
                out.write(f"{case['ordinal']}\n{case['lower']}\n{case['upper']}\n{case['identity']}\n")
    commands = [configuration(job,mode),
                [str(CMAKE), '--build', str(job/'build'), '--parallel', '4', '--target', *targets],
                [str(CTEST), '--test-dir', str(job/'build'), '--show-only=json-v1', '-R', pattern],
                [str(CTEST), '--test-dir', str(job/'build'), '--parallel', '1',
                 '--verbose', '--output-on-failure', '--no-tests=error', '-R', pattern]]
    publish(job / 'commands.json', commands)
    for number, command in enumerate(commands):
        if mode in ('actual-parent-affine','actual-parent-affine-integration') and number==1:
            import shutil
            shutil.copyfile(job/'rank-two-r188-failed-affine.txt',job/'build/rank-two-r188-failed-affine.txt')
        if mode in ('weighted-root-comparison','weighted-root-optimized') and number==1:
            import shutil
            shutil.copyfile(job/'rank-two-weighted-input.txt',job/'build/rank-two-weighted-profile-input.txt')
        if mode=='native-certification-profile' and number==1:
            import shutil
            shutil.copyfile(job/'rank-two-profile-input.txt',job/'build/rank-two-certification-profile-input.txt')
        if mode in ('native-fiber-profile','native-fiber-optimized-profile','native-fiber-factor-comparison','native-fiber-resultant-comparison') and number==1:
            import shutil
            shutil.copyfile(job/'rank-two-fiber-input.txt',job/'build/rank-two-fiber-profile-input.txt')
        print(json.dumps({'command': command}), flush=True)
        if number == 2:
            result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            with (job/'test-discovery.json').open('x') as out:
                out.write(result.stdout)
            print(result.stderr, file=sys.stderr, end='', flush=True)
            if result.returncode:
                return result.returncode
            discovered = json.loads(result.stdout)['tests']
            if mode=='strict-material-complete':
                expected={'spin-faster-rank-two-'+x for x in ('strict-material-mesh','strict-material-annulus',
                    'material-lifts','moving-mesh-edge','interior-mesh','finite-mesh-scope')}
                if len(discovered)!=len(expected) or {t['name'] for t in discovered}!=expected:
                    raise ValueError('TEST_DISCOVERY_SCOPE_MISMATCH')
            elif len(discovered) != len(targets):
                raise ValueError('TEST_DISCOVERY_SCOPE_MISMATCH')
        else:
            result = subprocess.run(command, cwd=ROOT)
            if result.returncode:
                return result.returncode
    binaries = {target: {'path': str(job/'build'/target),
                         'sha256': digest(job/'build'/target)} for target in targets}
    artifacts = {}
    if 'rank-two-complete-owner-material-regression' in targets:
        evidence=job/'build/rank-two-complete-owner-material-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise ValueError('MISSING_COMPLETE_OWNER_MATERIAL_EVIDENCE')
        artifacts['complete_owner_material']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-inventory-material-mesh-regression' in targets:
        evidence=job/'build/rank-two-inventory-material-mesh-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise ValueError('MISSING_INVENTORY_MATERIAL_EVIDENCE')
        artifacts['inventory_material_mesh']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-whole-interior-inventory-regression' in targets:
        evidence=job/'build/rank-two-whole-interior-inventory-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise ValueError('MISSING_WHOLE_INTERIOR_INVENTORY_EVIDENCE')
        artifacts['whole_interior_inventory']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-automatic-material-mesh-regression' in targets:
        evidence=job/'build/rank-two-automatic-material-mesh-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise ValueError('MISSING_AUTOMATIC_MATERIAL_EVIDENCE')
        artifacts['automatic_material_mesh']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-interior-inventory-regression' in targets:
        evidence=job/'build/rank-two-interior-inventory-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise ValueError('MISSING_INTERIOR_INVENTORY_EVIDENCE')
        artifacts['interior_inventory']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-mesh-scope-regression' in targets:
        evidence=job/'build/rank-two-finite-mesh-scope-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_FINITE_MESH_SCOPE_EVIDENCE')
        artifacts['finite_mesh_scope']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-strip-memo-regression' in targets:
        evidence=job/'build/rank-two-strip-memo-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_STRIP_MEMO_EVIDENCE')
        artifacts['strip_memo']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-refinement-reuse-regression' in targets:
        evidence=job/'build/rank-two-refinement-reuse-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_REFINEMENT_REUSE_EVIDENCE')
        artifacts['refinement_reuse']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-registry-reuse-regression' in targets:
        evidence=job/'build/rank-two-registry-reuse-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_REGISTRY_REUSE_EVIDENCE')
        artifacts['registry_reuse']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-complete-finite-coverage-regression' in targets:
        evidence=job/'build/rank-two-complete-finite-coverage-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_COMPLETE_FINITE_COVERAGE_EVIDENCE')
        artifacts['complete_finite_coverage']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-complete-finite-mesh-regression' in targets:
        evidence=job/'build/rank-two-complete-finite-mesh-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_COMPLETE_FINITE_MESH_EVIDENCE')
        artifacts['complete_finite_mesh']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-endpoint-timing-profile' in targets:
        evidence=job/'build/rank-two-endpoint-timing-profile-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_ACTUAL_ENDPOINT_TIMING_GEOMETRY')
        artifacts['profile_endpoint_binding']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-complete-finite-coverage-timing-profile' in targets:
        evidence=job/'build/rank-two-complete-finite-coverage-timing-profile-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_ACTUAL_WHOLE_TIMING_GEOMETRY')
        artifacts['profile_complete_finite_coverage']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-complete-endpoint-binding-regression' in targets:
        evidence=job/'build/rank-two-complete-endpoint-binding-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_COMPLETE_ENDPOINT_BINDING_EVIDENCE')
        artifacts['complete_endpoint_binding']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-complete-endpoint-import-regression' in targets:
        evidence=job/'build/rank-two-complete-endpoint-import-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_COMPLETE_ENDPOINT_IMPORT_EVIDENCE')
        artifacts['complete_endpoint_import']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-projection-endpoint-import-regression' in targets:
        evidence=job/'build/rank-two-projection-endpoint-import-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_PROJECTION_ENDPOINT_IMPORT_EVIDENCE')
        artifacts['projection_endpoint_import']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-endpoint-inventory-regression' in targets:
        evidence=job/'build/rank-two-finite-endpoint-inventory-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_FINITE_ENDPOINT_INVENTORY_EVIDENCE')
        artifacts['finite_endpoint_inventory']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-all-seam-endpoints-regression' in targets:
        evidence=job/'build/rank-two-all-seam-endpoints-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_ALL_SEAM_ENDPOINTS_EVIDENCE')
        artifacts['all_seam_endpoints']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-coverage-continuation-regression' in targets:
        evidence=job/'build/rank-two-finite-coverage-continuation-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_FINITE_COVERAGE_CONTINUATION_EVIDENCE')
        artifacts['finite_coverage_continuation']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-seam-endpoint-regression' in targets:
        evidence=job/'build/rank-two-seam-endpoint-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_SEAM_ENDPOINT_EVIDENCE')
        artifacts['seam_endpoint']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-seam-refinement-regression' in targets:
        evidence=job/'build/rank-two-seam-refinement-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_SEAM_REFINEMENT_EVIDENCE')
        artifacts['seam_refinement']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-seam-registry-regression' in targets:
        evidence=job/'build/rank-two-seam-registry-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_SEAM_REGISTRY_EVIDENCE')
        artifacts['seam_registry']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-shared-finite-seams-regression' in targets:
        evidence=job/'build/rank-two-shared-finite-seams-evidence.txt'
        if not evidence.is_file():
            raise RuntimeError('missing actual shared finite seam evidence')
        artifacts['shared_finite_seams']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-cell-topology-regression' in targets:
        evidence=job/'build/rank-two-finite-cell-topology-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:
            raise ValueError('MISSING_FINITE_CELL_TOPOLOGY_EVIDENCE')
        artifacts['finite_cell_topology']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-annulus-edge-repair-regression' in targets:
        evidence=job/('build/rank-two-annulus-all-edges-evidence.txt' if mode=='annulus-all-edges' else 'build/rank-two-annulus-asymmetric-edge-evidence.txt' if mode=='annulus-asymmetric-edge' else 'build/rank-two-annulus-edge-repair-evidence.txt')
        if not evidence.is_file() or not evidence.stat().st_size:raise RuntimeError('MISSING_ACTUAL_ANNULUS_EDGE_REPAIR_EVIDENCE')
        artifacts['annulus_edge_repair']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-strict-material-mesh-regression' in targets:
        evidence=job/('build/rank-two-strict-material-annulus-evidence.txt' if mode=='strict-material-annulus' else 'build/rank-two-strict-material-mesh-evidence.txt')
        if not evidence.is_file() or not evidence.stat().st_size:raise RuntimeError('MISSING_STRICT_MATERIAL_MESH_EVIDENCE')
        artifacts['strict_material_mesh']={'path':str(evidence),'sha256':digest(evidence)}
        if mode=='strict-material-complete':
            annulus=job/'build/rank-two-strict-material-annulus-evidence.txt'
            if not annulus.is_file() or not annulus.stat().st_size:raise RuntimeError('MISSING_STRICT_MATERIAL_ANNULUS_EVIDENCE')
            artifacts['strict_material_annulus']={'path':str(annulus),'sha256':digest(annulus)}
    if 'rank-two-material-lifts-regression' in targets:
        evidence=job/'build/rank-two-material-lifts-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise RuntimeError('MISSING_ACTUAL_MATERIAL_LIFT_EVIDENCE')
        artifacts['material_lifts']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-interior-mesh-regression' in targets:
        evidence=job/'build/rank-two-interior-mesh-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_INTERIOR_MESH_EVIDENCE')
        artifacts['interior_mesh']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-moving-mesh-edge-regression' in targets:
        evidence=job/'build/rank-two-moving-mesh-edge-evidence.txt'
        if not evidence.is_file() or not evidence.stat().st_size:raise RuntimeError('MISSING_MOVING_MESH_EDGE_EVIDENCE')
        artifacts['moving_mesh_edge']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-moving-interior-regression' in targets:
        evidence=job/'build/rank-two-moving-interior-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_MOVING_INTERIOR_EVIDENCE')
        artifacts['moving_interior']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-covered-finite-mesh-regression' in targets:
        evidence=job/'build/rank-two-covered-finite-mesh-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_COVERED_FINITE_MESH_EVIDENCE')
        artifacts['covered_finite_mesh']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-affine-dependency-regression' in targets:
        evidence=job/'build/rank-two-affine-dependency-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_AFFINE_DEPENDENCY_EVIDENCE')
        artifacts['affine_dependency']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-arc-coverage-regression' in targets:
        evidence=job/'build/rank-two-finite-arc-coverage-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_FINITE_ARC_COVERAGE_EVIDENCE')
        artifacts['finite_arc_coverage']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-predicate-regression' in targets:
        evidence=job/'build/rank-two-finite-predicate-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_FINITE_PREDICATE_EVIDENCE')
        artifacts['finite_predicate']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-triangulation-regression' in targets:
        evidence=job/'build/rank-two-finite-triangulation-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_FINITE_TRIANGULATION_EVIDENCE')
        artifacts['finite_triangulation']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-finite-partition-regression' in targets:
        evidence=job/'build/rank-two-finite-partition-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_FINITE_PARTITION_EVIDENCE')
        artifacts['finite_partition']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-mixed-intervals-regression' in targets:
        evidence=job/'build/rank-two-mixed-intervals-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_MIXED_INTERVALS_EVIDENCE')
        artifacts['mixed_intervals']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-mixed-arc-refinement-regression' in targets:
        evidence=job/'build/rank-two-mixed-arc-refinement-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_MIXED_ARC_REFINEMENT_EVIDENCE')
        artifacts['mixed_arc_refinement']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-registry-subdivision-regression' in targets:
        evidence=job/'build/rank-two-registry-subdivision-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_REGISTRY_SUBDIVISION_EVIDENCE')
        artifacts['registry_subdivision']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-neighborhood-registry-regression' in targets:
        evidence=job/'build/rank-two-neighborhood-registry-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_NEIGHBORHOOD_REGISTRY_EVIDENCE')
        artifacts['neighborhood_registry']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-endpoint-neighborhood-regression' in targets:
        evidence=job/'build/rank-two-endpoint-neighborhood-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size==0:
            raise ValueError('MISSING_ENDPOINT_NEIGHBORHOOD_EVIDENCE')
        artifacts['endpoint_neighborhood']={'path':str(evidence),'sha256':digest(evidence)}
    if 'rank-two-chord-sweep-regression' in targets:
        evidence = job/'build/rank-two-chord-sweep-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_CHORD_SWEEP_EVIDENCE')
        artifacts['chord_sweep'] = {'path': str(evidence), 'sha256': digest(evidence)}
        regular = Path(str(evidence)+'.regular')
        if not regular.is_file() or regular.stat().st_size == 0:
            raise ValueError('MISSING_SWEEP_REGULAR_PROOFS')
        artifacts['chord_sweep_regular'] = {'path': str(regular), 'sha256': digest(regular)}
    if 'rank-two-regular-registry-regression' in targets:
        evidence = job/'build/rank-two-regular-registry-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_REGULAR_REGISTRY_EVIDENCE')
        artifacts['regular_registry'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-registry-rotation-regression' in targets:
        evidence = job/'build/rank-two-registry-rotation-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_REGISTRY_ROTATION_EVIDENCE')
        artifacts['registry_rotation'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-chord-incidence-regression' in targets:
        evidence = job/'build/rank-two-chord-incidence-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_CHORD_INCIDENCE_EVIDENCE')
        artifacts['chord_incidence'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-arc-refinement-regression' in targets:
        evidence = job/'build/rank-two-arc-refinement-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_ARC_REFINEMENT_EVIDENCE')
        artifacts['arc_refinement'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-registry-edge-regression' in targets:
        evidence = job/'build/rank-two-registry-edge-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_REGISTRY_EDGE_EVIDENCE')
        artifacts['registry_edges'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-arc-registry-regression' in targets:
        name = 'rank-two-arc-registry-profile-evidence.txt' if mode == 'arc-registry-profile' else 'rank-two-arc-registry-evidence.txt'
        evidence = job/'build'/name
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_ARC_REGISTRY_EVIDENCE')
        artifacts['arc_registry'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-endpoint-fan-regression' in targets:
        evidence = job/'build/rank-two-endpoint-fan-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_ENDPOINT_FAN_EVIDENCE')
        artifacts['endpoint_fans'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-endpoint-arc-regression' in targets:
        evidence = job/'build/rank-two-endpoint-arc-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_ENDPOINT_ARC_EVIDENCE')
        artifacts['endpoint_arcs'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-arc-strip-regression' in targets:
        evidence = job/'build/rank-two-arc-strip-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_REGULAR_ARC_STRIP_EVIDENCE')
        artifacts['regular_strips'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-arc-sampling-regression' in targets:
        name = {'arc-sampling-profile': 'rank-two-arc-sampling-profile-evidence.txt',
                'arc-strip-mapping-profile': 'rank-two-arc-strip-mapping-profile-evidence.txt'
                }.get(mode, 'rank-two-arc-sampling-evidence.txt')
        evidence = job/'build'/name
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_ACTUAL_ARC_SAMPLING_EVIDENCE')
        artifacts['arc_samples'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-frozen-edges-regression' in targets:
        scope = mode.removeprefix('edge-') if mode.startswith('edge-state-') else {
            'edge-first': 'first', 'edge-all': 'all', 'edge-stages': 'stages',
            'edge-scopes': 'scopes'}[mode]
        evidence = job/f'build/rank-two-frozen-edges-{scope}.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_FROZEN_EDGE_CROSSCHECK')
        kind = {'stages': 'stage_diagnostic', 'scopes': 'scope_contract'}.get(scope, 'frozen_edges')
        artifacts[kind] = {'path': str(evidence), 'sha256': digest(evidence), 'scope': scope}
    if 'rank-two-arrangement-regression' in targets:
        evidence = job/'build/rank-two-global-arrangement-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_GLOBAL_ARRANGEMENT_EVIDENCE')
        artifacts['global_arrangements'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-normal-mesh-regression' in targets:
        evidence = job/'build/rank-two-normal-mesh-evidence.txt'
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_NORMAL_MESH_EVIDENCE')
        artifacts['normal_mesh'] = {'path': str(evidence), 'sha256': digest(evidence)}
    if 'rank-two-whole-arrangement-regression' in targets:
        name = 'rank-two-whole-arrangement-profile-evidence.txt' if mode == 'observer-profile' else 'rank-two-whole-arrangement-evidence.txt'
        evidence = job/'build'/name
        if not evidence.is_file() or evidence.stat().st_size == 0:
            raise ValueError('MISSING_WHOLE_ARRANGEMENT_EVIDENCE')
        artifacts['whole_arrangement'] = {'path': str(evidence), 'sha256': digest(evidence)}
    publish(job/'batch.json', {'schema': 'rank-two-source-build-batch-v1',
                              'mode': mode, 'binaries': binaries,
                              'artifacts': artifacts,
                              'paper_completion': 'NOT_GRANTED'})
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=MODES, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not args.output.is_absolute() or args.output != args.output.resolve():
        parser.error('output must be an absolute normalized new job directory')
    if args.worker:
        return worker(args.mode, args.output)
    sys.path.insert(0, str(ROOT / 'python'))
    from make_it_stand.verification.source_job import run_source_job
    seconds, gib, _, _ = MODES[args.mode]
    result = run_source_job(args.output,
        [str(PYTHON), str(Path(__file__).resolve()), '--worker', '--mode', args.mode,
         '--output', '{job}'], identities(args.mode),
        maximum_wall_ns=seconds * 10**9, maximum_memory_bytes=gib * 1024**3,
        maximum_disk_bytes=80 * 1024**3, minimum_start_free_bytes=100 * 1024**3,
        minimum_free_bytes=20 * 1024**3)
    print(json.dumps({'status': result.status, 'returncode': result.returncode,
                      'job': str(result.directory), 'elapsed_ns': result.elapsed_ns}), flush=True)
    return 0 if result.status == 'PROCESS_EXITED' and result.returncode == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
