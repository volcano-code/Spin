# Provenance and license boundary

Last audited: 2026-08-11 (Asia/Shanghai)

## Sources and authority

1. **Normative specification:** R. Prévost et al., *Make It Stand: Balancing
   Shapes for 3D Fabrication*, ACM TOG 32(4), 2013,
   DOI `10.1145/2461912.2461957`.
2. **Behavioral oracle:** the authors' beta source archive from
   `https://igl.ethz.ch/projects/make-it-stand/make-it-stand_src_beta.zip`.
3. **Independent implementation baseline:**
   `https://github.com/helloxinw/make-it-balance`, commit
   `71e08361ca255a1cf6b89de6e4c92ef33024f432`.

When the paper and beta implementation differ, the paper governs the local
algorithm and the discrepancy is recorded in `conformance_matrix.md`.

The PDF used for the R0 audit has SHA-256:

```text
95e33982139150026af5abc3e3eb062e9142b651acaefa2d221df93237c7d01c
```

## Official beta archive

The archive downloaded for the R0 audit has SHA-256:

```text
f516ac32c43680108a158f559e0e3caa2709017b2888fb35b0d9c428a9d8b834
```

It is stored outside this Git repository under the user's dedicated
`reference-only` area and its extracted files are read-only. No symlink or file
from that archive is included here.

The archive's `readme_src.txt` states that the source is distributed for
research purposes only and that commercial use is prohibited. Consequently:

- official code is inspected only as a behavioral oracle;
- official C++ source, configurations, precomputed `.vox`/`.bbw` files,
  binaries and model assets are not copied into this MPL-2.0 repository;
- future fixture comparisons must use an optional external path and skip
  cleanly when the restricted archive is absent;
- any clean-room reimplementation must cite its mathematical source and must
  not copy protected implementation text.

## Local repository and submodule

The parent repository retains the Mozilla Public License 2.0 in `LICENSE`.
Its SHA-256 at the R0 audit is:

```text
60222c28c1a7f6a92c7df98e5c5f4459e624e6e285e0b9b94467af5f6ab3343d
```

The libigl submodule is pinned to
`f06fb13014aadff7eb1d7423ac018620441a3657`. It contains its own licensing
notices and optional dependencies. A release must include the notices required
by the libigl components actually enabled in that release.

## Inherited sample data

The three tracked sample meshes were inherited from the GitHub baseline:

| File | SHA-256 |
|---|---|
| `data/bunny.off` | `dda85ca6dc03ed83a0b754611d5c0b1a3c8ba5234e690479ac4a04063e6a7d4d` |
| `data/humpty_outer.stl` | `28e8a31869c34903445263656f1214aee6a3a3974fcc890d6ad7b9a3a3e8f5ed` |
| `data/spheres_outer.stl` | `5e2f3107903644c9928a8abd6f0206ede452b98b8177effa36a2999f622cd684` |

Their inclusion in the baseline is recorded, but this audit did not establish
independent redistribution rights. Recheck those rights before publishing a
new binary/data release.
