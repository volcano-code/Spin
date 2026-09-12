# Build Graph assisted Task5w boundary selection

User explicitly requested Build Graph. Its MCP tools were not callable in this
session. Used the already installed local CLI (same graph package), without
installing dependencies or changing tool configuration:

`/home/lin/.local/share/code-review-graph/venv/bin/code-review-graph`

Initial status -> incremental update -> status: 410 files, 5,569 nodes, 59,303
edges; languages bash, c, cpp, python. Update reported 0 indexed files changed,
0 nodes and 0 edges changed. The graph metadata still named commit `59f30cc`,
while checkout HEAD was `60b31e0`. This is reported, not rewritten as current
scientific evidence. A first CLI invocation used the unsupported `--repo-root`
flag and failed with a usage message; the correct CLI flag is `--repo`.

Saved raw queries:

- task-5w-graph-before.json: graph/branch/commit metadata.
- task-5w-graph-adaptive-callees.json: exact file-qualified legacy adaptive endpoint,
  34 extracted callees, including generic SLQP and split/merge calls. Many C++
  members remain unresolved or ambiguous. This is a navigation hint, not proof
  of complete call-graph resolution.
- task-5w-graph-paper-importers.json: 0 resolved importers. Direct source inspection
  finds actual includes in SourceCageSlqp.h and SourceCageGeometry.h, so the empty
  graph result cannot establish absence of consumers.

Direct source confirmation: SourceAdaptiveSpinSolver.cpp constructs
AdaptiveFillFunctional, then calls solveSlqpLbfgsPaper, splits/merges a local
SourceAdaptiveDomain, and uses rest-source binary moment checks. It does not
rebind SourceCagePaperFunctional or call the actual geometry guarded solver.
SourceCageBinding retains an immutable domain snapshot plus full actual BBW.
SourceCageMomentFunctional and SourceCagePaperFunctional have partition-dependent
dimensions, derivatives and hashes. This supports a separate native transaction
reconstructing that chain; modifying the old fill-only result would silently
change legacy meaning.

Chosen verification boundary: new immutable source + remapped state -> fresh
actual geometry check -> fresh guarded solve. No old gradient, KKT, history or
passed-geometry fact is inherited. Rest partition conservation and represented
deformed moment conservation are explicitly distinct.

Graph refresh after implementation and final statistics will be recorded in the
controller verification report. Graph statistics are never paper acceptance.

## Incremental discovery limitation identified

While new source files were still untracked, an additional update again reported
zero parsed changes. Inspection of the installed package's incremental.py
get_changed_files (line838) confirms it uses `git diff --name-status -z BASE --`,
not untracked-file status. Existing modified CMake input is not a parsed language.
The CLI update does not expose an explicit changed-file list. No package/config
fix is required: refresh after the normal audited source commit, when new source
paths enter the diff. Do not stage unverified files solely to make graph counts
increase, and do not interpret this interim zero as evidence that new code exists
in the index.

After implementation commit `434be7d0a9a323048e058f18bb913cf9ee1a7bb5`, the normal
incremental update parsed 3 changed files (56 nodes, 867 edges reported by the
update). Final total: 413 files, 5,624 nodes and 60,014 edges. Graph metadata and
checkout commit match. Raw committed-update/status files preserve these values.
No new dependency installation, cloud embeddings or tool configuration change.
