#!/usr/bin/env bash
set -euo pipefail
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$(git -C "$task_dir" rev-parse --show-toplevel)"
task_build=build/task-5c-no-faster-20260909
cmake --build "$task_build" -j3 2>&1 | tee "$task_dir/task-5aa-controller-all-build.log"
ctest --test-dir "$task_build" --show-only=json-v1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)' \
  > "$task_dir/task-5aa-registered-tests.json"
/home/lin/anaconda3/bin/python -c '
import json,sys
tests=json.load(open(sys.argv[1]))["tests"]
assert len(tests)==48 and len({t["name"] for t in tests})==48
for path in sorted({t["command"][0] for t in tests}): print(path)
' "$task_dir/task-5aa-registered-tests.json" > "$task_dir/task-5aa-registered-program-paths.txt"
xargs -d '\n' sha256sum < "$task_dir/task-5aa-registered-program-paths.txt" \
  > "$task_dir/task-5aa-registered-test-binaries.sha256"
# The Python command is not the extension it loads: retain that binary as well.
awk '$2 ~ /\.so$/ {print $2}' "$task_dir/task-5z-registered-test-binaries.sha256" |
  xargs -r -d '\n' sha256sum >> "$task_dir/task-5aa-registered-test-binaries.sha256"
ctest --test-dir "$task_build" --output-on-failure -j1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)' 2>&1 |
  tee "$task_dir/task-5aa-controller-related-tests.log"
sha256sum --check "$task_dir/task-5aa-registered-test-binaries.sha256" \
  > "$task_dir/task-5aa-test-binaries-after.log"
