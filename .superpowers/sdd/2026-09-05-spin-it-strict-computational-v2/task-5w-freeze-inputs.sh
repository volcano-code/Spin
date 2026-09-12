#!/usr/bin/env bash
# Generate a fresh local input manifest, never a success/acceptance flag.
set -euo pipefail
task_script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_repo=$(git -C "$task_script_dir" rev-parse --show-toplevel)
cd "$task_repo"
task_manifest=${1:?Specify a new manifest output path}
test ! -e "$task_manifest"
task_list=$(mktemp build/task-5w-input-list-XXXXXX)
# Actual compiler-discovered archive inputs, including transitive header files.
ninja -C build/task-5c-no-faster-20260909 -t deps |
  awk '/^CMakeFiles\/make_it_stand_core.dir\// {active=1; next}
       /^[^ ]/ {active=0}
       active && /^    / {sub(/^    /, ""); print}' |
  sort -u > "$task_list"
test -s "$task_list"
xargs -d '\n' sha256sum < "$task_list" > "$task_manifest"
sha256sum CMakeLists.txt build/task-5c-no-faster-20260909/CMakeCache.txt \
  build/task-5c-no-faster-20260909/libmake_it_stand_core.a \
  build/task-5c-no-faster-20260909/lib/libhighs.a \
  /usr/bin/c++ /home/lin/anaconda3/bin/python \
  /home/lin/anaconda3/lib/libboost_system.a /home/lin/anaconda3/lib/libboost_thread.a \
  /usr/lib/x86_64-linux-gnu/libgmp.so /usr/lib/x86_64-linux-gnu/libmpfr.so \
  /usr/lib/x86_64-linux-gnu/libz.so tests/source_surface_test_utils.h \
  "$task_script_dir/task-5w-controller-export.cpp" \
  "$task_script_dir/task-5w-controller-export-helpers.h" \
  "$task_script_dir/task-5w-controller-check.py" \
  "$task_script_dir/task-5w-controller-negative.py" \
  "$task_script_dir/task-5u-controller-check.py" \
  "$task_script_dir/task-5w-controller-replay.sh" \
  "$task_script_dir/task-5w-freeze-inputs.sh" >> "$task_manifest"
echo "Captured input manifest: $task_manifest"
echo "Retained compiler dependency list: $task_list"
