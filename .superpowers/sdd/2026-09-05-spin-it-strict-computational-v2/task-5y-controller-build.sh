#!/usr/bin/env bash
set -euo pipefail
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_repo=$(git -C "$task_script_dir" rev-parse --show-toplevel)
cd "$task_repo"
task_run=$(mktemp -d build/task-5y-controller-build-XXXXXX)
task_core_before=$(sha256sum build/task-5c-no-faster-20260909/libmake_it_stand_core.a)
git rev-parse HEAD > "$task_run/source-commit.txt"
git diff --binary | sha256sum > "$task_run/tracked-diff-hash.txt"
/usr/bin/time -v -o "$task_run/build-resource.log" /usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math \
  -DMAKE_IT_STAND_HAS_SPIN_IT_2014=1 -DMAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=1 \
  -Iinclude -Itests -isystem libigl/external/eigen \
  "$task_script_dir/task-5y-controller-export.cpp" -o "$task_run/export" \
  build/task-5c-no-faster-20260909/libmake_it_stand_core.a \
  /home/lin/anaconda3/lib/libboost_system.a /usr/lib/x86_64-linux-gnu/libgmp.so \
  /usr/lib/x86_64-linux-gnu/libmpfr.so /home/lin/anaconda3/lib/libboost_thread.a \
  build/task-5c-no-faster-20260909/lib/libhighs.a /usr/lib/x86_64-linux-gnu/libz.so \
  > "$task_run/build.log" 2>&1
test "$task_core_before" = "$(sha256sum build/task-5c-no-faster-20260909/libmake_it_stand_core.a)"
sha256sum "$task_run/export" > "$task_run/executable.sha256"
echo "$task_run"
