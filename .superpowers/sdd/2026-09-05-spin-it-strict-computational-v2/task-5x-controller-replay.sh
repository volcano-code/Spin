#!/usr/bin/env bash
# Synthetic/known-solution controller only; no original-paper status promotion.
set -euo pipefail
task_script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_repo=$(git -C "$task_script_dir" rev-parse --show-toplevel)
cd "$task_repo"
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum --check "$task_script_dir/task-5x-controller-inputs.sha256" || {
    echo 'STALE_RESULT: frozen source/compiler/dependency/core/runtime differs' >&2
    exit 1
}
task_run_dir=$(mktemp -d build/task-5x-controller-replay-XXXXXX)
/usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math \
  -DMAKE_IT_STAND_HAS_SPIN_IT_2014=1 -DMAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=1 \
  -Iinclude -Itests -isystem libigl/external/eigen \
  "$task_script_dir/task-5x-controller-export.cpp" -o "$task_run_dir/export" \
  build/task-5c-no-faster-20260909/libmake_it_stand_core.a \
  /home/lin/anaconda3/lib/libboost_system.a /usr/lib/x86_64-linux-gnu/libgmp.so \
  /usr/lib/x86_64-linux-gnu/libmpfr.so /home/lin/anaconda3/lib/libboost_thread.a \
  build/task-5c-no-faster-20260909/lib/libhighs.a /usr/lib/x86_64-linux-gnu/libz.so
for task_repeat in 1 2; do
    /usr/bin/time -v timeout 3600 "$task_run_dir/export" \
      > "$task_run_dir/run$task_repeat.jsonl" 2> "$task_run_dir/resource$task_repeat.log"
done
cmp "$task_run_dir/run1.jsonl" "$task_run_dir/run2.jsonl"
/home/lin/anaconda3/bin/python "$task_script_dir/task-5x-controller-check.py" "$task_run_dir/run1.jsonl"
/home/lin/anaconda3/bin/python "$task_script_dir/task-5x-controller-negative.py" "$task_run_dir/run1.jsonl"
sha256sum "$task_run_dir/export" "$task_run_dir/run1.jsonl" "$task_run_dir/run2.jsonl"
echo "Evidence: $task_run_dir"
