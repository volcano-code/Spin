#!/usr/bin/env bash
# Run from the canonical project root after building the scoped native library.
set -euo pipefail
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_evidence=.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2
task_build=build/task-5c-no-faster-20260909
task_run=$(mktemp -d build/task-5u-controller-replay-XXXXXX)
/usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math \
  -DMAKE_IT_STAND_HAS_SPIN_IT_2014=1 -DMAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=1 \
  -Iinclude -Itests -isystem libigl/external/eigen \
  "$task_evidence/task-5u-controller-export.cpp" -o "$task_run/export" \
  "$task_build/libmake_it_stand_core.a" /home/lin/anaconda3/lib/libboost_system.a \
  /usr/lib/x86_64-linux-gnu/libgmp.so /usr/lib/x86_64-linux-gnu/libmpfr.so \
  /home/lin/anaconda3/lib/libboost_thread.a "$task_build/lib/libhighs.a" \
  /usr/lib/x86_64-linux-gnu/libz.so
sha256sum "$task_build/libmake_it_stand_core.a" "$task_run/export" \
  "$task_evidence/task-5u-controller-export.cpp" \
  "$task_evidence/task-5u-controller-check.py" \
  "$task_evidence/task-5u-controller-negative.py" > "$task_run/inputs.sha256"
for task_repeat in 1 2; do
  timeout 300 "$task_run/export" > "$task_run/run$task_repeat.jsonl" 2> "$task_run/native$task_repeat.log"
  python "$task_evidence/task-5u-controller-check.py" "$task_run/run$task_repeat.jsonl" | tee "$task_run/check$task_repeat.log"
done
cmp "$task_run/run1.jsonl" "$task_run/run2.jsonl"
python "$task_evidence/task-5u-controller-negative.py" "$task_run/run1.jsonl" | tee "$task_run/negative.log"
sha256sum "$task_run/run1.jsonl" "$task_run/run2.jsonl" > "$task_run/outputs.sha256"
printf '%s\n' "$task_run"
