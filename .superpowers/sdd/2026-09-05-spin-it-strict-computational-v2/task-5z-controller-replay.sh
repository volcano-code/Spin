#!/usr/bin/env bash
set -euo pipefail
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_repo=$(git -C "$task_script_dir" rev-parse --show-toplevel)
cd "$task_repo"
task_manifest=${1:?frozen manifest}
task_native=${2:?frozen actual executable}
task_run=$(mktemp -d build/task-5z-controller-replay-XXXXXX)
sha256sum --check "$task_manifest" > "$task_run/inputs-before.log"
for task_trial in 1 2; do
  /usr/bin/time -v -o "$task_run/resource$task_trial.log" timeout 3600 "$task_native" > "$task_run/run$task_trial.jsonl"
done
cmp "$task_run/run1.jsonl" "$task_run/run2.jsonl"
/home/lin/anaconda3/bin/python "$task_script_dir/task-5z-controller-check.py" "$task_run/run1.jsonl" | tee "$task_run/independent-check.log"
/home/lin/anaconda3/bin/python "$task_script_dir/task-5z-controller-negative.py" "$task_run/run1.jsonl" | tee "$task_run/negative-check.log"
sha256sum --check "$task_manifest" > "$task_run/inputs-after.log"
sha256sum "$task_run/run1.jsonl" "$task_run/run2.jsonl"
echo "Preserved final run $task_run"
