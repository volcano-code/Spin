#!/usr/bin/env bash
set -euo pipefail
task_script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_repo=$(git -C "$task_script_dir" rev-parse --show-toplevel)
cd "$task_repo"
task_manifest=${1:?new manifest path}
task_native=${2:?actual native executable}
task_export=${3:?actual preview export}
test ! -e "$task_manifest"
test -x "$task_native"
task_list=$(mktemp build/task-5z-input-list-XXXXXX)
ninja -C build/task-5c-no-faster-20260909 -t deps |
  awk '/^CMakeFiles\/make_it_stand_core.dir\// {active=1; next}
       /^[^ ]/ {active=0}
       active && /^    / {sub(/^    /, ""); print}' |
  sort -u > "$task_list"
test -s "$task_list"
xargs -d '\n' sha256sum < "$task_list" > "$task_manifest"
sha256sum CMakeLists.txt build/task-5c-no-faster-20260909/CMakeCache.txt \
  build/task-5c-no-faster-20260909/build.ninja \
  build/task-5c-no-faster-20260909/libmake_it_stand_core.a \
  build/task-5c-no-faster-20260909/lib/libhighs.a \
  /usr/bin/c++ /home/lin/anaconda3/bin/python \
  /home/lin/anaconda3/lib/libboost_system.a /home/lin/anaconda3/lib/libboost_thread.a \
  /usr/lib/x86_64-linux-gnu/libgmp.so /usr/lib/x86_64-linux-gnu/libmpfr.so \
  /usr/lib/x86_64-linux-gnu/libz.so tests/source_surface_test_utils.h \
  tests/spin_it_2014_source_cage_material_regions_regression.cpp \
  build/task-5c-no-faster-20260909/spin-it-2014-source-cage-material-regions-regression \
  "$task_native" "$task_script_dir/task-5z-controller-export.cpp" \
  "$task_script_dir/task-5z-controller-protocol.md" \
  "$task_script_dir/task-5z-controller-export-helpers.h" \
  "$task_script_dir/task-5y-controller-check.py" \
  "$task_script_dir/task-5w-controller-export-helpers.h" \
  "$task_script_dir/task-5z-controller-check.py" \
  "$task_script_dir/task-5z-controller-negative.py" \
  "$task_script_dir/task-5z-python-deps.py" \
  "$task_script_dir/task-5w-controller-check.py" \
  "$task_script_dir/task-5u-controller-check.py" \
  "$task_script_dir/task-5z-controller-replay.sh" \
  "$task_script_dir/task-5z-controller-build.sh" \
  "$task_script_dir/task-5z-freeze-inputs.sh" >> "$task_manifest"
task_python_list=$(mktemp build/task-5z-python-input-list-XXXXXX)
/home/lin/anaconda3/bin/python "$task_script_dir/task-5z-python-deps.py" "$task_export" > "$task_python_list"
xargs -d '\n' sha256sum < "$task_python_list" >> "$task_manifest"
{
  ldd "$task_native" /usr/bin/c++ /home/lin/anaconda3/bin/python
  while IFS= read -r task_module; do
    case "$task_module" in *.so) ldd "$task_module" ;; esac
  done < "$task_python_list"
} |
  awk '/=> \// {print $3} /^[[:space:]]*\// && $1 !~ /:$/ {print $1}' |
  sort -u | xargs -r -d '\n' sha256sum >> "$task_manifest"
echo "Captured $task_manifest; dependencies $task_list and $task_python_list"
