#!/usr/bin/env bash
# Copy the validated handoff packages into the compact geometry tree.
# Run after validate_svt_handoff_step2.py; then inspect the step-3 log below.
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_root="$project_root/ES_material"
target_root="$project_root/compact/tracking/silicon_disks"
log_path="$project_root/logs/stage_svt_handoff_step3.txt"

if [[ -e "$target_root/all_6rsu" || -e "$target_root/rsu_opt" ]]; then
  echo "Destination package already exists; inspect it before retrying." >&2
  exit 1
fi
if [[ ! -f "$source_root/rsu6/catalog.csv" || ! -f "$source_root/rsu_opt/catalog.csv" ]]; then
  echo "A source catalogue is missing." >&2
  exit 1
fi

mkdir -p "$target_root"
cp -a -- "$source_root/rsu6" "$target_root/all_6rsu"
cp -a -- "$source_root/rsu_opt" "$target_root/rsu_opt"

{
  printf 'SVT handoff staging, %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  for scenario in all_6rsu rsu_opt; do
    if [[ "$scenario" == all_6rsu ]]; then source_name=rsu6; else source_name=rsu_opt; fi
    source_dir="$source_root/$source_name"
    target_dir="$target_root/$scenario"
    source_count="$(find "$source_dir" -maxdepth 1 -type f | wc -l)"
    target_count="$(find "$target_dir" -maxdepth 1 -type f | wc -l)"
    [[ "$source_count" -eq 21 && "$target_count" -eq 21 ]]
    for source_file in "$source_dir"/*; do
      target_file="$target_dir/$(basename "$source_file")"
      cmp -s -- "$source_file" "$target_file"
    done
    printf '%s -> %s: %s files, byte-identical\n' "$source_dir" "$target_dir" "$target_count"
    sha256sum "$target_dir"/catalog.csv
  done
  printf 'RESULT PASS\n'
} > "$log_path"

printf 'PASS: staged both packages; details in %s\n' "$log_path"
