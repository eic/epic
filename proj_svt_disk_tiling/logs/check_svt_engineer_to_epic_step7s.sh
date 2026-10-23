#!/usr/bin/env bash
set -euo pipefail

# Read-only convention audit. Run from the parent epic repository.
repo_root=$(git rev-parse --show-toplevel)
cd "${repo_root}"

python3 proj_svt_disk_tiling/logs/audit_svt_engineer_to_epic_step7s.py \
  | tee proj_svt_disk_tiling/logs/check_svt_engineer_to_epic_step7s.txt
