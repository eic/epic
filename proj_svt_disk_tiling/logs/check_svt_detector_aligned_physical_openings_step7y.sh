#!/usr/bin/env bash
set -euo pipefail

# Read-only comparison of the current XML physical openings with the fixed
# EPIC beam-pipe envelope in the detector-aligned SVT layer frame.
repo_root=$(git rev-parse --show-toplevel)
python3 "${repo_root}/proj_svt_disk_tiling/logs/audit_svt_detector_aligned_physical_openings_step7y.py" \
  | tee "${repo_root}/proj_svt_disk_tiling/logs/check_svt_detector_aligned_physical_openings_step7y.txt"
