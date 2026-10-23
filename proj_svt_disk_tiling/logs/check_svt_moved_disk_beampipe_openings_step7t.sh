#!/usr/bin/env bash
set -euo pipefail

# Re-evaluate the fixed EPIC beam-pipe envelope at the two disk positions
# restored in Step 6d. The former branch positions are printed for comparison.
repo_root=$(git rev-parse --show-toplevel)
layout_root=$(dirname "${repo_root}")
helper="${repo_root}/proj_svt_disk_tiling/disk_layout/scripts/beampipe_envelope_circles.py"
output="${repo_root}/proj_svt_disk_tiling/logs/check_svt_moved_disk_beampipe_openings_step7t.txt"

{
  echo "Current Step-6d disk centers"
  for z_mm in -250 250 -450 450; do
    python3 "${helper}" "${z_mm}" --unit mm --repo-root "${layout_root}"
  done

  echo
  echo "Former branch disk centers"
  for z_mm in -255 255 -550 550; do
    python3 "${helper}" "${z_mm}" --unit mm --repo-root "${layout_root}"
  done

  echo
  echo "Current 10 mm layer faces farthest from the IP"
  for z_mm in -255 -455 -655 -855 -1025 255 455 705 955 1205; do
    python3 "${helper}" "${z_mm}" --unit mm --repo-root "${layout_root}"
  done
} > "${output}"

cat "${output}"
