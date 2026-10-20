#!/usr/bin/env bash
# Step 8j: build/export, then inspect visualization and ROOT overlaps.
# Run in the 26.09.0-stable EIC shell before the separate Geant4 check.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log_dir="$repo_root/proj_svt_disk_tiling/logs"
report="$log_dir/check_svt_upright_modules_visual_root_step8j_${run_stamp}.txt"
root_report="$log_dir/check_svt_upright_modules_root_step8j_${run_stamp}.txt"
display="$repo_root/proj_svt_disk_tiling/out/visualize_upright_modules_step8j_${run_stamp}.root"
mkdir -p "$repo_root/proj_svt_disk_tiling/out"
cd "$repo_root"
exec > >(tee "$report") 2>&1

printf 'SVT Step 8j visual/ROOT check: %s\n' "$run_stamp"
command -v dd_web_display >/dev/null || {
  printf '%s\n' 'ERROR: enter the 26.09.0-stable EIC shell first.' >&2
  exit 1
}
python3 -B proj_svt_disk_tiling/logs/validate_svt_upright_modules_step8j.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_fpc_registration_step8d.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_placement_transform_step5b.py
git diff --check

cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install -DEPIC_ECCE_LEGACY_COMPAT=OFF
cmake --build build -j4
cmake --install build
source "$repo_root/install/bin/thisepic.sh"

xml="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_tiles_new.xml"
dd_web_display -o "$display" --export "$xml"
checkOverlaps --option m --t 0.01 -c "$xml" 2>&1 | tee "$root_report"
if ! grep -Fq 'Number of illegal overlaps/extrusions : 0' "$root_report"; then
  printf '%s\n' 'FAIL: ROOT did not report zero illegal overlaps/extrusions.' >&2
  exit 1
fi
printf 'ROOT PASS. Visually review %s before running the Step 8j Geant4 script.\n' "$display"
