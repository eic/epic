#!/usr/bin/env bash
# Step 8j: Geant4 overlap check after the visual/ROOT script builds and installs.
# Deferred until the separate A-module FPC alignment fix is ready. Then run in
# the 26.09.0-stable EIC shell after reviewing the rebuilt geometry.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
report="$repo_root/proj_svt_disk_tiling/logs/check_svt_upright_modules_geant_step8j_${run_stamp}.txt"
geant_report="$repo_root/proj_svt_disk_tiling/logs/check_svt_upright_modules_geant_raw_step8j_${run_stamp}.txt"
cd "$repo_root"
exec > >(tee "$report") 2>&1

printf 'SVT Step 8j Geant4 overlap check: %s\n' "$run_stamp"
command -v dd_web_display >/dev/null || {
  printf '%s\n' 'ERROR: enter the 26.09.0-stable EIC shell first.' >&2
  exit 1
}
test -f "$repo_root/install/bin/thisepic.sh"
source "$repo_root/install/bin/thisepic.sh"
xml="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_tiles_new.xml"
test -f "$xml"
cmp "$repo_root/compact/tracking/silicon_disks_modules.xml" \
    "$DETECTOR_PATH/compact/tracking/silicon_disks_modules.xml"
printf 'Compact XML: %s\n' "$xml"
python scripts/checkOverlaps.py -c "$xml" 2>&1 | tee "$geant_report"

if ! grep -Fq 'Geometry overlaps check completed !' "$geant_report"; then
  printf '%s\n' 'FAIL: Geant4 did not report a completed overlap check.' >&2
  exit 1
fi
if grep -Eq 'G4Exception|Overlap is detected|^[[:space:]]*Overlap ov[0-9]+' "$geant_report"; then
  printf '%s\n' 'FAIL: Geant4 overlap warning; review the report.' >&2
  exit 1
fi
printf '%s\n' 'PASS: Geant4 overlap check completed without reported overlaps.'
