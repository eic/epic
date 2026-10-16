#!/usr/bin/env bash
# Run inside the same selected stable eic-shell after the Step 8d build, ROOT,
# and visual checks. This reuses the unchanged installed geometry.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_fpc_registration_geant_step8d.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8d bridge-FPC Geant4 overlap check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'RESULT: COMPLETE (review Geant4 output for overlap reports or exceptions)'
