#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment after the Step 8c
# build/export, ROOT, and visual checks have passed.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_corrugation_y_registration_geant_step8c.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8c corrected-corrugation Geant4 check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'RESULT: COMPLETE (review Geant4 output for overlap reports or exceptions)'
