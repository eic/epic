#!/usr/bin/env bash
# Continuation of Step 8a when the long Geant4 scan outlives an interactive
# orchestration session. Run inside the selected stable eic-shell environment.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_all_6rsu_geant_step8a.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8a all_6rsu Geant4 continuation: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'RESULT: COMPLETE (review Geant4 output for overlap reports or exceptions)'
