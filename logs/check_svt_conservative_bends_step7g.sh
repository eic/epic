#!/usr/bin/env bash
# Run inside ~/weic/eic-shell --version 26.09.0-stable.
# Check slide flat widths with conservative gaps replacing omitted 5 mm bends.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_conservative_bends_step7g.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 7g conservative-bend geometry check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf 'RESULT: COMPLETE (inspect conservative-bend component overlaps)\n'
