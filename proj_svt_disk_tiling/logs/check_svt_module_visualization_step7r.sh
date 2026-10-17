#!/usr/bin/env bash
# Run inside ~/weic/eic-shell --version 26.09.0-stable.
# Rebuild and export the tracking-only geometry for visual inspection.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
log_file="$repo_root/proj_svt_disk_tiling/logs/check_svt_module_visualization_step7r.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 7r module-visualization check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 proj_svt_disk_tiling/logs/validate_svt_module_visualization_step7r.py
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
printf 'RESULT: EXPORTED (visually confirm silicon/module daughters fill keep-clear gaps)\n'
