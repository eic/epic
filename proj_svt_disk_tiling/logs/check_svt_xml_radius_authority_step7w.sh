#!/usr/bin/env bash
# Run inside ~/weic/eic-shell --version 26.09.0-stable.
# Verify XML-owned outer radii, rebuild/export, and repeat both overlap checks.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
log_file="$repo_root/proj_svt_disk_tiling/logs/check_svt_xml_radius_authority_step7w.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 7w XML radius-authority check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 proj_svt_disk_tiling/logs/validate_svt_xml_radius_authority_step7w.py
python3 proj_svt_disk_tiling/logs/audit_svt_radial_boundaries_step6c.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'EXPECTATION: zero overlaps; constructed radii are numerically unchanged.'
printf '%s\n' 'RESULT: COMPLETE (review overlap output before committing)'
