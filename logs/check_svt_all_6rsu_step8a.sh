#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment.
# Validate the all_6rsu selection, rebuild/export, and run ROOT and Geant4
# overlap checks. Visual inspection remains the final review gate.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_all_6rsu_step8a.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8a all_6rsu check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 logs/validate_svt_all_6rsu_step8a.py
python3 logs/validate_svt_scenario_activation_step7d.py
python3 logs/validate_svt_tiling_delivery_layout_step7x.py
python3 logs/validate_svt_engineer_to_epic_step7v.py
python3 logs/validate_svt_physical_openings_step7z.py
python3 logs/validate_svt_xml_radius_authority_step7w.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'EXPECTATION: 2160 six-RSU modules and zero ROOT/Geant4 overlaps.'
printf '%s\n' 'RESULT: COMPLETE (visually inspect both sides and review overlap output before committing)'
