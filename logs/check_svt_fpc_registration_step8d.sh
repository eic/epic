#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment. This checks the
# corrected bridge-FPC/baseplate, epoxy-edge, and LEC/REC registration before
# any further corrugation-overlap change.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_fpc_registration_step8d.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8d bridge-FPC registration check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 logs/validate_svt_fpc_registration_step8d.py
python3 logs/validate_svt_module_params_step4a.py
python3 logs/validate_svt_module_builder_step4b.py
python3 logs/validate_svt_placement_transform_step5b.py
python3 logs/validate_svt_scenario_activation_step7d.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
printf '%s\n' 'VISUAL CHECK: FPCs touch the baseplate, avoid edge epoxy, and follow LEC/REC by 0.5 mm.'
printf '%s\n' 'RESULT: COMPLETE (review build, ROOT, and visualization before the Geant4 continuation)'
