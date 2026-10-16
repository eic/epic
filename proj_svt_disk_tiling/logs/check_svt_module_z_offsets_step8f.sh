#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment. This records the
# Step 8f rsu_opt translations through build, export, ROOT, and visual review.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
log_file="$repo_root/proj_svt_disk_tiling/logs/check_svt_module_z_offsets_step8f.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8f module-z-offset check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 proj_svt_disk_tiling/logs/validate_svt_module_z_offsets_step8f.py
python3 proj_svt_disk_tiling/logs/validate_svt_provisional_z_offset_step7j.py
python3 proj_svt_disk_tiling/logs/validate_svt_fpc_registration_step8d.py
python3 proj_svt_disk_tiling/logs/validate_svt_placement_transform_step5b.py
python3 proj_svt_disk_tiling/logs/validate_svt_scenario_activation_step7d.py
python3 proj_svt_disk_tiling/logs/audit_svt_fpc_pair_overlaps_step8d.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
printf '%s\n' 'VISUAL CHECK: inward modules move 30 um away; FPC/flat contact is removed; support gaps remain.'
printf '%s\n' 'RESULT: COMPLETE (review build, ROOT, and visualization before the Geant4 continuation)'
