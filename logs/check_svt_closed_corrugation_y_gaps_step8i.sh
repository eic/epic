#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment. Validate and build the
# Step 8i no-gap sharp-junction support, then export and run the ROOT check.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_closed_corrugation_y_gaps_step8i.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8i closed corrugation y-gap check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 logs/validate_svt_closed_corrugation_y_gaps_step8i.py
python3 logs/validate_svt_central_corrugation_web_step8h.py
python3 logs/validate_svt_restored_corrugation_step8g.py
python3 logs/validate_svt_module_z_offsets_step8f.py
python3 logs/validate_svt_corrugation_lattice_step7b.py
python3 logs/validate_svt_explicit_corrugation_builder_step7c.py
python3 logs/audit_svt_epoxy_facet_y_step8b.py
python3 logs/validate_svt_placement_transform_step5b.py
python3 logs/validate_svt_scenario_activation_step7d.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
printf '%s\n' 'VISUAL CHECK: every flat meets both neighboring webs with no y gap.'
printf '%s\n' 'VISUAL CHECK: module registration and the eight restored central webs remain correct.'
printf '%s\n' 'RESULT: COMPLETE (review ROOT and visualization before the required Geant4 continuation)'
