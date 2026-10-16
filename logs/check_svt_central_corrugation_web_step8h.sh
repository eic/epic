#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment. Validate and build the
# Step 8h central-web correction, then export and run the ROOT overlap check.
# The user explicitly waived a further Geant4 check for this focused change.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_central_corrugation_web_step8h.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8h central-corrugation-web check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
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
printf '%s\n' 'VISUAL CHECK: verify the central web is present on ED0--ED4 and HD0--HD2.'
printf '%s\n' 'VISUAL CHECK: HD3b/HD4 retain their y=0 central flat and both neighboring webs.'
printf '%s\n' 'RESULT: COMPLETE (no further Geant4 check requested for Step 8h)'
