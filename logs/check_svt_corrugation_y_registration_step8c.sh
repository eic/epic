#!/usr/bin/env bash
# Run inside the selected stable eic-shell environment.
# Validate the corrected facet phase, rebuild/export, and run the deterministic
# ROOT overlap check. Geant4 is intentionally deferred until visual review.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_corrugation_y_registration_step8c.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8c corrugation y-registration check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 logs/validate_svt_corrugation_lattice_step7b.py
python3 logs/audit_svt_epoxy_facet_y_step8b.py
python3 logs/validate_svt_explicit_corrugation_builder_step7c.py
python3 logs/validate_svt_scenario_activation_step7d.py
python3 logs/validate_svt_tiling_delivery_layout_step7x.py
python3 logs/validate_svt_all_6rsu_step8a.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
printf '%s\n' 'EXPECTATION: every edge-glue strip overlaps its neighboring same-surface facet in y.'
printf '%s\n' 'RESULT: COMPLETE (perform an edge-on visual review before running Geant4)'
