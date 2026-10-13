#!/usr/bin/env bash
# Run inside ~/weic/eic-shell --version 26.09.0-stable.
# Validate the versioned tiling layout, installed metadata, and rsu_opt export.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_tiling_delivery_layout_step7x.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 7x versioned-delivery check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
python3 logs/validate_svt_tiling_delivery_layout_step7x.py
python3 logs/validate_svt_xml_radius_authority_step7w.py
python3 logs/validate_svt_engineer_to_epic_step7v.py
git diff --check
./setup_local.sh 1 0
source "$repo_root/install/bin/thisepic.sh"

installed_delivery="$DETECTOR_PATH/compact/tracking/silicon_disks_tiling/20261001_ES"
for scenario in rsu_opt all_6rsu; do
    test -f "$installed_delivery/$scenario/catalog.csv"
    for disk in ED0 ED1 ED2 ED3 ED4 HD0 HD1 HD2 HD3b HD4; do
        test -f "$installed_delivery/$scenario/${disk}_modules.csv"
        test -f "$installed_delivery/$scenario/${disk}_metadata.txt"
    done
done
test -f "$installed_delivery/SVT_endcap_corrugation_rows_33_36.csv"

# Run away from the source root so successful loading proves that the installed
# package, including metadata, is complete.
mkdir -p "$repo_root/out"
cd "$repo_root/out"
detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
dd_web_display --export "$detector_config"
checkOverlaps --option m --t 0.01 -c "$detector_config"
printf '%s\n' 'EXPECTATION: installed package loads outside the source root; zero ROOT overlaps.'
printf '%s\n' 'RESULT: COMPLETE (review overlap output before committing)'
