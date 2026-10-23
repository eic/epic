#!/usr/bin/env bash
# Run after check_svt_closed_corrugation_y_gaps_step8i.sh and visual review.
# Reuse that exact installed rsu_opt geometry for the required Geant4 scan.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
log_file="$repo_root/proj_svt_disk_tiling/logs/check_svt_closed_corrugation_y_gaps_geant_step8i.txt"
build_log="$repo_root/proj_svt_disk_tiling/logs/check_svt_closed_corrugation_y_gaps_step8i.txt"
source_xml="$repo_root/compact/tracking/silicon_disks_modules.xml"
installed_xml="$repo_root/install/share/epic/compact/tracking/silicon_disks_modules.xml"
relative_lattice="compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
source_lattice="$repo_root/$relative_lattice"
installed_lattice="$repo_root/install/share/epic/$relative_lattice"

if [[ ! -f "$build_log" ]]; then
  printf '%s\n' 'ERROR: missing Step 8i build transcript.' >&2
  printf '%s\n' 'Run bash proj_svt_disk_tiling/logs/check_svt_closed_corrugation_y_gaps_step8i.sh first.' >&2
  exit 1
fi
for geometry_source in \
  "$source_xml" \
  "$source_lattice" \
  "$repo_root/src/SiEndcapModuleTracker_geo.cpp" \
  "$repo_root/proj_svt_disk_tiling/logs/generate_svt_corrugation_lattice_step7b.py"; do
  if [[ "$geometry_source" -nt "$build_log" ]]; then
    printf 'ERROR: geometry source is newer than the Step 8i build: %s\n' "$geometry_source" >&2
    printf '%s\n' 'Rerun bash proj_svt_disk_tiling/logs/check_svt_closed_corrugation_y_gaps_step8i.sh first.' >&2
    exit 1
  fi
done
if [[ ! -f "$installed_xml" ]] || ! cmp -s "$source_xml" "$installed_xml"; then
  printf '%s\n' 'ERROR: installed silicon-disks XML does not match source.' >&2
  printf '%s\n' 'Rerun bash proj_svt_disk_tiling/logs/check_svt_closed_corrugation_y_gaps_step8i.sh first.' >&2
  exit 1
fi
if [[ ! -f "$installed_lattice" ]] || ! cmp -s "$source_lattice" "$installed_lattice"; then
  printf '%s\n' 'ERROR: installed corrugation table does not match the Step 8i source.' >&2
  printf '%s\n' 'Rerun bash proj_svt_disk_tiling/logs/check_svt_closed_corrugation_y_gaps_step8i.sh first.' >&2
  exit 1
fi

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8i Geant4 continuation: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'RESULT: COMPLETE (review Geant4 output for support/module and support/support overlaps)'
