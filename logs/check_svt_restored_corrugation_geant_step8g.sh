#!/usr/bin/env bash
# Run after check_svt_restored_corrugation_step8g.sh and its visual review. This
# reuses that exact installed rsu_opt geometry for the independent Geant scan.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_file="$repo_root/logs/check_svt_restored_corrugation_geant_step8g.txt"
build_log="$repo_root/logs/check_svt_restored_corrugation_step8g.txt"
source_xml="$repo_root/compact/tracking/silicon_disks_modules.xml"
source_cpp="$repo_root/src/SiEndcapModuleTracker_geo.cpp"
installed_xml="$repo_root/install/share/epic/compact/tracking/silicon_disks_modules.xml"

if [[ ! -f "$build_log" ]]; then
  printf '%s\n' 'ERROR: missing Step 8g build transcript.' >&2
  printf '%s\n' 'Run bash logs/check_svt_restored_corrugation_step8g.sh first.' >&2
  exit 1
fi
if [[ "$source_xml" -nt "$build_log" || "$source_cpp" -nt "$build_log" ]]; then
  printf '%s\n' 'ERROR: Step 8g geometry sources are newer than the recorded build.' >&2
  printf '%s\n' 'Rerun bash logs/check_svt_restored_corrugation_step8g.sh first.' >&2
  exit 1
fi
if [[ ! -f "$installed_xml" ]] || ! cmp -s "$source_xml" "$installed_xml"; then
  printf '%s\n' 'ERROR: installed silicon-disks XML does not match the Step 8g source.' >&2
  printf '%s\n' 'Rerun bash logs/check_svt_restored_corrugation_step8g.sh first.' >&2
  exit 1
fi

exec > >(tee "$log_file") 2>&1
printf 'SVT step 8g Geant4 continuation: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
cd "$repo_root"
source "$repo_root/install/bin/thisepic.sh"

detector_config="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_modules.xml"
printf 'Detector configuration: %s\n' "$detector_config"
python scripts/checkOverlaps.py -c "$detector_config"
printf '%s\n' 'RESULT: COMPLETE (review Geant4 output for support/module overlaps)'
