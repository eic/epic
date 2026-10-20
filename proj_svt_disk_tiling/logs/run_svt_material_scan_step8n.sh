#!/usr/bin/env bash
# Scan Step 8l SVT endcap material versus eta at one fixed azimuth.
# Run inside the 26.09.0-stable EIC shell; this script does not rebuild EPIC.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
start_z_cm="${1:-0}"
phi_deg="${2:-20}"
if [[ $# -gt 2 || ! "$start_z_cm" =~ ^-?[0-9]+([.][0-9]+)?$ ||
      ! "$phi_deg" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
  printf 'Usage: bash %s [start_z_cm [phi_deg]]\n' "$0" >&2
  exit 2
fi
if ! command -v python3 >/dev/null || ! command -v root >/dev/null; then
  printf '%s\n' 'ERROR: enter the 26.09.0-stable EIC shell first.' >&2
  exit 1
fi
if [[ ! -f "$repo_root/install/bin/thisepic.sh" ]]; then
  printf '%s\n' 'ERROR: this EPIC checkout must be built and installed first.' >&2
  exit 1
fi
set +u  # Generated setup reads optional environment variables.
source "$repo_root/install/bin/thisepic.sh"
set -u
export DETECTOR_PATH="$repo_root/install/share/epic"

xml="$DETECTOR_PATH/epic_craterlake_tracking_only_disk_tiles_new.xml"
installed_module="$DETECTOR_PATH/compact/tracking/silicon_disks_modules.xml"
source_module="$repo_root/compact/tracking/silicon_disks_modules.xml"
if [[ ! -f "$xml" || ! -f "$installed_module" ]]; then
  printf '%s\n' 'ERROR: installed tracking-only disk-tiles geometry is missing.' >&2
  exit 1
fi
cmp "$source_module" "$installed_module" || {
  printf '%s\n' 'ERROR: installed SVT module XML differs from this checkout; rebuild/install first.' >&2
  exit 1
}
if ! grep -Fq 'SiEndcapTilingRSUCu_thickness' "$installed_module"; then
  printf '%s\n' 'ERROR: installed geometry does not contain the Step 8l RSU copper layer.' >&2
  exit 1
fi

run_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
output_dir="$repo_root/proj_svt_disk_tiling/matscan/step8n_${run_stamp}_phi${phi_deg}_z${start_z_cm}cm"
if [[ -e "$output_dir" ]]; then
  printf 'ERROR: output directory already exists: %s\n' "$output_dir" >&2
  exit 1
fi
mkdir -p "$output_dir"
cd "$output_dir"
exec > >(tee material_scan_stdio.txt) 2>&1

printf 'Step 8n SVT endcap material scan: %s\n' "$run_stamp"
printf 'HEAD=%s\nXML=%s\nStart=(0,0,%s) cm; phi=%s deg\n' \
  "$(git -C "$repo_root" rev-parse HEAD)" "$xml" "$start_z_cm" "$phi_deg"
printf 'Eta=[-4.2,4.2], 841 bins; stop radius=200 cm; epsilon=1e-3; quantity=X0\n'
sha256sum "$xml" "$installed_module" "$repo_root/scripts/subdetector_tests/material_scan.py"

# The script assigns material only to top-level DD4hep detectors. The current
# geometry exposes three SVT endcap assemblies, not six disk detector names.
# This is a single-phi slice, not a phi-averaged material budget.
MPLBACKEND=Agg python3 "$repo_root/scripts/subdetector_tests/material_scan.py" "$xml" \
  --path-r 200 \
  --start-point "0,0,$start_z_cm" \
  --eta-min -4.2 --eta-max 4.2 --eta-nbins 841 \
  --epsilon 1e-3 --phi "$phi_deg" --value-type X0 \
  --detectors InnerSiTrackerSubAssembly,MiddleSiEndcapSubAssembly,OuterSiEndcapSubAssembly

for result in material_scan_details.pdf material_scan_agg.csv material_scan.png; do
  if [[ ! -s "$result" ]]; then
    printf 'ERROR: missing or empty material-scan result: %s\n' "$result" >&2
    exit 1
  fi
done
IFS= read -r csv_header < material_scan_agg.csv
for detector in InnerSiTrackerSubAssembly MiddleSiEndcapSubAssembly OuterSiEndcapSubAssembly; do
  if [[ ",$csv_header," != *",$detector,"* ]]; then
    printf 'ERROR: material scan did not resolve detector %s; do not use the PDF.\n' "$detector" >&2
    exit 1
  fi
done
printf 'Completed scan: %s\n' "$output_dir"
