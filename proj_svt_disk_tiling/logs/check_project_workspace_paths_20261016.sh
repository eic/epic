#!/usr/bin/env bash
# Verify relocated project scripts without rebuilding or running overlap checks.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
exec > >(tee proj_svt_disk_tiling/logs/check_project_workspace_paths_20261016.txt) 2>&1

test -f proj_svt_disk_tiling/ES_material/rsu_opt/catalog.csv
test -f proj_svt_disk_tiling/disk_layout/scripts/beampipe_envelope_circles.py

for script in proj_svt_disk_tiling/logs/*.sh proj_svt_disk_tiling/tmp/*.sh; do
  bash -n "$script"
done

python3 -B - <<'PY'
import ast
from pathlib import Path

for directory in (Path("proj_svt_disk_tiling/logs"), Path("proj_svt_disk_tiling/tmp")):
    for script in directory.glob("*.py"):
        ast.parse(script.read_text(), filename=str(script))
PY

python3 -B proj_svt_disk_tiling/logs/validate_svt_closed_corrugation_y_gaps_step8i.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_tiling_delivery_layout_step7x.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_scenario_activation_step7d.py
git diff --check
git diff --cached --check
printf '%s\n' 'PASS: relocated project scripts and representative validators'
