#!/usr/bin/env bash
set -euo pipefail

# Read-only convention audit; no geometry build or source-data rewrite.
python3 proj_svt_disk_tiling/logs/validate_svt_coordinate_convention_step7o.py
git diff --check
