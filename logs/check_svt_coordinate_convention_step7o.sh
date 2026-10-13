#!/usr/bin/env bash
set -euo pipefail

# Read-only convention audit; no geometry build or source-data rewrite.
python3 logs/validate_svt_coordinate_convention_step7o.py
git diff --check
