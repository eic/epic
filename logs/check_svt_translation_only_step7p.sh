#!/usr/bin/env bash
set -euo pipefail

# Read-only transform audit. It does not compile or modify detector geometry.
python3 logs/audit_svt_translation_only_step7p.py
git diff --check
