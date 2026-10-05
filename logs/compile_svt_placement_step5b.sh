#!/usr/bin/env bash
# Run inside ~/weic/eic-shell --version 26.09.0-stable.
# This compiles the strict tiling-directory reader and coordinate transform.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$repo_root/build/rsu_opt_step5b"
log_file="$repo_root/logs/compile_svt_placement_step5b.txt"

exec > >(tee "$log_file") 2>&1
printf 'SVT step 5b compile check: %s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
printf 'Repository: %s\n' "$repo_root"
printf 'Build directory: %s\n' "$build_dir"
cmake -S "$repo_root" -B "$build_dir" \
  -DCMAKE_INSTALL_PREFIX="$repo_root/install"
cmake --build "$build_dir" --target epic --parallel 2
printf 'RESULT: PASS (epic target compiled)\n'
