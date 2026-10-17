#!/usr/bin/env bash
# Check the remote formatting merge before creating the merge commit.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
exec > >(tee proj_svt_disk_tiling/logs/check_svt_remote_merge_20261016.txt) 2>&1

python3 -B - <<'PY'
import re
import subprocess
from pathlib import Path

source_path = "src/SiEndcapModuleTracker_geo.cpp"
before = subprocess.check_output(
    ["git", "show", f"HEAD:{source_path}"], universal_newlines=True
)
after = Path(source_path).read_text()
token_pattern = re.compile(
    r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|'
    r"'(?:\\.|[^'\\])*'|[A-Za-z_]\w*|"
    r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[^\s]",
    re.DOTALL,
)

def code_tokens(source):
    return [
        token.group()
        for token in token_pattern.finditer(source)
        if not token.group().startswith(("//", "/*"))
    ]

local_tokens = code_tokens(before)
merged_tokens = code_tokens(after)
if local_tokens != merged_tokens:
    for position, (local, merged) in enumerate(zip(local_tokens, merged_tokens)):
        if local != merged:
            raise SystemExit(
                f"FAIL: C++ token {position} changed from {local!r} to {merged!r}"
            )
    raise SystemExit(
        f"FAIL: C++ token count changed from {len(local_tokens)} to {len(merged_tokens)}"
    )
print(f"PASS: merged C++ preserves all {len(local_tokens)} local code tokens")
PY

python3 -B proj_svt_disk_tiling/logs/validate_svt_engineer_to_epic_step7v.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_physical_openings_step7z.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_fpc_registration_step8d.py
python3 -B proj_svt_disk_tiling/logs/validate_svt_closed_corrugation_y_gaps_step8i.py
git diff --check
printf '%s\n' 'PASS: remote formatting merge preserved the local SVT implementation'
