#!/usr/bin/env python3
"""Validate the versioned SVT tiling layout and its installed-file contract."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TILING_ROOT = ROOT / "compact/tracking/silicon_disks_tiling"
DELIVERY = TILING_ROOT / "20261001_ES"
OLD_ROOT = ROOT / "compact/tracking/silicon_disks"
MODULES_XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CMAKE = ROOT / "CMakeLists.txt"
OLD_REFERENCE = "compact/tracking/" + "silicon_disks/"
PROVENANCE = (
    "From Ernst S. Module tiling 20261001 -\n"
    "https://conferences.lbl.gov/event/2600/"
)
SCENARIOS = {"rsu_opt": 2164, "all_6rsu": 2160}
DISKS = ("ED0", "ED1", "ED2", "ED3", "ED4", "HD0", "HD1", "HD2", "HD3b", "HD4")


assert not OLD_ROOT.exists(), OLD_ROOT
assert TILING_ROOT.is_dir() and DELIVERY.is_dir()

for readme in (TILING_ROOT / "README.md", DELIVERY / "README.md"):
    text = readme.read_text()
    assert PROVENANCE in text, readme
    assert "logs/" not in text, readme

for scenario, expected_count in SCENARIOS.items():
    package = DELIVERY / scenario
    assert (package / "catalog.csv").is_file()
    module_count = 0
    for disk in DISKS:
        placement = package / f"{disk}_modules.csv"
        metadata = package / f"{disk}_metadata.txt"
        assert placement.is_file(), placement
        assert metadata.is_file(), metadata
        with placement.open(newline="") as stream:
            module_count += sum(1 for _ in csv.DictReader(stream))
    assert module_count == expected_count, (scenario, module_count, expected_count)

corrugation = DELIVERY / "SVT_endcap_corrugation_rows_33_36.csv"
assert corrugation.is_file()

xml = MODULES_XML.read_text()
package_prefix = "compact/tracking/silicon_disks_tiling/20261001_ES/"
corrugation_path = (
    "compact/tracking/silicon_disks_tiling/20261001_ES/"
    "SVT_endcap_corrugation_rows_33_36.csv"
)
active_paths = {
    line.split('file="', 1)[1].split('"', 1)[0]
    for line in xml.splitlines()
    if "<module_placements" in line and 'format="tiling-directory"' in line
}
assert len(active_paths) == 1, active_paths
active_path = active_paths.pop()
assert active_path.startswith(package_prefix), active_path
active_scenario = active_path[len(package_prefix) :]
assert active_scenario in SCENARIOS, active_scenario
assert xml.count(f'file="{active_path}"') == 6
assert xml.count(f'rowwise_placement="{corrugation_path}"') == 10
assert OLD_REFERENCE not in xml

cmake = CMAKE.read_text()
assert 'PATTERN "*_metadata.txt"' in cmake
assert 'PATTERN "*.txt"' not in cmake

for pattern in ("*.py", "*.sh"):
    for path in (ROOT / "proj_svt_disk_tiling/logs").glob(pattern):
        assert OLD_REFERENCE not in path.read_text(), path

print("PASS: versioned 20261001_ES delivery contains both complete scenarios")
print(f"PASS: all 6 detectors atomically select {active_scenario}")
print("PASS: production XML has 10 shared corrugation references")
print("PASS: permanent READMEs contain provenance and no development-log links")
print("PASS: CMake installs only the required *_metadata.txt text files")
