#!/usr/bin/env python3
"""Validate the versioned SVT tiling layout and its installed-file contract."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
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
package_path = "compact/tracking/silicon_disks_tiling/20261001_ES/rsu_opt"
corrugation_path = (
    "compact/tracking/silicon_disks_tiling/20261001_ES/"
    "SVT_endcap_corrugation_rows_33_36.csv"
)
assert xml.count(f'file="{package_path}"') == 6
assert xml.count(f'rowwise_placement="{corrugation_path}"') == 10
assert OLD_REFERENCE not in xml

cmake = CMAKE.read_text()
assert 'PATTERN "*_metadata.txt"' in cmake
assert 'PATTERN "*.txt"' not in cmake

for pattern in ("*.py", "*.sh"):
    for path in (ROOT / "logs").glob(pattern):
        assert OLD_REFERENCE not in path.read_text(), path

print("PASS: versioned 20261001_ES delivery contains both complete scenarios")
print("PASS: production XML has 6 versioned package and 10 corrugation references")
print("PASS: permanent READMEs contain provenance and no development-log links")
print("PASS: CMake installs only the required *_metadata.txt text files")
