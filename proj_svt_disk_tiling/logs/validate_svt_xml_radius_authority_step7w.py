#!/usr/bin/env python3
"""Verify that XML owns SVT outer radii and metadata only checks them."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
MODULES = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
SCENARIOS = ("rsu_opt", "all_6rsu")
EXPECTED_RADII_MM = {
    "ED0": 230.0,
    "HD0": 230.0,
    "ED1": 405.0,
    "HD1": 405.0,
    "ED2": 430.0,
    "ED3": 430.0,
    "ED4": 430.0,
    "HD2": 430.0,
    "HD3b": 430.0,
    "HD4": 430.0,
}


def constant_expression(text, name):
    match = re.search(rf'<constant name="{name}"\s+value="([^"]+)"', text)
    assert match, name
    return match.group(1).strip()


def metadata_outer_radius(path):
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line.startswith("outer_radius_mm"):
            return float(line.split("=", 1)[1])
    raise AssertionError(f"{path}: missing outer_radius_mm")


definitions = DEFINITIONS.read_text()
modules = MODULES.read_text()
cpp = CPP.read_text()

source_constants = {
    "SiEndcapTilingInner_outer_radius": "230*mm",
    "SiEndcapTilingMiddle_outer_radius": "405*mm",
    "SiEndcapTilingOuter_outer_radius": "430*mm",
}
for name, expected in source_constants.items():
    assert constant_expression(definitions, name) == expected, name

definition_aliases = {
    "InnerTrackerEndcapP_rmax": "SiEndcapTilingInner_outer_radius",
    "InnerTrackerEndcapN_rmax": "InnerTrackerEndcapP_rmax",
    "TrackerEndcapPDisk1_rmax": "SiEndcapTilingMiddle_outer_radius",
    "TrackerEndcapNDisk1_rmax": "TrackerEndcapPDisk1_rmax",
}
for index in (2, 3, 4):
    definition_aliases[f"TrackerEndcapPDisk{index}_rmax"] = (
        "SiEndcapTilingOuter_outer_radius"
    )
    definition_aliases[f"TrackerEndcapNDisk{index}_rmax"] = (
        f"TrackerEndcapPDisk{index}_rmax"
    )
for name, expected in definition_aliases.items():
    assert constant_expression(definitions, name) == expected, name

for side in ("P", "N"):
    inner_layer = f"InnerTrackerEndcap{side}Layer1_rmax"
    inner_module = f"InnerTrackerEndcap{side}Mod1_rmax"
    assert constant_expression(modules, inner_module) == (
        f"InnerTrackerEndcap{side}_rmax"
    )
    assert constant_expression(modules, inner_layer) == f"{inner_module} + 1*um"
    for index in (1, 2, 3, 4):
        layer = f"TrackerEndcap{side}Layer{index}_rmax"
        module = f"TrackerEndcap{side}Mod{index}_rmax"
        assert constant_expression(modules, module) == (
            f"TrackerEndcap{side}Disk{index}_rmax"
        )
        assert constant_expression(modules, layer) == f"{module} + 1*um"

assert "disk.rmax = metadata.outer_radius" not in cpp
assert "expected_layer_outer_radius = metadata.outer_radius + 1.0 * um" in cpp
assert "disk.rmax - expected_layer_outer_radius" in cpp

for scenario in SCENARIOS:
    for disk, expected in EXPECTED_RADII_MM.items():
        actual = metadata_outer_radius(
            PACKAGE_ROOT / scenario / f"{disk}_metadata.txt"
        )
        assert actual == expected, (scenario, disk, actual, expected)

print("PASS: XML explicitly owns the 230/405/430 mm module-boundary radii")
print("PASS: every layer adds the explicit 0.001 mm containment allowance")
print("PASS: both metadata packages agree and C++ verifies without overriding XML")
