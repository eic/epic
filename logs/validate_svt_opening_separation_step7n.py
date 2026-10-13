#!/usr/bin/env python3
"""Validate source/physical SVT opening separation and XML source values."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
MODULES = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
PROVENANCE = (
    "From Ernst S. Module tiling 20261001 - "
    "https://conferences.lbl.gov/event/2600/"
)

definitions = DEFINITIONS.read_text()
modules = MODULES.read_text()
cpp = CPP.read_text()

assert PROVENANCE in definitions
assert modules.count(PROVENANCE) >= 7


def constant_expression(name):
    match = re.search(rf'<constant name="{name}"\s+value="([^"]+)"', definitions)
    assert match, name
    return match.group(1).strip()


# Supplied page-7 values, expressed in DD4hep layer-local x. Physical radii
# add the common 5 mm facility clearance without altering source metadata.
expected = {
    "InnerTrackerEndcapP_lepton_radius": "31.750*mm + Beampipe_bakeout_buffer",
    "InnerTrackerEndcapN_lepton_radius": "31.750*mm + Beampipe_bakeout_buffer",
    "TrackerEndcapPDisk2_lepton_center_x": "-3.1360*mm",
    "TrackerEndcapPDisk2_lepton_radius": "37.1360*mm + Beampipe_bakeout_buffer",
    "TrackerEndcapPDisk3_lepton_center_x": "-8.7360*mm",
    "TrackerEndcapPDisk3_lepton_radius": "42.7360*mm + Beampipe_bakeout_buffer",
    "TrackerEndcapPDisk4_lepton_center_x": "-14.3360*mm",
    "TrackerEndcapPDisk4_lepton_radius": "48.3360*mm + Beampipe_bakeout_buffer",
    "TrackerEndcapNDisk2_lepton_radius": "34.9000*mm + Beampipe_bakeout_buffer",
    "TrackerEndcapNDisk3_hadron_center_x": "-21.2544*mm",
    "TrackerEndcapNDisk3_hadron_radius": "15.9000*mm + Beampipe_bakeout_buffer",
    "TrackerEndcapNDisk4_hadron_center_x": "-25.5053*mm",
    "TrackerEndcapNDisk4_hadron_radius": "15.9000*mm + Beampipe_bakeout_buffer",
}
for name, expression in expected.items():
    assert constant_expression(name) == expression, name

assert "bool has_source_opening{false};" in cpp
assert "source_lepton_opening" in cpp
assert "source_hadron_opening" in cpp
assert "disk.has_source_opening = true;" in cpp
assert "disk.lepton_opening = metadata.opening_0" not in cpp
assert "disk.hadron_opening = metadata.opening_1" not in cpp
assert "disk.has_source_opening ? disk.source_lepton_opening : disk.lepton_opening" in cpp
assert "disk.has_source_opening ? disk.source_hadron_opening : disk.hadron_opening" in cpp

print("PASS: supplied page-7 openings are explicit, attributed XML parameters")
print("PASS: physical XML radii add the configurable 5 mm clearance")
print("PASS: metadata openings validate placement without replacing physical holes")
