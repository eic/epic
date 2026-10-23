#!/usr/bin/env python3
"""Validate atomic module-placement scenario/support activation in compact XML."""

import csv
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
XML_PATH = ROOT / "compact/tracking/silicon_disks_modules.xml"
DELIVERY = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
CORRUGATION = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
EXPECTED_DISKS = {"ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4"}
EXPECTED_LAYER_MAP = {
    "OuterTrackerEndcapN_disk4": "ED4",
    "OuterTrackerEndcapN_disk3": "ED3",
    "OuterTrackerEndcapN_disk2": "ED2",
    "MiddleTrackerEndcapN_disk1": "ED1",
    "InnerTrackerEndcapN_disk1": "ED0",
    "InnerTrackerEndcapP_disk1": "HD0",
    "MiddleTrackerEndcapP_disk1": "HD1",
    "OuterTrackerEndcapP_disk2": "HD2",
    "OuterTrackerEndcapP_disk3": "HD3b",
    "OuterTrackerEndcapP_disk4": "HD4",
}


root = ET.parse(XML_PATH).getroot()
placements = root.findall(".//module_placements")
active_placements = [item for item in placements if item.get("format") == "tiling-directory"]
assert len(active_placements) == 6
active_paths = {item.get("file") for item in active_placements}
assert len(active_paths) == 1, active_paths
active_path = active_paths.pop()
active_scenario = Path(active_path).name
assert active_scenario in ("rsu_opt", "all_6rsu"), active_scenario
assert active_path == (
    "compact/tracking/silicon_disks_tiling/20261001_ES/" + active_scenario
)
package = DELIVERY / active_scenario

frames = root.findall(".//frame[@type='corrugated']")
assert len(frames) == 10
assert {frame.get("rowwise_placement") for frame in frames} == {
    "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
}

constants = {item.get("name"): item.get("value") for item in root.findall(".//constant")}
assert constants["SiEndcapCorrugation_theta"] == "45*deg"
assert package.is_dir()
assert CORRUGATION.is_file()

with CORRUGATION.open(newline="") as stream:
    corrugation_disks = {row["disk"] for row in csv.DictReader(stream)}
assert corrugation_disks == EXPECTED_DISKS
source_text = (ROOT / "src/SiEndcapModuleTracker_geo.cpp").read_text()
for layer_key, source_key in EXPECTED_LAYER_MAP.items():
    assert f'{{"{layer_key}", "{source_key}"}}' in source_text
for disk in EXPECTED_DISKS:
    assert (package / f"{disk}_modules.csv").is_file()
    assert (package / f"{disk}_metadata.txt").is_file()
assert (package / "catalog.csv").is_file()

print(
    f"PASS: six detectors select {active_scenario} and ten layers select "
    "the shared 45-degree support"
)
