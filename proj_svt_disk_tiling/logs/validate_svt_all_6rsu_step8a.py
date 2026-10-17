#!/usr/bin/env python3
"""Validate that production atomically selects the complete all_6rsu package."""

import csv
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
XML_PATH = ROOT / "compact/tracking/silicon_disks_modules.xml"
PACKAGE = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/all_6rsu"
EXPECTED_PATH = "compact/tracking/silicon_disks_tiling/20261001_ES/all_6rsu"
EXPECTED_DISKS = (
    "ED0",
    "ED1",
    "ED2",
    "ED3",
    "ED4",
    "HD0",
    "HD1",
    "HD2",
    "HD3b",
    "HD4",
)


root = ET.parse(XML_PATH).getroot()
placements = [
    item
    for item in root.findall(".//module_placements")
    if item.get("format") == "tiling-directory"
]
assert len(placements) == 6
assert {item.get("file") for item in placements} == {EXPECTED_PATH}

with (PACKAGE / "catalog.csv").open(newline="") as stream:
    catalog = {row["type_id"]: row for row in csv.DictReader(stream)}

module_count = 0
used_types = set()
disk_counts = {}
for disk in EXPECTED_DISKS:
    with (PACKAGE / f"{disk}_modules.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert (PACKAGE / f"{disk}_metadata.txt").is_file()
    disk_counts[disk] = len(rows)
    module_count += len(rows)
    used_types.update(row["type_id"] for row in rows)

assert module_count == 2160, module_count
assert used_types
assert all(int(catalog[type_id]["rsu_count"]) == 6 for type_id in used_types)

print("disk,module_count")
for disk in EXPECTED_DISKS:
    print(f"{disk},{disk_counts[disk]}")
print("PASS: all six detector blocks select the all_6rsu package")
print("PASS: all ten placement files are present and contain 2160 modules total")
print(f"PASS: all {len(used_types)} used catalogue types contain exactly six RSUs")
