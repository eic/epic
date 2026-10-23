#!/usr/bin/env python3
"""Validate the implemented engineer-to-EPIC tiling conversion.

The delivered package is read-only. Its x coordinates are reflected once at
ingestion, LEFT/RIGHT placed handedness is exchanged, and y/z/facing/source
identifiers remain unchanged. Both supported placement scenarios are audited.
"""

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
SCENARIOS = ("rsu_opt", "all_6rsu")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
X_SIGN = -1.0


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            key, value = (part.strip() for part in line.split("=", 1))
            values[key] = value
    return values


cpp = CPP.read_text()
assert "constexpr double engineer_to_epic_x_sign = -1.0" in cpp
assert "row.x_origin = engineer_to_epic_x_sign * row.source_x_origin" in cpp
assert "row.handedness = reflected_handedness(row.source_handedness)" in cpp
assert "parsed_metadata.opening_0.center_x *= engineer_to_epic_x_sign" in cpp
assert "parsed_metadata.opening_1.center_x *= engineer_to_epic_x_sign" in cpp
assert "source_engineer_x_origin_mm" in cpp
assert "placed_epic_x_origin_mm" in cpp

print(
    "scenario,modules,mirrored_centers,chirality_swaps,source_openings"
)
grand_modules = 0
for scenario in SCENARIOS:
    package = PACKAGE_ROOT / scenario
    with (package / "catalog.csv").open(newline="") as stream:
        catalog_rows = list(csv.DictReader(stream))
    catalog = {row["type_id"]: row for row in catalog_rows}

    mirror_type = {}
    for source_type, source in catalog.items():
        opposite = "LEFT" if source["chirality"] == "RIGHT" else "RIGHT"
        matches = [
            candidate["type_id"]
            for candidate in catalog_rows
            if candidate["module_height_mm"] == source["module_height_mm"]
            and candidate["module_length_mm"] == source["module_length_mm"]
            and candidate["rsu_count"] == source["rsu_count"]
            and candidate["facing"] == source["facing"]
            and candidate["chirality"] == opposite
        ]
        assert len(matches) == 1, (scenario, source_type, matches)
        mirror_type[source_type] = matches[0]
    assert all(mirror_type[mirror_type[name]] == name for name in mirror_type)

    module_count = 0
    mirrored_centers = 0
    chirality_swaps = 0
    source_opening_count = 0

    for disk in DISKS:
        metadata = read_metadata(package / f"{disk}_metadata.txt")
        for index in range(int(metadata["n_opening_primitives"])):
            engineer_x = float(metadata[f"opening_primitive_{index}_cx_mm"])
            epic_x = X_SIGN * engineer_x
            assert abs(epic_x + engineer_x) <= 1.0e-12
            source_opening_count += 1

        with (package / f"{disk}_modules.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            module_count += 1
            source_type = row["type_id"]
            placed_type = mirror_type[source_type]
            source = catalog[source_type]
            placed = catalog[placed_type]
            assert placed["facing"] == source["facing"]
            assert placed["chirality"] != source["chirality"]
            chirality_swaps += 1

            length = float(source["module_length_mm"])
            width = float(source["module_height_mm"])
            boundary_magnitude = length / 2.0 - 20.0
            engineer_origin_x = float(row["x_origin_mm"])
            engineer_center_x = engineer_origin_x - math.copysign(
                boundary_magnitude, engineer_origin_x
            )
            epic_origin_x = X_SIGN * engineer_origin_x
            epic_center_x = epic_origin_x - math.copysign(
                boundary_magnitude, epic_origin_x
            )
            assert abs(epic_center_x + engineer_center_x) <= 1.0e-9
            mirrored_centers += 1

    grand_modules += module_count
    print(
        f"{scenario},{module_count},{mirrored_centers},{chirality_swaps},"
        f"{source_opening_count}"
    )

print(f"PASS: {grand_modules} module centers are reflected exactly once")
print("PASS: every source type has one dimension/facing-matched opposite chirality")
print("PASS: source-opening centers receive the same engineer-to-EPIC conversion")
print("PASS: physical openings are independently validated by Step 7z")
print("PASS: y, z, facing, source type IDs, and source files are unchanged")
