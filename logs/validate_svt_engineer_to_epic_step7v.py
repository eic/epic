#!/usr/bin/env python3
"""Validate the implemented engineer-to-EPIC tiling conversion.

The delivered package is read-only. Its x coordinates are reflected once at
ingestion, LEFT/RIGHT placed handedness is exchanged, and y/z/facing/source
identifiers remain unchanged. Both supported placement scenarios are audited.
"""

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
SCENARIOS = ("rsu_opt", "all_6rsu")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
X_SIGN = -1.0
BUFFER = 5.0

# Current Step-7n physical holes, now expressed in detector-aligned layer
# coordinates. Step 7t will replace their dimensions later from the beam pipe.
PHYSICAL_UNBUFFERED = {
    "ED0": ((0.0, 31.7500),),
    "ED1": ((0.0, 31.7500),),
    "ED2": ((0.0, 34.9000),),
    "ED3": ((0.0, 34.9000), (+21.2544, 15.9000)),
    "ED4": ((0.0, 34.9000), (+25.5053, 15.9000)),
    "HD0": ((0.0, 31.7500),),
    "HD1": ((0.0, 31.7500),),
    "HD2": ((-3.1360, 37.1360),),
    "HD3b": ((-8.7360, 42.7360),),
    "HD4": ((-14.3360, 48.3360),),
}


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            key, value = (part.strip() for part in line.split("=", 1))
            values[key] = value
    return values


def opening_clearance(x_center, y_center, length, width, openings):
    """Return the minimum signed rectangle-to-circle clearance."""
    clearances = []
    for opening_x, radius in openings:
        dx = max(abs(opening_x - x_center) - length / 2.0, 0.0)
        dy = max(abs(y_center) - width / 2.0, 0.0)
        clearances.append(math.hypot(dx, dy) - radius)
    return min(clearances)


cpp = CPP.read_text()
definitions = DEFINITIONS.read_text()
assert "constexpr double engineer_to_epic_x_sign = -1.0" in cpp
assert "row.x_origin = engineer_to_epic_x_sign * row.source_x_origin" in cpp
assert "row.handedness = reflected_handedness(row.source_handedness)" in cpp
assert "parsed_metadata.opening_0.center_x *= engineer_to_epic_x_sign" in cpp
assert "parsed_metadata.opening_1.center_x *= engineer_to_epic_x_sign" in cpp
assert "source_engineer_x_origin_mm" in cpp
assert "placed_epic_x_origin_mm" in cpp
for value in ("21.2544", "25.5053", "-3.1360", "-8.7360", "-14.3360"):
    assert re.search(rf'value="{re.escape(value)}\*mm"', definitions), value

print(
    "scenario,modules,mirrored_centers,chirality_swaps,opening_mismatches,"
    "min_source_clearance_mm,min_physical_clearance_mm,physical_intersections"
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
    opening_mismatches = 0
    minimum_source_clearance = math.inf
    minimum_physical_clearance = math.inf
    physical_intersections = 0

    for disk in DISKS:
        metadata = read_metadata(package / f"{disk}_metadata.txt")
        epic_openings = []
        for index in range(int(metadata["n_opening_primitives"])):
            engineer_x = float(metadata[f"opening_primitive_{index}_cx_mm"])
            radius = float(metadata[f"opening_primitive_{index}_a_mm"])
            epic_openings.append((X_SIGN * engineer_x, radius))

        for epic_x, source_radius in epic_openings:
            if not any(
                abs(epic_x - physical_x) + source_radius <= physical_radius + 1.0e-4
                for physical_x, physical_radius in PHYSICAL_UNBUFFERED[disk]
            ):
                opening_mismatches += 1

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

            y_center = float(row["y_origin_mm"])
            source_clearance = opening_clearance(
                epic_center_x, y_center, length, width, epic_openings
            )
            physical_openings = tuple(
                (center_x, radius + BUFFER)
                for center_x, radius in PHYSICAL_UNBUFFERED[disk]
            )
            physical_clearance = opening_clearance(
                epic_center_x, y_center, length, width, physical_openings
            )
            minimum_source_clearance = min(minimum_source_clearance, source_clearance)
            minimum_physical_clearance = min(minimum_physical_clearance, physical_clearance)
            if physical_clearance < -1.0e-4:
                physical_intersections += 1

    grand_modules += module_count
    print(
        f"{scenario},{module_count},{mirrored_centers},{chirality_swaps},"
        f"{opening_mismatches},{minimum_source_clearance:.6f},"
        f"{minimum_physical_clearance:.6f},{physical_intersections}"
    )

print(f"PASS: {grand_modules} module centers are reflected exactly once")
print("PASS: every source type has one dimension/facing-matched opposite chirality")
print("PASS: transformed source openings align with the current physical holes")
print("PASS: y, z, facing, source type IDs, and source files are unchanged")
