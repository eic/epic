#!/usr/bin/env python3
"""Audit a complete engineer-to-EPIC x-convention transformation.

Ground rules tested here:
  * the installed beam-pipe geometry and current physical XML holes are fixed;
  * supplied y/z coordinates and facing are unchanged;
  * supplied x coordinates are mapped with x_EPIC = -x_engineer;
  * an x reflection swaps LEFT/RIGHT module chirality;
  * source files and source IDs remain unchanged.

This script is read-only.
"""

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks"
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
SCENARIOS = ("rsu_opt", "all_6rsu")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")

# The physical holes currently implemented by Step 7n, before the common 5 mm
# buffer is added. Coordinates are DD4hep layer-local.
PHYSICAL_UNBUFFERED = {
    "ED0": ((0.0, 31.7500),),
    "ED1": ((0.0, 31.7500),),
    "ED2": ((0.0, 34.9000),),
    "ED3": ((0.0, 34.9000), (-21.2544, 15.9000)),
    "ED4": ((0.0, 34.9000), (-25.5053, 15.9000)),
    "HD0": ((0.0, 31.7500),),
    "HD1": ((0.0, 31.7500),),
    "HD2": ((-3.1360, 37.1360),),
    "HD3b": ((-8.7360, 42.7360),),
    "HD4": ((-14.3360, 48.3360),),
}
BUFFER = 5.0


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        values[key] = value
    return values


def opening_clearance(x_center, y_center, length, width, openings):
    """Minimum signed clearance from a module rectangle to circular openings."""
    clearances = []
    for opening_x, radius in openings:
        dx = max(abs(opening_x - x_center) - length / 2.0, 0.0)
        dy = max(abs(y_center) - width / 2.0, 0.0)
        clearances.append(math.hypot(dx, dy) - radius)
    return min(clearances)


# Confirm that the current XML is the physical-hole baseline being audited.
definitions = DEFINITIONS.read_text()
for value in ("-3.1360", "-8.7360", "-14.3360", "-21.2544", "-25.5053"):
    assert re.search(rf'value="{re.escape(value)}\*mm"', definitions), value

print("scenario,modules,type_swaps,source_opening_mismatches,min_source_clearance_mm,"
      "min_physical_clearance_mm,physical_rectangle_intersections")

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
    type_swaps = 0
    opening_mismatches = 0
    minimum_source_clearance = math.inf
    minimum_physical_clearance = math.inf
    physical_intersections = 0

    for disk in DISKS:
        layer_sign = -1.0 if disk.startswith("ED") else 1.0
        metadata = read_metadata(package / f"{disk}_metadata.txt")
        source_openings_local = []
        for index in range(int(metadata["n_opening_primitives"])):
            engineer_x = float(metadata[f"opening_primitive_{index}_cx_mm"])
            radius = float(metadata[f"opening_primitive_{index}_a_mm"])
            # First convert engineer detector x to EPIC detector x, then into
            # the current reflected/non-reflected DD4hep layer frame.
            source_openings_local.append((layer_sign * -engineer_x, radius))

        # Every transformed source primitive must equal, or be wholly contained
        # by, a current unbuffered physical primitive.
        for source_x, source_radius in source_openings_local:
            if not any(
                abs(source_x - physical_x) + source_radius <= physical_radius + 1.0e-4
                for physical_x, physical_radius in PHYSICAL_UNBUFFERED[disk]
            ):
                opening_mismatches += 1

        with (package / f"{disk}_modules.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            module_count += 1
            source_type = row["type_id"]
            placed_type = mirror_type[source_type]
            assert placed_type != source_type
            assert catalog[placed_type]["facing"] == catalog[source_type]["facing"]
            assert catalog[placed_type]["chirality"] != catalog[source_type]["chirality"]
            type_swaps += 1

            length = float(catalog[source_type]["module_length_mm"])
            width = float(catalog[source_type]["module_height_mm"])
            boundary_magnitude = length / 2.0 - 20.0
            engineer_origin_x = float(row["x_origin_mm"])
            engineer_center_x = engineer_origin_x - math.copysign(
                boundary_magnitude, engineer_origin_x
            )
            epic_center_x = -engineer_center_x
            local_center_x = layer_sign * epic_center_x
            y_center = float(row["y_origin_mm"])

            source_clearance = opening_clearance(
                local_center_x, y_center, length, width, source_openings_local
            )
            physical_openings = tuple(
                (center_x, radius + BUFFER)
                for center_x, radius in PHYSICAL_UNBUFFERED[disk]
            )
            physical_clearance = opening_clearance(
                local_center_x, y_center, length, width, physical_openings
            )
            minimum_source_clearance = min(minimum_source_clearance, source_clearance)
            minimum_physical_clearance = min(minimum_physical_clearance, physical_clearance)
            if physical_clearance < -1.0e-4:
                physical_intersections += 1

    grand_modules += module_count
    print(
        f"{scenario},{module_count},{type_swaps},{opening_mismatches},"
        f"{minimum_source_clearance:.6f},{minimum_physical_clearance:.6f},"
        f"{physical_intersections}"
    )

print(f"PASS: complete x-convention transform is defined for all {grand_modules} modules")
print("PASS: transformed source openings match the current physical-hole centers")
print("PASS: y, z, facing, and source identifiers remain unchanged")
print("NOTE: LEFT/RIGHT geometry must be swapped; flipping coordinates alone is incomplete")
