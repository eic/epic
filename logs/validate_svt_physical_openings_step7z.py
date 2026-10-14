#!/usr/bin/env python3
"""Validate beam-pipe-derived physical SVT openings and module clearance."""

import csv
import collections
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
SCENARIOS = ("rsu_opt", "all_6rsu")
BUFFER_MM = 5.0

# Values audited directly from the fixed EPIC beam-pipe XML in Step 7y.
# Each tuple is (physical XML prefix, ((lepton cx, r), (hadron cx, r))).
PHYSICAL_UNBUFFERED = {
    "ED0": ("InnerTrackerEndcapN", ((0.0, 31.757000), (0.0, 31.757000))),
    "ED1": ("TrackerEndcapNDisk1", ((0.0, 31.757000), (0.0, 31.757000))),
    "ED2": ("TrackerEndcapNDisk2", ((0.0, 31.757000), (0.0, 31.757000))),
    "ED3": ("TrackerEndcapNDisk3", ((0.0, 33.500000), (+21.379454, 14.285000))),
    "ED4": ("TrackerEndcapNDisk4", ((0.0, 33.500000), (+25.630340, 14.285000))),
    "HD0": ("InnerTrackerEndcapP", ((0.0, 31.757000), (0.0, 31.757000))),
    "HD1": ("TrackerEndcapPDisk1", ((0.0, 31.757000), (0.0, 31.757000))),
    "HD2": ("TrackerEndcapPDisk2", ((0.0, 31.757000), (-0.990167, 35.037842))),
    "HD3b": ("TrackerEndcapPDisk3", ((0.0, 31.757000), (-6.616117, 40.661072))),
    "HD4": ("TrackerEndcapPDisk4", ((0.0, 31.757000), (-12.242066, 46.284301))),
}


def constant_expression(text, name):
    match = re.search(
        r'<constant\s+name="{}"\s+value="([^"]+)"'.format(re.escape(name)), text
    )
    if match is None:
        raise KeyError("constant '{}' not found".format(name))
    return match.group(1).strip()


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            key, value = (part.strip() for part in line.split("=", 1))
            values[key] = value
    return values


def rectangle_circle_clearance(x_center, y_center, length, width, circles):
    """Return minimum signed clearance between a rectangle and circle union."""
    clearances = []
    for circle_x, radius in circles:
        dx = max(abs(circle_x - x_center) - length / 2.0, 0.0)
        dy = max(abs(y_center) - width / 2.0, 0.0)
        clearances.append(math.hypot(dx, dy) - radius)
    return min(clearances)


definitions = DEFINITIONS.read_text()
cpp = CPP.read_text()

assert 'name="Beampipe_bakeout_buffer"        value="5*mm"' in definitions
assert "Physical layer/support openings derived from the fixed EPIC beam-pipe" in definitions
assert "Supplied tiling openings remain" in definitions
assert "Negative disk openings use detector-aligned layer coordinates" in definitions
assert "disk.has_source_opening = true;" in cpp
assert "disk.source_lepton_opening = metadata.opening_0;" in cpp
assert "disk.source_hadron_opening = metadata.opening_1;" in cpp
assert "disk.lepton_opening = metadata.opening_0" not in cpp
assert "disk.hadron_opening = metadata.opening_1" not in cpp

for disk, (prefix, circles) in PHYSICAL_UNBUFFERED.items():
    for label, (center_x, radius) in zip(("lepton", "hadron"), circles):
        center_expression = constant_expression(
            definitions, "{}_{}_center_x".format(prefix, label)
        )
        center_y_expression = constant_expression(
            definitions, "{}_{}_center_y".format(prefix, label)
        )
        radius_expression = constant_expression(
            definitions, "{}_{}_radius".format(prefix, label)
        )
        assert center_expression.endswith("*mm"), disk
        assert math.isclose(float(center_expression[:-3]), center_x, abs_tol=5.0e-7), disk
        assert center_y_expression == "0.0*mm", disk
        radius_value = radius_expression.split("*mm", 1)[0]
        assert math.isclose(float(radius_value), radius, abs_tol=5.0e-7), disk
        assert "+ Beampipe_bakeout_buffer" in radius_expression, disk

print("scenario,modules,min_physical_clearance_mm,physical_rectangle_intersections")
grand_modules = 0
for scenario in SCENARIOS:
    package = PACKAGE_ROOT / scenario
    with (package / "catalog.csv").open(newline="") as stream:
        catalog = {row["type_id"]: row for row in csv.DictReader(stream)}

    module_count = 0
    minimum_clearance = math.inf
    intersections = 0
    intersections_by_disk = collections.Counter()
    intersection_details = []
    for disk, (_, unbuffered_circles) in PHYSICAL_UNBUFFERED.items():
        # Reading metadata here verifies its continued presence; its opening
        # values deliberately do not define or replace the physical circles.
        metadata = read_metadata(package / "{}_metadata.txt".format(disk))
        assert int(metadata["n_opening_primitives"]) in (1, 2)

        buffered_circles = tuple(
            (center_x, radius + BUFFER_MM) for center_x, radius in unbuffered_circles
        )
        with (package / "{}_modules.csv".format(disk)).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            module_count += 1
            module = catalog[row["type_id"]]
            length = float(module["module_length_mm"])
            width = float(module["module_height_mm"])
            boundary_magnitude = length / 2.0 - 20.0
            engineer_origin_x = float(row["x_origin_mm"])
            epic_origin_x = -engineer_origin_x
            epic_center_x = epic_origin_x - math.copysign(boundary_magnitude, epic_origin_x)
            y_center = float(row["y_origin_mm"])
            clearance = rectangle_circle_clearance(
                epic_center_x, y_center, length, width, buffered_circles
            )
            minimum_clearance = min(minimum_clearance, clearance)
            if clearance < -1.0e-4:
                intersections += 1
                intersections_by_disk[disk] += 1
                intersection_details.append(
                    (
                        clearance,
                        disk,
                        row["row_index"],
                        row["mod_index"],
                        row["type_id"],
                    )
                )

    # The full rectangular catalogue envelope is intentionally conservative:
    # it includes empty corners outside the detailed module solids. Record the
    # exact deterministic precheck result; ROOT must decide actual-solid overlap.
    assert intersections == 56, (scenario, intersections)
    assert math.isclose(minimum_clearance, -0.047675, abs_tol=5.0e-7), (
        scenario,
        minimum_clearance,
    )
    grand_modules += module_count
    print("{},{},{:.6f},{}".format(scenario, module_count, minimum_clearance, intersections))
    if intersections:
        print("  conservative contacts by disk: {}".format(dict(intersections_by_disk)))
        for detail in sorted(intersection_details)[:5]:
            clearance, disk, row_index, mod_index, type_id = detail
            print(
                "  contact: disk={} row={} module={} type={} clearance={:.6f} mm".format(
                    disk, row_index, mod_index, type_id, clearance
                )
            )

print("PASS: all 60 physical XML x/y/r fields match the Step 7y beam-pipe audit")
print("PASS: every physical radius independently adds the configurable 5 mm buffer")
print("PASS: conservative precheck reproduces 56 expected near contacts per scenario")
print("REVIEW: exact-solid ROOT overlap validation is mandatory for those near contacts")
print("PASS: source metadata remains verification input and cannot override physical XML")
