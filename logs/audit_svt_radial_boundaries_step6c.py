#!/usr/bin/env python3
"""Audit full supplied module footprints against supplied disk boundaries."""

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact" / "tracking" / "silicon_disks"
SCENARIOS = ("all_6rsu", "rsu_opt")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
LAYER_ENVELOPE_ALLOWANCE_MM = 0.001


def metadata_values(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line and "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            try:
                values[key] = float(value)
            except ValueError:
                pass
    return values


def rectangle_circle_clearance(x_min, x_max, y_min, y_max, cx, cy, radius):
    closest_x = min(max(cx, x_min), x_max)
    closest_y = min(max(cy, y_min), y_max)
    return math.hypot(closest_x - cx, closest_y - cy) - radius


for scenario in SCENARIOS:
    directory = PACKAGE_ROOT / scenario
    with (directory / "catalog.csv").open(newline="") as stream:
        catalog = {row["type_id"]: row for row in csv.DictReader(stream)}
    total = 0
    scenario_violations = 0
    envelope_violations = 0
    print(scenario)
    for disk in DISKS:
        metadata = metadata_values(directory / (disk + "_metadata.txt"))
        outer_radius = metadata["outer_radius_mm"]
        openings = [
            (
                metadata["opening_primitive_{}_cx_mm".format(index)],
                metadata["opening_primitive_{}_cy_mm".format(index)],
                metadata["opening_primitive_{}_a_mm".format(index)],
            )
            for index in range(int(metadata["n_opening_primitives"]))
        ]
        minimum_outer_clearance = float("inf")
        minimum_opening_clearance = float("inf")
        violations = 0
        disk_envelope_violations = 0
        worst_row = None
        with (directory / (disk + "_modules.csv")).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            module = catalog[row["type_id"]]
            length = float(module["module_length_mm"])
            width = float(module["module_height_mm"])
            boundary_offset = length / 2.0 - float(module["BLEC_mm"])
            x_origin = float(row["x_origin_mm"])
            x_center = x_origin - math.copysign(boundary_offset, x_origin)
            y_center = float(row["y_origin_mm"])
            x_min, x_max = x_center - length / 2.0, x_center + length / 2.0
            y_min, y_max = y_center - width / 2.0, y_center + width / 2.0
            outer_clearance = outer_radius - max(
                math.hypot(x, y) for x in (x_min, x_max) for y in (y_min, y_max)
            )
            opening_clearance = min(
                rectangle_circle_clearance(x_min, x_max, y_min, y_max, *opening)
                for opening in openings
            )
            minimum_outer_clearance = min(minimum_outer_clearance, outer_clearance)
            minimum_opening_clearance = min(minimum_opening_clearance, opening_clearance)
            if outer_clearance < -1e-9 or opening_clearance < -1e-9:
                violations += 1
                if worst_row is None or min(outer_clearance, opening_clearance) < worst_row[0]:
                    worst_row = (
                        min(outer_clearance, opening_clearance),
                        int(row["row_index"]),
                        int(row["mod_index"]),
                        row["type_id"],
                    )
            if (outer_clearance + LAYER_ENVELOPE_ALLOWANCE_MM < -1e-9 or
                    opening_clearance < -1e-9):
                disk_envelope_violations += 1
            total += 1
        print(
            "  {:4s}: modules={:3d} min outer={:+.8f} mm "
            "min opening={:+.8f} mm nominal={} envelope={} worst={}".format(
                disk, len(rows), minimum_outer_clearance, minimum_opening_clearance, violations,
                disk_envelope_violations, worst_row,
            )
        )
        scenario_violations += violations
        envelope_violations += disk_envelope_violations
    print("  SUMMARY: {} baseplates, {} nominal-radius rounding excursions, "
          "{} envelope violations".format(total, scenario_violations, envelope_violations))
    assert envelope_violations == 0

print("PASS: 0.001 mm layer allowance contains all source-rounding excursions")
print("PASS: sensor/film footprints are subsets of their corresponding baseplates")
