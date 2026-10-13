#!/usr/bin/env python3
"""Show where unclipped straight webs enter the slide keep-clear envelopes."""

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LATTICE = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
THICKNESS = 0.150
WEB_SPAN = 6.000 - THICKNESS
KEEP_WIDTH = {"A": 21.94, "B": 20.56}
KEEP_DEPTH = 1.00


def projections(polygon, axis):
    values = [x * axis[0] + z * axis[1] for x, z in polygon]
    return min(values), max(values)


def polygons_overlap(first, second):
    """Strict convex-polygon overlap test; touching faces are not overlaps."""
    for polygon in (first, second):
        for index, point in enumerate(polygon):
            next_point = polygon[(index + 1) % len(polygon)]
            edge = (next_point[0] - point[0], next_point[1] - point[1])
            axis = (-edge[1], edge[0])
            first_min, first_max = projections(first, axis)
            second_min, second_max = projections(second, axis)
            if first_max <= second_min + 1e-12 or second_max <= first_min + 1e-12:
                return False
    return True


by_disk = defaultdict(list)
with LATTICE.open(newline="") as stream:
    for row in csv.DictReader(stream):
        by_disk[row["disk"]].append(row)

total_web_violations = 0
for disk, rows in by_disk.items():
    rows.sort(key=lambda row: float(row["row_y_mm"]))
    reflected = disk.startswith("ED")
    flats = []
    for row in rows:
        global_surface = float(row["surface_offset_mm"])
        surface = -global_surface if reflected else global_surface
        sign = 1.0 if surface > 0.0 else -1.0
        center = surface - sign * THICKNESS / 2.0
        flats.append(
            {
                "y": float(row["row_y_mm"]),
                "surface": surface,
                "center": center,
                "width": float(row["flat_width_mm"]),
                "family": row["family"],
            }
        )

    webs = []
    for left, right in zip(flats, flats[1:]):
        left_edge = left["y"] + left["width"] / 2.0
        right_edge = right["y"] - right["width"] / 2.0
        center_y = (left_edge + right_edge) / 2.0
        center_z = (left["center"] + right["center"]) / 2.0
        z_span = right["center"] - left["center"]
        length = math.hypot(WEB_SPAN, z_span)
        tangent = (WEB_SPAN / length, z_span / length)
        normal = (-tangent[1], tangent[0])
        polygon = []
        for along, across in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            polygon.append(
                (
                    center_y + along * tangent[0] * length / 2.0
                    + across * normal[0] * THICKNESS / 2.0,
                    center_z + along * tangent[1] * length / 2.0
                    + across * normal[1] * THICKNESS / 2.0,
                )
            )
        webs.append(polygon)

    violations = 0
    for flat in flats:
        half_width = KEEP_WIDTH[flat["family"]] / 2.0
        surface = flat["surface"]
        channel_z = surface - math.copysign(KEEP_DEPTH, surface)
        z_low, z_high = sorted((surface, channel_z))
        envelope = [
            (flat["y"] - half_width, z_low),
            (flat["y"] + half_width, z_low),
            (flat["y"] + half_width, z_high),
            (flat["y"] - half_width, z_high),
        ]
        violations += sum(polygons_overlap(envelope, web) for web in webs)

    total_web_violations += violations
    print(f"{disk}: flats={len(flats)} web/keep-clear intersections={violations}")

assert total_web_violations > 0
print(
    f"EXPECTED: {total_web_violations} unclipped straight-web/keep-clear intersections; "
    "the C++ builder removes those support x sections over inward modules"
)
