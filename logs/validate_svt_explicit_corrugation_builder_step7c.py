#!/usr/bin/env python3
"""Validate inputs and formulas used by the explicit-flat support builder."""

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LATTICE = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
THICKNESS_MM = 0.150
HEIGHT_MM = 6.000
THETA_DEG = 45.000
CENTER_DELTA_MM = HEIGHT_MM - THICKNESS_MM
WEB_Y_SPAN_MM = CENTER_DELTA_MM / math.tan(math.radians(THETA_DEG))
FLAT_WIDTH_MM = {"A": 12.15, "B": 10.65}


by_disk = defaultdict(list)
with LATTICE.open(newline="") as stream:
    for row in csv.DictReader(stream):
        by_disk[row["disk"]].append(row)

for disk, rows in by_disk.items():
    rows.sort(key=lambda row: float(row["row_y_mm"]))
    reflected = disk.startswith("ED")
    flat_widths = []
    flat_centers = []
    for row in rows:
        pitch = float(row["pitch_mm"])
        flat_width = float(row["flat_width_mm"])
        assert flat_width == FLAT_WIDTH_MM[row["family"]]
        global_surface_offset = float(row["surface_offset_mm"])
        local_surface_offset = -global_surface_offset if reflected else global_surface_offset
        surface_sign = 1.0 if local_surface_offset > 0.0 else -1.0
        flat_widths.append(flat_width)
        flat_center = local_surface_offset - surface_sign * THICKNESS_MM / 2.0
        flat_centers.append(flat_center)
        recovered_global_surface = (
            -local_surface_offset if reflected else local_surface_offset
        )
        assert recovered_global_surface == global_surface_offset

    angles = []
    for index in range(len(rows) - 1):
        separation = float(rows[index + 1]["row_y_mm"]) - float(rows[index]["row_y_mm"])
        available_gap = separation - (flat_widths[index] + flat_widths[index + 1]) / 2.0
        assert abs(available_gap - WEB_Y_SPAN_MM) < 1e-9
        z_span = flat_centers[index + 1] - flat_centers[index]
        assert abs(abs(z_span) - CENTER_DELTA_MM) < 1e-9
        angles.append(abs(math.degrees(math.atan2(z_span, WEB_Y_SPAN_MM))))

    assert all(abs(angle - THETA_DEG) < 1e-9 for angle in angles)
    assert all(width > 0.0 for width in flat_widths)
    print(
        f"{disk}: flats={len(rows)} webs={len(rows) - 1} "
        f"flat_widths={sorted(set(round(width, 3) for width in flat_widths))} "
        f"web_angle={angles[0]:.3f} deg flat_web_gap=0.000 mm"
    )

print("PASS: enlarged simulation flats meet every 45-degree web without y gaps")
print(
    "NOTE: the sharp-junction approximation omits the 5 mm-radius bends and "
    "uses 10.65/12.15 mm B/A flats"
)
