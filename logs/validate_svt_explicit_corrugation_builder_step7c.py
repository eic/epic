#!/usr/bin/env python3
"""Validate inputs and formulas used by the explicit-flat support builder."""

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LATTICE = ROOT / "compact/tracking/silicon_disks/SVT_endcap_corrugation_rows_33_36.csv"
THICKNESS_MM = 0.150
HEIGHT_MM = 6.000
THETA_DEG = 45.000
CENTER_DELTA_MM = HEIGHT_MM - THICKNESS_MM
WEB_Y_SPAN_MM = CENTER_DELTA_MM / math.tan(math.radians(THETA_DEG))
FLAT_WIDTH_MM = {"A": 7.86, "B": 6.36}
KEEP_CLEAR_WIDTH_MM = {"A": 21.94, "B": 20.56}
KEEP_CLEAR_DEPTH_MM = 1.00
OMITTED_BEND_GAP_MM = 2.145


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
        assert abs(available_gap - 10.14) < 1e-9
        assert abs((available_gap - WEB_Y_SPAN_MM) / 2.0 - OMITTED_BEND_GAP_MM) < 1e-9
        z_span = flat_centers[index + 1] - flat_centers[index]
        assert abs(abs(z_span) - CENTER_DELTA_MM) < 1e-9
        angles.append(abs(math.degrees(math.atan2(z_span, WEB_Y_SPAN_MM))))

    assert all(abs(angle - THETA_DEG) < 1e-9 for angle in angles)
    assert all(width > 0.0 for width in flat_widths)
    print(
        f"{disk}: flats={len(rows)} webs={len(rows) - 1} "
        f"flat_widths={sorted(set(round(width, 3) for width in flat_widths))} "
        f"web_angle={angles[0]:.3f} deg bend_gap={OMITTED_BEND_GAP_MM:.3f} mm/side"
    )

print("PASS: slide flats and 45-degree webs preserve conservative omitted-bend gaps")
print(
    "NOTE: keep-clear envelopes "
    f"A={KEEP_CLEAR_WIDTH_MM['A']:.2f}x{KEEP_CLEAR_DEPTH_MM:.2f} mm and "
    f"B={KEEP_CLEAR_WIDTH_MM['B']:.2f}x{KEEP_CLEAR_DEPTH_MM:.2f} mm are recorded, "
    "and intersecting straight-piece x sections are omitted over inward modules"
)
