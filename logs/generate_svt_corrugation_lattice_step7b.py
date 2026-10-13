#!/usr/bin/env python3
"""Generate the documented 33/36 mm corrugation flat-center lattice.

The transition delivery labels every module row with an A or B baseplate
family.  A maps to a 36 mm same-surface corrugation pitch and B maps to a
33 mm pitch.  Consecutive rows alternate surfaces, so their within-family
separations are 18 and 16.5 mm.  At an A/B boundary, center the join by using
the mean separation, 17.25 mm.
"""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLACEMENTS = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/rsu_opt"
OUTPUT = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
PITCH_MM = {"A": 36.0, "B": 33.0}
FLAT_WIDTH_MM = {"A": 7.86, "B": 6.36}


def load_rows(disk):
    rows = defaultdict(list)
    with (PLACEMENTS / f"{disk}_modules.csv").open(newline="") as stream:
        for module in csv.DictReader(stream):
            rows[int(module["row_index"])].append(module)

    ordered = []
    for row_index, modules in sorted(rows.items()):
        families = {module["type_id"][0] for module in modules}
        z_values = {float(module["z_baseplate_mm"]) for module in modules}
        if len(families) != 1 or len(z_values) != 1:
            raise ValueError(f"{disk} row {row_index} mixes families or corrugation surfaces")
        ordered.append((row_index, families.pop(), z_values.pop()))
    return ordered


def centered_y_coordinates(rows):
    """Return decreasing y coordinates with the middle row/pair centered at zero."""
    count = len(rows)
    coordinates = [None] * count
    if count % 2:
        coordinates[count // 2] = 0.0
        left = right = count // 2
    else:
        left = count // 2 - 1
        right = count // 2
        center_step = (PITCH_MM[rows[left][1]] + PITCH_MM[rows[right][1]]) / 4.0
        coordinates[left] = center_step / 2.0
        coordinates[right] = -center_step / 2.0

    for position in range(left - 1, -1, -1):
        step = (PITCH_MM[rows[position][1]] + PITCH_MM[rows[position + 1][1]]) / 4.0
        coordinates[position] = coordinates[position + 1] + step
    for position in range(right + 1, count):
        step = (PITCH_MM[rows[position - 1][1]] + PITCH_MM[rows[position][1]]) / 4.0
        coordinates[position] = coordinates[position - 1] - step
    return coordinates


with OUTPUT.open("w", newline="") as stream:
    fieldnames = (
        "disk",
        "row_index",
        "row_y_mm",
        "surface_offset_mm",
        "family",
        "pitch_mm",
        "flat_width_mm",
    )
    writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    for disk in DISKS:
        rows = load_rows(disk)
        coordinates = centered_y_coordinates(rows)
        surface_values = {row[2] for row in rows}
        if len(surface_values) != 2:
            raise ValueError(f"{disk} does not contain exactly two corrugation surfaces")
        disk_center_z = (min(surface_values) + max(surface_values)) / 2.0
        for (row_index, family, surface_z), row_y in zip(rows, coordinates):
            writer.writerow(
                {
                    "disk": disk,
                    "row_index": row_index,
                    "row_y_mm": f"{row_y:.3f}",
                    "surface_offset_mm": f"{surface_z - disk_center_z:.3f}",
                    "family": family,
                    "pitch_mm": f"{PITCH_MM[family]:.3f}",
                    "flat_width_mm": f"{FLAT_WIDTH_MM[family]:.3f}",
                }
            )

print(f"wrote {OUTPUT.relative_to(ROOT)}")
