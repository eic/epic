#!/usr/bin/env python3
"""Audit external-epoxy overlap with corrugation facets in global y.

The delivered module y coordinate is the channel center.  Each module's two
2.5 mm edge glue lines must overlap the neighboring same-surface corrugation
facets, not a facet centered on the module row itself.
"""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DELIVERY = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
LATTICE = DELIVERY / "SVT_endcap_corrugation_rows_33_36.csv"
SCENARIOS = ("rsu_opt", "all_6rsu")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
MODULE_WIDTH_MM = {"A": 35.0, "B": 32.0}
GLUE_WIDTH_MM = 2.5


def interval_overlap(left, right):
    return max(0.0, min(left[1], right[1]) - max(left[0], right[0]))


def load_current_facets():
    facets = defaultdict(list)
    with LATTICE.open(newline="") as stream:
        for row in csv.DictReader(stream):
            facets[row["disk"]].append(
                {
                    "y": float(row["row_y_mm"]),
                    "surface": float(row["surface_offset_mm"]),
                    "family": row["family"],
                    "width": float(row["flat_width_mm"]),
                }
            )
    for rows in facets.values():
        rows.sort(key=lambda row: row["y"], reverse=True)
    return facets


def legacy_phase(corrected_rows):
    """Reconstruct the rejected table: no boundaries and same row surface."""
    return [dict(row, surface=-row["surface"]) for row in corrected_rows[1:-1]]


def glue_intervals(module_y, family):
    half_width = MODULE_WIDTH_MM[family] / 2.0
    return (
        (module_y - half_width, module_y - half_width + GLUE_WIDTH_MM),
        (module_y + half_width - GLUE_WIDTH_MM, module_y + half_width),
    )


def best_facet_overlap(glue, surface, facets):
    matches = []
    for facet in facets:
        if abs(facet["surface"] - surface) > 1.0e-9:
            continue
        half_width = facet["width"] / 2.0
        overlap = interval_overlap(glue, (facet["y"] - half_width, facet["y"] + half_width))
        if overlap > 0.0:
            matches.append(overlap)
    return max(matches, default=0.0)


current_by_disk = load_current_facets()
legacy_by_disk = {disk: legacy_phase(current_by_disk[disk]) for disk in DISKS}

for scenario in SCENARIOS:
    strip_count = 0
    legacy_contacts = 0
    production_contacts = 0
    production_minimum = GLUE_WIDTH_MM

    for disk in DISKS:
        placement_path = DELIVERY / scenario / f"{disk}_modules.csv"
        with placement_path.open(newline="") as stream:
            modules = list(csv.DictReader(stream))

        z_values = [float(module["z_baseplate_mm"]) for module in modules]
        disk_center_z = (min(z_values) + max(z_values)) / 2.0
        for module in modules:
            family = module["type_id"][0]
            module_y = float(module["y_origin_mm"])
            surface = float(module["z_baseplate_mm"]) - disk_center_z
            for glue in glue_intervals(module_y, family):
                strip_count += 1
                legacy_overlap = best_facet_overlap(glue, surface, legacy_by_disk[disk])
                production_overlap = best_facet_overlap(glue, surface, current_by_disk[disk])
                legacy_contacts += legacy_overlap > 1.0e-9
                production_contacts += production_overlap > 1.0e-9
                production_minimum = min(production_minimum, production_overlap)

    print(
        f"{scenario}: strips={strip_count} "
        f"legacy_y_contacts={legacy_contacts} "
        f"production_y_contacts={production_contacts} "
        f"production_min_overlap={production_minimum:.4f} mm"
    )
    assert legacy_contacts == 0
    assert production_contacts == strip_count
    assert production_minimum > 2.4

print("PASS: the rejected lattice misses every strip; production neighboring facets cover all strips")
