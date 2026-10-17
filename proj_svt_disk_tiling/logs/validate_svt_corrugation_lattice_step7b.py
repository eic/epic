#!/usr/bin/env python3
"""Validate the generated 33/36 mm corrugation facet-center lattice."""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LATTICE = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
SCENARIOS = ("rsu_opt", "all_6rsu")
PITCH_MM = {"A": 36.0, "B": 33.0}
MODULE_WIDTH_MM = {"A": 35.0, "B": 32.0}
GLUE_WIDTH_MM = 2.5


def interval_overlap(left, right):
    return max(0.0, min(left[1], right[1]) - max(left[0], right[0]))


by_disk = defaultdict(list)
with LATTICE.open(newline="") as stream:
    for row in csv.DictReader(stream):
        by_disk[row["disk"]].append(row)

for disk, rows in by_disk.items():
    for row, following in zip(rows, rows[1:]):
        index = int(row["row_index"])
        next_index = int(following["row_index"])
        assert next_index == index + 1
        expected_step = (PITCH_MM[row["family"]] + PITCH_MM[following["family"]]) / 4.0
        actual_step = float(row["row_y_mm"]) - float(following["row_y_mm"])
        assert abs(actual_step - expected_step) < 1e-9
        assert float(row["surface_offset_mm"]) == -float(following["surface_offset_mm"])

    coordinates = [float(row["row_y_mm"]) for row in rows]
    assert abs(max(coordinates) + min(coordinates)) <= 1.5
    assert {float(row["surface_offset_mm"]) for row in rows} == {-3.0, 3.0}

    transitions = sum(left["family"] != right["family"] for left, right in zip(rows, rows[1:]))
    print(
        f"{disk}: rows={len(rows)} y=[{min(coordinates):+.3f},{max(coordinates):+.3f}] "
        f"transitions={transitions}"
    )

for scenario in SCENARIOS:
    matched_rows = 0
    maximum_shift = 0.0
    maximum_b_shift = 0.0
    glue_strip_count = 0
    minimum_glue_overlap = GLUE_WIDTH_MM
    print(scenario)
    for disk, rows in by_disk.items():
        interior_by_index = {int(row["row_index"]): row for row in rows[1:-1]}
        placement_families = {}
        placement_y = defaultdict(set)
        placement_z = defaultdict(set)
        placement_modules = []
        with (PACKAGE_ROOT / scenario / f"{disk}_modules.csv").open(newline="") as stream:
            for module in csv.DictReader(stream):
                placement_modules.append(module)
                row_index = int(module["row_index"])
                placement_families.setdefault(row_index, module["type_id"][0])
                placement_y[row_index].add(float(module["y_origin_mm"]))
                placement_z[row_index].add(float(module["z_baseplate_mm"]))

        assert set(placement_families) <= set(interior_by_index)
        missing_rows = sorted(set(interior_by_index) - set(placement_families))
        expected_missing = [0, 47] if scenario == "all_6rsu" and disk == "ED4" else []
        assert missing_rows == expected_missing

        all_surface_z = {value for values in placement_z.values() for value in values}
        assert len(all_surface_z) == 2
        disk_center_z = (min(all_surface_z) + max(all_surface_z)) / 2.0
        for row_index, family in placement_families.items():
            lattice_row = interior_by_index[row_index]
            assert family == lattice_row["family"]
            assert len(placement_y[row_index]) == 1
            assert len(placement_z[row_index]) == 1
            source_y = next(iter(placement_y[row_index]))
            shift = abs(float(lattice_row["row_y_mm"]) - source_y)
            maximum_shift = max(maximum_shift, shift)
            if family == "B":
                maximum_b_shift = max(maximum_b_shift, shift)
            source_surface_offset = next(iter(placement_z[row_index])) - disk_center_z
            assert abs(source_surface_offset + float(lattice_row["surface_offset_mm"])) < 1e-9
            matched_rows += 1

        for module in placement_modules:
            family = module["type_id"][0]
            source_y = float(module["y_origin_mm"])
            source_surface_offset = float(module["z_baseplate_mm"]) - disk_center_z
            module_half_width = MODULE_WIDTH_MM[family] / 2.0
            glue_intervals = (
                (source_y - module_half_width, source_y - module_half_width + GLUE_WIDTH_MM),
                (source_y + module_half_width - GLUE_WIDTH_MM, source_y + module_half_width),
            )
            same_surface_facets = [
                row
                for row in rows
                if abs(float(row["surface_offset_mm"]) - source_surface_offset) < 1e-9
            ]
            for glue_interval in glue_intervals:
                overlaps = []
                for facet in same_surface_facets:
                    facet_half_width = float(facet["flat_width_mm"]) / 2.0
                    facet_y = float(facet["row_y_mm"])
                    overlaps.append(
                        interval_overlap(
                            glue_interval,
                            (facet_y - facet_half_width, facet_y + facet_half_width),
                        )
                    )
                best_overlap = max(overlaps, default=0.0)
                assert abs(best_overlap - GLUE_WIDTH_MM) < 1e-9
                glue_strip_count += 1
                minimum_glue_overlap = min(minimum_glue_overlap, best_overlap)

        print(f"  {disk}: matched={len(placement_families)} missing={missing_rows or 'none'}")

    assert maximum_b_shift < 0.05
    print(
        f"  matched_rows={matched_rows} max_module_shift={maximum_shift:.4f} "
        f"max_B_module_shift={maximum_b_shift:.4f} "
        f"glue_strips={glue_strip_count} min_glue_overlap={minimum_glue_overlap:.4f}"
    )

print("PASS: both packages agree with the pitches, joins, opposite-surface facets, and glue contact")
