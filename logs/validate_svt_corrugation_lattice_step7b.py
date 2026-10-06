#!/usr/bin/env python3
"""Validate the generated 33/36 mm corrugation flat-center lattice."""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LATTICE = ROOT / "compact/tracking/SVT_endcap_corrugation_rows_33_36.csv"
PLACEMENTS = ROOT / "compact/tracking/silicon_disks/rsu_opt"
PITCH_MM = {"A": 36.0, "B": 33.0}


by_disk = defaultdict(list)
with LATTICE.open(newline="") as stream:
    for row in csv.DictReader(stream):
        by_disk[row["disk"]].append(row)

for disk, rows in by_disk.items():
    placement_families = {}
    placement_y = defaultdict(set)
    with (PLACEMENTS / f"{disk}_modules.csv").open(newline="") as stream:
        for module in csv.DictReader(stream):
            row_index = int(module["row_index"])
            placement_families.setdefault(row_index, module["type_id"][0])
            placement_y[row_index].add(float(module["y_origin_mm"]))

    assert len(rows) == len(placement_families)
    for row, following in zip(rows, rows[1:]):
        index = int(row["row_index"])
        next_index = int(following["row_index"])
        assert next_index == index + 1
        assert row["family"] == placement_families[index]
        expected_step = (PITCH_MM[row["family"]] + PITCH_MM[following["family"]]) / 4.0
        actual_step = float(row["row_y_mm"]) - float(following["row_y_mm"])
        assert abs(actual_step - expected_step) < 1e-9
        assert float(row["surface_offset_mm"]) == -float(following["surface_offset_mm"])

    coordinates = [float(row["row_y_mm"]) for row in rows]
    assert abs(max(coordinates) + min(coordinates)) <= 1.5
    assert {float(row["surface_offset_mm"]) for row in rows} == {-3.0, 3.0}

    # B rows expose the construction lattice directly except for the documented
    # projective freedom at a transition.  The only observed B displacement is
    # 0.0391 mm at the first lower-y B row of HD3b and HD4.
    maximum_b_shift = 0.0
    for row in rows:
        index = int(row["row_index"])
        if row["family"] == "B":
            assert len(placement_y[index]) == 1
            source_y = next(iter(placement_y[index]))
            maximum_b_shift = max(maximum_b_shift, abs(float(row["row_y_mm"]) - source_y))
    assert maximum_b_shift < 0.05

    transitions = sum(left["family"] != right["family"] for left, right in zip(rows, rows[1:]))
    print(
        f"{disk}: rows={len(rows)} y=[{min(coordinates):+.3f},{max(coordinates):+.3f}] "
        f"transitions={transitions} max_B_module_shift={maximum_b_shift:.4f}"
    )

print("PASS: pitches, centered joins, A/B assignments, and source B-row anchors agree")
