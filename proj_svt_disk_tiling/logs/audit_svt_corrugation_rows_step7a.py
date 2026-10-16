#!/usr/bin/env python3
"""Measure corrugation information encoded by the delivered module rows."""

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = ROOT / "compact" / "tracking" / "silicon_disks"
SCENARIOS = ("all_6rsu", "rsu_opt")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
METADATA_Z = {
    "ED4": -1020.0,
    "ED3": -850.0,
    "ED2": -650.0,
    "ED1": -450.0,
    "ED0": -250.0,
    "HD0": 250.0,
    "HD1": 450.0,
    "HD2": 700.0,
    "HD3b": 950.0,
    "HD4": 1200.0,
}
RECOMMENDED_PITCH_MM = {
    disk: (33.0 if disk in ("ED0", "ED1", "HD0", "HD1") else 36.0)
    for disk in DISKS
}


def compress_indices(indices):
    ranges = []
    for value in sorted(indices):
        if not ranges or value != ranges[-1][1] + 1:
            ranges.append([value, value])
        else:
            ranges[-1][1] = value
    return [str(first) if first == last else "{}-{}".format(first, last) for first, last in ranges]


for scenario in SCENARIOS:
    directory = PACKAGE_ROOT / scenario
    print(scenario)
    for disk in DISKS:
        rows_by_index = defaultdict(list)
        with (directory / (disk + "_modules.csv")).open(newline="") as stream:
            for module in csv.DictReader(stream):
                rows_by_index[int(module["row_index"])].append(module)

        row_data = []
        type_groups = set()
        family_rows = defaultdict(list)
        for row_index, modules in sorted(rows_by_index.items()):
            y_values = {float(module["y_origin_mm"]) for module in modules}
            corrugation_z = {float(module["z_baseplate_mm"]) for module in modules}
            assert len(y_values) == 1
            assert len(corrugation_z) == 1
            type_groups.update(module["type_id"][:2] for module in modules)
            families = {module["type_id"][0] for module in modules}
            assert len(families) == 1, "one support row mixes A/B baseplate widths"
            family_rows[families.pop()].append(row_index)
            row_data.append((row_index, y_values.pop(), corrugation_z.pop()))

        y_steps = {
            round(abs(right[1] - left[1]), 6)
            for left, right in zip(row_data, row_data[1:])
        }
        same_surface_steps = {
            round(abs(right[1] - left[1]), 6)
            for left, right in zip(row_data, row_data[2:])
        }
        z_offsets = {round(row[2] - METADATA_Z[disk], 6) for row in row_data}
        inferred_pitch = 2.0 * next(iter(y_steps)) if len(y_steps) == 1 else None
        matches_recommendation = (
            inferred_pitch is not None and
            abs(inferred_pitch - RECOMMENDED_PITCH_MM[disk]) < 1e-9
        )
        print(
            "  {:4s}: rows={:2d} y=[{:+7.2f},{:+7.2f}] adjacent={} "
            "same-surface={} inferred pitch={} recommended={:.1f} match={} "
            "z offsets={} types={}".format(
                disk,
                len(row_data),
                min(row[1] for row in row_data),
                max(row[1] for row in row_data),
                sorted(y_steps),
                sorted(same_surface_steps),
                inferred_pitch,
                RECOMMENDED_PITCH_MM[disk],
                matches_recommendation,
                sorted(z_offsets),
                sorted(type_groups),
            )
        )
        print("        baseplate-family rows: {}".format(
            ", ".join("{}={}".format(family, compress_indices(indices))
                      for family, indices in sorted(family_rows.items()))
        ))

print("NOTE: module rows encode facet centers, not bend radii or transition geometry")
