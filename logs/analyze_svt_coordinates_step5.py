#!/usr/bin/env python3
"""Audit the supplied x-origin convention without changing placements."""

import csv
from collections import defaultdict
from decimal import Decimal
from math import hypot
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "compact/tracking/silicon_disks"


def read_catalog(path):
    with path.open() as stream:
        return {row["type_id"]: row for row in csv.DictReader(stream)}


def check_scenario(scenario, reverse_x=False):
    directory = ROOT / scenario
    catalog = read_catalog(directory / "catalog.csv")
    groups = defaultdict(list)
    count = 0
    outer_violations = 0
    examples = []
    zero_origins = 0
    for path in sorted(directory.glob("*_modules.csv")):
        metadata = {}
        with (directory / path.name.replace("_modules.csv", "_metadata.txt")).open() as stream:
            for line in stream:
                if "=" in line:
                    key, value = line.split("=", 1)
                    metadata[key.strip()] = value.split("#", 1)[0].strip()
        outer_radius = float(metadata["outer_radius_mm"])
        with path.open() as stream:
            for row in csv.DictReader(stream):
                design = catalog[row["type_id"]]
                length = Decimal(design["module_length_mm"])
                blec = Decimal(design["BLEC_mm"])
                x_origin = Decimal(row["x_origin_mm"])
                if x_origin == 0:
                    zero_origins += 1
                # Supplied origins lie at the outer end of each module's
                # RSU chain; the package extends toward x=0.
                x_direction = 1 if x_origin > 0 else -1
                if reverse_x:
                    x_direction *= -1
                x_center = x_origin - x_direction * (length / 2 - blec)
                x_min, x_max = x_center - length / 2, x_center + length / 2
                y_center = Decimal(row["y_origin_mm"])
                half_width = Decimal(design["module_height_mm"]) / 2
                if max(hypot(float(x), float(y)) for x in (x_min, x_max)
                       for y in (y_center - half_width, y_center + half_width)) > outer_radius + 0.0001:
                    outer_violations += 1
                    if len(examples) < 5:
                        examples.append((path.name, row["row_index"], row["type_id"],
                                         row["x_origin_mm"], str(x_center)))
                # Different sensor z planes may intentionally overlap in x.
                key = (row["disk_id"], row["row_index"], row["z_sensor_mm"])
                groups[key].append((x_min, x_max))
                count += 1

    overlaps = 0
    for intervals in groups.values():
        intervals.sort()
        for left, right in zip(intervals, intervals[1:]):
            if left[1] > right[0] + Decimal("0.0001"):
                overlaps += 1
    return count, overlaps, outer_violations, zero_origins, examples


def main():
    for scenario in ("all_6rsu", "rsu_opt"):
        count, overlaps, outer_violations, zero_origins, examples = check_scenario(scenario)
        _, reversed_overlaps, reversed_outer_violations, _, _ = check_scenario(scenario, reverse_x=True)
        print(f"{scenario}: {count} modules; {overlaps} same-row/same-sensor-z x overlaps "
              f"with radial-inward x, {reversed_overlaps} with reversed x; "
              f"outer-radius violations {outer_violations}/{reversed_outer_violations}; "
              f"zero x origins {zero_origins}")
        if examples:
            print("first violating rows:", examples)


if __name__ == "__main__":
    main()
