#!/usr/bin/env python3
"""Reconstruct projected epoxy/baseplate overlaps from the rsu_opt source files."""

import csv
import itertools
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "compact/tracking/silicon_disks/rsu_opt"
GLUE_WIDTH_MM = 2.5


with (PACKAGE / "catalog.csv").open(newline="") as stream:
    catalog = {row["type_id"]: row for row in csv.DictReader(stream)}


def overlap(first, second):
    return min(first[1], second[1]) - max(first[0], second[0])


def interval(center, width):
    return center - width / 2.0, center + width / 2.0


def module_record(disk, row):
    spec = catalog[row["type_id"]]
    length = float(spec["module_length_mm"])
    width = float(spec["module_height_mm"])
    x_origin = float(row["x_origin_mm"])
    inward_sign = 1.0 if x_origin > 0.0 else -1.0
    x_center = x_origin - inward_sign * (length / 2.0 - float(spec["BLEC_mm"]))
    y_center = float(row["y_origin_mm"])
    return {
        "disk": disk,
        "row": int(row["row_index"]),
        "module": int(row["mod_index"]),
        "type": row["type_id"],
        "facing": spec["facing"],
        "surface": float(row["z_baseplate_mm"]),
        "x": interval(x_center, length),
        "baseplate_y": interval(y_center, width),
        "epoxy_low_y": (y_center - width / 2.0, y_center - width / 2.0 + GLUE_WIDTH_MM),
        "epoxy_high_y": (y_center + width / 2.0 - GLUE_WIDTH_MM, y_center + width / 2.0),
    }


modules = []
for path in sorted(PACKAGE.glob("*_modules.csv")):
    disk = path.name[: -len("_modules.csv")]
    with path.open(newline="") as stream:
        modules.extend(module_record(disk, row) for row in csv.DictReader(stream))

epoxy_pairs = []
baseplate_pairs = []
for first, second in itertools.combinations(modules, 2):
    if first["disk"] != second["disk"] or first["surface"] != second["surface"]:
        continue
    x_overlap = overlap(first["x"], second["x"])
    if x_overlap <= 0.0:
        continue
    baseplate_y_overlap = overlap(first["baseplate_y"], second["baseplate_y"])
    if baseplate_y_overlap > 0.0:
        baseplate_pairs.append((first, second, x_overlap, baseplate_y_overlap))
    for edge in ("low", "high"):
        y_overlap = overlap(first[f"epoxy_{edge}_y"], second[f"epoxy_{edge}_y"])
        if y_overlap > 0.0:
            epoxy_pairs.append((first, second, edge, x_overlap, y_overlap))


def pair_label(first, second):
    return "+".join(sorted((first["facing"], second["facing"])))


print(f"modules={len(modules)}")
print(f"projected_epoxy_pairs={len(epoxy_pairs)}")
print("epoxy facing:", dict(sorted(Counter(pair_label(a, b) for a, b, *_ in epoxy_pairs).items())))
print("epoxy same row:", sum(a["row"] == b["row"] for a, b, *_ in epoxy_pairs))
print("epoxy by disk:", dict(sorted(Counter(a["disk"] for a, *_ in epoxy_pairs).items())))
print(
    "epoxy overlap ranges [mm]:",
    f"x={min(item[3] for item in epoxy_pairs):.4f}..{max(item[3] for item in epoxy_pairs):.4f}",
    f"y={min(item[4] for item in epoxy_pairs):.4f}..{max(item[4] for item in epoxy_pairs):.4f}",
)

print(f"projected_baseplate_pairs={len(baseplate_pairs)}")
print(
    "baseplate facing:",
    dict(sorted(Counter(pair_label(a, b) for a, b, *_ in baseplate_pairs).items())),
)
print("baseplate same row:", sum(a["row"] == b["row"] for a, b, *_ in baseplate_pairs))
print("baseplate by disk:", dict(sorted(Counter(a["disk"] for a, *_ in baseplate_pairs).items())))
print(
    "baseplate overlap ranges [mm]:",
    f"x={min(item[2] for item in baseplate_pairs):.4f}..{max(item[2] for item in baseplate_pairs):.4f}",
    f"y={min(item[3] for item in baseplate_pairs):.4f}..{max(item[3] for item in baseplate_pairs):.4f}",
)

print("\nfirst ten opposite-facing epoxy pairs:")
shown = 0
for first, second, edge, x_overlap, y_overlap in epoxy_pairs:
    if first["facing"] == second["facing"]:
        continue
    print(
        f"{first['disk']} row/module {first['row']}/{first['module']} {first['type']} + "
        f"{second['row']}/{second['module']} {second['type']} edge={edge} "
        f"overlap={x_overlap:.4f}x{y_overlap:.4f} mm"
    )
    shown += 1
    if shown == 10:
        break
