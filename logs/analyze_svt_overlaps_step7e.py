#!/usr/bin/env python3
"""Classify concise overlap records from the Step 7d geometry log."""

import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs/check_svt_rsu_opt_geometry_step7d.txt"
PATTERN = re.compile(
    r"^ = Overlap (ov\d+): (.+?) overlapping (.+?) ovlp=([0-9.eE+-]+)$"
)


def leaf(name):
    return name.rsplit("/", 1)[-1]


def is_module(name):
    return re.match(r"^[AB][56][01]{2}_\d+$", leaf(name)) is not None


def frame_kind(name):
    value = leaf(name)
    if "_frame_flat" in value:
        return "flat"
    if "_frame_web" in value:
        return "web"
    if "_frame_uniform" in value:
        return "fallback"
    return None


def classify(first, second):
    first_frame = frame_kind(first)
    second_frame = frame_kind(second)
    first_module = is_module(first)
    second_module = is_module(second)
    if first_module and second_module:
        return "module-module"
    if first_frame and second_module:
        return f"module-{first_frame}"
    if second_frame and first_module:
        return f"module-{second_frame}"
    if first_frame and second_frame:
        return "support-support"
    if "BeamPipe" in first or "BeamPipe" in second or "beam" in first.lower() or "beam" in second.lower():
        return "beampipe-other"
    return "other"


records = []
with LOG.open(errors="replace") as stream:
    for line in stream:
        match = PATTERN.match(line.rstrip())
        if match:
            overlap_id, first, second, depth = match.groups()
            # scripts/checkOverlaps.py reports overlap depths in centimetres.
            depth_mm = 10.0 * float(depth)
            records.append((overlap_id, first, second, depth_mm, classify(first, second)))

if not records:
    raise SystemExit("no concise overlap records found")

print(f"records={len(records)} first={records[0][0]} last={records[-1][0]}")
by_class = defaultdict(list)
for record in records:
    by_class[record[4]].append(record)

for category, items in sorted(by_class.items(), key=lambda item: (-len(item[1]), item[0])):
    depths = Counter(round(item[3], 9) for item in items)
    layers = Counter(item[1].split("/", 1)[0] for item in items if "/" in item[1])
    maximum = max(items, key=lambda item: item[3])
    print(f"\n{category}: count={len(items)} max={maximum[3]:.9g} mm")
    print("  depths:", ", ".join(f"{depth:g}:{count}" for depth, count in depths.most_common(8)))
    if layers:
        print("  layers:", ", ".join(f"{name}:{count}" for name, count in layers.most_common()))
    for overlap_id, first, second, depth, _ in items[:5]:
        print(f"  {overlap_id} {depth:.9g} mm: {leaf(first)} <> {leaf(second)}")

fallback_count = sum(frame_kind(first) == "fallback" or frame_kind(second) == "fallback"
                     for _, first, second, _, _ in records)
print(f"\nfallback_support_records={fallback_count}")
