#!/usr/bin/env python3
"""Classify unique Geant4 component overlaps from the Step 7f log."""

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LOG = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "logs/check_svt_component_geometry_step7f.txt"
START = re.compile(r"Overlap is detected for volume (.+?)(?: \([^)]*\))?(?: with (.+?)(?: \([^)]*\))?)?$")
ENCAPSULATED = re.compile(r"apparently fully encapsulating volume (.+?)(?: \([^)]*\))? at the same level!")
DEPTH = re.compile(r" by ([0-9.eE+-]+) (um|mm)")
MODULE = re.compile(r"!([AB][56][01]{2})_\d+#.*!(component[^#:]+)")
DISK = re.compile(r"((?:Inner|Middle|Outer)TrackerEndcap[PN]_disk\d+)")


def frame_kind(name):
    if "_frame_flat" in name:
        return "flat"
    if "_frame_web" in name:
        return "web"
    if "_frame_uniform" in name:
        return "fallback"
    return None


def module_component(name):
    match = MODULE.search(name)
    if not match:
        return None
    module_type, component = match.groups()
    if component.startswith("component0_external_epoxy"):
        kind = "external-epoxy"
    elif component.startswith("component1"):
        kind = "baseplate"
    elif component.startswith("component2"):
        kind = "film"
    elif component.startswith("component3"):
        if "_sensor_" in component:
            kind = "sensor-active"
        elif "_bias_" in component:
            kind = "sensor-bias"
        elif "_backbone_" in component:
            kind = "sensor-backbone"
        elif "_powerswitch_" in component:
            kind = "sensor-powerswitch"
        else:
            kind = "sensor-other"
    elif component.startswith("component4"):
        kind = "bridge-fpc"
    elif component.startswith("component5"):
        kind = "ancasic"
    else:
        kind = "module-other"
    facing = "outward" if module_type[2] == "0" else "inward"
    return module_type, facing, kind


def classify(first, second):
    first_frame, second_frame = frame_kind(first), frame_kind(second)
    first_component, second_component = module_component(first), module_component(second)
    if first_frame and second_component:
        return f"support-{first_frame}/{second_component[1]}-{second_component[2]}"
    if second_frame and first_component:
        return f"support-{second_frame}/{first_component[1]}-{first_component[2]}"
    if first_component and second_component:
        kinds = sorted((f"{first_component[1]}-{first_component[2]}",
                        f"{second_component[1]}-{second_component[2]}"))
        return "module-component/" + "+".join(kinds)
    if first_frame and second_frame:
        return "support-support"
    if "BeamPipe" in first or "BeamPipe" in second:
        return "beampipe-layer"
    return "other"


lines = LOG.read_text(errors="replace").splitlines()
records = []
index = 0
while index < len(lines):
    match = START.search(lines[index].strip())
    if not match:
        index += 1
        continue
    first, second = match.groups()
    depth_mm = None
    cursor = index + 1
    while cursor < min(index + 8, len(lines)):
        encapsulated = ENCAPSULATED.search(lines[cursor].strip())
        if encapsulated:
            second = encapsulated.group(1)
        depth = DEPTH.search(lines[cursor])
        if depth:
            value, unit = depth.groups()
            depth_mm = float(value) / 1000.0 if unit == "um" else float(value)
        if "G4Exception-END" in lines[cursor]:
            break
        cursor += 1
    if second:
        records.append((first, second, depth_mm, classify(first, second)))
    index = cursor + 1

unique = {}
for record in records:
    key = tuple(sorted((record[0], record[1])))
    old = unique.get(key)
    if old is None or (record[2] is not None and (old[2] is None or record[2] > old[2])):
        unique[key] = record

print(f"reports={len(records)} unique_pairs={len(unique)}")
by_class = defaultdict(list)
for record in unique.values():
    by_class[record[3]].append(record)

for category, items in sorted(by_class.items(), key=lambda item: (-len(item[1]), item[0])):
    measured = [item[2] for item in items if item[2] is not None]
    maximum = max(measured) if measured else None
    print(f"{category}: count={len(items)} max_mm={maximum if maximum is not None else 'encapsulation'}")

print("\nmodule types in support overlaps:")
types = Counter()
for first, second, _, category in unique.values():
    if category.startswith("support-"):
        component = module_component(first) or module_component(second)
        types[component[0]] += 1
print(", ".join(f"{name}:{count}" for name, count in sorted(types.items())))

print("\nsupport overlaps by detector disk:")
support_disks = Counter()
for first, second, _, category in unique.values():
    if not category.startswith("support-"):
        continue
    disk = DISK.search(first) or DISK.search(second)
    support_disks[disk.group(1) if disk else "unknown"] += 1
for name, count in sorted(support_disks.items()):
    print(f"{name}: {count}")

fallback_count = sum(frame_kind(first) == "fallback" or frame_kind(second) == "fallback"
                     for first, second, _, _ in unique.values())
print(f"fallback_support_pairs={fallback_count}")
