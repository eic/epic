#!/usr/bin/env python3
"""Validate the current provisional inward/outward module translations."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"

text = XML.read_text()
clearance_match = re.search(
    r'name="SiEndcapTilingInwardModule_z_clearance"\s+value="([0-9.]+)\*mm"', text
)
relative_match = re.search(
    r'name="SiEndcapTilingOppositeFacing_z_offset"\s+value="([0-9.]+)\*mm"', text
)
assert clearance_match and relative_match
clearance = float(clearance_match.group(1))
relative_offset = float(relative_match.group(1))
outward_offset = clearance + relative_offset

epoxy = 0.080
baseplate = 0.150
assert abs(clearance - 0.030) < 1e-12
assert abs(relative_offset - (epoxy + baseplate)) < 1e-12

# Coordinate d is positive from the corrugation reference into the channel.
unshifted_epoxy = (-epoxy, 0.0)
unshifted_baseplate = (-(epoxy + baseplate), -epoxy)
inward_epoxy = tuple(value - clearance for value in unshifted_epoxy)
inward_baseplate = tuple(value - clearance for value in unshifted_baseplate)
outward_epoxy = tuple(value - outward_offset for value in unshifted_epoxy)
outward_baseplate = tuple(value - outward_offset for value in unshifted_baseplate)


def overlap(first, second):
    return min(first[1], second[1]) - max(first[0], second[0])


for first in (inward_epoxy, inward_baseplate):
    for second in (outward_epoxy, outward_baseplate):
        assert overlap(first, second) <= 1e-12

cpp = CPP.read_text()
assert 'description.constant<double>("SiEndcapTilingInwardModule_z_clearance")' in cpp
assert 'description.constant<double>("SiEndcapTilingOppositeFacing_z_offset")' in cpp
assert 'module_template.outward_facing' in cpp
assert 'applied_module_normal_z_offset_mm' in cpp
assert 'placed_sensor_reference_z_mm' in cpp

print(f"inward clearance={clearance:.3f} mm away from the channel")
print(
    f"outward offset={outward_offset:.3f} mm = common clearance {clearance:.3f} "
    f"+ preserved relative offset {relative_offset:.3f} mm"
)
print(f"inward epoxy/baseplate:  {inward_epoxy}, {inward_baseplate} mm relative to surface")
print(f"outward epoxy/baseplate: {outward_epoxy}, {outward_baseplate} mm relative to surface")
print("PASS: translated backing stacks touch at one boundary and do not overlap")
print("NOTE: placed sensor references shift by -0.030 mm inward and +0.260 mm outward")
