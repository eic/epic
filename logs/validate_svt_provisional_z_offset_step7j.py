#!/usr/bin/env python3
"""Validate the provisional outward-module z offset against backing-stack slabs."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"

text = XML.read_text()
match = re.search(
    r'name="SiEndcapTilingOutwardModule_z_offset"\s+value="([0-9.]+)\*mm"', text
)
assert match
offset = float(match.group(1))

epoxy = 0.080
baseplate = 0.150
assert abs(offset - (epoxy + baseplate)) < 1e-12

# Coordinate d is positive from the corrugation reference into the channel.
inward_epoxy = (-epoxy, 0.0)
inward_baseplate = (-(epoxy + baseplate), -epoxy)
outward_epoxy = (inward_epoxy[0] - offset, inward_epoxy[1] - offset)
outward_baseplate = (inward_baseplate[0] - offset, inward_baseplate[1] - offset)


def overlap(first, second):
    return min(first[1], second[1]) - max(first[0], second[0])


for first in (inward_epoxy, inward_baseplate):
    for second in (outward_epoxy, outward_baseplate):
        assert overlap(first, second) <= 1e-12

cpp = CPP.read_text()
assert 'description.constant<double>("SiEndcapTilingOutwardModule_z_offset")' in cpp
assert 'module_template.outward_facing' in cpp
assert 'applied_outward_z_offset_mm' in cpp
assert 'placed_sensor_reference_z_mm' in cpp

print(f"offset={offset:.3f} mm = epoxy {epoxy:.3f} + baseplate {baseplate:.3f} mm")
print(f"inward epoxy/baseplate:  {inward_epoxy}, {inward_baseplate} mm relative to surface")
print(f"outward epoxy/baseplate: {outward_epoxy}, {outward_baseplate} mm relative to surface")
print("PASS: provisional translated backing stacks do not overlap")
print("NOTE: outward placed sensor references differ from source by 0.230 mm")
