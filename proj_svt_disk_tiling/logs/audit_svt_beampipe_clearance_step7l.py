#!/usr/bin/env python3
"""Audit SVT physical beampipe openings without changing detector geometry.

The branch-local XML openings were generated for the former 50 mm endcap
envelope.  This script compares their pipe-axis centers with the far face of
the current 10 mm tiled layer.  All values are in mm and radians.
"""

import math


LAYER_HALF_THICKNESS = 5.0
FORMER_LAYER_HALF_THICKNESS = 25.0
BAKEOUT_BUFFER = 5.0

# Current central_beampipe.xml / compact/definitions.xml parameters.
UPSTREAM_AXIS_ANGLE = 0.025
DOWNSTREAM_AXIS_ANGLE = 0.045 / 2.0
DOWNSTREAM_AXIS_INTERSECTION = 661.0

# Disk centers from the delivered placement package and current SVT XML.
DISKS = {
    "ED3": -850.0,
    "ED4": -1020.0,
    "HD2": 700.0,
    "HD3b": 950.0,
    "HD4": 1200.0,
}

# Historical pre-Step-7n XML openings before the common 5 mm buffer. They are
# retained here to make the Step 7l provenance finding reproducible.
PRE_STEP7N_XML_OPENINGS = {
    "ED3": ((0.0, 33.500), (-21.880, 14.285)),
    "ED4": ((0.0, 33.500), (-26.130, 14.285)),
    "HD2": ((0.0, 31.757), (-1.440, 35.488)),
    "HD3b": ((0.0, 31.757), (-7.066, 41.111)),
    "HD4": ((0.0, 31.757), (-12.692, 46.734)),
}

# Source-layout openings are reproduced here only for comparison.  They remain
# source data and are not proposed as physical layer-envelope dimensions.
METADATA_OPENINGS = {
    "ED3": ((0.0, 34.900), (-21.2544, 15.900)),
    "ED4": ((0.0, 34.900), (-25.5053, 15.900)),
    "HD2": ((3.1360, 37.1360),),
    "HD3b": ((8.7360, 42.7360),),
    "HD4": ((14.3360, 48.3360),),
}


def far_face(z_center):
    """Return the layer face farther from the interaction point."""
    return z_center + math.copysign(LAYER_HALF_THICKNESS, z_center)


def expected_displaced_center(disk, half_thickness):
    """Expected displaced-pipe center in the SVT layer-local x coordinate."""
    z_center = DISKS[disk]
    z_abs = abs(z_center) + half_thickness
    if disk.startswith("ED"):
        return -z_abs * math.tan(UPSTREAM_AXIS_ANGLE)
    return -(z_abs - DOWNSTREAM_AXIS_INTERSECTION) * math.tan(DOWNSTREAM_AXIS_ANGLE)


print(
    "disk  z_far  metadata displaced (cx,r)  pre-7n XML (cx,r+5)      "
    "expected cx@5  former cx@25  excess projection"
)
for disk in DISKS:
    metadata = METADATA_OPENINGS[disk][-1]
    xml = PRE_STEP7N_XML_OPENINGS[disk][-1]
    expected = expected_displaced_center(disk, LAYER_HALF_THICKNESS)
    former = expected_displaced_center(disk, FORMER_LAYER_HALF_THICKNESS)
    excess = abs(former - expected)
    print(
        f"{disk:4s} {far_face(DISKS[disk]):+7.1f}  "
        f"({metadata[0]:+8.4f},{metadata[1]:7.4f})       "
        f"({xml[0]:+7.3f},{xml[1] + BAKEOUT_BUFFER:6.3f})        "
        f"{expected:+9.4f}    {former:+9.4f}       {excess:6.4f}"
    )

    # Rounded XML centers agree with the former 25 mm half-layer construction.
    assert abs(xml[0] - former) < 5.0e-4
    # The old envelope is conservative by less than 0.51 mm of axis projection.
    assert 0.44 < excess < 0.51

print()
print("Historical pre-Step-7n centered physical openings after buffer:")
for disk in DISKS:
    center, radius = PRE_STEP7N_XML_OPENINGS[disk][0]
    print(f"  {disk}: cx={center:+.3f} r={radius + BAKEOUT_BUFFER:.3f}")

print()
print("PASS: pre-Step-7n XML openings followed the physical pipe direction and retained")
print("      a conservative former-envelope projection; metadata must not replace")
print("      them when constructing the layer mother or clipping carbon support")
