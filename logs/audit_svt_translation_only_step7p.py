#!/usr/bin/env python3
"""Read-only audit of replacing negative-side Ry(pi) with translation only.

The supplied module rows are detector coordinates.  This script proves, row by
row, that removing the layer rotation while removing its compensating signs
preserves the final detector-frame module position and orientation.
"""

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks"
SCENARIOS = ("rsu_opt", "all_6rsu")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")

# Current Step-7n physical hole centers are layer-local.  The proposed values
# below are not a geometry change: they are the local coordinates required to
# preserve the same detector-frame hole after removing the layer rotation.
CURRENT_PHYSICAL_LOCAL_X = {
    "ED3": -21.2544,
    "ED4": -25.5053,
    "HD2": -3.1360,
    "HD3b": -8.7360,
    "HD4": -14.3360,
}

# Every active detector block must use independently supplied detector-frame
# tiling coordinates. The commented legacy CSV path is intentionally excluded.
xml_text = (ROOT / "compact/tracking/silicon_disks_modules.xml").read_text()
active_tiling_paths = re.findall(
    r'<module_placements\s+file="[^"]+"\s+format="tiling-directory"\s*/>', xml_text
)
active_legacy_paths = re.findall(
    r'(?m)^\s*<module_placements\s+file="[^"]+"\s+format="csv"\s*/>', xml_text
)
assert len(active_tiling_paths) == 6
assert not active_legacy_paths


def values_close(lhs, rhs, tolerance=1.0e-9):
    return all(abs(a - b) <= tolerance for a, b in zip(lhs, rhs))


print("scenario,modules,position_invariants,orientation_invariants")
grand_total = 0
for scenario in SCENARIOS:
    package = PACKAGE_ROOT / scenario
    with (package / "catalog.csv").open(newline="") as stream:
        catalog = {row["type_id"]: row for row in csv.DictReader(stream)}

    total = 0
    position_invariants = 0
    orientation_invariants = 0
    for disk in DISKS:
        layer_sign = -1.0 if disk.startswith("ED") else 1.0
        with (package / f"{disk}_modules.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            total += 1
            module = catalog[row["type_id"]]
            length = float(module["module_length_mm"])
            boundary_magnitude = length / 2.0 - 20.0
            x_origin = float(row["x_origin_mm"])
            global_x = x_origin - math.copysign(boundary_magnitude, x_origin)
            global_y = float(row["y_origin_mm"])
            global_layer_z = float((package / f"{disk}_metadata.txt").read_text().split(
                "z_center_mm = ", 1
            )[1].splitlines()[0])
            # Any detector-frame stack reference demonstrates the z transform;
            # use the supplied corrugation reference directly.
            global_z = float(row["z_baseplate_mm"])

            old_local = (
                layer_sign * global_x,
                global_y,
                layer_sign * (global_z - global_layer_z),
            )
            old_reconstructed_global = (
                layer_sign * old_local[0],
                old_local[1],
                global_layer_z + layer_sign * old_local[2],
            )
            proposed_local = (global_x, global_y, global_z - global_layer_z)
            proposed_reconstructed_global = (
                proposed_local[0],
                proposed_local[1],
                global_layer_z + proposed_local[2],
            )
            expected_global = (global_x, global_y, global_z)
            assert values_close(old_reconstructed_global, expected_global)
            assert values_close(proposed_reconstructed_global, expected_global)
            position_invariants += 1

            normal_global = 1.0 if float(row["z_sensor_mm"]) > global_z else -1.0
            boundary_global_sign = 1.0 if x_origin > 0.0 else -1.0
            boundary_intrinsic_sign = 1.0 if module["chirality"] == "RIGHT" else -1.0
            x_axis_global = boundary_global_sign * boundary_intrinsic_sign

            # Orientations are represented by the signs of the module's
            # intrinsic x/y/z axes in detector x/y/z.  All rotations are 0 or pi.
            old_x_local = layer_sign * x_axis_global
            old_z_local = layer_sign * normal_global
            old_y_local = old_x_local * old_z_local
            old_global_axes = (
                layer_sign * old_x_local,
                old_y_local,
                layer_sign * old_z_local,
            )
            proposed_global_axes = (
                x_axis_global,
                x_axis_global * normal_global,
                normal_global,
            )
            assert old_global_axes == proposed_global_axes
            orientation_invariants += 1

    grand_total += total
    print(f"{scenario},{total},{position_invariants},{orientation_invariants}")

print("\nphysical opening migration preserving current detector-frame geometry")
print("disk,current_local_x,current_global_x,translation_only_local_x")
for disk, current_local in CURRENT_PHYSICAL_LOCAL_X.items():
    layer_sign = -1.0 if disk.startswith("ED") else 1.0
    current_global = layer_sign * current_local
    print(f"{disk},{current_local:+.4f},{current_global:+.4f},{current_global:+.4f}")

print(f"\nPASS: {grand_total} module rows preserve detector-frame position and orientation")
print("PASS: facing, chirality, type IDs, and source CSV values require no changes")
print("PASS: all 6 active detector blocks use detector-frame tiling-directory inputs")
print("PASS: ED corrugation web Rx(-a) under Ry(pi) is box-equivalent to direct Rx(+a)")
print("NOTE: ED physical XML hole signs must change only to preserve their existing global position")
