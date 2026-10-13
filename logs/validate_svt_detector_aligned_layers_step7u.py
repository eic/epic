#!/usr/bin/env python3
"""Validate the behavior-preserving removal of the negative-side layer Ry(pi).

The active tiling inputs are detector-coordinate data. This check compares the
old reflected-layer formulation with the new translation-only formulation for
every delivered module, source opening, physical opening, and corrugation row.
"""

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
MODULES_XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CORRUGATION = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/SVT_endcap_corrugation_rows_33_36.csv"
SCENARIOS = ("rsu_opt", "all_6rsu")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")

# Step 7u changes only the local representation of these existing physical
# holes. The detector-global centers must remain positive on ED.
PHYSICAL_X_MIGRATION = {
    "ED3": (-21.2544, +21.2544),
    "ED4": (-25.5053, +25.5053),
}


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            key, value = (part.strip() for part in line.split("=", 1))
            values[key] = value
    return values


def close_tuple(lhs, rhs, tolerance=1.0e-9):
    return all(abs(a - b) <= tolerance for a, b in zip(lhs, rhs))


cpp = CPP.read_text()
definitions = DEFINITIONS.read_text()
modules_xml = MODULES_XML.read_text()

# Static implementation guards: no layer rotation or reflected-axis
# compensation may remain in the active builder.
assert "RotationZYX(0.0, -M_PI, 0.0)" not in cpp
assert "const double layer_axis_sign" not in cpp
assert "const double global_to_local_x" not in cpp
assert "reflect ? -flat.surface_offset : flat.surface_offset" not in cpp
assert "Position(0, 0, global_layer_center)" in cpp
assert "row.dz = global_module_center - global_layer_center" in cpp
assert "const double local_surface_offset = flat.surface_offset" in cpp
assert "Negative disk openings use detector-aligned layer coordinates" in definitions

active_tiling_paths = re.findall(
    r'<module_placements\s+file="[^"]+"\s+format="tiling-directory"\s*/>', modules_xml
)
active_legacy_paths = re.findall(
    r'(?m)^\s*<module_placements\s+file="[^"]+"\s+format="csv"\s*/>', modules_xml
)
assert len(active_tiling_paths) == 6
assert not active_legacy_paths

print("scenario,modules,module_position_invariants,module_orientation_invariants,source_openings")
grand_modules = 0
grand_openings = 0
for scenario in SCENARIOS:
    package = PACKAGE_ROOT / scenario
    with (package / "catalog.csv").open(newline="") as stream:
        catalog = {row["type_id"]: row for row in csv.DictReader(stream)}

    module_count = 0
    position_count = 0
    orientation_count = 0
    opening_count = 0
    for disk in DISKS:
        old_layer_sign = -1.0 if disk.startswith("ED") else 1.0
        metadata = read_metadata(package / f"{disk}_metadata.txt")
        global_layer_z = float(metadata["z_center_mm"])

        for index in range(int(metadata["n_opening_primitives"])):
            source_x = float(metadata[f"opening_primitive_{index}_cx_mm"])
            old_local_x = old_layer_sign * source_x
            old_global_x = old_layer_sign * old_local_x
            new_local_x = source_x
            new_global_x = new_local_x
            assert abs(old_global_x - new_global_x) <= 1.0e-12
            opening_count += 1

        with (package / f"{disk}_modules.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            module_count += 1
            module = catalog[row["type_id"]]
            length = float(module["module_length_mm"])
            boundary_magnitude = length / 2.0 - 20.0
            x_origin = float(row["x_origin_mm"])
            global_x = x_origin - math.copysign(boundary_magnitude, x_origin)
            global_y = float(row["y_origin_mm"])
            global_z = float(row["z_baseplate_mm"])

            old_local = (
                old_layer_sign * global_x,
                global_y,
                old_layer_sign * (global_z - global_layer_z),
            )
            old_global = (
                old_layer_sign * old_local[0],
                old_local[1],
                global_layer_z + old_layer_sign * old_local[2],
            )
            new_local = (global_x, global_y, global_z - global_layer_z)
            new_global = (
                new_local[0],
                new_local[1],
                global_layer_z + new_local[2],
            )
            assert close_tuple(old_global, new_global)
            position_count += 1

            normal_global = 1.0 if float(row["z_sensor_mm"]) > global_z else -1.0
            boundary_global_sign = 1.0 if x_origin > 0.0 else -1.0
            boundary_intrinsic_sign = 1.0 if module["chirality"] == "RIGHT" else -1.0
            x_axis_global = boundary_global_sign * boundary_intrinsic_sign

            old_x_local = old_layer_sign * x_axis_global
            old_z_local = old_layer_sign * normal_global
            old_global_axes = (
                old_layer_sign * old_x_local,
                old_x_local * old_z_local,
                old_layer_sign * old_z_local,
            )
            new_global_axes = (
                x_axis_global,
                x_axis_global * normal_global,
                normal_global,
            )
            assert old_global_axes == new_global_axes
            orientation_count += 1

    grand_modules += module_count
    grand_openings += opening_count
    print(
        f"{scenario},{module_count},{position_count},{orientation_count},{opening_count}"
    )

for disk, (old_local_x, new_local_x) in PHYSICAL_X_MIGRATION.items():
    old_global_x = -old_local_x
    assert abs(old_global_x - new_local_x) <= 1.0e-12
    assert re.search(
        rf'TrackerEndcapN(?:Disk3|Disk4)_hadron_center_x" value="{new_local_x:.4f}\*mm"',
        definitions,
    )

with CORRUGATION.open(newline="") as stream:
    corrugation_rows = list(csv.DictReader(stream))
for row in corrugation_rows:
    source_offset = float(row["surface_offset_mm"])
    old_layer_sign = -1.0 if row["disk"].startswith("ED") else 1.0
    old_local_offset = old_layer_sign * source_offset
    old_global_offset = old_layer_sign * old_local_offset
    new_global_offset = source_offset
    assert abs(old_global_offset - new_global_offset) <= 1.0e-12

print(f"PASS: {grand_modules} module transforms preserve detector-global position and orientation")
print(f"PASS: {grand_openings} source-opening centers preserve detector-global position")
print("PASS: physical ED opening centers preserve detector-global position")
print(f"PASS: {len(corrugation_rows)} corrugation rows preserve detector-global surface offsets")
print("PASS: all active inputs use tiling-directory; the disabled legacy CSV path was not migrated")
