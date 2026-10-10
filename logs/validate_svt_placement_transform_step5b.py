#!/usr/bin/env python3
"""Reconstruct every supplied module origin and sensor-z anchor."""

import csv
from collections import Counter
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "compact/tracking/silicon_disks"
DISKS = (
    ("ED4", "OuterTrackerEndcapN_disk4", -1),
    ("ED3", "OuterTrackerEndcapN_disk3", -1),
    ("ED2", "OuterTrackerEndcapN_disk2", -1),
    ("ED1", "MiddleTrackerEndcapN_disk1", -1),
    ("ED0", "InnerTrackerEndcapN_disk1", -1),
    ("HD0", "InnerTrackerEndcapP_disk1", 1),
    ("HD1", "MiddleTrackerEndcapP_disk1", 1),
    ("HD2", "OuterTrackerEndcapP_disk2", 1),
    ("HD3b", "OuterTrackerEndcapP_disk3", 1),
    ("HD4", "OuterTrackerEndcapP_disk4", 1),
)
EXPECTED_COUNTS = {"all_6rsu": 2160, "rsu_opt": 2164}

EPOXY = Decimal("0.080")
BASEPLATE = Decimal("0.150")
FILM = Decimal("0.060")
SENSOR = Decimal("0.050")
BRIDGE_FPC = Decimal("0.110")
ANCASIC = Decimal("0.380")


def metadata(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.split("#", 1)[0].strip()
    return values


def local_z_references(outward):
    serial_module = BASEPLATE + FILM + SENSOR + BRIDGE_FPC + ANCASIC
    total = serial_module + (EPOXY if outward else 0)
    bottom = -total / 2
    sensor_reference = bottom + (EPOXY if outward else 0) + BASEPLATE + FILM + SENSOR
    corrugation_reference = bottom if outward else bottom + BASEPLATE + EPOXY
    return total, corrugation_reference, sensor_reference


def main():
    for scenario, expected_count in EXPECTED_COUNTS.items():
        directory = ROOT / scenario
        with (directory / "catalog.csv").open() as stream:
            catalog = {row["type_id"]: row for row in csv.DictReader(stream)}
        keys = set()
        types = Counter()
        rotations = Counter()
        count = 0

        for disk_id, (stem, disk_key, layer_axis_sign) in enumerate(DISKS):
            meta = metadata(directory / f"{stem}_metadata.txt")
            assert int(meta["disk_id"]) == disk_id
            layer_center_global = Decimal(meta["z_center_mm"])
            assert (layer_center_global > 0) == (layer_axis_sign > 0)
            with (directory / f"{stem}_modules.csv").open() as stream:
                for row in csv.DictReader(stream):
                    assert int(row["disk_id"]) == disk_id
                    key = (disk_id, int(row["row_index"]), int(row["mod_index"]))
                    assert key not in keys
                    keys.add(key)
                    design = catalog[row["type_id"]]
                    outward = design["facing"] == "OUTWARD"
                    assert design["chirality"] in ("LEFT", "RIGHT")
                    length = Decimal(design["module_length_mm"])
                    boundary_magnitude = length / 2 - Decimal(design["BLEC_mm"])
                    boundary_local = boundary_magnitude * (
                        1 if design["chirality"] == "RIGHT" else -1
                    )
                    x_origin_global = Decimal(row["x_origin_mm"])
                    radial_sign = 1 if x_origin_global > 0 else -1
                    x_center_global = x_origin_global - radial_sign * boundary_magnitude
                    x_center_layer = layer_axis_sign * x_center_global
                    boundary_layer = layer_axis_sign * x_origin_global - x_center_layer

                    z_corrugation_global = Decimal(row["z_baseplate_mm"])
                    z_sensor_global = Decimal(row["z_sensor_mm"])
                    normal_global = 1 if z_sensor_global > z_corrugation_global else -1
                    normal_layer = layer_axis_sign * normal_global
                    x_axis_layer = boundary_layer / boundary_local
                    assert abs(x_axis_layer) == 1
                    rotation_y = 0 if normal_layer > 0 else 180
                    rotation_z = 0 if x_axis_layer * normal_layer > 0 else 180

                    total, corrugation_local, sensor_local = local_z_references(outward)
                    expected_delta = sensor_local - corrugation_local
                    assert expected_delta == (Decimal("0.340") if outward else Decimal("0.030"))
                    assert abs(z_sensor_global - z_corrugation_global) == expected_delta
                    module_center_global = z_corrugation_global - normal_global * corrugation_local
                    module_center_layer = layer_axis_sign * (
                        module_center_global - layer_center_global
                    )

                    reconstructed_x = x_center_global + (layer_axis_sign * x_axis_layer) * boundary_local
                    reconstructed_corrugation = module_center_global + normal_global * corrugation_local
                    reconstructed_sensor = module_center_global + normal_global * sensor_local
                    assert reconstructed_x == x_origin_global
                    assert reconstructed_corrugation == z_corrugation_global
                    assert reconstructed_sensor == z_sensor_global
                    assert total == (Decimal("0.830") if outward else Decimal("0.750"))
                    assert abs(module_center_layer) < Decimal("10")

                    types[row["type_id"]] += 1
                    rotations[(rotation_z, rotation_y)] += 1
                    count += 1

        assert count == expected_count
        assert set(types) == set(catalog)
        print(f"{scenario}: PASS, {count} modules, {len(types)} types")
        print(f"  disk map: {', '.join(f'{i}:{item[1]}' for i, item in enumerate(DISKS))}")
        print(f"  (rotation_z, rotation_y) degrees: {dict(sorted(rotations.items()))}")
    print("All LEC-boundary midpoints and both z reference surfaces reconstruct exactly.")


if __name__ == "__main__":
    main()
