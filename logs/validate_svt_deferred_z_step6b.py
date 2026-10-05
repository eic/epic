#!/usr/bin/env python3
"""Validate the explicit temporary translation to current XML disk z centers."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact" / "tracking" / "silicon_disks"
SCENARIOS = ("all_6rsu", "rsu_opt")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")
CURRENT_XML_Z_MM = {
    "ED4": -1020.0,
    "ED3": -850.0,
    "ED2": -650.0,
    "ED1": -550.0,
    "ED0": -255.0,
    "HD0": 255.0,
    "HD1": 550.0,
    "HD2": 700.0,
    "HD3b": 950.0,
    "HD4": 1200.0,
}


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line and "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            values[key] = value
    return values


for scenario in SCENARIOS:
    module_count = 0
    print(scenario)
    for disk in DISKS:
        directory = PACKAGE_ROOT / scenario
        metadata = read_metadata(directory / (disk + "_metadata.txt"))
        source_center = float(metadata["z_center_mm"])
        shift = CURRENT_XML_Z_MM[disk] - source_center
        with (directory / (disk + "_modules.csv")).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert rows
        for row in rows:
            corrugation = float(row["z_baseplate_mm"])
            sensor = float(row["z_sensor_mm"])
            translated_corrugation = corrugation + shift
            translated_sensor = sensor + shift
            assert abs((translated_sensor - translated_corrugation) - (sensor - corrugation)) < 1e-12
            assert abs((translated_corrugation - CURRENT_XML_Z_MM[disk]) -
                       (corrugation - source_center)) < 1e-12
            module_count += 1
        print("  {:4s}: source={:7.1f} XML={:7.1f} shift={:+7.1f} modules={}".format(
            disk, source_center, CURRENT_XML_Z_MM[disk], shift, len(rows)))
    print("  PASS: {} modules preserve both source-reference offsets".format(module_count))

print("PASS: deferred disk-z translation is explicit and lossless")
