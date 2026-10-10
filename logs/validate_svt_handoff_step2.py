#!/usr/bin/env python3
"""Read-only source-package checks for step 2; write a reproducible log.

Run: python3 logs/validate_svt_handoff_step2.py
Afterward, review logs/validate_svt_handoff_step2.txt and update the step log.
"""

import csv
import hashlib
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "ES_material"
OUTPUT = ROOT / "logs" / "validate_svt_handoff_step2.txt"
DISKS = {
    "ED4": (0, -1020), "ED3": (1, -850), "ED2": (2, -650),
    "ED1": (3, -450), "ED0": (4, -250), "HD0": (5, 250),
    "HD1": (6, 450), "HD2": (7, 700), "HD3b": (8, 950),
    "HD4": (9, 1200),
}
PLACEMENT_COLUMNS = {
    "disk_id", "row_index", "mod_index", "type_id", "x_origin_mm",
    "y_origin_mm", "z_baseplate_mm", "z_sensor_mm",
}


def table(path, columns):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = columns - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        return list(reader)


def metadata(path):
    """Read key/value metadata without changing the delivered file."""
    values = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def check_package(name, lines, issues):
    """Validate IDs, type references, metadata totals, and source hashes."""
    folder = SOURCE / name
    expected = {"catalog.csv"}
    for disk in DISKS:
        expected.update((f"{disk}_metadata.txt", f"{disk}_modules.csv"))
    actual = {file.name for file in folder.iterdir() if file.is_file()}
    if actual != expected:
        issues.append(f"{name}: missing={sorted(expected - actual)}, extra={sorted(actual - expected)}")
    catalog_rows = table(folder / "catalog.csv", {
        "type_id", "rsu_count", "module_height_mm", "module_length_mm",
        "active_height_mm", "active_length_mm", "sensor_height_mm",
        "sensor_length_mm", "BLEC_mm", "BREC_mm", "LEC_mm", "REC_mm",
        "PER_mm", "facing", "chirality",
    })
    catalog = {row["type_id"]: row for row in catalog_rows}
    if len(catalog) != len(catalog_rows):
        issues.append(f"{name}: duplicate catalogue type ID")
    for kind in catalog_rows:
        code = kind["type_id"]
        if len(code) != 4 or code[0] not in "AB" or code[1] not in "56" or code[2] not in "01" or code[3] not in "01":
            issues.append(f"{name}: malformed type ID {code}")
            continue
        if int(kind["rsu_count"]) != int(code[1]) or float(kind["module_height_mm"]) != {"A": 35.0, "B": 32.0}[code[0]]:
            issues.append(f"{name}: type count/width mismatch for {code}")
        if kind["facing"] != {"0": "OUTWARD", "1": "INWARD"}[code[2]] or kind["chirality"] != {"0": "RIGHT", "1": "LEFT"}[code[3]]:
            issues.append(f"{name}: facing/chirality mismatch for {code}")
        active_length = float(kind["active_length_mm"])
        if abs(active_length - int(code[1]) * 21.666) > 0.0001:
            issues.append(f"{name}: RSU-chain length mismatch for {code}")
        if abs(float(kind["sensor_length_mm"]) - active_length - float(kind["LEC_mm"]) - float(kind["REC_mm"])) > 0.0001:
            issues.append(f"{name}: sensor length mismatch for {code}")
        if abs(float(kind["module_length_mm"]) - active_length - float(kind["BLEC_mm"]) - float(kind["BREC_mm"])) > 0.0001:
            issues.append(f"{name}: baseplate length mismatch for {code}")
        if abs(float(kind["active_height_mm"]) - float(kind["sensor_height_mm"]) + 2 * float(kind["PER_mm"])) > 0.0001:
            issues.append(f"{name}: nominal active height mismatch for {code}")
    counts = Counter()
    area_total = 0.0
    versions = set()
    lines.append(f"PACKAGE {name}: {folder}")
    lines.append("disk count 5RSU 6RSU nominal_active_mm2")
    for disk, (disk_id, z_center) in sorted(DISKS.items(), key=lambda pair: pair[1][0]):
        rows = table(folder / f"{disk}_modules.csv", PLACEMENT_COLUMNS)
        info = metadata(folder / f"{disk}_metadata.txt")
        versions.add((info["catalog_version"], info["code_version"], info["generated_at"]))
        keys = set()
        disk_counts = Counter()
        area = 0.0
        for row in rows:
            key = (int(row["disk_id"]), int(row["row_index"]), int(row["mod_index"]))
            if key in keys or key[0] != disk_id:
                issues.append(f"{name}/{disk}: duplicate or wrong-disk key {key}")
            keys.add(key)
            kind = catalog.get(row["type_id"])
            if kind is None:
                issues.append(f"{name}/{disk}: unknown type {row['type_id']}")
                continue
            disk_counts[int(kind["rsu_count"])] += 1
            area += float(kind["active_height_mm"]) * float(kind["active_length_mm"])
            if abs(abs(float(row["z_baseplate_mm"]) - z_center) - 3.0) > 0.0001:
                issues.append(f"{name}/{disk}: unexpected baseplate z at {key}")
            offset = abs(float(row["z_sensor_mm"]) - float(row["z_baseplate_mm"]))
            if min(abs(offset - 0.03), abs(offset - 0.34)) > 0.0001:
                issues.append(f"{name}/{disk}: unexpected sensor z offset at {key}")
        if int(info["disk_id"]) != disk_id or float(info["z_center_mm"]) != z_center:
            issues.append(f"{name}/{disk}: metadata disk ID or z mismatch")
        if int(info["total_modules"]) != len(rows):
            issues.append(f"{name}/{disk}: metadata module count mismatch")
        if info["outer_radius_target"] != "module" or abs(float(info["outer_radius_mm"]) - float(info["module_r_out_mm"])) > 0.0001:
            issues.append(f"{name}/{disk}: metadata module outer-radius mismatch")
        if int(info["n_opening_primitives"]) != 2:
            issues.append(f"{name}/{disk}: unexpected opening count")
        for index in (0, 1):
            for suffix in ("cx_mm", "cy_mm", "a_mm", "b_mm"):
                float(info[f"opening_primitive_{index}_{suffix}"])
        if abs(float(info["total_active_area_mm2"]) - area) > 0.001:
            issues.append(f"{name}/{disk}: metadata nominal area mismatch")
        counts.update(disk_counts)
        area_total += area
        lines.append(f"{disk:4} {len(rows):5} {disk_counts[5]:4} {disk_counts[6]:4} {area:.4f}")
    if len(versions) != 1:
        issues.append(f"{name}: inconsistent metadata versions {sorted(versions)}")
    lines.append(f"TOTAL {sum(counts.values())} modules; 5RSU={counts[5]}, 6RSU={counts[6]}, nominal_active_mm2={area_total:.6f}")
    lines.append(f"VERSIONS {sorted(versions)}")
    lines.append("SHA256 (source files, relative to ES_material):")
    for filename in sorted(expected):
        file = folder / filename
        if file.is_file():
            digest = hashlib.sha256(file.read_bytes()).hexdigest()
            lines.append(f"{digest}  {name}/{filename}")
    return counts


def main():
    lines = ["SVT disk handoff validation; source files read-only"]
    issues = []
    six = check_package("rsu6", lines, issues)
    mixed = check_package("rsu_opt", lines, issues)
    if six != Counter({6: 2160}):
        issues.append(f"rsu6: unexpected overall types {dict(six)}")
    if mixed != Counter({5: 265, 6: 1899}):
        issues.append(f"rsu_opt: unexpected overall types {dict(mixed)}")
    lines.append(f"RESULT {'PASS' if not issues else 'FAIL'}; issues={len(issues)}")
    lines.extend(f"ISSUE {issue}" for issue in issues)
    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"{'PASS' if not issues else 'FAIL'}: {OUTPUT}; issues={len(issues)}")
    raise SystemExit(bool(issues))


if __name__ == "__main__":
    main()
