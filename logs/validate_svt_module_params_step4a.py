#!/usr/bin/env python3
"""Check staged tiling dimensions against both delivered module catalogues."""

import csv
from decimal import Decimal
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
SCENARIOS = ("all_6rsu", "rsu_opt")


def number(constants, name, unit):
    value = constants[name]
    suffix = "*" + unit
    if not value.endswith(suffix):
        raise ValueError(f"{name}: expected a {unit} literal, found {value}")
    return Decimal(value[: -len(suffix)])


def main():
    tree = ET.parse(XML)
    constants = {
        item.attrib["name"]: item.attrib["value"]
        for item in tree.iter("constant")
    }
    rsu_length = number(constants, "SiEndcapRSU_length", "mm")
    lec = number(constants, "SiEndcapLEC_length", "mm")
    rec = number(constants, "SiEndcapREC_length", "mm")
    blec = number(constants, "SiEndcapTilingBLEC_length", "mm")
    brec = number(constants, "SiEndcapTilingBREC_length", "mm")
    widths = {
        "A": number(constants, "SiEndcapTilingBaseplateA_width", "mm"),
        "B": number(constants, "SiEndcapTilingBaseplateB_width", "mm"),
    }
    for name, expected in (("Sensor", 50), ("Film", 60),
                           ("BaseplateCF", 150), ("ExternalEpoxy", 80)):
        actual = number(constants, f"SiEndcapTiling{name}_thickness", "um")
        assert actual == expected, (name, actual, expected)
    for n in (5, 6):
        prefix = f"SiEndcapTiling{n}RSU"
        assert constants[f"{prefix}_package_length"] == (
            f"{prefix}_count*SiEndcapRSU_length + "
            "SiEndcapTilingBLEC_length + SiEndcapTilingBREC_length"
        )
        assert constants[f"{prefix}_sensor_length"] == (
            f"{prefix}_count*SiEndcapRSU_length + "
            "SiEndcapLEC_length + SiEndcapREC_length"
        )

    for scenario in SCENARIOS:
        count = 0
        with (ROOT / "compact/tracking/silicon_disks" / scenario / "catalog.csv").open() as stream:
            for row in csv.DictReader(stream):
                n = int(row["rsu_count"])
                assert n == int(constants[f"SiEndcapTiling{n}RSU_count"])
                package = n * rsu_length + blec + brec
                sensor = n * rsu_length + lec + rec
                assert package == Decimal(row["module_length_mm"]), row["type_id"]
                assert sensor == Decimal(row["sensor_length_mm"]), row["type_id"]
                assert widths[row["type_id"][0]] == Decimal(row["module_height_mm"]), row["type_id"]
                assert blec == Decimal(row["BLEC_mm"]), row["type_id"]
                assert brec == Decimal(row["BREC_mm"]), row["type_id"]
                assert lec == Decimal(row["LEC_mm"]), row["type_id"]
                assert rec == Decimal(row["REC_mm"]), row["type_id"]
                count += 1
        print(f"{scenario}: PASS, {count} catalogue types")
    print("XML parse and slide stack: PASS")
    print(f"5 RSU: package {5 * rsu_length + blec + brec} mm, sensor {5 * rsu_length + lec + rec} mm")
    print(f"6 RSU: package {6 * rsu_length + blec + brec} mm, sensor {6 * rsu_length + lec + rec} mm")


if __name__ == "__main__":
    main()
