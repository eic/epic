#!/usr/bin/env python3
"""Analytical checks for the supplied five/six-RSU module dimensions."""

import csv
from decimal import Decimal
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"


def literal(constants, name, unit):
    value = constants[name]
    suffix = "*" + unit
    assert value.endswith(suffix), (name, value)
    return Decimal(value[: -len(suffix)])


def main():
    constants = {
        item.attrib["name"]: item.attrib["value"]
        for item in ET.parse(XML).iter("constant")
    }
    rsu = literal(constants, "SiEndcapRSU_length", "mm")
    blec = literal(constants, "SiEndcapTilingBLEC_length", "mm")
    brec = literal(constants, "SiEndcapTilingBREC_length", "mm")
    lec = literal(constants, "SiEndcapLEC_length", "mm")
    rec = literal(constants, "SiEndcapREC_length", "mm")
    glue_width = literal(constants, "SiEndcapTilingEpoxyGlueLine_width", "mm")
    widths = {
        "A": literal(constants, "SiEndcapTilingBaseplateA_width", "mm"),
        "B": literal(constants, "SiEndcapTilingBaseplateB_width", "mm"),
    }
    epoxy = literal(constants, "SiEndcapTilingExternalEpoxy_thickness", "um")
    carbon = literal(constants, "SiEndcapTilingBaseplateCF_thickness", "um")
    film = literal(constants, "SiEndcapTilingFilm_thickness", "um")
    sensor = literal(constants, "SiEndcapTilingSensor_thickness", "um")
    assert constants["SiEndcapTilingLEC_thickness"] == "SiEndcapTilingSensor_thickness"
    assert constants["SiEndcapTilingREC_thickness"] == "SiEndcapTilingSensor_thickness"
    assert epoxy + carbon + film + sensor == 340
    baseplate_to_sensor_center = carbon / 2 + film + sensor / 2
    assert baseplate_to_sensor_center == 160
    sensitive_y = (literal(constants, "SiEndcapRSU_width", "mm") -
                   2 * (literal(constants, "SiEndcapRSU_periphery_width", "mm") +
                        literal(constants, "SiEndcapRSU_bias_width", "mm")))
    sensitive_x = (rsu - 2 * literal(constants, "SiEndcapRSU_backbone_width", "mm") -
                   6 * literal(constants, "SiEndcapRSU_powerswitch_width", "mm"))
    assert (sensitive_x, sensitive_y) == (Decimal("21.426"), Decimal("18.394"))

    seen = set()
    placement_z_differences = set()
    for scenario in ("all_6rsu", "rsu_opt"):
        directory = ROOT / "compact/tracking/silicon_disks" / scenario
        path = directory / "catalog.csv"
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            name = row["type_id"]
            seen.add(name)
            n = int(row["rsu_count"])
            assert name[:2] in ("A6", "B5", "B6")
            assert name[2:] in ("00", "01", "10", "11")
            width = widths[name[0]]
            length = n * rsu + blec + brec
            chain = n * rsu
            assert length == Decimal(row["module_length_mm"])
            assert width == Decimal(row["module_height_mm"])
            assert chain + lec + rec == Decimal(row["sensor_length_mm"])
            assert 0 < 2 * glue_width < width
            assert sensor + film + carbon + epoxy == 340
            # Chirality swaps which end of the baseplate carries BLEC/LEC.
            lec_on_left = row["chirality"] == "LEFT"
            left_span = blec if lec_on_left else brec
            rsu_start = -length / 2 + left_span
            rsu_end = rsu_start + chain
            assert -length / 2 < rsu_start < rsu_end < length / 2
            left_chip = lec if lec_on_left else rec
            right_chip = rec if lec_on_left else lec
            assert rsu_start - left_chip >= -length / 2
            assert rsu_end + right_chip <= length / 2
        print(f"{scenario}: PASS, {len(rows)} types")
        for placement_file in directory.glob("*_modules.csv"):
            with placement_file.open() as stream:
                for placement in csv.DictReader(stream):
                    placement_z_differences.add(
                        Decimal(placement["z_sensor_mm"]) -
                        Decimal(placement["z_baseplate_mm"])
                    )
    assert seen == {group + suffix for group in ("A6", "B5", "B6")
                    for suffix in ("00", "01", "10", "11")}
    print("Module dimensions, end clearances, and epoxy widths: PASS")
    print("Slide stack (epoxy/carbon/film/sensor): 80/150/60/50 um")
    print("Baseplate-center to sensor-center separation: 160 um (before facing transform)")
    print("Detailed sensitive area (5/6 RSU, mm^2):",
          5 * sensitive_x * sensitive_y, 6 * sensitive_x * sensitive_y)
    print("CSV z_sensor - z_baseplate (mm):", sorted(placement_z_differences))


if __name__ == "__main__":
    main()
