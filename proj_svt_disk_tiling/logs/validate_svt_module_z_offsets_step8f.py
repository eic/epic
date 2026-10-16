#!/usr/bin/env python3
"""Validate the Step 8f module translations."""

from decimal import Decimal
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"


def literal_mm(constants, name):
    value = constants[name]
    assert value.endswith("*mm"), (name, value)
    return Decimal(value[:-3])


def shifted(interval, offset):
    return tuple(value + offset for value in interval)


def overlap(first, second):
    return min(first[1], second[1]) - max(first[0], second[0])


def main():
    root = ET.parse(XML).getroot()
    constants = {
        node.attrib["name"]: node.attrib["value"]
        for node in root.iter("constant")
    }
    inward_clearance = literal_mm(
        constants, "SiEndcapTilingInwardModule_z_clearance"
    )
    facing_offset = literal_mm(
        constants, "SiEndcapTilingOppositeFacing_z_offset"
    )
    assert inward_clearance == Decimal("0.030")
    assert facing_offset == Decimal("0.230")

    inward_applied = -inward_clearance
    outward_applied = inward_clearance + facing_offset
    assert outward_applied == Decimal("0.260")
    assert outward_applied + inward_applied == facing_offset

    # Coordinate d is positive from the corrugation surface into the channel.
    inward = {
        "baseplate": shifted((Decimal("-0.230"), Decimal("-0.080")), inward_applied),
        "external_epoxy": shifted((Decimal("-0.080"), Decimal("0.000")), inward_applied),
        "film": shifted((Decimal("-0.080"), Decimal("-0.020")), inward_applied),
        "sensor": shifted((Decimal("-0.020"), Decimal("0.030")), inward_applied),
        "fpc_kapton": shifted((Decimal("-0.080"), Decimal("0.000")), inward_applied),
        "fpc_aluminum": shifted((Decimal("0.000"), Decimal("0.030")), inward_applied),
        "ancasic_glue": shifted((Decimal("0.030"), Decimal("0.110")), inward_applied),
        "ancasic": shifted((Decimal("0.110"), Decimal("0.410")), inward_applied),
    }
    assert inward["fpc_kapton"] == (Decimal("-0.110"), Decimal("-0.030"))
    assert inward["fpc_aluminum"] == (Decimal("-0.030"), Decimal("0.000"))
    assert inward["sensor"] == (Decimal("-0.050"), Decimal("0.000"))
    assert inward["ancasic"] == (Decimal("0.080"), Decimal("0.380"))
    carbon_flat = (Decimal("0.000"), Decimal("0.150"))
    assert overlap(inward["fpc_kapton"], carbon_flat) < 0
    assert overlap(inward["fpc_aluminum"], carbon_flat) == 0

    # Express the outward backing stack in the same +channel coordinate. Its
    # positive outward-normal displacement therefore appears with a minus sign.
    outward_backing = {
        "external_epoxy": shifted(
            (Decimal("-0.080"), Decimal("0.000")), -outward_applied
        ),
        "baseplate": shifted(
            (Decimal("-0.230"), Decimal("-0.080")), -outward_applied
        ),
    }
    assert outward_backing["baseplate"] == (
        Decimal("-0.490"), Decimal("-0.340")
    )
    assert outward_backing["external_epoxy"] == (
        Decimal("-0.340"), Decimal("-0.260")
    )
    assert inward["baseplate"][0] == outward_backing["external_epoxy"][1]
    for inward_piece in (inward["baseplate"], inward["external_epoxy"]):
        for outward_piece in outward_backing.values():
            assert overlap(inward_piece, outward_piece) <= 0

    selected = {
        node.attrib["file"]
        for node in root.iter("module_placements")
        if node.attrib.get("format") == "tiling-directory"
    }
    assert selected == {
        "compact/tracking/silicon_disks_tiling/20261001_ES/rsu_opt"
    }

    cpp = CPP.read_text()
    assert 'description.constant<double>("SiEndcapTilingInwardModule_z_clearance")' in cpp
    assert 'description.constant<double>("SiEndcapTilingOppositeFacing_z_offset")' in cpp
    assert 'applied_module_normal_z_offset_mm' in cpp
    assert 'applied_outward_z_offset_mm' not in cpp

    print(
        f"applied offsets along module normal: inward={inward_applied} mm, "
        f"outward={outward_applied} mm; relative={facing_offset} mm"
    )
    print(f"inward component intervals from corrugation surface: {inward}")
    print(f"outward backing intervals in +channel coordinate: {outward_backing}")
    print("PASS: FPC and sensor end at the carbon boundary")
    print("NOTE: inward edge epoxy now ends 0.030 mm before the corrugation surface")


if __name__ == "__main__":
    main()
