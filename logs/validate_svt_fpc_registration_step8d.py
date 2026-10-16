#!/usr/bin/env python3
"""Validate the Step 8d bridge-FPC registration on the tiling baseplate."""

from decimal import Decimal
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"


def literal(constants, name, unit):
    value = constants[name]
    suffix = "*" + unit
    assert value.endswith(suffix), (name, value)
    return Decimal(value[: -len(suffix)])


def overlap(first, second):
    return min(first[1], second[1]) - max(first[0], second[0])


def fpc_y_interval(module_width, target_length, top_offset, glue_width):
    top_edge = module_width / 2 + top_offset
    lower_glue_inner_edge = -module_width / 2 + glue_width
    effective_length = min(target_length, top_edge - lower_glue_inner_edge)
    return (top_edge - effective_length, top_edge)


def main():
    constants = {
        item.attrib["name"]: item.attrib["value"]
        for item in ET.parse(XML).iter("constant")
    }
    widths = {
        "A": literal(constants, "SiEndcapTilingBaseplateA_width", "mm"),
        "B": literal(constants, "SiEndcapTilingBaseplateB_width", "mm"),
    }
    glue_width = literal(constants, "SiEndcapTilingEpoxyGlueLine_width", "mm")
    left_target = literal(constants, "SiEndcapLeftBridgeFPC_length", "mm")
    left_offset = literal(constants, "SiEndcapLeftBridgeFPC_y_offset", "mm")
    right_target = literal(constants, "SiEndcapRightBridgeFPC_length", "mm")
    right_offset = literal(constants, "SiEndcapRightBridgeFPC_y_offset", "mm")
    gap = literal(constants, "SiEndcapBridgeFPC_LEC_REC_gap", "mm")
    assert gap == Decimal("0.5")

    for family, module_width in widths.items():
        lower_epoxy = (-module_width / 2, -module_width / 2 + glue_width)
        upper_epoxy = (module_width / 2 - glue_width, module_width / 2)
        intervals = {
            "left": fpc_y_interval(module_width, left_target, left_offset, glue_width),
            "right": fpc_y_interval(module_width, right_target, right_offset, glue_width),
        }
        for interval in intervals.values():
            assert overlap(interval, lower_epoxy) <= 0
            assert overlap(interval, upper_epoxy) <= 0
        print(
            f"{family}: left_y={intervals['left']} length={intervals['left'][1] - intervals['left'][0]} mm; "
            f"right_y={intervals['right']} length={intervals['right'][1] - intervals['right'][0]} mm"
        )

    lec = literal(constants, "SiEndcapLEC_length", "mm")
    rec = literal(constants, "SiEndcapREC_length", "mm")
    blec = literal(constants, "SiEndcapTilingBLEC_length", "mm")
    brec = literal(constants, "SiEndcapTilingBREC_length", "mm")
    left_fpc_x = literal(constants, "SiEndcapLeftBridgeFPC_width", "mm")
    right_fpc_x = literal(constants, "SiEndcapRightBridgeFPC_width", "mm")
    assert lec + gap + left_fpc_x <= blec
    assert rec + gap + right_fpc_x <= brec

    epoxy = literal(constants, "SiEndcapTilingExternalEpoxy_thickness", "um")
    baseplate = literal(constants, "SiEndcapTilingBaseplateCF_thickness", "um")
    film = literal(constants, "SiEndcapTilingFilm_thickness", "um")
    sensor = literal(constants, "SiEndcapTilingSensor_thickness", "um")
    kapton = literal(constants, "SiEndcapBridgeFPC_Kapton_thickness", "um")
    aluminum = literal(constants, "SiEndcapBridgeFPC_Aluminum_thickness", "um")
    ancasic_glue = literal(constants, "SiEndcapAncASICGlue_thickness", "um")
    ancasic = literal(constants, "SiEndcapAncASIC_thickness", "um")
    fpc = kapton + aluminum
    electronics = fpc + ancasic_glue + ancasic
    assert fpc == 110
    assert electronics == 490

    totals = {}
    relative_intervals = {}
    for facing in ("inward", "outward"):
        total = baseplate + max(film + sensor, electronics)
        if facing == "outward":
            total += epoxy
        bottom = -total / 2
        baseplate_bottom = bottom + (epoxy if facing == "outward" else 0)
        baseplate_top = baseplate_bottom + baseplate
        fpc_interval = (baseplate_top, baseplate_top + fpc)
        ancasic_glue_interval = (fpc_interval[1], fpc_interval[1] + ancasic_glue)
        ancasic_interval = (ancasic_glue_interval[1], ancasic_glue_interval[1] + ancasic)
        corrugation_reference = (
            bottom if facing == "outward" else bottom + baseplate + epoxy
        )
        relative_fpc = tuple(value - corrugation_reference for value in fpc_interval)
        relative_ancasic_glue = tuple(
            value - corrugation_reference for value in ancasic_glue_interval
        )
        relative_ancasic = tuple(
            value - corrugation_reference for value in ancasic_interval
        )
        assert fpc_interval[0] == baseplate_top
        assert ancasic_interval[1] == total / 2
        totals[facing] = total
        relative_intervals[facing] = (
            relative_fpc,
            relative_ancasic_glue,
            relative_ancasic,
        )
        print(
            f"{facing}: total={total} um prototype_relative_to_corrugation_before_offset: "
            f"FPC={relative_fpc} um AncASIC_glue={relative_ancasic_glue} um "
            f"AncASIC={relative_ancasic} um"
        )
    assert totals == {"inward": Decimal("640"), "outward": Decimal("720")}
    assert relative_intervals == {
        "inward": ((Decimal("-80"), Decimal("30")),
                    (Decimal("30"), Decimal("110")),
                    (Decimal("110"), Decimal("410"))),
        "outward": ((Decimal("230"), Decimal("340")),
                     (Decimal("340"), Decimal("420")),
                     (Decimal("420"), Decimal("720"))),
    }

    cpp = CPP.read_text()
    assert cpp.count("SiEndcapBridgeFPC_LEC_REC_gap") == 2
    assert cpp.count("bridge_fpc.stacked_in_z = false") == 1
    assert cpp.count("ancasic.stacked_in_z = false") == 1
    assert cpp.count("lower_epoxy_inner_edge") == 4
    assert cpp.count("const double baseplate_top_z") == 1
    print("PASS: FPCs touch the baseplate, stop at the epoxy, and follow LEC/REC by 0.5 mm")


if __name__ == "__main__":
    main()
