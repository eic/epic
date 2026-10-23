#!/usr/bin/env python3
"""Validate Step 8g's continuous, module-independent corrugation support."""

from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"


def main():
    xml_root = ET.parse(XML).getroot()
    constants = {
        node.attrib["name"]: node.attrib["value"]
        for node in xml_root.iter("constant")
    }

    # Step 8f's module translations remain the sole clearance adjustment.
    assert constants["SiEndcapTilingInwardModule_z_clearance"] == "0.030*mm"
    assert constants["SiEndcapTilingOppositeFacing_z_offset"] == "0.230*mm"

    # The virtual keep-clear dimensions only drove the removed ad-hoc cuts.
    assert not any(name.startswith("SiEndcapTilingKeepClear") for name in constants)

    selected = {
        node.attrib["file"]
        for node in xml_root.iter("module_placements")
        if node.attrib.get("format") == "tiling-directory"
    }
    assert selected == {
        "compact/tracking/silicon_disks_tiling/20261001_ES/rsu_opt"
    }

    cpp = CPP.read_text()
    for removed_token in (
        "CorrugationKeepClear",
        "keep_clears",
        "keep_clears.push_back",
        "SiEndcapTilingKeepClear",
    ):
        assert removed_token not in cpp

    signature = (
        "int place_corrugated_frame(Detector& description, Volume& layer_vol, "
        "const DiskBoundary& disk,\n                           "
        "const CorrugatedFrameConfig& config)"
    )
    call = "place_corrugated_frame(description, layer_vol, disk, frame_config)"
    assert signature in cpp
    assert call in cpp

    # Boundary/opening clipping remains: only module-dependent clipping is gone.
    assert "contained_x_intervals(disk, footprint)" in cpp
    assert "build_disk_solid(layer_name, disk)" in cpp

    print("active placement: rsu_opt")
    print("module-normal offsets: inward=-0.030 mm, outward=+0.260 mm")
    print("PASS: corrugation has no module-dependent x cutouts")
    print("PASS: disk-boundary and physical-opening clipping remain enabled")


if __name__ == "__main__":
    main()
