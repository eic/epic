#!/usr/bin/env python3
"""Check the endcap RSU layer stack, sensor reference, and Cu clearance."""

from decimal import Decimal
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"


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
    inactive = literal(constants, "SiEndcapTilingRSUInactiveSi_thickness", "um")
    active = literal(constants, "SiEndcapTilingRSUActiveSi_thickness", "um")
    copper = literal(constants, "SiEndcapTilingRSUCu_thickness", "um")
    silicon = literal(constants, "SiEndcapTilingSensor_thickness", "um")
    film = literal(constants, "SiEndcapTilingFilm_thickness", "um")
    epoxy = literal(constants, "SiEndcapTilingExternalEpoxy_thickness", "um")
    carbon = literal(constants, "SiEndcapTilingBaseplateCF_thickness", "um")
    clearance = literal(constants, "SiEndcapTilingInwardModule_z_clearance", "mm") * 1000
    opposite = literal(constants, "SiEndcapTilingOppositeFacing_z_offset", "mm") * 1000

    assert (inactive, active, copper) == (36, 14, Decimal("2.5"))
    assert inactive + active == silicon == 50
    assert silicon + copper == Decimal("52.5")
    source_inward_delta = silicon + film - epoxy
    source_outward_delta = silicon + film + carbon + epoxy
    assert (source_inward_delta, source_outward_delta) == (30, 340)
    assert clearance == source_inward_delta + copper == Decimal("32.5")
    assert opposite == 230

    # Distances from the inward-facing corrugation surface after translation.
    silicon_top = source_inward_delta - clearance
    inactive_interval = (silicon_top - silicon, silicon_top - active)
    active_interval = (inactive_interval[1], silicon_top)
    copper_interval = (silicon_top, silicon_top + copper)
    assert inactive_interval == (Decimal("-52.5"), Decimal("-16.5"))
    assert active_interval == (Decimal("-16.5"), Decimal("-2.5"))
    assert copper_interval == (Decimal("-2.5"), Decimal("0"))

    cpp = CPP.read_text()
    for token in (
        'Box inactive_solid(rsu_pitch_x / 2.0, comp_y / 2.0, inactive_si / 2.0)',
        'Box active_solid(rsu_pitch_x / 2.0, comp_y / 2.0, active_si / 2.0)',
        'Box copper_solid(rsu_pitch_x / 2.0, comp_y / 2.0, copper / 2.0)',
        'Volume(base_name + "_active_si", active_solid, material)',
        'active_layer_volume.placeVolume(',
        'prototype.active_layers.push_back(active_pv)',
        'prototype.sensitive_active_layer_indices.push_back(repeat_idx)',
        'active_elements.at(prototype.sensitive_active_layer_indices.at(idx))',
        'description.visAttributes("InvisibleWithDaughters")',
        'description.visAttributes("SVTInactiveSiVis")',
        'description.visAttributes("SVTRSUCopperVis")',
        'description.material("Copper")',
        'attach_sensitive(box_volume, pv, box_z - z_position)',
        'attach_sensitive(box_volume, pv, active_center_z - z_position)',
        'return bottom + epoxy + description.constant<double>("SiEndcapTilingBaseplateCF_thickness")',
    ):
        assert token in cpp, token
    assert 'const double sensor_branch = description.constant<double>("SiEndcapTilingFilm_thickness") +' in cpp
    for name in ("display.xml", "display_detailed.xml", "display_geoviewer.xml"):
        display = (ROOT / "compact" / name).read_text()
        assert '<vis name="SVTInactiveSiVis"' in display
        assert '<vis name="SVTRSUCopperVis"' in display
        assert 'ref="AnlTeal"' in display
    print("RSU layers (um from inward corrugation surface):")
    print(f"  inactive Si: {inactive_interval}; active Si: {active_interval}; Cu: {copper_interval}")
    print("PASS: 50 um Si + 2.5 um Cu; CSV Si-top reference retained; Cu reaches carbon without overlap")
    print(f"Module-normal offsets: inward=-{clearance / 1000} mm, outward=+{(clearance + opposite) / 1000} mm")


if __name__ == "__main__":
    main()
