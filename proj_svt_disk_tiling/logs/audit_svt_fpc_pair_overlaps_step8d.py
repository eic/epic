#!/usr/bin/env python3
"""Reconstruct bridge-FPC versus LEC/REC film overlaps for Step 8d.

The DD4hep module prototypes are assemblies, so Geant names collisions using
assembly copy numbers that do not identify the source CSV row.  This audit
applies the same engineer-to-EPIC x reflection, handedness reflection, and
module rotations as SiEndcapModuleTracker_geo.cpp, then reports the original
disk/row/module keys for every FPC/end-film intersection.
"""

import argparse
import csv
import itertools
from collections import namedtuple
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
XML = ROOT / "compact/tracking/silicon_disks_modules.xml"

RSU_LENGTH = 21.666
FILM_THICKNESS = 0.060
SENSOR_THICKNESS = 0.050
EXTERNAL_EPOXY_THICKNESS = 0.080
BASEPLATE_THICKNESS = 0.150
FPC_KAPTON_THICKNESS = 0.080
FPC_ALUMINUM_THICKNESS = 0.030
ANCASIC_GLUE_THICKNESS = 0.080
ANCASIC_THICKNESS = 0.300
INWARD_Z_CLEARANCE = 0.030
OPPOSITE_FACING_Z_OFFSET = 0.230
EPOXY_GLUE_LINE_WIDTH = 2.5
FPC_GAP = 0.5
LEFT_FPC_X = 10.0
LEFT_FPC_Y = 24.5
LEFT_FPC_TOP_OFFSET = -7.5
RIGHT_FPC_X = 4.0
RIGHT_FPC_Y = 25.78
RIGHT_FPC_TOP_OFFSET = -6.22


Box = namedtuple("Box", "disk row module type_id facing name x y z")


def source_key(box):
    return box.disk, box.row, box.module


def interval(center, width):
    return center - width / 2.0, center + width / 2.0


def transformed_interval(module_center, axis_sign, local_center, width):
    return interval(module_center + axis_sign * local_center, width)


def overlap(first, second):
    return min(first[1], second[1]) - max(first[0], second[0])


def active_package():
    paths = {
        node.attrib["file"]
        for node in ET.parse(XML).iter("module_placements")
        if node.attrib.get("format") == "tiling-directory"
    }
    if len(paths) != 1:
        raise RuntimeError(f"expected one active tiling package, found {sorted(paths)}")
    return ROOT / paths.pop()


def assert_xml_dimensions():
    constants = {
        node.attrib["name"]: node.attrib["value"]
        for node in ET.parse(XML).iter("constant")
    }
    expected = {
        "SiEndcapRSU_length": "21.666*mm",
        "SiEndcapTilingInwardModule_z_clearance": "0.030*mm",
        "SiEndcapTilingOppositeFacing_z_offset": "0.230*mm",
        "SiEndcapTilingSensor_thickness": "50*um",
        "SiEndcapTilingFilm_thickness": "60*um",
        "SiEndcapTilingBaseplateCF_thickness": "150*um",
        "SiEndcapTilingExternalEpoxy_thickness": "80*um",
        "SiEndcapTilingEpoxyGlueLine_width": "2.5*mm",
        "SiEndcapBridgeFPC_Kapton_thickness": "80*um",
        "SiEndcapBridgeFPC_Aluminum_thickness": "30*um",
        "SiEndcapLeftBridgeFPC_length": "24.5*mm",
        "SiEndcapLeftBridgeFPC_width": "10.0*mm",
        "SiEndcapLeftBridgeFPC_y_offset": "-7.5*mm",
        "SiEndcapRightBridgeFPC_length": "25.78*mm",
        "SiEndcapRightBridgeFPC_width": "4.0*mm",
        "SiEndcapRightBridgeFPC_y_offset": "-6.22*mm",
        "SiEndcapBridgeFPC_LEC_REC_gap": "0.5*mm",
        "SiEndcapAncASICGlue_thickness": "80*um",
        "SiEndcapAncASIC_thickness": "300*um",
    }
    differences = {
        name: (constants.get(name), value)
        for name, value in expected.items()
        if constants.get(name) != value
    }
    if differences:
        raise RuntimeError(f"audit dimensions disagree with XML: {differences}")


def fpc_y_geometry(module_width, target_length, top_offset):
    top_edge = module_width / 2.0 + top_offset
    lower_epoxy_inner_edge = -module_width / 2.0 + EPOXY_GLUE_LINE_WIDTH
    box_y = min(target_length, max(0.0, top_edge - lower_epoxy_inner_edge))
    return box_y, top_edge - box_y / 2.0


def local_z_geometry(facing):
    outward = facing == "OUTWARD"
    total = BASEPLATE_THICKNESS + max(
        FILM_THICKNESS + SENSOR_THICKNESS,
        FPC_KAPTON_THICKNESS + FPC_ALUMINUM_THICKNESS
        + ANCASIC_GLUE_THICKNESS + ANCASIC_THICKNESS,
    )
    if outward:
        total += EXTERNAL_EPOXY_THICKNESS
    bottom = -total / 2.0
    baseplate_top = bottom + (EXTERNAL_EPOXY_THICKNESS if outward else 0.0)
    baseplate_top += BASEPLATE_THICKNESS
    corrugation_reference = (
        bottom
        if outward
        else bottom + BASEPLATE_THICKNESS + EXTERNAL_EPOXY_THICKNESS
    )
    film = (baseplate_top, baseplate_top + FILM_THICKNESS)
    kapton = (baseplate_top, baseplate_top + FPC_KAPTON_THICKNESS)
    aluminum = (
        kapton[1],
        kapton[1] + FPC_ALUMINUM_THICKNESS,
    )
    return corrugation_reference, film, kapton, aluminum


def global_z_interval(module_center, normal_sign, local_interval):
    values = [module_center + normal_sign * value for value in local_interval]
    return min(values), max(values)


def build_boxes(package):
    with (package / "catalog.csv").open(newline="") as stream:
        catalog = {row["type_id"]: row for row in csv.DictReader(stream)}

    film_boxes = []
    fpc_boxes = []
    for path in sorted(package.glob("*_modules.csv")):
        disk = path.name[: -len("_modules.csv")]
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for source in rows:
            spec = catalog[source["type_id"]]
            type_id = source["type_id"]
            facing = spec["facing"]
            module_width = float(spec["module_height_mm"])
            module_length = float(spec["module_length_mm"])
            rsu_count = int(spec["rsu_count"])
            chain_length = rsu_count * RSU_LENGTH
            lec_length = float(spec["LEC_mm"])
            rec_length = float(spec["REC_mm"])
            lec_span = float(spec["BLEC_mm"])
            rec_span = float(spec["BREC_mm"])

            # The source chirality is reflected together with engineer x.
            handedness = "left" if spec["chirality"] == "RIGHT" else "right"
            lec_after_rsu = handedness == "right"
            end_span = rec_span if handedness == "right" else lec_span
            sensor_x = -module_length / 2.0 + end_span + chain_length / 2.0
            boundary_local = sensor_x + (
                chain_length / 2.0 if handedness == "right" else -chain_length / 2.0
            )

            source_x = float(source["x_origin_mm"])
            epic_x_origin = -source_x
            module_x = epic_x_origin - (
                abs(boundary_local) if epic_x_origin > 0.0 else -abs(boundary_local)
            )
            local_x_sign = (epic_x_origin - module_x) / boundary_local

            z_surface = float(source["z_baseplate_mm"])
            z_sensor = float(source["z_sensor_mm"])
            normal_sign = 1.0 if z_sensor > z_surface else -1.0
            local_y_sign = local_x_sign * normal_sign
            corrugation_reference, film_z, kapton_z, aluminum_z = local_z_geometry(facing)
            module_z = z_surface - normal_sign * corrugation_reference
            applied_z_offset = (
                INWARD_Z_CLEARANCE + OPPOSITE_FACING_Z_OFFSET
                if facing == "OUTWARD"
                else -INWARD_Z_CLEARANCE
            )
            module_z += normal_sign * applied_z_offset

            common = {
                "disk": disk,
                "row": int(source["row_index"]),
                "module": int(source["mod_index"]),
                "type_id": type_id,
                "facing": facing,
            }
            module_y = float(source["y_origin_mm"])

            lec_local_x = sensor_x + (
                chain_length / 2.0 + lec_length / 2.0
                if lec_after_rsu
                else -chain_length / 2.0 - lec_length / 2.0
            )
            rec_local_x = sensor_x + (
                -chain_length / 2.0 - rec_length / 2.0
                if lec_after_rsu
                else chain_length / 2.0 + rec_length / 2.0
            )
            for name, local_x, box_x in (
                ("lec_glue", lec_local_x, lec_length),
                ("rec_glue", rec_local_x, rec_length),
            ):
                film_boxes.append(
                    Box(
                        name=name,
                        x=transformed_interval(module_x, local_x_sign, local_x, box_x),
                        y=interval(module_y, 19.564),
                        z=global_z_interval(module_z, normal_sign, film_z),
                        **common,
                    )
                )

            for side, after_rsu, side_span, endcap_length, target_x, target_y, top_offset in (
                ("left", lec_after_rsu, lec_span, lec_length,
                 LEFT_FPC_X, LEFT_FPC_Y, LEFT_FPC_TOP_OFFSET),
                ("right", not lec_after_rsu, rec_span, rec_length,
                 RIGHT_FPC_X, RIGHT_FPC_Y, RIGHT_FPC_TOP_OFFSET),
            ):
                box_x = min(target_x, max(0.0, side_span - endcap_length - FPC_GAP))
                local_x = (
                    chain_length / 2.0 + endcap_length + FPC_GAP + box_x / 2.0
                    if after_rsu
                    else -chain_length / 2.0 - endcap_length - FPC_GAP - box_x / 2.0
                )
                local_x += sensor_x
                box_y, local_y = fpc_y_geometry(module_width, target_y, top_offset)
                for material, local_z in (
                    ("kapton", kapton_z),
                    ("aluminum", aluminum_z),
                ):
                    fpc_boxes.append(
                        Box(
                            name=f"{side}_fpc_{material}",
                            x=transformed_interval(module_x, local_x_sign, local_x, box_x),
                            y=transformed_interval(module_y, local_y_sign, local_y, box_y),
                            z=global_z_interval(module_z, normal_sign, local_z),
                            **common,
                        )
                    )
    return fpc_boxes, film_boxes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("package", nargs="?", type=Path)
    args = parser.parse_args()
    assert_xml_dimensions()
    package = args.package.resolve() if args.package else active_package()

    fpc_boxes, film_boxes = build_boxes(package)
    overlaps = []
    films_by_disk = {}
    for disk, boxes in itertools.groupby(sorted(film_boxes, key=lambda box: box.disk), key=lambda box: box.disk):
        films_by_disk[disk] = list(boxes)
    for fpc in fpc_boxes:
        for film in films_by_disk[fpc.disk]:
            if source_key(fpc) == source_key(film):
                continue
            dimensions = tuple(overlap(getattr(fpc, axis), getattr(film, axis)) for axis in "xyz")
            if all(value > 1.0e-6 for value in dimensions):
                overlaps.append((fpc, film, dimensions))

    print(f"package={package.relative_to(ROOT)}")
    print(f"fpc_boxes={len(fpc_boxes)} film_end_boxes={len(film_boxes)}")
    print(f"fpc/end-film overlaps={len(overlaps)}")
    physical_pairs = {
        tuple(sorted((source_key(fpc), source_key(film))))
        for fpc, film, _ in overlaps
    }
    print(f"unique module pairs={len(physical_pairs)}")
    for fpc, film, dimensions in overlaps:
        print(
            f"{fpc.disk}: row/module {fpc.row}/{fpc.module} {fpc.type_id} {fpc.name} vs "
            f"{film.row}/{film.module} {film.type_id} {film.name}; "
            f"overlap_xyz_mm={dimensions[0]:.6f},{dimensions[1]:.6f},{dimensions[2]:.6f}"
        )


if __name__ == "__main__":
    main()
