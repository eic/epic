#!/usr/bin/env python3
"""Audit SVT physical openings in the current detector-aligned layer frame.

The fixed EPIC beam-pipe XML is the geometry authority.  The existing
``disk_layout`` helper is used only to evaluate that beam-pipe envelope; no
module placement from ``disk_layout`` enters this calculation.

For each 10 mm SVT layer, evaluate the two beam-pipe outer-wall circles at the
layer face farthest from the interaction point.  Compare those unbuffered
circles with the physical opening constants currently in
``definitions_craterlake.xml``.  The common 5 mm bakeout allowance is reported
separately and is not folded into the beam-pipe dimensions below.
"""

import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
BEAMPIPE_HELPER = ROOT / "proj_svt_disk_tiling/disk_layout/scripts/beampipe_envelope_circles.py"

# disk, signed side, layer-center constant, physical-opening constant prefix
DISKS = (
    ("ED0", -1.0, "InnerTrackerEndcapN_zmin", "InnerTrackerEndcapN"),
    ("ED1", -1.0, "TrackerEndcapNDisk1_zmin", "TrackerEndcapNDisk1"),
    ("ED2", -1.0, "TrackerEndcapNDisk2_zmin", "TrackerEndcapNDisk2"),
    ("ED3", -1.0, "TrackerEndcapNDisk3_zmin", "TrackerEndcapNDisk3"),
    ("ED4", -1.0, "TrackerEndcapNDisk4_zmin", "TrackerEndcapNDisk4"),
    ("HD0", +1.0, "InnerTrackerEndcapP_zmin", "InnerTrackerEndcapP"),
    ("HD1", +1.0, "TrackerEndcapPDisk1_zmin", "TrackerEndcapPDisk1"),
    ("HD2", +1.0, "TrackerEndcapPDisk2_zmin", "TrackerEndcapPDisk2"),
    ("HD3b", +1.0, "TrackerEndcapPDisk3_zmin", "TrackerEndcapPDisk3"),
    ("HD4", +1.0, "TrackerEndcapPDisk4_zmin", "TrackerEndcapPDisk4"),
)


def load_helper(path):
    spec = importlib.util.spec_from_file_location("beampipe_envelope_circles", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


helper = load_helper(BEAMPIPE_HELPER)
definitions_text = DEFINITIONS.read_text()
cpp_text = CPP.read_text()
raw_constants = helper.load_constants(DEFINITIONS)
constants = helper.ConstantResolver(raw_constants)
beam_pipe = helper.load_geometry(ROOT.parent)

# These guards ensure the expected values below are in the same detector-
# aligned x frame used by the current physical XML openings.
assert "RotationZYX(0.0, -M_PI, 0.0)" not in cpp_text
assert "Both layer frames are detector-aligned" in cpp_text
assert "Negative disk openings use detector-aligned layer coordinates" in definitions_text

layer_thickness_mm = constants["SiTrackerEndcapLayer_thickness"]
bakeout_buffer_mm = constants["Beampipe_bakeout_buffer"]
assert math.isclose(layer_thickness_mm, 10.0, abs_tol=1.0e-12)
assert math.isclose(bakeout_buffer_mm, 5.0, abs_tol=1.0e-12)

print("Physical-opening audit in detector-aligned EPIC coordinates")
print("beam-pipe source: compact/central_beampipe.xml plus compact/definitions.xml")
print("layer positions/openings: compact/tracking/definitions_craterlake.xml")
print("layer thickness: {:.3f} mm".format(layer_thickness_mm))
print("bakeout allowance (excluded from table radii): {:.3f} mm".format(bakeout_buffer_mm))
print()
print(
    "disk,outer_face_z_mm,circle,expected_cx_mm,expected_r_mm,"
    "current_cx_mm,current_unbuffered_r_mm,delta_cx_mm,delta_r_mm"
)

changed_fields = 0
max_center_delta = 0.0
max_y_delta = 0.0
max_radius_delta = 0.0
results = []

for disk, side, z_constant, prefix in DISKS:
    center_magnitude_mm = constants[z_constant]
    center_z_mm = side * center_magnitude_mm
    outer_face_z_mm = center_z_mm + side * layer_thickness_mm / 2.0
    expected_circles = helper.circles_at_z(outer_face_z_mm, beam_pipe)
    assert len(expected_circles) == 2

    disk_results = []
    for xml_label, expected in zip(("lepton", "hadron"), expected_circles):
        center_name = "{}_{}_center_x".format(prefix, xml_label)
        center_y_name = "{}_{}_center_y".format(prefix, xml_label)
        radius_name = "{}_{}_radius".format(prefix, xml_label)
        current_center_mm = constants[center_name]
        current_y_mm = constants[center_y_name]
        current_buffered_radius_mm = constants[radius_name]
        current_radius_mm = current_buffered_radius_mm - bakeout_buffer_mm
        delta_center_mm = expected.x_mm - current_center_mm
        delta_y_mm = expected.y_mm - current_y_mm
        delta_radius_mm = expected.r_mm - current_radius_mm

        if abs(delta_center_mm) > 5.0e-7:
            changed_fields += 1
        if abs(delta_y_mm) > 5.0e-7:
            changed_fields += 1
        if abs(delta_radius_mm) > 5.0e-7:
            changed_fields += 1
        max_center_delta = max(max_center_delta, abs(delta_center_mm))
        max_y_delta = max(max_y_delta, abs(delta_y_mm))
        max_radius_delta = max(max_radius_delta, abs(delta_radius_mm))

        print(
            "{},{:.3f},{},{:.6f},{:.6f},{:.6f},{:.6f},{:+.6f},{:+.6f}".format(
                disk,
                outer_face_z_mm,
                xml_label,
                expected.x_mm,
                expected.r_mm,
                current_center_mm,
                current_radius_mm,
                delta_center_mm,
                delta_radius_mm,
            )
        )
        disk_results.append((xml_label, expected.x_mm, expected.r_mm))

    results.append((disk, outer_face_z_mm, disk_results))

result_by_disk = {disk: circles for disk, _, circles in results}
assert result_by_disk["ED3"][1][1] > 0.0
assert result_by_disk["ED4"][1][1] > 0.0
for disk in ("HD2", "HD3b", "HD4"):
    assert result_by_disk[disk][1][1] < 0.0

print()
print("Expected unbuffered physical circles for the XML update")
print("disk | outer face z [mm] | lepton (cx, r) [mm] | hadron (cx, r) [mm]")
print("--- | ---: | ---: | ---:")
for disk, outer_face_z_mm, disk_results in results:
    lepton = disk_results[0]
    hadron = disk_results[1]
    print(
        "{} | {:.3f} | ({:.6f}, {:.6f}) | ({:.6f}, {:.6f})".format(
            disk,
            outer_face_z_mm,
            lepton[1],
            lepton[2],
            hadron[1],
            hadron[2],
        )
    )

print()
print("changed XML x/y/r fields: {} of {}".format(changed_fields, len(DISKS) * 6))
print("maximum absolute center correction: {:.6f} mm".format(max_center_delta))
print("maximum absolute y correction: {:.6f} mm".format(max_y_delta))
print("maximum absolute radius correction: {:.6f} mm".format(max_radius_delta))
print("PASS: all ten disks were evaluated in detector-aligned EPIC coordinates")
print("PASS: ED displaced-circle signs are positive; HD displaced-circle signs are negative")
print("PASS: source tiling metadata and disk_layout module placements were not used")
