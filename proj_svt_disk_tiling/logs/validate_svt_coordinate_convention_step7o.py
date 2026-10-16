#!/usr/bin/env python3
"""Compare provider opening signs with physical XML openings in one frame.

Page 2 defines detector +z as hadron-going, +y upward, and +x toward the
inside of the ring.  The delivered module/metadata coordinates are already in
that detector frame and must not be reflected to fit another subsystem.
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES/rsu_opt"
DEFINITIONS = ROOT / "compact/tracking/definitions_craterlake.xml"

DISKS = {
    "ED3": ("TrackerEndcapNDisk3_hadron_center_x", -1.0),
    "ED4": ("TrackerEndcapNDisk4_hadron_center_x", -1.0),
    "HD2": ("TrackerEndcapPDisk2_lepton_center_x", +1.0),
    "HD3b": ("TrackerEndcapPDisk3_lepton_center_x", +1.0),
    "HD4": ("TrackerEndcapPDisk4_lepton_center_x", +1.0),
}


def metadata_value(disk, key):
    text = (PACKAGE / f"{disk}_metadata.txt").read_text()
    match = re.search(rf"^{key}\s*=\s*([-+0-9.]+)", text, re.MULTILINE)
    assert match, (disk, key)
    return float(match.group(1))


def xml_mm(name):
    text = DEFINITIONS.read_text()
    match = re.search(rf'<constant name="{name}"\s+value="([-+0-9.]+)\*mm', text)
    assert match, name
    return float(match.group(1))


print("disk,provider_global_x_mm,expected_layer_local_x_mm,current_xml_local_x_mm,status")
failures = 0
for disk, (constant, layer_axis_sign) in DISKS.items():
    primitive = 1 if disk.startswith("ED") else 0
    provider_global = metadata_value(disk, f"opening_primitive_{primitive}_cx_mm")
    expected_local = layer_axis_sign * provider_global
    current_local = xml_mm(constant)
    status = "MATCH" if abs(expected_local - current_local) < 1.0e-6 else "OPPOSITE"
    failures += status == "OPPOSITE"
    print(f"{disk},{provider_global:+.4f},{expected_local:+.4f},{current_local:+.4f},{status}")

print("\nConvention calculation:")
print("  positive-z layers: x_local = x_detector")
print("  negative-z layers: x_local = -x_detector (existing Ry(pi) placement)")
print("  therefore every displaced provider opening is at positive layer-local x")
print(f"RESULT: {failures} physical XML centers are opposite the provider openings")
print("NOTE: this diagnoses a cross-subsystem disagreement; it does not authorize")
print("      reflecting the delivered modules or moving a physical beampipe hole")
assert failures == 5, "update this audit if either subsystem convention changes"
