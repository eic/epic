#!/usr/bin/env python3
"""Audit supplied SVT disk metadata against the current XML boundaries.

This script is intentionally read-only.  The current XML values below are the
resolved numerical values from compact/tracking/definitions_craterlake.xml;
keeping them explicit makes the comparison easy to audit without requiring a
DD4hep runtime environment.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "compact" / "tracking" / "silicon_disks"
SCENARIOS = ("all_6rsu", "rsu_opt")
DISKS = ("ED4", "ED3", "ED2", "ED1", "ED0", "HD0", "HD1", "HD2", "HD3b", "HD4")

# Resolved from definitions_craterlake.xml.  Opening radii include the XML's
# Beampipe_bakeout_buffer = 5 mm.  Coordinates are the XML layer-local values.
CURRENT_XML = {
    "ED4": (-1020.0, 421.4, ((0.0, 0.0, 38.5), (-26.130, 0.0, 19.285))),
    "ED3": (-850.0, 421.4, ((0.0, 0.0, 38.5), (-21.880, 0.0, 19.285))),
    "ED2": (-650.0, 421.4, ((0.0, 0.0, 36.757), (0.0, 0.0, 36.757))),
    "ED1": (-450.0, 415.0, ((0.0, 0.0, 36.757), (0.0, 0.0, 36.757))),
    "ED0": (-250.0, 240.0, ((0.0, 0.0, 36.757), (0.0, 0.0, 36.757))),
    "HD0": (250.0, 240.0, ((0.0, 0.0, 36.757), (0.0, 0.0, 36.757))),
    "HD1": (450.0, 415.0, ((0.0, 0.0, 36.757), (0.0, 0.0, 36.757))),
    "HD2": (700.0, 421.4, ((0.0, 0.0, 36.757), (-1.440, 0.0, 40.488))),
    "HD3b": (950.0, 421.4, ((0.0, 0.0, 36.757), (-7.066, 0.0, 46.111))),
    "HD4": (1200.0, 421.4, ((0.0, 0.0, 36.757), (-12.692, 0.0, 51.734))),
}


def read_metadata(path):
    values = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        try:
            values[key] = float(value)
        except ValueError:
            continue
    return values


def supplied_tuple(values):
    assert int(values["n_opening_primitives"]) == 2
    openings = []
    for index in range(2):
        a = values[f"opening_primitive_{index}_a_mm"]
        b = values[f"opening_primitive_{index}_b_mm"]
        assert a == b, "Step 6a currently expects the supplied primitives to be circles"
        openings.append(
            (
                values[f"opening_primitive_{index}_cx_mm"],
                values[f"opening_primitive_{index}_cy_mm"],
                a,
            )
        )
    return values["z_center_mm"], values["outer_radius_mm"], tuple(openings)


metadata = {
    scenario: {
        disk: supplied_tuple(read_metadata(PACKAGE_ROOT / scenario / f"{disk}_metadata.txt"))
        for disk in DISKS
    }
    for scenario in SCENARIOS
}

assert metadata["all_6rsu"] == metadata["rsu_opt"], (
    "The two scenarios do not have identical disk-boundary metadata"
)

print("PASS: all_6rsu and rsu_opt disk-boundary metadata are identical")
print("All dimensions are mm. XML opening radii include its 5 mm bakeout buffer.")
print(
    "disk   z_sup   z_xml      dz   r_sup   r_xml      dr   "
    "supplied openings (cx,cy,r)              current XML openings (cx,cy,r)"
)
for disk in DISKS:
    z_sup, r_sup, openings_sup = metadata["rsu_opt"][disk]
    z_xml, r_xml, openings_xml = CURRENT_XML[disk]
    print(
        f"{disk:4s} {z_sup:7.1f} {z_xml:7.1f} {z_sup-z_xml:7.1f} "
        f"{r_sup:7.1f} {r_xml:7.1f} {r_sup-r_xml:7.1f}   "
        f"{str(openings_sup):42s} {openings_xml}"
    )

print("PASS: Step 6a boundary audit completed without modifying source data")
