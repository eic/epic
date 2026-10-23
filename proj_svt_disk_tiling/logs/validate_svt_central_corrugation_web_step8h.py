#!/usr/bin/env python3
"""Validate placement eligibility of every CSV-defined central support web."""

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CPP = ROOT / "src/SiEndcapModuleTracker_geo.cpp"
LATTICE = (
    ROOT
    / "compact/tracking/silicon_disks_tiling/20261001_ES"
    / "SVT_endcap_corrugation_rows_33_36.csv"
)
THICKNESS_MM = 0.150
HEIGHT_MM = 6.000
WEB_Y_SPAN_MM = HEIGHT_MM - THICKNESS_MM  # theta = 45 degrees


def flat_center_z(row):
    surface_z = float(row["surface_offset_mm"])
    return surface_z - math.copysign(THICKNESS_MM / 2.0, surface_z)


def main():
    rows_by_disk = defaultdict(list)
    with LATTICE.open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows_by_disk[row["disk"]].append(row)

    crossing_pairs = {}
    centered_flats = {}
    for disk, rows in rows_by_disk.items():
        rows.sort(key=lambda row: float(row["row_y_mm"]))
        centered_flats[disk] = [
            int(row["row_index"])
            for row in rows
            if abs(float(row["row_y_mm"])) < 1.0e-12
        ]
        for lower, upper in zip(rows, rows[1:]):
            lower_edge = (
                float(lower["row_y_mm"]) + float(lower["flat_width_mm"]) / 2.0
            )
            upper_edge = (
                float(upper["row_y_mm"]) - float(upper["flat_width_mm"]) / 2.0
            )
            web_center_y = (lower_edge + upper_edge) / 2.0
            web_z_span = flat_center_z(upper) - flat_center_z(lower)
            web_angle = math.atan2(web_z_span, WEB_Y_SPAN_MM)
            web_length = math.hypot(WEB_Y_SPAN_MM, web_z_span)
            projected_half_y = (
                abs(math.cos(web_angle)) * web_length / 2.0
                + abs(math.sin(web_angle)) * THICKNESS_MM / 2.0
            )
            if web_center_y - projected_half_y < 0.0 < web_center_y + projected_half_y:
                assert disk not in crossing_pairs
                crossing_pairs[disk] = (
                    int(lower["row_index"]),
                    int(upper["row_index"]),
                    web_center_y,
                    projected_half_y,
                )

    assert set(crossing_pairs) == {
        "ED0",
        "ED1",
        "ED2",
        "ED3",
        "ED4",
        "HD0",
        "HD1",
        "HD2",
    }
    assert centered_flats["HD3b"] == [23]
    assert centered_flats["HD4"] == [23]
    assert all(not centered_flats[disk] for disk in crossing_pairs)

    cpp = CPP.read_text()
    assert "Sloped pieces are not split" not in cpp
    assert "skip them instead of crossing y=0" not in cpp
    assert "const bool is_axis_aligned" in cpp
    assert "const bool has_direct_y_span" in cpp
    assert (
        "y_max > 0.0 && is_axis_aligned && has_direct_y_span" in cpp
    )

    for disk, (lower, upper, center, half_span) in sorted(crossing_pairs.items()):
        print(
            f"{disk}: central web rows {lower}--{upper}, "
            f"y_center={center:+.3f} mm, projected y=[{center-half_span:+.3f},"
            f"{center+half_span:+.3f}] mm"
        )
    print("HD3b/HD4: row 23 flat is centered at y=0; adjacent webs do not cross y=0")
    print("PASS: all eight CSV-defined central webs remain eligible for placement")


if __name__ == "__main__":
    main()
