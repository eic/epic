#!/usr/bin/env python3
"""Validate the Step 8i no-gap sharp-junction corrugation table."""

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
GENERATOR = ROOT / "proj_svt_disk_tiling/logs/generate_svt_corrugation_lattice_step7b.py"
EXPECTED_WIDTH_MM = {"A": 12.15, "B": 10.65}
CARBON_THICKNESS_MM = 0.150
CORRUGATION_HEIGHT_MM = 6.000
WEB_ANGLE_DEG = 45.000
WEB_Y_SPAN_MM = (CORRUGATION_HEIGHT_MM - CARBON_THICKNESS_MM) / math.tan(
    math.radians(WEB_ANGLE_DEG)
)
DIMENSION_TOLERANCE_MM = 1.0e-9


def main():
    by_disk = defaultdict(list)
    with LATTICE.open(newline="") as stream:
        for row in csv.DictReader(stream):
            width = float(row["flat_width_mm"])
            assert abs(width - EXPECTED_WIDTH_MM[row["family"]]) < 1.0e-12
            by_disk[row["disk"]].append(row)

    flat_count = 0
    web_count = 0
    negative_roundoff_count = 0
    minimum_slack_mm = 0.0
    for disk, rows in sorted(by_disk.items()):
        rows.sort(key=lambda row: float(row["row_y_mm"]))
        gaps = []
        for lower, upper in zip(rows, rows[1:]):
            # Match the C++ arithmetic order exactly. Mathematically identical
            # joins can acquire a tiny negative slack from binary subtraction.
            lower_edge = float(lower["row_y_mm"]) + float(lower["flat_width_mm"]) / 2.0
            upper_edge = float(upper["row_y_mm"]) - float(upper["flat_width_mm"]) / 2.0
            edge_separation = upper_edge - lower_edge
            slack_mm = edge_separation - WEB_Y_SPAN_MM
            minimum_slack_mm = min(minimum_slack_mm, slack_mm)
            negative_roundoff_count += slack_mm < 0.0
            assert slack_mm + DIMENSION_TOLERANCE_MM >= 0.0
            gap_per_end = (edge_separation - WEB_Y_SPAN_MM) / 2.0
            assert abs(gap_per_end) < 1.0e-9
            gaps.append(gap_per_end)
        flat_count += len(rows)
        web_count += len(rows) - 1
        print(
            f"{disk}: flats={len(rows)} webs={len(rows)-1} "
            f"widths={sorted({float(row['flat_width_mm']) for row in rows})} "
            f"flat_web_gap={max(abs(gap) for gap in gaps):.3f} mm"
        )

    assert flat_count == 446
    assert web_count == 436
    assert negative_roundoff_count > 0
    generator = GENERATOR.read_text()
    assert 'FLAT_WIDTH_MM = {"A": 12.15, "B": 10.65}' in generator
    cpp = CPP.read_text()
    assert "const double dimension_tolerance = 1.0e-9 * mm;" in cpp
    assert "available_y_gap + dimension_tolerance < web_y_span" in cpp
    print(f"totals: flats={flat_count} webs={web_count}")
    print(
        f"floating-point audit: negative_slacks={negative_roundoff_count} "
        f"minimum={minimum_slack_mm:.3e} mm tolerance={DIMENSION_TOLERANCE_MM:.1e} mm"
    )
    print("PASS: every straight flat meets its adjacent web with zero y gap")


if __name__ == "__main__":
    main()
