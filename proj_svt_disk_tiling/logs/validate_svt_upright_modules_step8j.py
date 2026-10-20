#!/usr/bin/env python3
"""Step 8j: check negative-Y electronics for both tiling scenarios."""

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "compact/tracking/silicon_disks_tiling/20261001_ES"
CPP = (ROOT / "src/SiEndcapModuleTracker_geo.cpp").read_text()

assert "row.tiling_coordinates && row.y_origin <= 0.0" in CPP
assert 'cache_key += "_upright_electronics_y"' in CPP
assert "electronics_y_sign * (component.y_offset + pos_y)" in CPP
assert CPP.count("electronics_y_sign * (component.y_offset + pos_y)") == 2

counts = Counter()
for scenario in ("rsu_opt", "all_6rsu"):
    for path in sorted((PACKAGE / scenario).glob("*_modules.csv")):
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                y_origin = float(row["y_origin_mm"])
                x_origin = float(row["x_origin_mm"])
                normal_z = (
                    1 if float(row["z_sensor_mm"]) > float(row["z_baseplate_mm"]) else -1
                )
                # The engineer-to-EPIC x reflection exchanges placed handedness.
                placed_right = row["type_id"].endswith("1")
                epic_x_sign = 1 if -x_origin > 0 else -1
                local_boundary_sign = 1 if placed_right else -1
                local_x_global = epic_x_sign * local_boundary_sign
                local_y_global = local_x_global * normal_z
                electronics_y_sign = -1 if y_origin <= 0 else 1
                assert local_y_global * electronics_y_sign == 1, (
                    scenario,
                    path.name,
                    row["row_index"],
                    row["mod_index"],
                )
                counts[(scenario, "upright" if electronics_y_sign < 0 else "standard")] += 1

for (scenario, variant), count in sorted(counts.items()):
    print(f"{scenario}: {variant} module prototypes used by {count} placements")
print(f"PASS: {sum(counts.values())} placements keep FPC/AncASIC +y toward global +Y")
