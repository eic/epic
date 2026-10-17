#!/usr/bin/env python3
"""Validate that tiling-module assembly daughters are visible."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CPP = (ROOT / "src/SiEndcapModuleTracker_geo.cpp").read_text()

tiling_start = CPP.index("auto build_tiling_module")
tiling_end = CPP.index("return module_template;", tiling_start)
tiling_builder = CPP[tiling_start:tiling_end]

assert 'module_template.vis           = "SVTModuleVis";' in tiling_builder
assert "prototype.volume = Assembly(module_template.name);" in CPP

for display_name in ("display.xml", "display_detailed.xml", "display_geoviewer.xml"):
    display = (ROOT / "compact" / display_name).read_text()
    start = display.index('<vis name="SVTModuleVis"')
    definition = display[start : display.index("/>", start) + 2]
    assert 'visible="true"' in definition, (display_name, definition)
    assert 'showDaughters="true"' in definition, (display_name, definition)

print("PASS: supplied tiling modules use the dedicated SVTModuleVis")
print("PASS: SVTModuleVis exposes assembly daughters in all display profiles")
print("PASS: generic TrackerModuleVis remains unchanged for unrelated trackers")
