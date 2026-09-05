# DIRC Simplified Bar-Sensitive Implementation

## Summary

Successfully refactored the DIRC detector from a complex 4-surface design (entry/mid/exit planes) to a simpler, more robust 2-surface design:
- **surface=0**: MCP optical photon hits (existing)
- **surface=1**: Quartz bar charged particle hits (new, per-step)

## Changes Made

### 1. `src/DIRC_geo.cpp` (~55 lines removed, 2 lines added)

**Removed:**
- `entry_plane_vol`, `mid_plane_vol`, `exit_plane_vol` volume definitions
- All plane placement code (attempted 40× placements inside bar loop)
- Related printout debug statements

**Added:**
- `bar_vol.setSensitiveDetector(sens)` — marks bars as sensitive detectors
- `bar_vol.setLimitSet(desc.limitSet("dirc_bar_limits"))` — 0.5mm step limit for fine sampling
- `.addPhysVolID("surface", 1)` to bar placements in loop

**Result:** Cleaner, more maintainable geometry with single responsibility per volume.

### 2. `compact/pid/dirc.xml`

**Replaced:**
- `dirc_plane_limits` (0.005mm step) → `dirc_bar_limits` (0.5mm step)

**Simplified readout:**
- Removed 4-surface MultiSegmentation (entry/mid/exit)
- Now 2-surface MultiSegmentation: mcp (surface=0, 3mm grid) + bar (surface=1, 0.5mm grid)
- Updated ID field: `surface:2` → `surface:1` (1 bit instead of 2 bits)

**Updated documentation:**
- Clarified that bars now generate per-Geant4-step hits
- Added note about using `CollectSingleDeposits=True` in npsim

## Expected Behavior

When npsim runs with the `CollectSingleDeposits=True` option:

```bash
npsim --compactFile epic_dirc_only.xml -G -N10 \
  --gun.particle e- --gun.distribution uniform \
  --action.tracker "Geant4TrackerWeightedAction:CollectSingleDeposits=True" \
  --outputFile output.edm4hep.root
```

**Expected hits:**
- **surface=0** (MCP): Many optical photon hits (existing behavior, ~1 per photon)
- **surface=1** (bar): Many charged particle hits (10–100 per electron traversal, depending on bar geometry and incident angle)

Step size is controlled by `dirc_bar_limits` (0.5mm), so fine-grained position sampling along track.

## Verification

Analysis script created: `scripts/analyze_dirc_bar_hits.py`

```bash
python3 scripts/analyze_dirc_bar_hits.py output.edm4hep.root
```

Outputs:
- Total hits by surface
- Hits per event distribution
- Success/warning indicators

## Build Status

✅ **Build successful** (all 80+ detector configurations compiled)
✅ **Code committed** (commit 4df7e7e4d)
✅ **Ready for simulation**

## Next Steps

1. Run npsim with `CollectSingleDeposits=True` flag
2. Use analysis script to verify surface=1 hits appear
3. Optionally integrate into full detector chain (currently dirc_only configuration)

## Environment Notes

Current environment issue: npsim.py requires PyROOT, but free-threaded Python doesn't have it.
Workaround: Run npsim in an environment with ROOT available (e.g., eic-shell with proper setup).

## Advantages Over Previous Design

| Aspect | Entry/Mid/Exit Planes | Bar-Sensitive |
|--------|-----|---|
| Code complexity | High (3 volumes × 40 bars) | Low (1 modification) |
| Geometry robustness | Fragile (loop-placement bug) | Robust (standard pattern) |
| Optical purity | Risk of interfaces | None (single material) |
| Position sampling | 3 discrete points | Continuous (every 0.5mm) |
| Maintenance burden | High (3 surfaces to manage) | Low (2 surfaces) |
