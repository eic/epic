# SVT Endcap Corrugated Module Implementation Log

Append substantial project changes here. This file complements git history and is intended for debugging, reproducibility, and future handoff.

Entry template:

```text
## YYYY-MM-DD HH:MM UTC - short title

Files changed:
- path

Intent:
- What this change is trying to accomplish.

Implementation:
- What changed technically.

Validation:
- Commands/scripts run and outcome.

Known limitations / next step:
- What remains unresolved.
```

## 2026-05-25 UTC - Project documentation scaffold

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Start a careful, step-by-step implementation record for the SVT endcap corrugated module geometry update.
- Capture current understanding of the design requirements, implementation phases, validation gates, and logging convention.

Implementation:
- Added an agent-friendly project guide.
- Added this append-only implementation log.

Validation:
- Documentation-only change. No geometry build or overlap checks run.

Known limitations / next step:
- The worktree currently contains pre-existing prototype geometry changes. Decide whether to keep them as a starting point or discard/reset them before implementing the new `EIC_LAS_6RSU_CORR` path.

## 2026-05-25 UTC - Initial design decisions recorded

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Record first implementation choices from user review of the open questions.

Implementation:
- Set first implementation package/support width to `32.00 mm`.
- Decided to ignore carbon fiber bottom-edge cutouts initially.
- Decided corrugated module handedness should be explicit in the CSV, using `left` / `right`, rather than inferred from placement.
- Recorded that 5-RSU corrugated dimensions are unknown; assume the same end dimensions as 6-RSU later, but focus first on 6-RSU.
- Recorded that the current prototype can be retained only if it remains separate and cleanly removable.

Validation:
- Documentation-only change. No geometry build or overlap checks run.

Known limitations / next step:
- Need to choose placeholder DD4hep materials for LEC, REC, and FPC.
- Need to implement the `EIC_LAS_6RSU_CORR` path cleanly, keeping the current prototype separable.

## 2026-05-25 UTC - Component width decision recorded

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Record that the implementation should use the design widths for different module components rather than forcing all layers to one common width.

Implementation:
- Updated the project plan to use a `32.00 mm` carbon fiber/support envelope and a `19.56 mm` RSU/sensor band with `6.22 mm` margins.
- Noted that glue/adhesive width should follow the component it bonds when clear, using simple rectangular approximations at first.

Validation:
- Documentation-only change. No geometry build or overlap checks run.

Known limitations / next step:
- During implementation, decide explicitly whether each glue/adhesive layer follows the full support envelope or the RSU/sensor band.

## 2026-05-25 UTC - Clean baseline preference recorded

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Record the decision to begin implementation from a cleaner geometry/code baseline.

Implementation:
- Updated Phase 0 to prefer discarding or setting aside the current dirty source/XML prototype changes before implementing `EIC_LAS_6RSU_CORR`.
- Clarified that documentation, drawing PDFs, and project notes should not be reset as part of cleanup.

Validation:
- Documentation-only change. No geometry build or overlap checks run.

Known limitations / next step:
- Before cleanup, list the exact dirty files and decide which source/XML changes to discard, which local artifacts to keep, and whether to preserve a patch/reference copy of the current prototype.

## 2026-05-25 23:33 UTC - Phase 0 worktree cleanup

Files changed:
- `tmp/phase0_cleanup_2026-05-25.sh`
- `tmp/file_actions.log`
- `tmp/retained/2026-05-25/phase0_cleanup/prototype_geometry_changes.patch`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Start the corrugated module implementation from a cleaner geometry/code baseline while preserving the current prototype for reference.

Implementation:
- Wrote and ran `tmp/phase0_cleanup_2026-05-25.sh`.
- Saved the dirty tracked prototype/display diff to `tmp/retained/2026-05-25/phase0_cleanup/prototype_geometry_changes.patch`.
- Restored these tracked files to `HEAD`:
  - `compact/display_detailed.xml`
  - `compact/display_geoviewer.xml`
  - `compact/tracking/silicon_disks_modules.xml`
  - `src/SiEndcapModuleTracker_geo.cpp`
- Moved untracked local artifacts into retention:
  - `configurations/craterlake_svt_only.yml`
  - `configurations/craterlake_tracking_only_disk_tiles_new.yml`
  - `setup_local.sh`
  - `install/`
- Logged the moved files in `tmp/file_actions.log`.

Validation:
- `git status --short` now shows no modified tracked source/XML/display files.
- Remaining untracked items are project docs/PDF references and the `tmp/` retention area.
- The retained prototype patch has 199 lines.

Known limitations / next step:
- `setup_local.sh` and `install/` were moved to retention and can be restored from `tmp/retained/2026-05-25/phase0_cleanup/` if needed for validation.
- Next implementation step should begin from restored tracked geometry files and add only the new parameterized `EIC_LAS_6RSU_CORR` path.

## 2026-05-25 UTC - Debug helpers restored

Files changed:
- `setup_local.sh`
- `configurations/craterlake_tracking_only_disk_tiles_new.yml`
- `tmp/file_actions.log`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Restore local helper files that streamline debugging and geometry checks.

Implementation:
- Restored `setup_local.sh` from `tmp/retained/2026-05-25/phase0_cleanup/setup_local.sh`.
- Restored `configurations/craterlake_tracking_only_disk_tiles_new.yml` from `tmp/retained/2026-05-25/phase0_cleanup/configurations/craterlake_tracking_only_disk_tiles_new.yml`.
- Logged both restore actions in `tmp/file_actions.log`.

Validation:
- Documentation/file-restore only. No geometry build or overlap checks run.

Known limitations / next step:
- `install/` remains retained and was not restored.

## 2026-05-25 UTC - Phase 1 XML parameter scaffold

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Add the XML constants needed for the future `EIC_LAS_6RSU_CORR` implementation without changing current module behavior.

Implementation:
- Added RSU dimensions:
  - `SiEndcapRSU_width`
  - `SiEndcapRSU_length`
  - `SiEndcapRSU_periphery_width`
  - `SiEndcapRSU_bias_width`
  - `SiEndcapRSU_backbone_width`
- Added six-RSU corrugated package dimensions:
  - `SiEndcapModule_width_corrugated = 32.0*mm`
  - `SiEndcapModule6RSU_package_length = 152.02*mm`
  - end extensions and sensor-strip margins.
- Kept strip decomposition, CF cutout dimensions, and y-margin values as drawing/documentation notes until they are wired into geometry.
- Added new corrugated-module thickness constants for sensor, adhesive, and CF.
- Recorded placeholder material choices in the project guide:
  - FPC/base flex: `Kapton`
  - optional metal traces/readout strips: `Copper`
  - LEC/REC chip-like placeholders: `Silicon`
  - adhesive/glue: `SVT_Endcap_Glue`

Validation:
- Ran `xmllint --noout compact/tracking/silicon_disks_modules.xml`; it passed.
- No DD4hep build/export was run in this phase because no code consumes the new constants yet.

Known limitations / next step:
- Constants are intentionally unused until the C++ module builder is extended in later phases.
- The exact LEC/REC geometry remains approximate.

## 2026-05-26 UTC - Phase 1 verified after main merge

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Resume after latest main-branch merge and confirm the Phase 1 parameter scaffold is still present and valid.

Implementation:
- Inspected `compact/tracking/silicon_disks_modules.xml`.
- Confirmed the corrugated RSU/module constants are present and remain the only tracked geometry diff for Phase 1.

Validation:
- Ran `xmllint --noout compact/tracking/silicon_disks_modules.xml`; it passed.

Known limitations / next step:
- No C++ code consumes these constants yet. Proceed to Phase 2 when ready.

## 2026-05-26 UTC - Phase 2 component offsets

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Enable future corrugated module components to be placed away from the module center while preserving current module behavior.

Implementation:
- Added `x_offset` and `y_offset` fields to `ComponentTemplate`, both defaulting to `0.0`.
- Changed component placement in `build_module_prototype(...)` from `Position(0, 0, z_position)` to `Position(component.x_offset, component.y_offset, z_position)`.
- Existing aggregate initializers for legacy module components do not specify offsets, so they keep the default zero offsets.

Validation:
- Ran `xmllint --noout compact/tracking/silicon_disks_modules.xml`; it passed.
- Ran `git diff --check`; it passed.
- Attempted `cmake --build build -j4`; it failed because the existing build directory references missing `/opt/local/bin/gmake`.
- Attempted `make -C build -j4`; it failed because the existing build directory references a missing absolute CMake path under `/opt/software/...`.

Known limitations / next step:
- No successful compile check was possible with the current stale build directory/environment.
- Next phase can add `EIC_LAS_6RSU_CORR` using these offsets; before a full build/export, regenerate or restore the local build/install environment.

## 2026-05-26 UTC - Phase 3 initial 6-RSU corrugated module type

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`
- `tmp/validate_phase3_build.sh`

Intent:
- Add the first clean `EIC_LAS_6RSU_CORR` module type alongside the legacy modules, without changing the placement CSV.

Implementation:
- Added `build_corrugated_6rsu(...)` inside `builtin_module_templates(...)`.
- Added `EIC_LAS_6RSU_CORR` to the built-in module map.
- The new module currently uses:
  - package length `SiEndcapModule6RSU_package_length = 152.02 mm`
  - package/support width `SiEndcapModule_width_corrugated = 32.0 mm`
  - full-envelope carbon-fiber support layer
  - full-envelope adhesive layer
  - sensitive silicon strip with length `6 * SiEndcapRSU_length` and width `SiEndcapRSU_width`
  - x offset from the 11 mm left extension and 4.5 mm sensor-left margin
- Deferred localized FPC/readout/end boxes until the builder can represent same-layer side-by-side passive pieces without adding each local box thickness to the module stack.

Validation:
- Ran `xmllint --noout compact/tracking/silicon_disks_modules.xml`; it passed.
- Ran `git diff --check`; it passed.
- Added and ran `tmp/validate_phase3_build.sh`.
- First configure attempt exposed stale build-cache paths.
- Updated the script to run `cmake --fresh`.
- Fresh configure then failed because DD4hep is not available in this shell:
  `Could not find a package configuration file provided by "DD4hep"`.

Known limitations / next step:
- No successful compile/export check yet; run inside the proper EPIC/DD4hep environment.
- The new module is not used by the CSV yet, so detector output should remain unchanged until rows are switched or a test row/config is introduced.
- Passive FPC/readout/end features are intentionally deferred.

## 2026-05-26 UTC - Test CSV for corrugated 6-RSU visualization

Files changed:
- `compact/tracking/SVT_endcap_modules_corrugation_6rsu_corr_test.csv`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Provide a small opt-in placement CSV for visualizing the new `EIC_LAS_6RSU_CORR` prototype without modifying the production placement CSV.

Implementation:
- Added four sparse enabled rows on `InnerTrackerEndcapN_disk1` using `EIC_LAS_6RSU_CORR`.
- Added two disabled legacy comparison rows using `EIC_LAS_6RSU` and `EIC_LAS_5RSU`.
- Left the production XML unchanged; it still points at `compact/tracking/SVT_endcap_modules_corrugation_all.csv`.

Validation:
- Checked the CSV columns with `awk`; all data rows have 8 fields.
- No geometry export was run yet.

Known limitations / next step:
- To visualize, temporarily point the relevant `<module_placements file=...>` entries in `compact/tracking/silicon_disks_modules.xml` at `compact/tracking/SVT_endcap_modules_corrugation_6rsu_corr_test.csv`, or make a separate test compact XML/config.
- The test rows do not include handedness yet; handedness support is a later phase.

## 2026-05-26 UTC - Repaired installed materials-map symlink

Files changed:
- `tmp/fix_materials_map_symlink_2026-05-26.sh`
- `tmp/file_actions.log`
- `install/share/epic/calibrations/543a6557cfe39ac6`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Fix `dd_web_display` failure caused by a self-referential installed calibration symlink:
  `filesystem error: status: Too many levels of symbolic links [calibrations/materials-map.cbor]`.

Implementation:
- Confirmed `install/share/epic/calibrations/materials-map.cbor` pointed to `543a6557cfe39ac6`.
- Confirmed `install/share/epic/calibrations/543a6557cfe39ac6` was a symlink to itself through the absolute install path.
- Wrote and ran `tmp/fix_materials_map_symlink_2026-05-26.sh`.
- Moved the bad self-link to `tmp/retained/2026-05-26/materials_map_symlink_fix/543a6557cfe39ac6.self_symlink`.
- Restored the real 29 MB payload from the Phase 0 retained install tree.

Validation:
- Verified `test -f install/share/epic/calibrations/materials-map.cbor` succeeds.
- Verified `install/share/epic/calibrations/543a6557cfe39ac6` is now a regular file.
- Could not rerun `dd_web_display` in this shell because it is not on `PATH`.

Known limitations / next step:
- Re-run the visualization command in the configured EPIC/DD4hep environment.
- There is another self-referential installed calibration symlink, `9508b7435023068a`, but it was not part of the reported material-map failure and was left unchanged.

## 2026-05-26 UTC - Temporary XML test CSV switch noted

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Preserve the temporary XML change that points `InnerTrackerEndcapN` at the corrugated 6-RSU test CSV for easier visualization.

Implementation:
- User changed the `InnerTrackerEndcapN` `<module_placements>` entry from the production CSV to:
  `compact/tracking/SVT_endcap_modules_corrugation_6rsu_corr_test.csv`.
- The original production CSV line is currently kept inside an XML comment.
- This is intentionally retained for testing convenience.

Validation:
- User ran `dd_web_display`; after repairing the material-map symlink, export proceeded with only an expected warning about no valid modules for `InnerTrackerEndcapP_disk1`.

Known limitations / next step:
- Before finalizing or committing production geometry, remind the user to revert this temporary XML test switch back to `compact/tracking/SVT_endcap_modules_corrugation_all.csv`.

## 2026-05-26 UTC - Phase 4 explicit handedness parsing

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `compact/tracking/SVT_endcap_modules_corrugation_6rsu_corr_test.csv`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Add explicit left/right handedness support to the placement CSV path for corrugated modules.

Implementation:
- Added `handedness` to `ModuleRow`.
- Added `parse_handedness(...)`, accepting only empty, `left`, or `right`.
- Corrugated module rows, identified by module names ending in `_CORR`, now require explicit `left` or `right` handedness.
- Non-corrugated legacy rows remain backward compatible with missing/empty handedness.
- Module prototype cache keys now include handedness when present, preparing for future handed prototype variants.
- Updated the corrugated 6-RSU test CSV to include a `handedness` column with left/right examples.

Validation:
- Ran `awk` field-count check on `compact/tracking/SVT_endcap_modules_corrugation_6rsu_corr_test.csv`; data rows have 9 fields.
- Ran `xmllint --noout compact/tracking/silicon_disks_modules.xml`; it passed.
- Ran `git diff --check`; it passed.
- No DD4hep export was run from this shell.

Known limitations / next step:
- Left/right handedness still does not change the local geometry.
- The temporary XML test switch in `compact/tracking/silicon_disks_modules.xml` is intentionally retained for visualization and should be reverted before finalizing production geometry.
- Before opening a pull request, remove the temporary test CSV and the working markdown files added for this implementation unless any are intentionally promoted to permanent documentation.

## 2026-05-27 UTC - Phase 5 RSU four-region active pattern

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`
- `tmp/validate_phase5_rsu_four_region_2026-05-27.sh`

Intent:
- Reimplement Phase 5 on top of the Phase 4 checkpoint commit so the substantial RSU subdivision change is isolated in a fresh working-tree diff.
- Replace the single sensitive silicon strip in `EIC_LAS_6RSU_CORR` with a vertex-barrel-style approximation where each RSU has four active regions and inactive silicon gaps.

Implementation:
- Extended `ComponentTemplate` with `x_repeat` and `rsu_four_region_pattern`.
- Marked the corrugated 6-RSU silicon component as six repeated RSU pitches using the four-region pattern.
- In `build_module_prototype(...)`, split each RSU pitch into two halves along local x and two halves along local y.
- For each RSU, create four sensitive silicon rectangles.
- Model backbone, bias, and periphery/readout regions as passive silicon boxes using `TrackerServiceVis`.
- Follow the vertex-barrel `upper`/`lower` convention along local y: bias strips are placed at the internal boundary between the two y halves, while periphery/readout strips are placed at the outer edges.
- Dimensions are taken from XML constants:
  `SiEndcapRSU_length`, `SiEndcapRSU_width`, `SiEndcapRSU_backbone_width`,
  `SiEndcapRSU_bias_width`, and `SiEndcapRSU_periphery_width`.

Validation:
- Created and ran `tmp/validate_phase5_rsu_four_region_2026-05-27.sh`.
- `git diff --check` passed.
- `xmllint --noout compact/tracking/silicon_disks_modules.xml` passed.
- CSV field-count check confirmed active test rows have 9 fields with explicit handedness.

Known limitations / next step:
- This has not yet been compiled or exported in the DD4hep environment from this shell.
- Left/right handedness swaps the sensor-band start side, but detailed asymmetric FPC/readout/end features are still not modeled.
- The temporary XML test switch in `compact/tracking/silicon_disks_modules.xml` is intentionally retained for visualization and should be reverted before finalizing production geometry.

## 2026-05-27 UTC - Fixed handed corrugated prototype selection

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Fix the issue where `left` and `right` handed corrugated rows had separate cache keys but still built identical module prototypes.

Implementation:
- Added `corrugated_6rsu_sensor_x_offset(...)` to calculate the sensor-band local x offset from the package length, RSU chain length, and handed end-extension/sensor-margin constants.
- Added `with_corrugated_handedness(...)`, which copies the selected module template and updates sensitive component offsets before prototype construction.
- `left` now uses `SiEndcapModule6RSU_left_extension + SiEndcapModule6RSU_sensor_left_margin`.
- `right` now uses `SiEndcapModule6RSU_right_extension + SiEndcapModule6RSU_sensor_right_margin`.

Validation:
- Re-ran `tmp/validate_phase5_rsu_four_region_2026-05-27.sh`; it passed.
- Ran `git diff --check`; it passed.
- Ran `xmllint --noout compact/tracking/silicon_disks_modules.xml`; it passed.
- Sanity-calculated the current sensor-band offsets:
  `left = +4.488 mm`, `right = -4.512 mm`.

Known limitations / next step:
- This fixes the modeled sensor-band handedness, but detailed asymmetric FPC/readout/end pieces are still deferred.

## 2026-05-28 UTC - SVT-specific visualization attributes

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `compact/display.xml`
- `compact/display_detailed.xml`
- `compact/display_geoviewer.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Make the corrugated RSU module easier to inspect by using SVT-specific visualization attributes and separating glue color from inactive RSU silicon.

Implementation:
- Changed module support components from `TrackerSupportVis` to `SVTSupportVis`.
- Changed active silicon components from `TrackerLayerVis` to `SVTSensorVis`.
- Changed inactive RSU backbone/bias/periphery silicon boxes from `TrackerServiceVis` to `SVTReadoutVis`.
- Added `SVTGlueVis` and assigned glue/adhesive components to it.
- Added the SVT visualization definitions to `display_detailed.xml` and `display_geoviewer.xml` in addition to the existing main `display.xml` definitions.

Validation:
- Ran `git diff --check`; it passed.
- Ran `xmllint --noout` on `compact/display.xml`, `compact/display_detailed.xml`, `compact/display_geoviewer.xml`, and `compact/tracking/silicon_disks_modules.xml`; it passed.

Known limitations / next step:
- This has not yet been visually checked with `dd_web_display` in the configured DD4hep environment.

## 2026-05-28 UTC - Updated remaining implementation roadmap

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Update the project documentation with the agreed remaining workflow after successful visual inspection and overlap checks.

Roadmap at the time, superseded by the later Phase 6/7/8 reorder entries:
- Migrate the real placement CSV to `EIC_LAS_6RSU_CORR` with explicit `left/right` handedness.
- Rebuild/export, carefully inspect handedness, and rerun both overlap checks.
- Run a more thorough `eicrecon`/ACTS validation to catch reconstruction-level issues.
- Add simple, parameterized LEC/REC details.
- Repeat geometry, overlap, handedness, and `eicrecon` validation after LEC/REC details are added.
- Cleanup temporary XML switches, the test CSV, and working markdown files before PR unless any are intentionally promoted to permanent documentation.

Validation:
- Documentation-only update; no geometry or XML validation required.

## 2026-06-05 19:13 UTC - RSU 12-tile detail decision

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Record the next RSU-detail implementation choice before changing geometry code.

Implementation plan:
- Use `1RSU_top_half.png` as the current RSU architecture reference.
- Keep the existing six-RSU corrugated module and two-by-two RSU-half structure.
- Refine each RSU x-half from one active rectangle into three active tile columns.
- Add passive power-switch strips using the existing `SiEndcapRSU_powerswitch_width` constant.
- Preserve passive backbone, biasing, and readout/periphery silicon regions.
- Keep only active tile regions sensitive.

Expected geometry change:
- Each RSU changes from `2 x-halves * 2 y-halves = 4` sensitive regions to
  `2 x-halves * 3 tile columns * 2 y-halves = 12` sensitive regions.
- A 6-RSU module therefore changes from 24 sensitive regions to 72 sensitive regions.

Validation to run:
- Check XML parsing and whitespace.
- Check C++ diff formatting.
- Increase the `sensor` ID bitfield if needed for 72 sensitive regions per module.
- Rebuild/export and rerun overlap checks in the configured DD4hep environment before relying on the new geometry.

## 2026-06-05 19:13 UTC - Implement RSU 12-tile approximation

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `compact/tracking/silicon_disks_modules.xml`
- `RSU_DRAWING_SUMMARY.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`
- `tmp/validate_rsu_twelve_tile_2026-06-05.sh`

Intent:
- Refine the corrugated RSU sensitive-region model to better match `1RSU_top_half.png`.

Implementation:
- Renamed the RSU component flag from `rsu_four_region_pattern` to `rsu_twelve_tile_pattern`.
- Split each local-x half-RSU into three active tile columns.
- Added passive silicon power-switch strips after each tile column using `SiEndcapRSU_powerswitch_width`.
- Preserved passive silicon backbone, bias, and periphery/readout strips.
- Kept only active tile rectangles sensitive.
- Expanded `TrackerEndcapHits` from `sensor:5` to `sensor:7` to cover 72 sensitive regions in each 6-RSU corrugated module.
- Updated project/drawing documentation to describe the 12-tile approximation and current dimensions.

Validation:
- Created and ran `tmp/validate_rsu_twelve_tile_2026-06-05.sh`.
- `git diff --check` passed.
- `xmllint --noout compact/tracking/silicon_disks_modules.xml` passed.
- Confirmed the old `rsu_four_region_pattern` name is absent.
- Confirmed the new powerswitch constant, `rsu_twelve_tile_pattern`, and `sensor:7` hooks are present.

Known limitations / next step:
- This has not yet been compiled, exported with `dd_web_display`, or overlap-checked in the configured DD4hep environment.
- Visual inspection should confirm that the 12 active tile regions per RSU and passive power-switch strips are clear before moving to LEC/REC details.

## 2026-06-05 UTC - Updated remaining phase order

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Record the revised implementation order after completing the RSU 12-tile approximation.

Roadmap:
- Phase 6: add simple passive LEC/REC detail.
- Phase 7: make corrugation geometry flexible by allowing row-wise `h`, `d`, and `theta` values for each disk through a dedicated corrugation CSV.
- Phase 8: add FPC and AncASIC detail.
- Phase 9: update the real placement CSV for corrugated-informed geometry and change the placement reference point from the bottom-left corner to the midpoint at the RSU-LEC boundary.
- Phase 10: run `eicrecon`/ACTS compatibility and performance checks.
- Phase 11: cleanup temporary files and prepare the pull request.

Validation:
- Documentation-only update. No geometry build or overlap checks run.

Known limitations / next step:
- Before implementing Phase 6, confirm practical first-pass dimensions and placement conventions for the LEC/REC boxes.

## 2026-06-05 UTC - Phase 6 LEC/REC XML scaffold

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Add XML-driven dimensions for first-pass passive LEC/REC geometry before changing the C++ module builder.

Implementation:
- Added LEC/REC lengths, widths, and thicknesses.
- Set `SiEndcapLEC_length = 4.5 mm` and `SiEndcapREC_length = 1.5 mm`, matching the `136.02 = 4.50 + 6 * 21.67 + 1.50` strip-length decomposition.
- Set LEC/REC widths to `SiEndcapRSU_width` for the first approximation.
- Set LEC/REC thicknesses to `50*um` as chip-like silicon placeholders.
- Removed the temporary FEC aliases and kept REC naming in the XML.
- Removed zero-valued LEC/REC gap constants because they are not needed until a nonzero physical separation is specified.
- Left FPC geometry for the next phase, as planned.

Validation:
- XML-only change; run `xmllint` and `git diff --check` before implementation continues.

Known limitations / next step:
- The LEC/REC dimensions are first-pass interpretations of the strip-length margins. The next C++ step should place simple passive boxes and visually verify whether this convention matches the drawing and handed module behavior.

## 2026-06-05 UTC - Normalize LEC/REC naming

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Remove ambiguous FEC naming introduced during the Phase 6 scaffold and keep the implementation aligned with the LEC/REC terminology used in the design notes.

Implementation:
- Removed temporary `SiEndcapFEC_*` aliases.
- Kept `SiEndcapLEC_length = 4.5 mm` and `SiEndcapREC_length = 1.5 mm`.
- Removed unused legacy LEC/REC FPC width constants after the bridge FPC constants were introduced.
- Removed zero-valued LEC/REC gap constants; these can be reintroduced later if a real physical gap is specified.
- Updated the current roadmap wording from LEC/FEC to LEC/REC.

Validation:
- XML and whitespace checks should be rerun after this cleanup.

## 2026-06-05 UTC - Phase 6 LEC/REC passive boxes

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Implement the first-pass LEC/REC geometry using the XML constants added for Phase 6.

Implementation:
- Added `LocalBoxTemplate` for passive boxes that are placed inside a stack component without adding to the total module thickness.
- Added LEC and REC passive silicon boxes to the corrugated 6-RSU sensor/electronics layer.
- Used `SVTReadoutVis` for the LEC/REC boxes.
- Kept LEC/REC non-sensitive; no sensor IDs or surfaces are created for them.
- Used handed placement:
  - `left`: LEC before the RSU chain, REC after it.
  - `right`: REC before the RSU chain, LEC after it.

Validation:
- Created and ran `tmp/validate_phase6_lec_rec_2026-06-05.sh`.
- `git diff --check` passed.
- `xmllint --noout compact/tracking/silicon_disks_modules.xml` passed.
- Confirmed the LEC/REC XML constants and C++ local-box hooks are present.
- Confirmed no `SiEndcapFEC_*` or zero-gap constants remain in the XML/C++ implementation.

Known limitations / next step:
- This is still a simple rectangular approximation. It needs DD4hep export, visual handedness inspection, and overlap checks before moving to FPC/AncASIC detail.

## 2026-06-06 UTC - Make LEC/REC visible in display

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `compact/display.xml`
- `compact/display_detailed.xml`
- `compact/display_geoviewer.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Address visual inspection feedback that LEC/REC boxes were not visible in the DD4hep display.

Implementation:
- Changed the corrugated 6-RSU module volume from `TrackerModuleVis` to `SVTModuleVis`.
- Added `SVTModuleVis` with `showDaughters="true"` in the standard, detailed, and geoviewer display XMLs.
- Added `SVTElectronicsVis` in the same display XMLs.
- Changed LEC/REC local boxes from `SVTReadoutVis` to `SVTElectronicsVis` so they are visually distinct from inactive RSU readout/periphery strips.

Validation:
- Updated and ran `tmp/validate_phase6_lec_rec_2026-06-05.sh`.
- `git diff --check` passed.
- `xmllint --noout` passed for:
  - `compact/tracking/silicon_disks_modules.xml`
  - `compact/display.xml`
  - `compact/display_detailed.xml`
  - `compact/display_geoviewer.xml`
- Confirmed `SVTModuleVis`, `SVTElectronicsVis`, and LEC/REC local-box hooks are present.

Known limitations / next step:
- If LEC/REC are still absent after rebuilding/exporting, confirm the geometry was rebuilt from this source and inspect the exported ROOT geometry for `component2_lec` and `component2_rec` volumes.

## 2026-06-07 UTC - Move LEC/REC into RSU pattern placement

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `tmp/validate_phase6_lec_rec_2026-06-05.sh`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Fix visual/export feedback that component 2 contained backbone, powerswitch, and periphery volumes but no LEC/REC volumes.

Implementation:
- Removed the generic `LocalBoxTemplate` / `local_boxes` scaffold.
- Added explicit `rsu_end_electronics` and `lec_after_rsu` flags to `ComponentTemplate`.
- Placed `component2_lec` and `component2_rec` directly inside the `rsu_twelve_tile_pattern` branch using the same `place_box(...)` helper that creates the visible RSU backbone, powerswitch, bias, periphery, and sensor boxes.
- Used `SiEndcapLEC_thickness` and `SiEndcapREC_thickness` explicitly for those boxes.
- Preserved handedness:
  - `left`: LEC before the RSU chain, REC after it.
  - `right`: REC before the RSU chain, LEC after it.
- Kept `SVTElectronicsVis` for LEC/REC.

Validation:
- Updated and ran `tmp/validate_phase6_lec_rec_2026-06-05.sh`.
- `git diff --check` passed.
- `xmllint --noout` passed for the tracking and display XML files.
- Confirmed the direct `component%d_lec` / `component%d_rec` placement hooks are present.
- Confirmed the removed `LocalBoxTemplate` / `local_boxes` scaffold is absent.

Known limitations / next step:
- Rebuild/export and confirm `component2_lec` and `component2_rec` are present in the visualization or exported geometry tree.

## 2026-06-07 UTC - Phase 7a corrugation row CSV scaffold

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `compact/tracking/silicon_disks_modules.xml`
- `compact/tracking/SVT_endcap_corrugation_rows_uniform.csv`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Make corrugated support construction row-configurable so users can provide row-wise pitch, angle, and height values.

Implementation:
- Added `rowwise_placement` support to corrugated XML `<frame>` blocks.
- Added `CorrugationRow` and a dedicated CSV loader for corrugation rows.
- Supported exact disk keys plus `*` / `all` rows for shared tooling tables. Exact disk rows override shared rows for that disk.
- Supported CSV columns:
  `disk,row_y_mm,h_mm,d_mm,theta_deg,enabled,comment`.
- Also allowed `half_pitch_mm` or full `pitch_mm` in place of `d_mm`.
- Interpreted `row_y_mm` as the positive-y lower-flat center for one corrugation cell.
- Mirrored nonzero rows to negative y; `row_y_mm = 0` is placed only once.
- Kept XML `h`, `d`, and `theta` as fallback behavior if no row CSV is configured or no rows match a disk.
- Added `compact/tracking/SVT_endcap_corrugation_rows_uniform.csv`, an explicit-row table using the current uniform values `h=6.0 mm`, `d=17.385 mm`, and `theta=35 deg`.
- Pointed all existing corrugated frame blocks in `silicon_disks_modules.xml` at the uniform reference CSV.

Validation:
- Created and ran `tmp/validate_phase8_corrugation_rows_2026-06-07.sh`.
- `git diff --check` passed.
- `xmllint --noout compact/tracking/silicon_disks_modules.xml` passed.
- Confirmed the uniform corrugation row CSV has seven fields on each row.
- Confirmed source hooks for `rowwise_placement`, `CorrugationRow`, `load_corrugation_rows`, mirrored row placement, and `pitch_mm` fallback are present.
- Confirmed all corrugated frame blocks in `silicon_disks_modules.xml` reference `compact/tracking/SVT_endcap_corrugation_rows_uniform.csv`.

Known limitations / next step:
- Phase 7a is intended to reproduce the previous uniform corrugation while exercising the CSV path. A follow-up nonuniform test CSV should vary at least one row and be checked visually before using this to drive production placement generation.

## 2026-06-07 UTC - Rename corrugation CSV XML hook

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`
- `tmp/validate_phase8_corrugation_rows_2026-06-07.sh`

Intent:
- Use clearer naming for the corrugation row CSV XML attribute.

Implementation:
- Renamed the corrugated-frame XML attribute from `row_file` to `rowwise_placement`.
- Renamed the corresponding C++ configuration field to `rowwise_placement`.
- Updated the Phase 7 validation script to check the new name.

Validation:
- Run the Phase 7 validation script after this entry.

## 2026-06-07 UTC - Normalize phase numbering

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Avoid future confusion after deciding to implement flexible corrugation before FPC/AncASIC detail.

Implementation:
- Renumbered flexible corrugation geometry as Phase 7.
- Renumbered FPC / AncASIC detail as Phase 8.
- Kept production placement CSV migration as Phase 9, reconstruction/ACTS validation as Phase 10, and cleanup/PR preparation as Phase 11.
- Updated the Phase 7a corrugation row CSV scaffold log entry.
- Marked the older May roadmap entry as superseded instead of leaving conflicting phase numbers.

Validation:
- Documentation-only update; run `git diff --check`.

## 2026-06-08 UTC - Phase 8 FPC planning notes

Files changed:
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `RSU_DRAWING_SUMMARY.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Record the FPC design-report dimensions before implementing Phase 8 geometry.
- Keep the main FPC explicitly deferred because its length is row-dependent.

Implementation notes:
- Phase 8 should start with simple, parameterized bridge FPC geometry.
- The left bridge FPC can use a first-pass rectangular approximation of `27.432 mm x 10.0 mm`.
- The right bridge FPC can use a first-pass rectangular approximation of `19.7612 mm x 4.0 mm`.
- The bridge FPC stack can be approximated as `80 um` Kapton plus effective aluminum layers from the reported fill factors.
- The main FPC should be treated as a later row-level geometry feature, not a fixed module-local box.

Main FPC to-do:
- Preserve page-2 reference values for later implementation: connector width `12.5476 mm`, body width `9.1186 mm`, sensor-to-main-FPC center distance about `17.385-17.392 mm`, RSU overlap `1 mm`, and top/bottom fill factors `0.45/0.68`.
- Decide whether main FPC length should be derived from the production placement CSV or a dedicated row-level FPC CSV.
- Implement and validate main FPC only after the row placement/reference-point migration is stable.

Validation:
- Documentation-only update; run `git diff --check`.

## 2026-06-08 UTC - Phase 8a bridge FPC XML scaffold

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Add explicit XML constants for first-pass bridge FPC geometry before touching the DD4hep placement code.
- Keep main FPC implementation deferred because its length is row-dependent.

Implementation:
- Added shared bridge FPC Kapton and nominal aluminum thickness constants.
- Added left bridge FPC rectangular-approximation dimensions: `27.432 mm x 10.0 mm`.
- Added left bridge FPC effective aluminum thicknesses from the report fill factors: `0.69 * 15 um` and `0.94 * 15 um`.
- Added right bridge FPC rectangular-approximation dimensions: `19.7612 mm x 4.0 mm`.
- Added right bridge FPC effective aluminum thicknesses from the report fill factors: `0.42 * 15 um` and `0.0 um`.
- Bridge FPC thickness is represented by explicit Kapton and effective aluminum layer constants.

Validation:
- Run `git diff --check`.
- Next implementation step should wire these constants into the corrugated module builder as passive, non-sensitive local boxes.

## 2026-06-08 UTC - Phase 8b bridge FPC passive boxes

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Place the first-pass bridge FPC material in the corrugated 6-RSU module prototype.
- Keep bridge FPCs passive and separate from sensitive RSU tile volumes.

Implementation:
- Added a dedicated bridge-FPC component pattern to the module-template builder.
- Placed bridge FPCs as a separate passive layer above the sensor layer so the FPC stack does not inflate the active sensor thickness.
- Used report `height` as the local-y span and report `width` as the local-x span.
- Placed the left bridge FPC on the LEC side and the right bridge FPC on the REC side.
- Mirrored the FPC side assignment with the existing module handedness convention.
- Built each bridge FPC as layered boxes:
  - bottom effective aluminum, skipped when the thickness is zero
  - Kapton
  - top effective aluminum
- Used `Kapton` and `Aluminum` materials with `SVTReadoutVis`.

Assumptions:
- The bridge FPC rectangular equivalents are material approximations, not exact outlines.
- The main FPC remains deferred because its length depends on row-level module placement.

Validation:
- Run `git diff --check`.
- Next validation should rebuild/export and visually confirm the bridge FPC volumes move with left/right handedness.

Local validation note:
- `git diff --check` passed in the plain login-shell environment.
- A build/install validation script was added at `tmp/validate_phase8_bridge_fpc_2026-06-08.sh`.
- Running the build script from the plain shell was blocked because the configured EPIC toolchain file `/opt/local/etc/cmake/find_package_resolve_symlinks.cmake` was not available.
- Re-run the script from the usual `eic-shell`/container environment before visual inspection.

## 2026-06-08 UTC - Phase 8b bridge FPC x-placement refinement

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Fix the visual-export warning where a bridge FPC was centered in a short `6.524 mm` end gap.
- Keep bridge FPC boxes inside the module parent envelope while allowing them to overlap inward over the RSU footprint in their own passive z-layer.

Implementation:
- Replaced the side-gap-centering helper with an edge-aware x-placement helper.
- If a bridge FPC fits inside the local end gap, it remains centered in that gap.
- If it is wider than the local end gap, it is placed flush with the package edge and allowed to extend inward over the RSU footprint.
- Removed the warning for this expected oversized-end-gap case.

Validation:
- Run `git diff --check`.
- Rebuild/export again and confirm the previous bridge-FPC width warning is gone.

## 2026-06-10 UTC - Phase 8b bridge FPC handedness and clearance fix

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Fix visual-inspection feedback that the bridge FPC x offsets were not mirroring correctly.
- Enforce the intended `0.5 mm` gap between each endcap piece and its bridge FPC, and another `0.5 mm` gap between each bridge FPC and the module edge.

Root cause:
- The bridge-FPC component used the default left-handed RSU-chain `x_offset`.
- `with_corrugated_handedness(...)` updated the sensitive RSU component offset, but did not update the bridge-FPC component offset.
- As a result, the LEC/REC boxes and bridge FPC boxes could be referenced to different RSU-chain centers after handedness was applied.
- The previous edge-aware helper also did not explicitly place FPCs between the LEC/REC outer edge and the module edge with the requested clearances.

Implementation:
- Updated the bridge-FPC component `x_offset` with the same handed RSU-chain offset used by the sensitive RSU component.
- Replaced the bridge-FPC x placement helper with one that uses:
  - the handed RSU-chain edge,
  - the LEC or REC length,
  - the module edge,
  - and `SiEndcapModuleEndClearance`.
- Left bridge FPC is placed between the LEC outer edge and the nearest module edge.
- Right bridge FPC is placed between the REC outer edge and the nearest module edge.
- The implemented FPC widths are chosen to preserve both `0.5 mm` gaps.

Validation:
- Run `git diff --check`.
- Rebuild/export and visually confirm:
  - left/right FPCs move with handedness,
  - each FPC is between its corresponding endcap and module edge,
  - the `0.5 mm` endcap-FPC and FPC-edge gaps are visible.

## 2026-06-10 UTC - Phase 8b remove handedness-dependent FPC width

Files changed:
- `src/SiEndcapModuleTracker_geo.cpp`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Remove the small handedness-dependent bridge-FPC x span caused by recomputing the available end space from rounded RSU/endcap dimensions.

Root cause:
- The FPC placement helper used the actual RSU-chain edge plus LEC/REC length to infer the remaining end space.
- The constants are slightly over-constrained: `6 * 21.666 + 4.5 + 1.5 = 135.996 mm`, while the drawing-level strip length is `136.02 mm`.
- That `0.024 mm` difference made the inferred LEC-side and REC-side clear spans differ by handedness after mirroring.

Implementation:
- Bridge FPC x placement now uses the nominal side spans directly:
  - LEC side: `SiEndcapModule6RSU_left_extension`
  - REC side: `SiEndcapModule6RSU_right_extension`
- This makes the placed FPC x span independent of handedness.
- With `SiEndcapModuleEndClearance = 0.5 mm`, the LEC-side available bridge-FPC x span is `10.0 mm`.
- The implemented left bridge FPC x span is therefore set to `10.0 mm`, so the dimensions remain exactly mirrored under handedness.

Validation:
- Run `git diff --check`.
- Rebuild/export and confirm the FPC dimensions and offsets are mirrored exactly.

## 2026-06-10 UTC - Remove unused legacy FPC width constants

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`
- `SVT_ENDCAP_CORRUGATED_RSU_HANDOFF.md`

Intent:
- Remove obsolete LEC/REC FPC width names now that bridge FPC dimensions are represented by explicit left/right bridge constants.

Implementation:
- Deleted the obsolete legacy LEC/REC FPC width constants from the active XML.
- Updated old documentation/log references so they do not preserve the dead constant names as active implementation guidance.

Validation:
- Run `git diff --check`.
- Search for the obsolete legacy names to confirm no references remain.

## 2026-06-10 UTC - Audit and remove unused corrugated XML constants

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`
- `RSU_DRAWING_SUMMARY.md`
- `SVT_ENDCAP_CORRUGATED_RSU_HANDOFF.md`

Intent:
- Remove `SiEndcap...` constants that were only scaffolding and had no active C++ or XML-attribute references.
- Keep drawing-derived values in documentation when they are useful for future implementation, but stop presenting them as active XML variables.

Removed active XML constants:
- Generic FPC thickness, superseded by explicit bridge FPC Kapton/aluminum layer constants.
- Corrugated inner/support cutout dimensions that are not yet modeled.
- Drawing-level strip length, since current placement derives RSU-chain and end features from the active package/extension/margin constants.
- Left/right CF bottom-cut lengths, since the cutouts are intentionally ignored at this stage.
- RSU y margin, since current code directly uses the RSU width and package width.

Validation:
- Run an active-reference audit for all `SiEndcap...` constants in `compact/tracking/silicon_disks_modules.xml`.
- Run `git diff --check`.

## 2026-06-10 UTC - Update corrugated sensor and endcap thickness convention

Files changed:
- `compact/tracking/silicon_disks_modules.xml`
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`

Intent:
- Align the corrugated sensor/electronics layer thickness with the latest disk material-budget estimate.
- Avoid maintaining independent LEC/REC thickness constants that can drift from the layer they occupy.

Implementation:
- Updated `SiEndcapSensor_thickness` to `65*um` as an effective thickness from Nikki's disk material budget dated 2026-06-09.
- Set `SiEndcapLEC_thickness` and `SiEndcapREC_thickness` to `SiEndcapSensor_thickness`.
- This keeps LEC/REC in the same local z layer as the RSU sensor/electronics layer and makes their thickness follow future sensor-layer updates automatically.
- Updated `SiEndcapCorrugationCF_thickness` to `0.15*mm`.

Validation:
- Run `git diff --check`.

## 2026-06-10 UTC - Phase 8 bridge FPC y-placement follow-up

To-do:
- Revisit the left/LEC-side bridge FPC y placement.
- The current implementation centers the rectangular left bridge FPC in local y, so the whole FPC remains contained within the module parent volume and its center aligns with the RSU band.
- In the real design, the top edge of the left bridge FPC should align with the top edge of the RSU.
- Because the left bridge FPC y span is larger than the RSU/module y span, that alignment will force the lower part of the left bridge FPC to extend below the current module envelope.
- Before implementing this, decide whether the module parent volume should be enlarged in y for this passive feature, or whether the FPC should be represented outside the module-local prototype in a higher-level placement.
