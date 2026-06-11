# SVT Endcap Corrugated Module Geometry Project

This document is the working project guide for updating the SVT endcap module geometry in the `epic` detector codebase. It is written for humans and future coding agents. Keep it current as the design and implementation evolve.

Related notes:

- `SVT_ENDCAP_CORRUGATED_RSU_HANDOFF.md`: earlier handoff/status note.
- `RSU_DRAWING_SUMMARY.md`: dimensions and design observations extracted from `rsu.pdf`.
- `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`: append-only implementation log.

## Goal

Replace the current simplified SVT endcap module model with a clearer, parameterized corrugated/endcap-packaged module model that captures tracking-relevant geometry without over-modeling every mechanical detail.

The implementation should:

- Preserve existing simple module types during early development.
- Add new corrugated module types incrementally.
- Keep dimensions in XML constants where practical, because the mechanical design is still changing.
- Keep placement data in the CSV, including handedness if that is the cleanest representation.
- Mark only active silicon sensor regions as sensitive.
- Represent passive silicon/FPC/CF/adhesive/readout geometry only to the level needed for material and tracking studies.

## Current Geometry State

Main implementation:

```text
src/SiEndcapModuleTracker_geo.cpp
```

Main XML:

```text
compact/tracking/silicon_disks_modules.xml
```

Placement CSV:

```text
compact/tracking/SVT_endcap_modules_corrugation_all.csv
```

Current built-in module names:

```text
EIC_LAS_6RSU
EIC_LAS_5RSU
```

Current simplified nominal footprints:

```text
EIC_LAS_6RSU: 130.0 mm x 30.0 mm
EIC_LAS_5RSU: 105.0 mm x 30.0 mm
```

Current stack:

```text
CarbonFiber sheet, passive
Glue, passive
Silicon sensor, sensitive
```

Current corrugated prototype:

- `EIC_LAS_6RSU_CORR` is implemented separately from the legacy `EIC_LAS_6RSU`.
- `ComponentTemplate` supports local component offsets, `x_repeat`, and `rsu_twelve_tile_pattern`.
- `ComponentTemplate` can also carry passive local boxes that share a stack layer without adding to the total module thickness.
- The corrugated 6-RSU sensor band is split into six RSU pitches.
- Each RSU pitch is approximated as twelve active silicon tile rectangles with inactive silicon bias, power-switch, periphery/readout, and backbone regions.
- First-pass passive LEC/REC boxes are placed beside the RSU chain. In `left` handed modules, LEC is before the RSU chain and REC is after it; in `right` handed modules this is mirrored.
- Left/right handedness is parsed from CSV and included in the prototype cache key. The current implementation swaps the sensor-band start side using the left/right end-extension and sensor-margin constants. More detailed asymmetric FPC/readout/end features are still deferred.

## Target Geometry Concept

Add new module names first, rather than replacing the existing ones:

```text
EIC_LAS_6RSU_CORR
EIC_LAS_5RSU_CORR
```

Preferred handedness nomenclature:

```text
left
right
```

Do not use top/bottom as the public naming convention for the module variant. Page 4 of `rsu_all.pdf` shows mirrored module handedness. The RSU strip and asymmetric FPC/readout/end features move with the handedness.

The CSV may gain an explicit handedness column if that remains the cleanest approach:

```csv
disk,module,x_min_mm,y_min_mm,dz_mm,facing,handedness,enabled,comment
```

Rationale:

- Handedness is placement information.
- The CSV is already the placement source of truth.
- Explicit values avoid fragile inference from `x`, `y`, `dz`, or corrugation row.

## Design Inputs Captured So Far

From `rsu.pdf` and `rsu_all.pdf`:

```text
RSU pitch/length along module x: 21.67 mm
RSU active/sensor width:        19.56 mm
Six-RSU strip length:          136.02 mm
Full package length:           152.02 mm
Implementation package width:    32.00 mm
Drawing package width marker:    30.00 mm
Inner width marker:             29.00 mm
RSU y margins:                   6.22 mm on each side
```

Useful decompositions:

```text
136.02 = 4.50 + 6 * 21.67 + 1.50
152.02 = 11.00 + 136.02 + 5.00
32.00  = 6.22 + 19.56 + 6.22
```

Carbon fiber support cutouts:

```text
bottom-left cutout:   3.00 x 11.00 mm^2
bottom-right cutout:  3.00 x 5.00 mm^2
```

For the first implementation, use a simplified rectangular `32 mm` package/support envelope and ignore the bottom-edge cutouts.

Use distinct component widths rather than forcing all layers to the same rectangle:

```text
carbon fiber/support envelope: 32.00 mm wide
RSU/sensor active band:        19.56 mm wide
RSU y margins:                  6.22 mm on each side
```

Glue/adhesive width should follow the component it bonds when that is clear. For the first implementation, use simple rectangular glue/adhesive boxes and keep the choice explicit in the module-building helper.

The original mechanical design is a work in progress, so constants and helper logic should be easy to update.

## Modeling Philosophy

Do model:

- Full package/envelope length and width at a practical approximation level.
- Left/right handedness.
- Different local box sizes for support, glue/adhesive, and sensor bands when the design specifies them.
- Separate passive support/readout/FPC regions where they matter for material or geometry.
- RSU active regions and inactive silicon gaps at a useful granularity.
- Twelve active tile regions per RSU using the current RSU architecture approximation.

Do not over-model at first:

- Tiny screw holes, pins, tabs, or small local features unless they matter for overlaps or material.
- Carbon fiber bottom-edge cutouts.
- Exact curved FPC shapes unless needed later.
- Every individual tile, if four active rectangular regions per RSU is sufficient.

Sensitive-volume rule:

- Active silicon only gets `setSensitiveDetector`, `PhysVolID("sensor", ...)`, and `VolPlane`.
- Passive silicon gaps may use material `Silicon`, but must not be sensitive.
- CF, glue, adhesive, FPC, LEC/REC/readout structures are passive.

## Implementation Best Practices

- Go step by step and validate after each meaningful change.
- Prefer XML constants for dimensions and thicknesses.
- Keep C++ helpers simple and readable.
- Avoid burying design dimensions directly inside placement logic.
- Keep the CSV as placement data. Module construction belongs in C++/XML constants.
- Keep existing `EIC_LAS_6RSU` and `EIC_LAS_5RSU` until the new corrugated module is validated.
- Prefer small, inspectable commits or change sets.
- Never discard dirty worktree changes without explicit user approval.
- When running multi-step validation or simulation commands, first write a script for the record, then run the script.

## Logging Convention

Every substantial change should be recorded in:

```text
SVT_ENDCAP_CORRUGATED_MODULE_LOG.md
```

Each entry should include:

```text
date/time
files changed
summary of intent
implementation details
validation run
known limitations / next step
```

This log complements git history. It is intended to make debugging and reproducibility easier while the geometry design is moving.

## Suggested Stepwise Plan

### Phase 0: Baseline Decision

Decide whether to keep the current dirty prototype changes or reset to a cleaner baseline before implementing `EIC_LAS_6RSU_CORR`.

Decision:

- Keep documentation files.
- Start from a cleaner geometry/code baseline before implementing `EIC_LAS_6RSU_CORR`.
- Discard or set aside the current dirty source/XML prototype changes before new implementation work begins.
- Retain the current prototype only as a reference if it can be preserved separately and removed cleanly later.
- Do not reset documentation, drawing PDFs, or project notes as part of this cleanup.
- Before cleanup, list the exact files/actions and log the baseline operation.

### Phase 1: Parameter Scaffold

Add XML constants for the corrugated module dimensions and thicknesses.

Active constants:

```xml
SiEndcapModule6RSU_package_length
SiEndcapModule_width_corrugated = 32.0*mm
SiEndcapSensor_thickness
SiEndcapAdhesive_thickness
SiEndcapModuleCF_thickness
```

Drawing-only values such as CF cutouts, RSU y margins, and strip decomposition
should stay in documentation until they are actually wired into the geometry.

Validation:

- Confirm XML parsing/build succeeds.
- Confirm old modules still build.

### Phase 2: Component Placement Offsets

Extend internal component definitions with local offsets:

```cpp
double x_offset{0.0};
double y_offset{0.0};
```

Use these offsets in `build_module_prototype(...)`.

Validation:

- Existing simple modules should remain geometrically unchanged when offsets default to zero.
- Build/export a compact geometry and compare old module count and placement behavior.

### Phase 3: New Corrugated Module Type

Add `EIC_LAS_6RSU_CORR` alongside existing modules.

Initial approximation:

- CF support as a simple `32 mm` wide envelope rectangle.
- Adhesive/glue layer using an explicit design-driven width choice.
- RSU active silicon strip with six RSU pitches and `19.56 mm` sensor-band width.
- Passive FPC/readout/end pieces as simple boxes, deferred until same-layer side-by-side components are represented without artificial thickness stacking.
- No CSV switch yet, or only a tiny test CSV/config if needed.

Validation:

- Build geometry.
- Export with `dd_web_display`.
- Run overlap checks.
- Visually inspect module envelope and sensitive areas.

### Phase 4: Left/Right Handedness

Add explicit handedness support.

Preferred CSV option:

```csv
handedness
```

Allowed values:

```text
left
right
```

Backward-compatible behavior:

- Missing handedness defaults to a conservative value only for old module names.
- New corrugated module rows should specify handedness explicitly.
- Do not infer handedness from corrugation row or position.
- For the first geometry implementation, handedness should at least swap the sensor-band start side:
  `left` uses the left extension and left sensor margin, while `right` uses the right extension and right sensor margin.

Validation:

- Place one left and one right module in a small test geometry.
- Confirm mirrored local `y` offsets and asymmetric end features.
- Confirm sensitive surfaces remain attached correctly.

### Phase 5: RSU Active Pattern

Initial implementation used a vertex-barrel-style approximation:

```text
12 physical tiles -> 4 active rectangular regions
inactive gaps/regions -> passive Silicon
```

This has now been refined using `1RSU_top_half.png`:

```text
2 local-x halves * 3 tile columns * 2 local-y halves
= 12 active rectangular tile regions per RSU
```

Implementation notes:

- Use the RSU architecture dimensions as the current approximation:
  `tile_x = (RSU_length/2 - backbone_width - 3 * powerswitch_width) / 3`
  and `active_y = RSU_width/2 - bias_width - periphery_width`.
- Keep inactive silicon regions passive, using the same silicon material but without the sensitive detector.
- Put bias strips at the internal boundary between the two local-y halves, periphery/readout strips at the outer edges, backbone strips at the start of each x-half, and power-switch strips after each tile column.
- Keep all dimensions driven by XML constants so later design updates do not require C++ rewrites.

Validation:

- Confirm each RSU produces twelve sensitive volumes.
- Confirm inactive silicon is present but not sensitive.
- Confirm sensor IDs and `VolPlane` entries are sensible.

### Phase 6: LEC / REC Detail Pass

Add first-pass detail for the local/end electronics before changing the production placement CSV.

Scope:

- Add simple, parameterized LEC and REC approximations.
- Place LEC/REC as passive local boxes in the sensor/electronics layer so they do not inflate the full module stack thickness.
- Keep these features passive.
- Use placeholder materials already agreed for now:
  `Kapton` for FPC/base flex, `Copper` for optional metal/readout traces, and `Silicon` for chip-like placeholders.
- Avoid CAD-level detail unless it is needed for material budget, envelope, overlap, or tracking performance.

Validation:

- Confirm geometry export succeeds.
- Visually inspect left/right handedness and asymmetric end features.
- Confirm active RSU tile regions remain the only sensitive silicon volumes.
- Run overlap checks before building further electronics detail on top.

### Phase 7: Flexible Corrugation Geometry

Make the corrugated support geometry row-configurable instead of using one fixed corrugation shape per disk/layer.

Target design:

- Add a dedicated corrugation CSV.
- Allow `h`, `d`, and `theta` to be specified for each row for a given disk.
- Keep the XML defaults as fallbacks for simple configurations and early debugging.
- Keep the corrugation CSV separate from the module placement CSV so support-shape parameters and module placement data remain auditable independently.
- Use positive-y row definitions and mirror nonzero rows to negative y.
- Allow `disk="*"` rows for a shared tooling table across all disks, while also supporting exact disk-key rows for maximum flexibility.

Initial CSV schema:

```csv
disk,row_y_mm,h_mm,d_mm,theta_deg,enabled,comment
```

Column meanings:

- `disk`: exact disk key, or `*` / `all` for rows shared by all disks. Exact disk rows override shared rows for that disk.
- `row_y_mm`: positive-y lower-flat center for one corrugation cell. Nonzero rows are mirrored to negative y.
- `h_mm`: corrugation height.
- `d_mm`: half-pitch, matching the existing XML `d` convention. The parser also accepts `half_pitch_mm` or full `pitch_mm`.
- `theta_deg`: web angle in degrees.
- `enabled`: optional row on/off switch.
- `comment`: free-form provenance/debug text.

Phase 7a implementation:

- Add row CSV parsing and mirrored row placement.
- Add `compact/tracking/SVT_endcap_corrugation_rows_uniform.csv` as an explicit-row reference table that reproduces the current uniform corrugation dimensions.
- Point the existing frame XML blocks at the uniform reference CSV.

Validation:

- Confirm row-wise `h`, `d`, and `theta` values are parsed and applied to the intended disk rows.
- Confirm corrugated support pieces remain inside their disk volumes.
- Run geometry export and overlap checks with at least one non-uniform test configuration.

### Phase 8: FPC / AncASIC Detail Pass

Add the next layer of passive electronics detail after the LEC/REC model and flexible corrugation scaffold are stable.

Scope:

- Add simple, parameterized bridge FPC geometry first.
- Defer main FPC implementation until a row-level model is defined. The main FPC length depends on the number of modules in a row, so it should not be hard-coded as a fixed module-local feature.
- Add simple, parameterized AncASIC placeholders.
- Keep materials and dimensions XML-driven where practical.
- Preserve handedness behavior: all asymmetric electronics features should move consistently with the module handedness.

Current FPC assumptions from `FPC Design Report 20260331.pdf`:

- Left bridge FPC first-pass rectangular approximation: `27.432 mm x 10.0 mm`.
- Bridge FPC stack approximation: `80 um` Kapton plus `30 um` Aluminum.
- Right bridge FPC first-pass rectangular approximation: `19.7612 mm x 4.0 mm`.
- Main FPC reference values from page 2 should be logged for later implementation: connector width `12.5476 mm`, body width `9.1186 mm`, sensor-to-main-FPC center distance about `17.385-17.392 mm`, RSU overlap `1 mm`, and top/bottom fill factors `0.45/0.68`.

Main FPC to-do:

- Decide whether the main FPC belongs in the disk/row builder rather than the module prototype builder.
- Add a row-level CSV or derive main FPC length from the production placement CSV after the corrugated reference point migration.
- Validate row-dependent main FPC lengths with visual inspection and overlap checks before using it in production geometry.

Validation:

- Repeat geometry export, visual handedness inspection, and overlap checks.
- Confirm the added passive features do not change sensitive-volume placement or sensor IDs unexpectedly.

### Phase 9: Production Placement CSV Migration

Update the real placement CSV to use corrugated-informed geometry and the new reference point convention.

Initial target:

```text
EIC_LAS_6RSU_CORR
```

Later, add `EIC_LAS_5RSU_CORR` only after the 6-RSU implementation is stable and the 5-RSU dimensions/handedness assumptions are clear.

Rules:

- Corrugated module rows must always provide explicit `left` or `right` handedness.
- Do not infer handedness from row, disk side, `x/y`, or `dz`.
- Update the placement reference point from the current bottom-left corner convention to the midpoint at the RSU-LEC boundary.
- Regenerate or migrate placements using the corrugated module envelope and the row-wise corrugation geometry, rather than simply renaming legacy `EIC_LAS_6RSU` rows.
- Keep the temporary test CSV available only until the production CSV path is validated.

Validation:

- Full geometry export from the production CSV path.
- Careful visual inspection of left/right handedness.
- Confirm the new placement reference point lands at the RSU-LEC boundary midpoint.
- Confirm the RSU sensor band, electronics, glue, support, and active sensors remain visually distinguishable.
- Run both overlap checks again.
- Confirm the expanded `sensor` bitfield remains sufficient for all sensitive volumes.

### Phase 10: Reconstruction / ACTS Validation

Run a more thorough reconstruction validation after the updated production placement and corrugated geometry migration.

Goal:

- Catch ACTS surface, volume manager, readout-ID, or reconstruction issues that are not visible in DD4hep geometry export or overlap checks.
- Check reconstruction compatibility and performance against an agreed reference workflow.

Validation:

- Run an agreed `eicrecon` workflow using the updated detector geometry.
- Save the run script and relevant logs under `tmp/` or another project-local validation area.
- Record command, detector XML/config, input events, and result summary in the implementation log.
- Compare compatibility/performance against the previous reference geometry or agreed baseline.

### Phase 11: Cleanup and PR Preparation

After electronics detail, flexible corrugation geometry, production placement, overlap checks, and reconstruction validation pass:

- Revert temporary XML test switches.
- Remove the temporary corrugated test CSV.
- Remove or trim working markdown files unless any are intentionally promoted to permanent documentation.
- Inspect the final diff for unrelated changes.
- Prepare a concise PR summary with validation results, assumptions, known approximations, and deferred details.

## Open Questions

- What placeholder materials should represent LEC, REC, and FPC in DD4hep?
- Does the 5-RSU corrugated package have distinct end dimensions? For now, assume the same end dimensions as the 6-RSU package, but focus implementation on 6-RSU first.
- Should the existing current prototype be retained as a stepping stone or discarded before a cleaner implementation? Current decision: retain it only if it stays separate and can be removed cleanly later.

Resolved for first implementation:

- Use `32.00 mm` as the package/support width in code.
- Use distinct component widths where the design specifies them, especially `32.00 mm` support width versus `19.56 mm` RSU/sensor width.
- Ignore CF bottom-edge cutouts initially.
- Require explicit `left`/`right` handedness in the CSV for corrugated modules.
- Start implementation from a cleaner geometry/code baseline rather than layering on the current dirty prototype changes.
- Use placeholder materials already available in the repo:
  - FPC/base flex: `Kapton`
  - optional metal traces/readout strips: `Copper`
  - LEC/REC chip-like placeholders: `Silicon`
  - adhesive/glue: `SVT_Endcap_Glue`

## Validation Workflow Notes

For EIC geometry checks, use the local EPIC workflow documented by the user:

```text
cd ~/eic_dir
./eic-shell --version 26.02.0-stable
cd /global/homes/p/pdatta/weic/epic
source setup_local.sh
```

When running a multi-step validation, first create a script under a project-local scripts or tmp area, then run it. Record the script path and result in the implementation log.

Useful validation commands inside the configured environment:

```text
dd_web_display --export ${DETECTOR_PATH}/{detector_config}.xml
checkOverlaps --option m --t 0.01 -c ${DETECTOR_PATH}/{detector_config}.xml
python scripts/checkOverlaps.py -c ${DETECTOR_PATH}/{detector_config}.xml
```

For reconstruction validation, use the agreed local `eicrecon` workflow and first write a small run script for reproducibility. Record the script path, detector XML/config, input sample, and result summary in `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`.

Use the narrowest validation that exercises the current change before running broader checks.

## Pre-PR Cleanup

Once the implementation is complete and before opening a pull request:

- Revert the temporary XML test switch back to the production placement CSV.
- Remove `compact/tracking/SVT_endcap_modules_corrugation_6rsu_corr_test.csv`.
- Remove the working markdown files added for this implementation unless any are intentionally promoted to permanent documentation:
  `RSU_DRAWING_SUMMARY.md`,
  `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`,
  and `SVT_ENDCAP_CORRUGATED_MODULE_PROJECT.md`.

## Agent Checklist Before Editing

- Read this file.
- Read `SVT_ENDCAP_CORRUGATED_MODULE_LOG.md`.
- Check `git status --short`.
- Inspect relevant current diffs before touching source files.
- Identify whether changes are documentation, geometry construction, XML constants, CSV placement, or validation scripts.
- Log substantial changes after implementation and validation.
- Do not reset, delete, or overwrite user changes without explicit approval.
