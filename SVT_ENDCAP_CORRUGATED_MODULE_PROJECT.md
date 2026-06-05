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
- `ComponentTemplate` supports local component offsets, `x_repeat`, and `rsu_four_region_pattern`.
- The corrugated 6-RSU sensor band is split into six RSU pitches.
- Each RSU pitch is approximated as four active silicon rectangles with inactive silicon bias, periphery/readout, and backbone regions.
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
- Four active regions per RSU if using the vertex-barrel-style 12-tile approximation.

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

Likely constants:

```xml
SiEndcapModule6RSU_package_length
SiEndcapModule6RSU_strip_length
SiEndcapModule_width_corrugated = 32.0*mm
SiEndcapModule_inner_width_corrugated
SiEndcapRSU_y_margin
SiEndcapModule_support_bottom_cut_width
SiEndcapModule6RSU_left_cut_length
SiEndcapModule6RSU_right_cut_length
SiEndcapSensor_thickness
SiEndcapFPC_thickness
SiEndcapAdhesive_thickness
SiEndcapModuleCF_thickness
```

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

### Phase 5: RSU Four-Region Active Pattern

Update each RSU to use the vertex-barrel-style approximation:

```text
12 physical tiles -> 4 active rectangular regions
inactive gaps/regions -> passive Silicon
```

Implementation notes:

- Use the vertex-barrel dimensions as the current approximation:
  `active_x = RSU_length/2 - backbone_width`
  and `active_y = RSU_width/2 - bias_width - periphery_width`.
- Keep inactive silicon regions passive, using the same silicon material but without the sensitive detector.
- Put bias strips at the internal boundary between the two local-y halves, and periphery/readout strips at the outer edges, following the vertex-barrel `upper`/`lower` convention.
- Keep all dimensions driven by XML constants so later design updates do not require C++ rewrites.

Validation:

- Confirm each RSU produces four sensitive volumes.
- Confirm inactive silicon is present but not sensitive.
- Confirm sensor IDs and `VolPlane` entries are sensible.

### Phase 6: Production CSV Migration

Update the real placement CSV to use corrugated module names and explicit handedness.

Initial target:

```text
EIC_LAS_6RSU_CORR
```

Later, add `EIC_LAS_5RSU_CORR` only after the 6-RSU implementation is stable and the 5-RSU dimensions/handedness assumptions are clear.

Rules:

- Corrugated module rows must always provide explicit `left` or `right` handedness.
- Do not infer handedness from row, disk side, `x/y`, or `dz`.
- Keep the temporary test CSV available only until the production CSV path is validated.

Validation:

- Full geometry export from the production CSV path.
- Careful visual inspection of left/right handedness.
- Confirm the RSU sensor band moves to the intended handed side.
- Confirm the inactive RSU regions, glue, support, and active sensors remain visually distinguishable.
- Run both overlap checks again.
- Confirm the expanded `sensor` bitfield remains sufficient for all sensitive volumes.

### Phase 7: Reconstruction / ACTS Validation

Run a more thorough reconstruction validation after the production CSV migration.

Goal:

- Catch ACTS surface, volume manager, readout-ID, or reconstruction issues that are not visible in DD4hep geometry export or overlap checks.

Validation:

- Run an agreed `eicrecon` workflow using the updated detector geometry.
- Save the run script and relevant logs under `tmp/` or another project-local validation area.
- Record command, detector XML/config, input events, and result summary in the implementation log.

### Phase 8: FEC / LEC Detail Pass

Add more detail to the front-end/readout features after the baseline corrugated module and production placement path are validated.

Scope:

- Add simple, parameterized FEC/LEC approximations.
- Use placeholder materials already agreed for now:
  `Kapton` for FPC/base flex, `Copper` for optional metal/readout traces, and `Silicon` for chip-like placeholders.
- Keep these features passive.
- Avoid CAD-level detail unless it is needed for material budget, envelope, overlap, or tracking performance.

Validation:

- Repeat full geometry export.
- Re-check left/right handedness.
- Re-run both overlap checks.
- Re-run the `eicrecon`/ACTS validation.
- Compare against the pre-FEC/LEC validation results.

### Phase 9: Cleanup and PR Preparation

After production placement, geometry validation, overlap checks, and reconstruction validation pass:

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
