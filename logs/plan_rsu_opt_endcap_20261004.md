# Plan: implement the supplied `rsu_opt` SVT endcap design

Scope: use `ES_material/rsu_opt/` as the preferred layout and keep the supplied `ES_material/rsu6/` layout selectable. Treat the delivered catalogues, module coordinates, and metadata as source data. Do not use placements or optimizers under `disk_layout/`. Record each implementation step in `logs/rsu_opt_endcap_20261004.md`, including files changed, formulas, commands, log paths, numerical checks, and remaining discrepancies.

## 1. Audit dimensions before geometry changes

Compare every current XML/C++ dimension and stack component with the PDF's Corrugation, Tiling, and Summary slides and both catalogues. Record matches, differences, and missing components in the log. This comparison is recorded in `logs/rsu_opt_endcap_20261004.md`. The user selected the slide's 50 µm sensor, 60 µm film adhesive, 150 µm baseplate, and separate 80 µm epoxy for both module lengths. Retain the current code's inactive bias, backbone, and powerswitch strips: 18.394 mm sensitive y width and 21.426 mm sensitive x length per RSU. These differ from the catalogue's 18.514 mm and 21.666 mm per RSU tiling metrics, producing 1.7487% less modeled sensitive area for the supplied mixed layout. Preserve both area definitions in validation. Check the PDF's 33/36 mm pitch, 45° angle, 6 mm depth, 0.15 mm skin, 5 mm bend radius, and 20.56 × 1.00 mm keep-clear envelope against the present frame.

## 2. Validate and preserve delivered data

Check both catalogues and all ten per-disk files for schema, type references, unique `(disk_id,row_index,mod_index)`, z centers, dimensions, and module/active-area totals against metadata. Record source checksums and generation versions. Expected mixed-design total: 2,164 modules, including 265 five-RSU and 1,899 six-RSU. Keep `ES_material/` read-only.

## 3. Package each scenario independently

Use `compact/tracking/silicon_disks/all_6rsu/` and `compact/tracking/silicon_disks/rsu_opt/`. Each directory holds its own `catalog.csv`, ten `*_modules.csv`, ten `*_metadata.txt`, and a short provenance README. Copy the source data byte-for-byte. If the plugin needs one file, generate a deterministic combined CSV inside each scenario directory while retaining the originals.

## 4. Build configurable five/six-RSU modules

Replace the current hard-coded six-RSU corrugated prototype with one builder parameterized by catalogue type: RSU count, 32/35 mm baseplate width, length, active/sensor dimensions, facing, and chirality. Put editable mechanical parameters in XML or the catalogue, with clear ownership: RSU pitch, sensor, baseplate, film and epoxy layers, BLEC/BREC, LEC/REC, bridge FPC, AncASIC, clearances and offsets. Initially share end electronics between five and six RSUs as confirmed by the user. The delivered mixed catalogue has B5xx, B6xx, and A6xx; do not invent A5xx placements. Validate derived lengths and clearances against catalogue values.

Review gates for this step: (4a) stage and check catalogue-aligned XML dimensions without changing the active prototype; (4b) build the detailed five/six-RSU solids and their distinct material interfaces, then validate representative types and z references. Stop for user review after each gate.

## 5. Place modules in the delivered coordinate convention

Read `disk_id,row_index,mod_index,type_id,x_origin_mm,y_origin_mm,z_baseplate_mm,z_sensor_mm` directly. Map disk IDs 0–9 to the ten current layers. Convert the LEC–first-RSU boundary origin into the solid transform, including facing and chirality. Make the modeled baseplate and sensor reference surfaces agree with both supplied absolute z coordinates, including negative-side reflection. Preserve `(disk,row,module,type)` identifiers for auditability and report unknown or malformed rows as errors.

Review gates: (5a) audit x/y and both z columns against the staged packages and resolve the physical z reference surfaces; (5b) add the strict CSV-directory reader and coordinate transform. Do not select a scenario until the transform can reproduce both z anchors.

## 6. Apply per-disk boundaries

Use each metadata file for z center, outer radius, and the two opening primitives. Check the entire baseplate and sensor footprint against those boundaries and confirm the layer z envelope contains the full component stack. Keep the metadata's module-radius target distinct from its sensor-radius summaries.

## 7. Rebuild the corrugated support for the delivered layout

ED0/ED1/HD0/HD1 use 33 mm pitch. The remaining disks use central 36 mm regions and outer 33 mm regions, with transitions on slopes as stated in the PDF. Build disk-specific support rows and transitions with the approved dimensional parameters. Module y origins may move within the keep-clear envelope, so derive facet locations from the source algorithm or documented support-row convention, not directly from `y_origin_mm`. Resolve any missing bend/transition specification before labeling this step exact.

## 8. Expose a scenario switch

Provide two named configurations or a compact XML include selection pointing to the corresponding self-contained input package. Use the same geometry plugin for both. Document exact commands and verify changing scenario requires no C++ edit or rebuild.

## 9. Validate and record results

Run narrow catalogue, count, transform, surface-z, pitch-transition, and area checks first. Compare the catalogue-based area to metadata and the bias/readout-strip-sensitive area to the modeled silicon separately. Then compile and export both geometries in the stable 26.09.0 container and run `checkOverlaps --option m --t 0.01` and `scripts/checkOverlaps.py`. Inspect representative 5/6-RSU, A/B, chirality/facing, and transition cases. Compare all per-disk counts and source-reported active areas to metadata. Run a small recorded simulation/reconstruction probe if a remaining geometry or hit-coverage question warrants it. Write the script before any multistep run; store output in named parent-repo logs.

## 10. Complete handoff

Update the log and README with the final stack decision, source checksums, configuration names, commands, validation results, and any remaining approximations. Acceptance requires both scenarios to load, every supplied module to be placed, reference coordinates to match within CSV rounding, and the support to follow the delivered 33/36 mm regions.
