# Ernst S. SVT module tiling, 2026-11-10

From Ernst S. Module tiling 20261110 (`transition/` variant: 33/36 mm pitch,
modules may shift for projectivity) - https://conferences.lbl.gov/event/2600/

Byte-for-byte copies of `ES_material/20261110/transition/{rsu_opt,rsu6}/`
(staged by `proj_svt_disk_tiling/logs/stage_svt_tiling_20261110_step11a.sh`):

- `rsu_opt/`: 2,164 modules (12 catalogue types, five- and six-RSU).
- `rsu6/`: 2,150 six-RSU modules (8 catalogue types).

Each scenario has `catalog.csv`, `disk_geometry_used.csv`, and per disk
`<disk>_modules.csv`, `<disk>_rows.csv`, `<disk>_metadata.txt`. SVG plots,
xlsx files, text reports and CAD exports remain only in `ES_material/`.

Production XML selects `rsu_opt` (Step 11b). Switching scenarios requires changing the
scenario component of all six `module_placements` paths.

Differences from `20261001_ES/` (handled by the loader since Step 11b, 2026-11-11):

- Disk 8 files are named `HD3_*` (was `HD3b_*`).
- Catalogue `chirality` is `NEGATIVE`/`POSITIVE` (was `RIGHT`/`LEFT`;
  A600 RIGHT -> NEGATIVE, A601 LEFT -> POSITIVE, same pattern for all types).
- OUTWARD modules' `z_baseplate_mm` and `z_sensor_mm` are shifted by 0.23 mm:
  `z_baseplate_mm - facet_z_mm` (row file) is now +/-0.23 mm, was 0 in
  20261001 (INWARD unchanged; x, y, types unchanged for rsu_opt). The shift
  is +0.23 mm along the module normal (away from the carbon) for every
  outward module, i.e. the delivery now includes the outward-vs-inward z
  offset. The builder's `SiEndcapTilingOppositeFacing_z_offset` is therefore
  0 (0.230 mm is needed only for 20261001).
- New `<disk>_rows.csv` (`disk_id,row_index,row_type,x_outer_mm,row_y_mm,
  facet_z_mm,half_pitch_mm`; `row_type` BLOCKED/UNBLOCKED, half pitch 16.5/18 mm)
  is the delivered per-row corrugation table, intended to replace the
  reconstructed `20261001_ES/SVT_endcap_corrugation_rows_33_36.csv`.
  NOT yet used: the corrugation is still built from that 20261001 table.
