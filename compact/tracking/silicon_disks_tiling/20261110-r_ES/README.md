# Ernst S. SVT module tiling, 2026-11-10, "-r" package (systematic studies)

Byte-for-byte copies of `ES_material/20261110-r/transition/{rsu_opt,rsu6}/`
(staged by `proj_svt_disk_tiling/logs/stage_svt_tiling_20261110r_step11f.sh`; the provider's
`flip_facing_build.log` in each source scenario says every disk was "flipped").

Same file layout as `20261110_ES/`. Relative to `20261110_ES/` (same scenario), per
(disk_id,row_index,mod_index): x, y and the row unchanged; every module type is replaced by
the type with the opposite facing AND opposite chirality (A600<->A611, A601<->A610, same for
B5xx/B6xx); `rows.csv` identical; the z columns follow the same convention (inward
z_baseplate = facet_z, outward +0.230 mm along the module normal).

Not selected by production XML (which uses `20261110_ES/rsu_opt`). Load it through an
overlay: `proj_svt_disk_tiling/logs/make_svt_overlay.sh 20261110-r_ES <rsu_opt|rsu6>`.
Intended for systematic studies only.
