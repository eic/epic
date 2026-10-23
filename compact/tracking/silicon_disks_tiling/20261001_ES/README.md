# Ernst S. SVT module tiling, 2026-10-01

From Ernst S. Module tiling 20261001 -
https://conferences.lbl.gov/event/2600/

This delivery provides two independently selectable placement scenarios:

- `rsu_opt/`: 2,164 modules—265 five-RSU and 1,899 six-RSU modules.
- `all_6rsu/`: 2,160 six-RSU modules.

Production XML currently selects `rsu_opt`. Switching scenarios requires
changing the scenario component of all six `module_placements` paths between
`all_6rsu` and `rsu_opt`; the shared corrugation table does not change.

Each scenario contains its own catalogue, ten disk placement CSVs, and ten
metadata files. Disk IDs 0–9 run from ED4 at negative z through HD4 at positive
z. The source module origin is the midpoint of the LEC–first-RSU boundary.
`z_baseplate_mm` is the corrugation reference surface despite its historical
column name; `z_sensor_mm` is the source sensor reference surface.

The source coordinates use the engineer convention with positive x toward the
center of the ring. The geometry loader applies the documented conversion
`x_EPIC = -x_engineer`, exchanges LEFT/RIGHT placed chirality, and preserves y,
z, and INWARD/OUTWARD facing. Source values and type identifiers remain intact
for provenance.

The catalogue's `active_*` values are nominal tiling metrics. Detailed DD4hep
sensitive silicon additionally excludes bias, backbone, and powerswitch strips,
so its modeled active area is smaller.

`SVT_endcap_corrugation_rows_33_36.csv` is shared by both scenarios. It is a
documented reconstruction from the delivery's A/B row families: B rows use the
33 mm profile, A rows use the 36 mm profile, and the centered transition uses a
17.25 mm adjacent-facet separation. Delivered module y coordinates are channel
centers. The facet at the same y is on the opposite corrugation surface, while
the module's two 2.5 mm edge-glue lines attach to the neighboring same-surface
facets. The table includes one boundary facet beyond each end of the module-row
lattice. It was not supplied as an authoritative corrugation table. The
simulation uses straight carbon flats/webs with explicit inward-module cutouts
rather than exact 5 mm-radius bends.

The current geometry applies two clearly named provisional normal translations:
inward modules move 0.030 mm away from the corrugation channel, and the
opposing-face relative offset remains 0.230 mm. Outward modules therefore
receive a total 0.260 mm outward translation. These translations remain
pending confirmation of the intended three-dimensional stack separation.
