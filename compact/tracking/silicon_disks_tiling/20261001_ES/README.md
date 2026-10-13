# Ernst S. SVT module tiling, 2026-10-01

From Ernst S. Module tiling 20261001 -
https://conferences.lbl.gov/event/2600/

This delivery provides two independently selectable placement scenarios:

- `rsu_opt/`: 2,164 modules—265 five-RSU and 1,899 six-RSU modules.
- `all_6rsu/`: 2,160 six-RSU modules.

Production XML currently selects `rsu_opt`. Selecting the alternative requires
changing the scenario component of all six `module_placements` paths from
`rsu_opt` to `all_6rsu`; the shared corrugation table does not change.

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
17.25 mm adjacent-flat separation. It was not supplied as an authoritative
corrugation table. The simulation uses straight carbon flats/webs with explicit
inward-module cutouts rather than exact 5 mm-radius bends.

The current geometry also applies a clearly named provisional 0.230 mm outward
module z offset pending confirmation of the intended three-dimensional
inward/outward stack separation.
