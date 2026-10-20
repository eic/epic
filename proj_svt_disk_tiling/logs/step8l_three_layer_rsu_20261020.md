# Step 8l: three-layer endcap RSU

The vertex-barrel implementation in `compact/tracking/vertex_barrel.xml`
specifies a full-area 36 um inactive-Si substrate, a 14 um Si layer with
twelve sensitive tiles and passive silicon gaps, then a full-area 2.5 um Cu
layer. Apply that same layer order and thickness to each physical endcap RSU
in the 5- and 6-RSU tiling modules. The existing 50 um endcap silicon remains
36 + 14 um, and Cu is additional: the physical RSU is 52.5 um thick. This
implements the user's explicit material-stack choice; it does not change the
source CSVs, module footprints, LEC/REC, film, baseplate, or external epoxy.

The endcap builder keeps the existing 12-tile x/y partition and IDs. Inactive
Si and Cu each cover the complete 19.564 x 21.666 mm RSU. As in the barrel,
the 14 um active-Si layer is one full-area, nonsensitive mother volume per
RSU. The existing sensitive tiles and passive-silicon backbone, power,
periphery, and bias pieces are its daughters. Daughter material replaces
mother material in its region; it is not counted twice. The DD4hep sensitive
surfaces now use
the center of that 14 um active layer, rather than the center of the original
50 um tile boxes. The modeled sensitive x/y area is unchanged. The legacy
`EIC_LAS_6RSU_CORR` prototype is left untouched.

The DD4hep DetElement hierarchy now follows the physical tiling hierarchy:
module -> per-RSU active-Si mother -> sensitive tile. Each tile keeps its
sensor ID and surface, attached beneath the correct mother placement so its
global transform can be reconstructed. Legacy direct-child surfaces remain
unchanged.

The mother adds one geometry node per physical RSU (12,719 placements in
`rsu_opt`) and one navigation level for the active tiles. It adds no material
and does not change tile IDs or sensitive area, but geometry export and
overlap checks may take longer.

For visualization, the active-Si mother is invisible with visible daughters;
the sensitive active tiles stay gold, inactive Si is teal, passive active-Si
bands remain dark green, and full-area Cu is translucent dark green so the
tiles can be inspected from that side. These visualization attributes are
defined consistently in all three compact display variants.

The placement CSV's `z_sensor_mm` still denotes the top of the 50 um Si, not
the top of Cu, so the source/model reference differences stay 30 um inward
and 340 um outward. Increase the common module-normal clearance from 30 to
32.5 um: inward translation is -0.0325 mm and outward is +0.2625 mm, keeping
the 0.230 mm relative facing offset. After translation, the inward Si top is
2.5 um before the corrugation carbon surface and Cu occupies that final
2.5 um; its outer face touches carbon without positive-volume overlap.
FPC aluminum also ends 2.5 um before carbon. The module envelope remains
controlled by the taller FPC/AncASIC branch (640/720 um inward/outward).

Numerical checks: `validate_svt_three_layer_rsu_step8l.py`, the existing
module-dimension, Z-offset, corrugation, FPC, and placement-transform
validators must pass. The older current-geometry Z-offset validators have
been updated to 32.5 um. The Step 8l visual/ROOT and separate Geant4 scripts
record the expensive checks.

Validation on 2026-10-20: all static validators passed, both Step 8l shell
scripts passed `bash -n`, and `git diff --check` passed. The EIC-shell build,
install, and geometry export completed; the user reviewed the visualization
and found it good. The ROOT overlap report started at 15:47:30 UTC found zero
illegal overlaps/extrusions. The Geant4 run started at 15:58:06 UTC completed
and reported
`Geometry overlaps check completed !` and the script's final
`PASS: Geant4 overlap check completed without reported overlaps.` The earlier
15:30:57 UTC Geant4 run did not record a completion line; it is not counted as
a pass. Reports are retained beside this note under the corresponding
`check_svt_three_layer_rsu_*_step8l_<timestamp>.txt` names.

Because silicon and copper material distribution changes, the production
ACTS material map must be regenerated and validated before using this
geometry for reconstruction, even if overlap checks pass.
