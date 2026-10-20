# Step 8j: negative-Y electronics orientation

The disk visualization showed the asymmetric module electronics upside down
below global Y=0. The placement transform maps local +y to global -Y for
negative-Y rows and Y=0 rows. In `rsu_opt` this covers 1,087 of 2,164
modules; in `all_6rsu` it covers 1,085 of 2,160. The rotation preserves the
provided LEC-boundary x midpoint and sensor-facing z direction, so it remains
unchanged.

The tiling builder now uses a separate cached prototype for those rows. Only
bridge-FPC and AncASIC local-y positions are reversed inside that prototype;
after the unchanged module rotation, these asymmetric electronics have the
same global-Y orientation as on positive-Y rows. Silicon tiles and IDs,
baseplate, film, external epoxy, x/z coordinates, source type IDs, and module
rotation remain unchanged. This does not address the separate A-module FPC
x-alignment issue.

Validation:

- `validate_svt_upright_modules_step8j.py` passed all 4,324 placements in both
  scenarios. Existing FPC-registration and placement-transform validators,
  shell syntax checks, and `git diff --check` passed.
- The EIC-shell build and install completed using
  `check_svt_upright_modules_visual_root_step8j.sh`. The exported geometry
  passed the user's visual review.
- The ROOT overlap check reported **zero illegal overlaps/extrusions** at
  0.01 mm tolerance. The recorded reports are
  `check_svt_upright_modules_visual_root_step8j_20261020T140524Z.txt` and
  `check_svt_upright_modules_root_step8j_20261020T140524Z.txt`.
- The separate `check_svt_upright_modules_geant_step8j.sh` has **not** been
  run. At the user's request, defer the expensive Geant volume-overlap check
  until after the A-module FPC alignment is fixed. Do not describe Step 8j as
  Geant-validated before that check completes.
