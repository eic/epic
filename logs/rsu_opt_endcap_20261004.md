# `rsu_opt` SVT endcap implementation log

Date: 2026-10-04. Scope: parent `epic` repository only. Source data: `ES_material/rsu6/` and `ES_material/rsu_opt/`. Plan: `logs/plan_rsu_opt_endcap_20261004.md`.

## Step 1: dimension audit and selected conventions

| Dimension | Supplied design | Current implementation | Finding |
|---|---:|---:|---|
| Small/large corrugation pitch | 33 / 36 mm | 34.77 mm uniform | Different; large pitch and transitions absent |
| Corrugation angle | 45° | 35° | Different |
| Corrugation depth | 6 mm | 6 mm | Matches |
| Carbon skin | 0.15 mm | 0.15 mm | Matches |
| Mid-surface bend radius | 5 mm | Straight box support | Missing |
| Keep-clear envelope | 20.56 × 1.00 mm | No explicit envelope | Missing |
| Sensor width | 19.56 mm on slide; 19.564 mm in catalogue | 19.564 mm | Matches catalogue |
| Active width | 18.514 mm in catalogue | 18.394 mm from `19.564 - 2 × (0.525 + 0.060)` | Different by 0.120 mm; retain the code's bias-strip convention for sensitive silicon |
| RSU chain / active length, six RSUs | 129.996 mm | `6 × 21.666 = 129.996 mm` nominal chain; active boxes are shortened by x readout strips | Nominal length matches, effective sensitive area needs a separate check |
| RSU chain / active length, five RSUs | 108.330 mm | No detailed five-RSU chain | Missing |
| Sensor length, five / six RSUs | 114.330 / 135.996 mm | Six-RSU chain plus 4.5 mm LEC and 1.5 mm REC; no five-RSU detailed version | Six-RSU nominal span matches; five-RSU missing |
| Sensor periphery | 0.525 mm | 0.525 mm | Matches |
| Sensor thickness | 0.05 mm on slide | 0.065 mm | Different; use 0.05 mm for both 5/6 RSUs in the new design |
| Film adhesive | 0.06 mm on slide | No separate film layer | Missing |
| Baseplate thickness | 0.15 mm | 0.15 mm | Matches |
| External epoxy adhesive | 0.08 mm on slide | One 0.08 mm module adhesive layer | Same nominal thickness; placement/meaning differs |
| Baseplate widths | 32 / 35 mm | 32 mm only | 35 mm missing |
| Five-RSU module length | 138.330 mm in mixed catalogue | Legacy flat 105 mm; no detailed corrugated five-RSU module | Missing |
| Six-RSU module length | 159.996 mm in catalogue | 152.02 mm | Different |
| BLEC / BREC | 20 / 10 mm | 11 / 5 mm end extensions in current prototype | Different conventions/dimensions; do not equate without coordinate audit |
| LEC / REC | 4.5 / 1.5 mm | 4.5 / 1.5 mm | Matches |
| 33 mm zero-skin flat width | 6.36 mm on Corrugation slide | About 8.82 mm from present straight-segment formula at 34.77 mm and 35° | Different construction; verify again when curved 45° support is built |
| Adhesion to corrugation facet | 2.5 mm glue line on Tiling slide | No dedicated glue-line constraint | Missing |
| Sensor keep-clear margin | About 0.5 mm each side in y, 1 mm in z on Tiling slide | No explicit keep-clear solid/check | Missing |

The mixed catalogue has B5xx, B6xx, and A6xx. It has no A5xx rows. The PDF's 2,164 total agrees with the mixed-design metadata: 265 five-RSU and 1,899 six-RSU. ED0/ED1/HD0/HD1 are uniform 33 mm regions; other disks have 36 mm central and 33 mm outer regions.

**Decision (user, 2026-10-04):** Keep the current C++ y-sensitive convention, including the 0.060 mm bias strip in each y half. Each half-width subtracts one 0.525 mm periphery and one 0.060 mm bias strip, so sensitive y width is 18.394 mm.

**Decision (user, 2026-10-04):** Keep the current C++ passive backbone and powerswitch strips in x. The PDF's module-catalogue slide lists `active_length_mm = 108.330` for five RSUs and `129.996` for six RSUs, exactly `rsu_count × 21.666 mm`; it does not subtract these detailed readout strips. C++ leaves sensitive x length per RSU of `21.666 - 2 × 0.060 - 6 × 0.020 = 21.426 mm`. With the confirmed y-bias convention, modeled sensitive area is 1,970.54922 / 2,364.659064 mm² per five/six-RSU module, versus catalogue area of 2,005.62162 / 2,406.745944 mm². For the delivered 265 five-RSU and 1,899 six-RSU modules, the catalogue sum is 5,101,900.276956 mm² and the modeled sensitive-silicon sum is 5,012,683.105836 mm²: 89,217.171120 mm² less, or **1.7487% below the catalogue sum**. This is the expected roughly 2% difference from retaining the detailed inactive strips in both x and y. Preserve the source catalogue and metadata areas as tiling metrics; report DD4hep sensitive-silicon area separately, without modifying source CSVs.

**Decision (user, 2026-10-04):** Use the Corrugation slide's 50 µm sensor, 60 µm film adhesive, 150 µm baseplate, and distinct 80 µm epoxy adhesive for both five and six RSUs. The current single 80 µm module glue component must be split/relocated according to these physical interfaces; changing its value alone would not reproduce the slide. Keep these four thicknesses user-configurable. The slide's ±2 mm module-face labels and the CSV's baseplate z values at disk center ±3 mm require a reference-surface transform check during placement implementation; do not assume they describe the same surface.

The mechanical audit is ready to guide implementation. No source data or geometry changed in this step.

## Step 2: validate and preserve the delivered data

Recorded checker: `logs/validate_svt_handoff_step2.py`.
Command: `python3 logs/validate_svt_handoff_step2.py`.
Detailed report and SHA-256 for all 42 source files: `logs/validate_svt_handoff_step2.txt`.
Result: **PASS, zero issues**.

| Package | Modules | Five-RSU | Six-RSU | Catalogue nominal active area |
|---|---:|---:|---:|---:|
| `ES_material/rsu6/` | 2,160 | 0 | 2,160 | 5,198,571.239040 mm² |
| `ES_material/rsu_opt/` | 2,164 | 265 | 1,899 | 5,101,900.276956 mm² |

For each package the checker found the complete catalogue, ten placement CSVs,
and ten metadata files; no unexpected files. It checked placement headers,
unique `(disk_id,row_index,mod_index)` keys within each disk, type references,
the type-code facing/chirality/count/width meanings, module/sensor/active length
identities, disk IDs and z centers, baseplate and sensor z offsets, module counts,
nominal active-area totals, outer-radius targets, two-opening metadata, and
consistent provenance fields. Every per-disk nominal area matched metadata to
within 0.001 mm². Both packages report catalog version `2026-09-25`, code
version `5601d4b-dirty`, and generation time `2026-10-01T16:46:05Z`.

The metadata areas use the catalogue's nominal active dimensions. The
separately documented 1.7487% smaller detailed sensitive-silicon area is an
intentional model convention, not a source-data failure. No delivered source
file, geometry XML, or C++ file was modified during step 2.

## Step 3: stage the selectable source packages

Recorded script: `logs/stage_svt_handoff_step3.sh`.
Command: `bash logs/stage_svt_handoff_step3.sh`.
Verification report: `logs/stage_svt_handoff_step3.txt`.
Result: **PASS**. Copied the 21 files from `ES_material/rsu6/` into
`compact/tracking/silicon_disks/all_6rsu/`, and the 21 files from
`ES_material/rsu_opt/` into `compact/tracking/silicon_disks/rsu_opt/`.
The script compared every source and destination file byte for byte and
recorded destination catalogue SHA-256 hashes. Added
`compact/tracking/silicon_disks/README.md` for provenance and layout usage.
The original `ES_material` packages remain in place and unchanged. Neither
geometry XML nor C++ selects these files yet; that is planned after the module
model and coordinate conventions are implemented.

## Subsequent steps

Append each completed plan step here with date, files changed, equations and reference surfaces, commands or run script, named log/output paths, numerical checks, and unresolved issues.

## Step 4a: stage configurable module dimensions (2026-10-04)

Files changed: `compact/tracking/silicon_disks_modules.xml` (new tiling-only
constants), `logs/validate_svt_module_params_step4a.py` (read-only checker),
this log, and `logs/plan_rsu_opt_endcap_20261004.md` (review gates).

The new XML block specifies RSU counts 5/6, A/B baseplate widths 35/32 mm,
BLEC/BREC lengths 20/10 mm, and the slide-4 stack: 50 µm silicon, 60 µm
film adhesive, 150 µm carbon baseplate, and distinct 80 µm external epoxy.
The existing XML RSU length (21.666 mm), LEC/REC lengths (4.5/1.5 mm),
and editable bridge-FPC/AncASIC dimensions remain shared. The derived
envelopes are `N × 21.666 + 20 + 10 = 138.330/159.996 mm` and sensor
lengths are `N × 21.666 + 4.5 + 1.5 = 114.330/135.996 mm` for N=5/6.
The 35 mm A baseplate applies only to A6 types in the supplied catalogues;
no A5 placement is implied.

Commands: `python3 logs/validate_svt_module_params_step4a.py` and
`git diff --check`. Both passed. The checker parsed the XML, checked the
four slide stack thicknesses and length expressions, and matched all 8
all-six-RSU and 12 mixed-catalogue types (counts, widths, package/sensor
lengths, BLEC/BREC and LEC/REC). No catalogue CSV was changed.

This gate intentionally does **not** change the active C++ six-RSU prototype:
it still has its earlier 152.02 mm package, 65 µm sensor and single 80 µm
glue layer. The new parameters are not consumed until step 4b. That step
must give the film and external epoxy separate solids and resolve their
z-reference surfaces; changing one old glue value would be incorrect.
Existing LEC/REC thickness aliases also still point to the prototype's
65 µm sensor and must be wired to the new tiling sensor thickness in 4b.
No placement package is selected yet, and no geometry build or overlap
check was run at this parameter-only gate.

## Step 4b: detailed tiling-module builder (2026-10-04; review gate)

Files changed: `src/SiEndcapModuleTracker_geo.cpp`,
`compact/tracking/silicon_disks_modules.xml`,
`logs/validate_svt_module_builder_step4b.py`, and this log.
The existing `EIC_LAS_6RSU_CORR` prototype is retained. The new builder
registers exactly the supplied A6, B5, and B6 families, each with the four
facing/chirality suffixes. It uses the XML-configured RSU count/widths,
BLEC/BREC, LEC/REC, shared bridge FPC and AncASIC dimensions, and 5/6-RSU
package lengths. It does not invent an A5 family.

For the new modules, local +z is from external epoxy toward silicon. The
component order is two 2.5 mm-wide epoxy glue lines under the 32/35 mm
carbon baseplate (80 µm), carbon baseplate (150 µm), silicon-side film
(60 µm), and silicon including the RSU pattern and LEC/REC (50 µm). The
bridge FPC (80 µm Kapton + 30 µm aluminium) and AncASIC (80 µm existing
glue + 300 µm chip) retain their prior editable placeholder dimensions.
They extend the enclosing air volume to 830 µm total z thickness, although
these end components are localized rather than full-area layers. The film
is placed under the RSU chain and LEC/REC, not repurposed as an FPC glue;
the slide does not specify a new bridge-FPC glue layer.

The RSU chain begins at `-package_length/2 + BLEC` for left chirality and
at `-package_length/2 + BREC` for right chirality, with LEC/REC swapped
accordingly. The existing detailed silicon subdivision is unchanged:
18.394 mm sensitive y and 21.426 mm sensitive x per RSU, yielding
1,970.549220 / 2,364.659064 mm² for 5/6-RSU modules. The physical
baseplate-center to sensor-center separation in local +z is
`150/2 + 60 + 50/2 = 160 µm` before a facing transform.

Commands: `python3 logs/validate_svt_module_builder_step4b.py`,
`python3 logs/validate_svt_module_params_step4a.py`, and `git diff --check`.
All passed. The new read-only checker compared all 8 all-six-RSU and 12
mixed-catalogue types with derived footprints, chip clearances, epoxy
widths, material-stack thicknesses, and sensitive-area formulas. No source
CSV or placement reader was changed. The current XML still selects the
legacy placement file; the new catalogue types are registered but not yet
placed by the existing CSV schema.

Important verification limit: no C++ compile or DD4hep export was run at
this gate. The configured compiler `/opt/local/bin/clang++` is absent on
this host, and the prescribed `~/eic_dir/eic-shell` launcher is absent.
The four unique source values of `z_sensor_mm - z_baseplate_mm` are
−0.340, −0.030, +0.030, and +0.340 mm; these are **not** the local 0.160 mm
center-to-center separation and must not be equated without the CSV's
reference-surface/facing transform in step 5. Runtime build and overlap
checks remain required before this implementation is considered validated.

**Compilation handoff update (user, 2026-10-04):** The stable EIC image is
`26.09.0-stable` under `/cvmfs/singularity.opensciencegrid.org/eicweb`.
The user's launcher is `~/weic/eic-shell`, not `~/eic_dir/eic-shell`.
The user offered to run the compile check; the agent was not authorized to
launch the shell. `logs/compile_svt_module_step4b.sh` records a focused
out-of-tree CMake build of the `epic` target and writes its output to
`logs/compile_svt_module_step4b.txt`. Run it **inside** the stable shell.
At the time of handoff, compile status was pending; no pass was claimed.
This update changed no geometry, source CSV, or placement.

**Compilation result (user run, 2026-10-05 01:47 UTC): PASS.** The saved
`logs/compile_svt_module_step4b.txt` shows CMake configuration with Clang
20.1.8, compilation of `src/SiEndcapModuleTracker_geo.cpp`, successful
linking of `lib/libepic.so`, `[100%] Built target epic`, and the script's
`RESULT: PASS`. The isolated build directory was `build/rsu_opt_step4b/`.
CMake warned that `EPIC_ECCE_LEGACY_COMPAT` was not used by this project;
this did not affect compilation. This closes the C++ compile check for
step 4b. DD4hep geometry construction/export, actual five-RSU placement,
reference-surface z alignment, and overlap checks remain unverified and are
separate later gates.

## Step 5a: placement-coordinate audit (review checkpoint)

Files changed: `logs/analyze_svt_coordinates_step5.py` (read-only audit),
this log, and `logs/plan_rsu_opt_endcap_20261004.md` (split step 5 into
audited review gates and update the stable image to 26.09.0). No geometry,
placement CSV, catalogue, or metadata file changed. Command:
`python3 logs/analyze_svt_coordinates_step5.py`.

The supplied `x_origin_mm` is the LEC–first-RSU boundary, not a module
center. Across both packages the module envelope extends toward `x = 0`:
`x_center_global = x_origin_mm - sign(x_origin_mm) ×
(module_length_mm/2 − BLEC_mm)`. The offset magnitude is 49.165 mm for
five RSUs or 59.998 mm for six RSUs. This radial-inward direction gives
**zero** module-corner outer-radius violations for all 2,160 `all_6rsu`
and 2,164 `rsu_opt` placements, using the supplied per-disk outer-radius
metadata. Reversing the direction produces 1,528 and 1,518 violations,
respectively. No source x origin is zero. This is strong evidence for the
global x footprint, but does not alone determine local-y chirality or z.
The simple same-row/same-sensor-z projected-x check found 10 intervals
overlapping for `all_6rsu` and zero for `rsu_opt`; it is **not** a DD4hep
volume-overlap test and does not apply the beampipe openings.

The staged filenames map disk IDs to current XML keys as follows:
0–4 are ED4, ED3, ED2, ED1, ED0 (negative-z outer, outer, outer,
middle, inner layers), and 5–9 are HD0, HD1, HD2, HD3b, HD4
(positive-z inner, middle, outer, outer, outer layers). The HD3 file is
named `HD3b_modules.csv`; it must not be guessed as `HD3_modules.csv`.

The z mapping is not yet specified tightly enough to place solids. The
four source `z_sensor_mm − z_baseplate_mm` values are ±0.030 and ±0.340 mm,
whereas the modeled carbon-center to silicon-center displacement is
0.160 mm along the module normal. The slide labels module back/front
faces at ±2 mm and the CSV baseplate values are at disk center ±3 mm.
Those numbers must be related through named physical surfaces (and
possibly facing-dependent definitions); treating either CSV z field as
the corresponding solid center would silently misplace silicon.
An explicit source definition of both z reference surfaces has been
requested from the user. The strict C++ directory parser, full 3D
transform, scenario selection, and geometry export are deferred until
that convention is resolved. This checkpoint is a measured partial
result, not completion of step 5.

**Provider clarification (received 2026-10-05):** Slide 2 defines the
module coordinate system and each `*_modules.csv` row gives that module
origin in detector coordinates. The CSV field `z_baseplate_mm` is a
misleading name: physically it is the module's corrugation/facet reference
surface, not the center of the 150 µm module baseplate and not the
corrugated CF frame mid-surface. Preserve the source header for provenance,
but map it internally to `z_corrugation_surface`.

The provider also confirmed the slide-3/4 stack convention. For an
outward-facing module, the sensor reference surface differs from the
corrugation surface by the complete stack:
`80 µm epoxy + 150 µm baseplate + 60 µm film + 50 µm sensor = 340 µm`.
For an inward-facing module, the sensor intrudes into the corrugation
channel by the sensor thickness plus the adhesive difference; the
baseplate does not enter:
`50 µm sensor + 60 µm film − 80 µm epoxy = 30 µm`.
This exactly explains the four signed source differences ±0.340 and
±0.030 mm for both placement packages. The previously computed 0.160 mm
baseplate-center-to-sensor-center separation is still correct within the
solid model, but is not the CSV reference-surface separation.

Step 5b may therefore implement the transform using the midpoint of the
LEC–first-RSU boundary for x/y and the corrugation reference surface for z,
then validate the modeled sensor reference surface against `z_sensor_mm`.
The global sign is taken directly from the supplied pair of z coordinates,
so it remains correct on either disk side and corrugation facet.

## Step 5b: strict package reader and coordinate transform

Files changed after commit `9afa33cf5`: `src/SiEndcapModuleTracker_geo.cpp`,
`logs/validate_svt_placement_transform_step5b.py`,
`logs/compile_svt_placement_step5b.sh`, and this log. Production XML and
all supplied CSV/metadata files remain unchanged.

The plugin now accepts `format="tiling-directory"` with the existing
`file` attribute pointing to either self-contained scenario directory.
It reads and validates `catalog.csv` plus the ten exact disk files,
including `HD3b_modules.csv`. Missing headers/files, catalogue/template
dimension disagreements, unknown types, duplicate `(disk,row,module)`
keys, ambiguous zero directions, and inconsistent z offsets are fatal.
The legacy single-CSV reader remains available and unchanged in behavior.

Disk IDs map to the existing XML layer keys as documented in step 5a.
Each row retains `disk_id`, `row_index`, `mod_index`, and `type_id` in its
DD4hep module name and `VariantParameters`. The supplied x/y point is used
as the **midpoint** of the LEC–first-RSU boundary. The module extends
toward x=0; rotations of 0 or 180 degrees about local z/y handle radial
direction, chirality, sensor normal, and the negative-side layer reflection.

The z transform renames `z_baseplate_mm` internally to
`z_corrugation_surface`. It anchors the physical corrugation interface and
checks the modeled sensor reference against `z_sensor_mm`. Following the
provider clarification, outward epoxy is serially below the baseplate and
the local sensor-reference displacement is 0.340 mm. For inward modules,
the two epoxy edge lines lie above the baseplate edge rather than adding a
serial layer; film plus sensor protrude 0.030 mm past the epoxy contact
surface. Consequently the enclosing module thickness is 0.830 mm outward
and 0.750 mm inward with the current bridge-FPC/AncASIC placeholders.

Command: `python3 logs/validate_svt_placement_transform_step5b.py`.
Result: **PASS** for all 2,160 `all_6rsu` modules and all 2,164 `rsu_opt`
modules. Every LEC-boundary midpoint, corrugation z reference, and sensor z
reference reconstructs exactly in detector coordinates. All four required
rotation combinations are exercised in both packages. `git diff --check`
and the compile-script shell syntax check also pass.

C++ compilation: **PASS** (user run, 2026-10-05 19:38 UTC). The recorded
`logs/compile_svt_placement_step5b.txt` shows Clang 20.1.8 compiling
`src/SiEndcapModuleTracker_geo.cpp`, linking `lib/libepic.so`, and
`[100%] Built target epic`. The command is preserved in
`logs/compile_svt_placement_step5b.sh` for repetition.

This closes the reader/transform compile gate. The production XML still
selects the legacy CSV, so neither new scenario is active yet. Geometry
export, per-disk metadata boundaries, complete placement counts in DD4hep,
and overlap checks remain later review gates.

## Step 6a: supplied-versus-current disk-boundary audit

This is a read-only review checkpoint before changing geometry behavior.
Files added/changed: `logs/audit_svt_disk_boundaries_step6a.py` and this log.
Command: `python3 logs/audit_svt_disk_boundaries_step6a.py`. The script reads
the staged metadata without modifying it and compares it with resolved values
from `compact/tracking/definitions_craterlake.xml`.

The `all_6rsu` and `rsu_opt` packages have identical z centers, outer radii,
and opening primitives on every disk. Scenario selection therefore changes
module types and placements, not disk boundaries. The supplied disk centers
are ED4/ED3/ED2/ED1/ED0 = -1020/-850/-650/-450/-250 mm and
HD0/HD1/HD2/HD3b/HD4 = 250/450/700/950/1200 mm. Relative to the current XML,
ED1 and HD1 move 100 mm toward the interaction point and ED0/HD0 move 5 mm;
the other six centers agree.

The supplied outer radii are 230 mm for ED0/HD0, 405 mm for ED1/HD1, and
430 mm for all outer disks. The current XML envelope radii resolve to 240,
415, and 421.4 mm respectively. Thus the supplied inner and middle disks are
10 mm smaller, while the supplied outer disks are 8.6 mm larger. The XML also
adds a separate 0.001 mm layer-envelope allowance; that implementation detail
is not included in the comparison above.

The opening definitions differ materially. The current XML opening radii
include `Beampipe_bakeout_buffer = 5 mm`, while metadata reports the source
primitives directly. On ED0/ED1/HD0/HD1, metadata uses duplicated concentric
31.750 mm circles, compared with the XML's 36.757 mm buffered circles. ED2
introduces an offset second 15.900 mm circle that does not exist in the
current XML. ED3 and ED4 retain two-circle openings but change their radii
and centers. HD2/HD3b/HD4 each supply two identical circles centered at
+3.136/+8.736/+14.336 mm with radii 37.136/42.736/48.336 mm, rather than the
current two distinct XML circles. These are topology/coordinate changes, not
rounding effects.

The metadata coordinates are detector-coordinate source values, whereas the
existing `<beampipe_opening>` constants are consumed in a layer-local frame.
Step 6b must explicitly transform source opening centers through the same
negative-side layer reflection used for module placements. It must also make
an explicit policy decision about the XML-only 5 mm bakeout buffer: silently
adding it to source metadata would no longer implement the supplied boundary
exactly. No such choice or geometry modification is made in this audit.

## Step 6b: defer supplied disk-z movement

User decision: consume the delivered CSVs unchanged but retain the current XML
disk z centers until the proposed positions have been checked for conflicts
with other detector systems. This is an explicit temporary compatibility
translation, not a reinterpretation or edit of the source coordinates.

For each tiling-directory row, the plugin now reads the corresponding
`*_metadata.txt` and calculates
`disk_z_shift = current_XML_global_layer_center - metadata_z_center_mm`.
It adds this same shift to the source corrugation and sensor z references
before constructing the module transform. Therefore their signed 0.030 or
0.340 mm separation and the complete material-stack placement are unchanged.
The shifts are ED4/ED3/ED2 = 0, ED1 = -100 mm, ED0 = -5 mm, HD0 = +5 mm,
HD1 = +100 mm, and HD2/HD3b/HD4 = 0.

Every placed tiling module records `source_disk_center_z_mm`,
`source_corrugation_surface_z_mm`, `source_sensor_reference_z_mm`, and
`applied_disk_z_shift_mm` as `VariantParameters`. This preserves both the
provider coordinates and the temporary translation in the constructed
geometry. The production XML still selects the legacy placement file.

Validation command: `python3 logs/validate_svt_deferred_z_step6b.py`. It checks
both packages and all 4,324 rows, confirming that translation preserves the
source sensor-to-corrugation separation and every module's offset from its
source disk center. Moving to the supplied disk z positions later requires
only making the XML layer centers agree with metadata; the derived shifts then
become zero without changing a CSV or a C++ shift table.

Compile handoff: `logs/compile_svt_deferred_z_step6b.sh` records a focused
out-of-tree build and writes `logs/compile_svt_deferred_z_step6b.txt`. Run the
script inside `~/weic/eic-shell --version 26.09.0-stable`.

**Compilation result (user run, 2026-10-05): PASS.** The saved output shows
`src/SiEndcapModuleTracker_geo.cpp` compiled successfully, `lib/libepic.so`
linked, `[100%] Built target epic`, and the script's final `RESULT: PASS`.
