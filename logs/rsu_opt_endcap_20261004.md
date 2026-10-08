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

## Step 6e: remove the temporary disk-z translation

Commit `e9706b826` closed the reviewed upstream z-alignment checkpoint. Since
all ten XML disk centers now equal the source metadata, the user requested
removal of the temporary translation to avoid a misleading configuration path.

The module row no longer stores a source disk center or applied z shift. CSV
`z_baseplate_mm` (internally the corrugation reference) and `z_sensor_mm` are
used directly in the module transform, and the `source_disk_center_z_mm` and
`applied_disk_z_shift_mm` module parameters have been removed. The two useful
source reference values remain recorded on each module.

Instead of translating a mismatch, tiling-directory construction now requires
the reflected/global XML layer center to equal metadata `z_center_mm` within
0.00005 mm and throws a named disk error otherwise. This makes any future
layout/XML divergence explicit. The historical Step 6b validator was updated
in place to require exact center agreement and direct use of all source z
references; its filename is retained so earlier log commands remain valid.

Validation commands:
`python3 logs/validate_svt_deferred_z_step6b.py`,
`python3 logs/validate_svt_placement_transform_step5b.py`, and
`python3 logs/audit_svt_radial_boundaries_step6c.py`.
Compile handoff: run `bash logs/compile_svt_direct_z_step6e.sh` inside
`~/weic/eic-shell --version 26.09.0-stable`. It writes
`logs/compile_svt_direct_z_step6e.txt`.

**Compilation result (user run, 2026-10-05): PASS.** The saved output shows
`src/SiEndcapModuleTracker_geo.cpp` compiled successfully, `lib/libepic.so`
linked, `[100%] Built target epic`, and the script's final `RESULT: PASS`.

## Step 6c: supplied radial boundaries

Commit `e474f8c7c` closed the reviewed Step 6a/6b checkpoint. Step 6c changes
the tiling-directory path only: legacy CSV geometry continues using its XML
boundaries, and production XML still selects that legacy path.

The strict package reader now requires and parses each metadata file's module
outer-radius target and two opening primitives. It accepts the delivered
`a_mm == b_mm` circular primitives and rejects missing fields, non-circles,
non-positive dimensions, wrong disk IDs, and unsupported radius targets. For
an active tiling directory, these values replace the layer's XML radial
boundary and opening definitions. The supplied centers are detector-global;
their x coordinates are sign-flipped for reflected negative-z layer solids so
the placed subtraction returns to the source global coordinate. No 5 mm XML
bakeout buffer is added: the placements already have a minimum opening
clearance of 4.99995 mm in the four-decimal source data.

The nominal `outer_radius_mm` explicitly targets module corners. Because the
CSV coordinates are rounded to 0.0001 mm, tangent corners extend only
0.0000286--0.0000558 mm beyond the nominal radius. The constructed layer
solid retains the existing 0.001 mm envelope allowance, i.e.
`solid_rmax = metadata_outer_radius + 0.001 mm`; this is not folded into or
reported as the design radius. Layer `VariantParameters` record the source z
center, nominal outer radius, and applied allowance.

Command: `python3 logs/audit_svt_radial_boundaries_step6c.py`. Result: PASS
for all 2,160 `all_6rsu` and 2,164 `rsu_opt` full baseplate footprints. There
are 384/386 nominal-radius rounding excursions respectively, but zero exceed
the 0.001 mm layer allowance and no footprint intersects an opening. Sensor
and film footprints are strict subsets of the checked full baseplates. The
deferred-z regression check also remains PASS for every row.

Compile handoff: run `bash logs/compile_svt_radial_boundaries_step6c.sh`
inside `~/weic/eic-shell --version 26.09.0-stable`. It writes
`logs/compile_svt_radial_boundaries_step6c.txt`.

**Compilation result (user run, 2026-10-05): PASS.** The saved output shows
`src/SiEndcapModuleTracker_geo.cpp` compiled successfully, `lib/libepic.so`
linked, `[100%] Built target epic`, and the script's final `RESULT: PASS`.
DD4hep construction/export remains a separate later validation gate.

## Step 6d: restore upstream-aligned disk z positions

Before proceeding to support geometry, `origin/main` was fetched and compared
with this branch. The branch is already based on the fetched main commit
`f9617538a` (51 commits ahead, 0 behind), so merge and normal rebase operations
would be no-ops. Two later branch changes had nevertheless moved the inner
disk center from 250 to 255 mm and the middle disk center from 450 to 550 mm.

This focused forward fix restores the current upstream/main values:
`InnerTrackerEndcapP_zmin = 25.0 cm` and
`TrackerEndcapPDisk1_zmin = 45.0 cm`. Their negative-side aliases therefore
restore ED0/HD0 to -250/+250 mm and ED1/HD1 to -450/+450 mm. All ten XML disk
centers now equal the supplied metadata centers exactly. No radial boundary,
opening, layer-thickness, module, or source-data value changes in this step.

The generic C++ translation remains deliberately in place, but its derived
`applied_disk_z_shift_mm` is now zero on every disk. The rerunnable Step 6a and
6b audit tables were updated to the current XML state, and the z validator now
requires every derived shift to be zero. Earlier log sections remain as the
historical record of why the temporary shift mechanism was introduced; this
section supersedes their nonzero-shift status.

Validation commands:
`python3 logs/audit_svt_disk_boundaries_step6a.py`,
`python3 logs/validate_svt_deferred_z_step6b.py`, and
`python3 logs/audit_svt_radial_boundaries_step6c.py`.
Compile handoff: run `bash logs/compile_svt_upstream_z_alignment_step6d.sh`
inside `~/weic/eic-shell --version 26.09.0-stable`; output is written to
`logs/compile_svt_upstream_z_alignment_step6d.txt`.

**Compilation result (user run, 2026-10-05): PASS.** The saved output shows
`src/SiEndcapModuleTracker_geo.cpp` compiled successfully, `lib/libepic.so`
linked, `[100%] Built target epic`, and the script's final `RESULT: PASS`.

## Step 7a: corrugated-support source audit

This is a read-only checkpoint. Files added/changed:
`logs/audit_svt_corrugation_rows_step7a.py` and this log. No geometry, XML,
catalogue, placement, or metadata value changed. The PDF corrugation and
beampipe-pitch slides were rendered locally for visual inspection because the
login node lacked `pdftotext`.

The slide defines two physical profiles: 33/36 mm pitch, 6 mm surface-to-
surface depth, 45-degree webs, 5 mm mid-surface bend radius, 0.15 mm skin,
and zero-skin flat widths 6.36/7.86 mm. Their inward-module keep-clear
envelopes are 20.56 x 1.00 and 21.94 x 1.00 mm respectively. The pitch slide
calls uniform 33 mm reasonable for ED0/ED1/HD0/HD1 and says 36 mm is a
reasonable starting point for the remaining disks; it does not itself define
disk-specific transition coordinates.

The current support is not either slide profile. XML supplies a uniform
34.77 mm pitch, 35-degree web, 6 mm depth, and 0.15 mm skin. C++ constructs
straight flat/web boxes with no 5 mm bends or explicit keep-clear profile.
The shared row file repeats this same 34.77 mm profile on every disk and is
therefore not a source for the delivered design.

Command: `python3 logs/audit_svt_corrugation_rows_step7a.py`. Both placement
packages encode the same support-row families and corrugation surfaces; only
module counts/types differ. ED0, ED1, HD0, and HD1 are unambiguous uniform
33 mm lattices: adjacent opposite-surface module rows are separated by
16.5 mm, same-surface rows by 33 mm, and their corrugation references alternate
exactly at disk center +/-3 mm.

The delivered outer disks are not uniform 36 mm module-row layouts. They use
central 35 mm-wide `A` baseplate rows and outer 32 mm-wide `B` rows. In
`rsu_opt`, the row-index regions are:

| Disk | Central A rows | Outer B rows |
|---|---:|---:|
| ED4 | 11-35 | 0-10, 36-47 |
| ED3 | 15-33 | 0-14, 34-47 |
| ED2 | 17-31 | 0-16, 32-47 |
| HD2 | 17-31 | 0-16, 32-47 |
| HD3b | 12-34 | 0-11, 35-46 |
| HD4 | 8-38 | 0-7, 39-46 |

This confirms central 36 mm and outer 33 mm *nominal support families* in the
provided design, despite the slide's simpler recommendation language. The
measured module y spacings are not exact support pitches in these disks:
same-surface separations range roughly 33-37 mm around the A/B transitions.
This agrees with the documented source-algorithm behavior that module y origins
may move within the keep-clear envelope. Consequently module centers are valid
for identifying row family, face, and ordering, but must not be copied as
corrugation bend/flat centerlines.

What is fully specified for Step 7b: the uniform 33 mm profile and its absolute
row lattice on ED0/ED1/HD0/HD1. What remains absent from the delivered CSV and
metadata: explicit outer-disk corrugation centerline coordinates through the
33-to-36 mm transition, including which curved segment absorbs the phase
change. Step 7b can therefore implement and validate the four uniform disks
without inventing transition geometry. Exact outer-disk construction should
wait for a provider row/transition definition or be labeled as a documented
reconstruction rather than source-exact.

**Additional source audit (`ES_material/20261001/`, user request):** The user
confirmed that `A` rows map to the 36 mm support profile and `B` rows to the
33 mm profile. `ES_material/20261001/README.docx` identifies three studies:
baseline = 34.77 mm with centered modules, transition = 33/36 mm with modules
allowed to shift according to projectivity, and projective = a hypothetical
unbuildable maximum-acceptance layout. It explicitly instructs users to take
`transition/rsu6` or `transition/rsu_opt`; the staged preferred package is the
transition `rsu_opt` design.

The transition delivery contains catalogues, module coordinates, metadata,
tiling SVGs, and row/module-gap and acceptance reports. It contains no source
code, configuration, support-row table, corrugation-center field, phase anchor,
or 33/36 join rule. Inspection of the SVG XML confirms that its rectangles and
row labels are generated from module/baseplate coordinates; it does not draw or
encode the underlying corrugation centerline. `row_gap_report.txt` likewise
reports module-row separations. Thus the 20261001 material confirms the pitch
families and confirms that modules may shift away from channel centers, but it
does not supply the two missing construction definitions: an absolute support
phase per outer disk and the geometry of each unlike-pitch half-cell/join.

## Step 7b: reconstruct the supplied 33/36 mm flat-center lattice

The reviewed A/B row table plus the agreed centered join convention closes the
two apparent gaps above. Type A denotes the 35 mm baseplate on the 36 mm
same-surface corrugation pitch; type B denotes the 32 mm baseplate on the 33 mm
pitch. Since consecutive row indices alternate corrugation surfaces, adjacent
flat-center separations are 18.0 mm within A, 16.5 mm within B, and the centered
A/B join is `(18.0 + 16.5)/2 = 17.25 mm`.

`logs/generate_svt_corrugation_lattice_step7b.py` reads the staged
`rsu_opt` placement package only to obtain each disk's ordered row indices,
A/B family, and alternating source surface. It anchors an odd row count at
`y=0`; for an even row count it centers the middle pair about `y=0`, then
propagates the three allowed separations outward. It writes the derived,
human-readable flat-center table
`compact/tracking/silicon_disks/SVT_endcap_corrugation_rows_33_36.csv`. The source placement
CSVs remain unchanged.

ED4 demonstrates that this is the intended absolute phase rather than an
arbitrary choice. Its central A rows 23/24 land at `+9/-9 mm`; propagating the
A lattice and centered joins puts B rows 10 and 36 at `+242.25/-224.25 mm`
and the outer rows at `+407.25/-405.75 mm`, exactly matching the delivered
unshifted B-module rows. The same check holds throughout the delivery except
for a 0.0391 mm module displacement at the lower A/B transition of HD3b and
HD4, consistent with the transition study's explicit allowance for projective
module motion within a channel.

Validation command:
`python3 logs/validate_svt_corrugation_lattice_step7b.py`. It requires every
adjacent separation to follow the A/B formula, every family to match the
source row, centered phase, alternating source surfaces, and B-module residuals
below 0.05 mm. This checkpoint defines row flat centers only; it deliberately
does not yet connect the table to the C++ support builder, because that builder's
older CSV describes whole symmetric cells rather than explicit alternating
flat centers. The next small implementation step is to teach the builder this
explicit-row representation and validate its generated support pieces.

## Step 7c: explicit-flat corrugation builder

Commit `fd2eb3c5b` closes the reviewed Step 7a/7b audit and reconstructed
lattice checkpoint. Step 7c adds an independent explicit-flat input path to
`src/SiEndcapModuleTracker_geo.cpp`; the legacy symmetric whole-cell loader and
builder remain available without behavior changes.

The new loader recognizes the columns `disk`, `row_index`, `row_y_mm`,
`surface_offset_mm`, and `pitch_mm`. It loads both positive- and negative-y
flat centers for the requested disk and does not mirror them. The builder
places each flat once and connects each adjacent pair with one straight carbon
web. Module placement is not involved, so the supplied projective module
offsets remain unchanged.

The explicit offsets are physical corrugation surfaces at +/-3 mm. For the
0.15 mm carbon skin, flat-volume centers are therefore at +/-2.925 mm. The
sharp-junction simulation approximation omits the slide's 5 mm bend radius but
preserves the 45-degree carbon mid-surface web: its z and y spans are both
`6.0 - 0.15 = 5.85 mm`. Consequently the connected flat widths are 10.65 mm
for the 33 mm profile and 12.15 mm for the 36 mm profile. At a centered
17.25 mm A/B join, half of those two flat widths plus the 5.85 mm web span
again equals the supplied flat-center separation exactly.

`python3 logs/validate_svt_explicit_corrugation_builder_step7c.py` passes for
all ten disks. It verifies one web per adjacent pair, positive flat widths,
continuous joins, alternating +/-2.925 mm flat-volume centers, and 45-degree
webs. The recorded compilation command is
`bash logs/compile_svt_explicit_corrugation_step7c.sh`; it must be run inside
the user's configured EIC shell. Production XML still points to the legacy
34.77 mm whole-cell table, so this checkpoint does not yet activate or alter
the detector geometry.

The derived table was subsequently relocated to
`compact/tracking/silicon_disks/SVT_endcap_corrugation_rows_33_36.csv` so the
corrugation input remains grouped with the two silicon-disk placement
packages. The generator and both validators use the organized path.

**Compilation result (user run, 2026-10-06): PASS.** The saved output in
`logs/compile_svt_explicit_corrugation_step7c.txt` shows
`src/SiEndcapModuleTracker_geo.cpp` compiled, `lib/libepic.so` linked,
`[100%] Built target epic`, and `RESULT: PASS`.

**Follow-up TODO:** request an authoritative corrugation-row table with future
placement deliveries, including unpopulated support rows. Also promote the
current reconstruction generator from `logs/` to a maintained `scripts/`
helper with explicit package, output, pitch-map, join, and centering inputs.
Until a provider table is available, keep the checked-in table identified as
a documented reconstruction and compare the two tables when it arrives.

## Step 7d: atomically activate `rsu_opt` modules and support

The production fragment now selects
`compact/tracking/silicon_disks/rsu_opt` with
`format="tiling-directory"` in all six silicon-endcap detector blocks. All ten
corrugated layer blocks select the organized explicit-flat table, and the
corrugation angle changes from the legacy 35 degrees to the slide value of
45 degrees. Module placements and their projective offsets remain byte-for-byte
source values. Switching to the supplied all-six-RSU scenario requires changing
the shared package name in the six detector declarations to `all_6rsu`; the
metadata boundaries and reconstructed support are common to both deliveries.

The activation review caught a coordinate issue that was invisible for the
old z-symmetric support. CSV `surface_offset_mm` is detector-global, whereas
ED detectors are built in a reflected local assembly. The explicit builder now
negates the surface offset before placing ED carbon volumes, so the final
global corrugation face remains the supplied +/-3 mm face. HD offsets pass
through unchanged. This transformation affects only support coordinates;
module transforms already handled reflection independently.

Validation command:
`python3 logs/validate_svt_scenario_activation_step7d.py`. Result: PASS. It
requires exactly six `rsu_opt` package selectors, ten explicit support-table
selectors, the 45-degree constant, all ten placement/metadata file pairs, the
catalogue, and complete corrugation disk coverage. The Step 7b and Step 7c
validators also remain PASS, including the reflected/global surface recovery.
The Step 7b validator now checks both selectable packages permanently:
`rsu_opt` matches all 426 populated rows, while `all_6rsu` matches all 424
populated rows and deliberately leaves only ED4 support rows 0 and 47 empty.
There are no A/B-family or +/-3 mm surface mismatches. Both packages preserve
the same maximum 0.7891 mm projective module shift and 0.0391 mm B-row shift.

The full geometry gate is recorded in
`logs/check_svt_rsu_opt_geometry_step7d.sh`. Run it inside the configured EIC
shell; it configures/builds/installs through `setup_local.sh`, exports
`epic_craterlake_tracking_only_disk_modules.xml`, and runs both overlap checks.
The production activation remains uncommitted until that output is reviewed.

**First geometry run (2026-10-06): invalid support fallback.** The input CSV
was installed, but the support loader searched its `disk` column with XML layer
keys such as `OuterTrackerEndcapN_disk4`; the table contains source keys such
as `ED4`. Every explicit lookup therefore returned empty, emitted missing
legacy `h_mm`/`theta_deg` warnings, and constructed `frame_uniform...` fallback
pieces. The resulting support-related overlap list does not test the intended
33/36 mm geometry.

The loader now converts all ten XML layer keys to the same ED/HD source keys
used by the tiling-package loader before reading explicit support rows. The
activation validator permanently checks all ten mappings. The geometry script
now reports `RESULT: COMPLETE` rather than `PASS`, because successful overlap
commands may still report physical intersections. A second full geometry run
is required; its log should contain explicit `frame_flat...`/`frame_web...`
pieces and no fallback warnings before its overlaps are interpreted.

## Step 7e: classify the valid explicit-support overlap run

The second Step 7d run loaded explicit flats/webs with zero fallback records.
`python3 logs/analyze_svt_overlaps_step7e.py` parses only the concise
`scripts/checkOverlaps.py` records, avoiding Geant4's repeated bidirectional
warnings. Its reported overlap values are centimetres and are converted to mm.
The 4,454 records divide into only three classes:

- 1,675 module-mother/module-mother overlaps, all exactly 0.230 mm;
- 1,059 module-mother/flat overlaps, 0.370--0.374886 mm;
- 1,720 module-mother/web overlaps, 0.358033--0.358747 mm.

Every concise record names a rectangular air module mother (`A...` or `B...`).
None names an actual sensor, baseplate, adhesive, bridge-FPC, AncASIC, or other
module child solid. A direct source-coordinate reconstruction reproduces the
0.230 mm module-mother overlaps along z: that depth is exactly the modeled
150 um baseplate plus 80 um external-epoxy reference inset. The approximately
0.37 mm support intersections likewise reflect the empty rectangular mother
envelope spanning across the physical corrugation contact. These enclosing air
boxes overlap even though the real component shapes need not do so, making the
current hierarchy invalid for navigation and hiding the component-level test.

## Step 7f: remove artificial module air envelopes

The tiling module prototype is now a reusable DD4hep `Assembly` containing the
same component solids and sensitive surfaces, instead of a full package-sized
air `Box`. This follows existing repository module-assembly patterns and does
not alter any component size, material, local coordinate, global transform,
source metadata, or sensitive-surface attachment. It removes only the
artificial solid whose overlaps dominated Step 7e, allowing the next geometry
check to test the actual components against each other and the support.

Lightweight Step 7c/7d validators and `git diff --check` pass. The full rerun is
recorded in `logs/check_svt_component_geometry_step7f.sh` and writes a separate
`logs/check_svt_component_geometry_step7f.txt`, preserving the Step 7d evidence.

**Step 7f geometry result (2026-10-06): component overlaps remain.** The air
mother-box overlaps disappear, exposing 1,092 unique component-level pairs
from 1,096 Geant4 reports. Command:
`python3 logs/analyze_svt_component_overlaps_step7f.py`.

The support intersections occur exclusively on inward-facing module types:

- 388 web/AncASIC pairs, with measured penetration up to 0.149982 mm;
- 271 flat/sensor-bias pairs;
- 135 flat/powerswitch-strip pairs;
- 88 flat/backbone-strip pairs;
- 26 web/**active-silicon** pairs;
- 10 web/bridge-FPC pairs, plus eight other web/inactive-silicon pairs.

There are also 147 inward/outward external-epoxy intersections and ten
inward/outward baseplate intersections between neighboring module assemblies.
Nine independent beampipe/layer-envelope reports remain, with the largest
reported penetration 4.90662 mm. No legacy corrugation fallback is present.

The support result invalidates the connected sharp-corner approximation as
currently constructed. The slide's nominal zero-skin straight-flat widths are
6.36/7.86 mm, whereas forcing straight 45-degree webs to meet flats without
the omitted 5 mm bends enlarged them to 10.65/12.15 mm. That additional carbon
occupies part of the documented inward-module keep-clear region. Because active
silicon and electronics are affected, do not suppress these reports or move
source modules heuristically. The next support revision must use the slide's
keep-clear geometry—either explicit rounded bends or a documented conservative
gap/cutout approximation—then rerun component overlaps. The module-component
and beampipe categories should be addressed separately after support clearance.

## Step 7g: conservative approximation of the omitted 5 mm bends

The explicit table now records the slide's nominal zero-skin straight-flat
width per row: 6.36 mm for the 33 mm/B profile and 7.86 mm for the 36 mm/A
profile. The loader requires and uses `flat_width_mm` directly rather than
enlarging flats to force a sharp connection.

For every B/B, A/A, and centered A/B interval, the space between facing flat
edges is 10.14 mm. A 45-degree straight carbon mid-surface web spanning the
5.85 mm separation between the 0.15 mm skin midplanes uses 5.85 mm in y. It is
centered in the available interval, leaving 2.145 mm on each side for the
omitted curved bend. The finite rotated skin projection reduces the empty
solid-to-solid gap slightly but leaves it positive. This deliberately
underestimates bend material rather than putting invented carbon inside the
documented module-clearance region.

The generator, checked-in CSV, loader, and analytical validator all carry the
explicit flat widths. Step 7b verifies both placement packages unchanged;
Step 7c verifies 6.36/7.86 mm flats, 45-degree webs, constant 10.14 mm
flat-edge intervals, and 2.145 mm omitted-bend gaps for all disks and A/B
joins. All lightweight checks pass. The full rerun is recorded separately in
`logs/check_svt_conservative_bends_step7g.sh` and writes
`logs/check_svt_conservative_bends_step7g.txt`.

**Step 7g geometry result (2026-10-06): the conservative bend gaps do not
resolve the component overlaps.** The rerun contains 1,100 unique pairs from
1,103 Geant4 reports, compared with 1,092 unique pairs from 1,096 reports in
Step 7f. The principal classes are 404 web/inward-AncASIC, 340
flat/inward-sensor-bias, 99 flat/inward-powerswitch, 45
flat/inward-backbone, and 18 web/inward-active-silicon pairs. There are also
162 inward/outward external-epoxy and nine inward/outward baseplate pairs,
plus ten beampipe/layer-envelope reports. No fallback support was built.

The support intersections are distributed across all ten disks (54 pairs per
HD0 side and approximately 99--103 per other disk side), so they are not an
isolated A/B-transition defect. They are consistent with the supplied inward
module convention: for example, HD0 row 1 has a corrugation reference at
247.000 mm and an inward sensor reference at 247.030 mm. The provider states
that inward modules intrude into the corrugation channel by the sensor and
adhesive delta. The current support, however, models every flat and web as an
uninterrupted box across the complete disk x-span. It therefore leaves solid
carbon in the intended intrusion/clearance region. Changing flat width in y
cannot address that missing x-dependent clearance.

Do not move the source modules, change their z references, or suppress these
reports. The next support step must represent module openings/clearances (or
another explicitly documented material-removal approximation) in the flats
and webs. The module-to-module and beampipe/layer-envelope reports remain
independent follow-up categories.

## Step 7h: placement-aware inward keep-clear gaps

The 6 mm corrugation depth is outer-face to outer-face and already includes the
0.15 mm corrugated carbon. The correct flat intervals are therefore
`-3.000..-2.850 mm` and `+2.850..+3.000 mm`, with midplanes at +/-2.925 mm and
a 5.85 mm midplane separation. The exploratory opposite-side interpretation
was not retained.

The slide defines virtual keep-clear regions only for inward-facing modules:
20.56 x 1.00 mm for B/33 mm and 21.94 x 1.00 mm for A/36 mm. In the slide,
detector-local x is out of the page and the quoted width is along local y. The
rectangle extends from the applicable outer corrugation face by 1 mm toward
the channel. It is a required empty region, not a physical air volume; inward
module material is intended to occupy it.

The straight-flat/web simulation approximation now uses each inward module's
actual transformed package x interval and source y/surface coordinates. If a
finite-thickness flat or oblique support piece enters that module's virtual y-z
rectangle, the piece's disk-clipped x intervals are split at both module edges
and the intersecting x section is omitted. Outward modules create no gaps.
This conservatively shortens both flat and oblique carbon over the module
footprint instead of attempting the slide's 5 mm-radius contour. The A/B widths
and 1 mm depth are named XML constants consumed directly by the C++ builder.
No source module coordinate, corrugation row, or physical module component is
changed.

The recorded full rerun remains `logs/check_svt_keep_clear_step7h.sh`. Its
component overlap result must confirm that support/inward-module intersections
are removed; module/module and beampipe/layer-envelope categories are separate.

**Step 7h geometry result (2026-10-07): support clearance passes.** The build,
export, and both overlap commands completed. Parsing the Geant4 report gives
197 reports and 190 unique pairs. There are **zero** flat/module, web/module,
support/support, or fallback-support pairs. This removes every support-related
category seen in Steps 7f/7g, including active silicon, inactive sensor strips,
AncASICs, and bridge FPCs.

The remaining unique pairs are 176 inward/outward external-epoxy pairs, four
inward/outward baseplate pairs, and ten beampipe/layer-envelope reports (maximum
reported penetration 4.25735 mm). The module-pair reports are independent of
the support clipping; Geant4's sampled report count can vary between runs, so
the pair geometry—not the small count change from Step 7g—is the next item to
audit. The beampipe reports remain a separate layer-solid/opening issue.

## Step 7i: source-coordinate module-pair audit

`python3 logs/audit_svt_module_pair_overlaps_step7i.py` reconstructs the exact
axis-aligned package, baseplate, and two 2.5 mm edge-glue footprints from the
rsu_opt catalogue and source placements. It uses the already validated
LEC-boundary transform; no geometry or source data is changed.

The result shows that the Geant4 report substantially undersamples thin partial
intersections. There are 1,676 projected full-package baseplate pairs and 3,352
edge-epoxy pairs (low and high edge for each baseplate pair). Every pair is in
the same source row, has one inward- and one outward-facing module, and shares
the same corrugation reference surface. The pattern occurs on all ten disks.
Projected x-overlap lengths range from 24.0000 to 124.8780 mm; baseplates share
their full 32 or 35 mm y width, while each epoxy pair shares the full 2.5 mm
edge-strip width. Representative ED0 row 1 B511/B501 boxes overlap by
71.0688 mm in x on both glue edges.

Therefore the 176 epoxy and four baseplate pairs seen by Geant4 are not isolated
defects and must not be repaired by merging only the sampled epoxy volumes. At
least one current mechanical assumption is incomplete: the carbon baseplate
may not span the full catalogue package length, inward/outward modules may use
a different mechanical stack relationship, or the delivered placement optimizes
sensor coverage without enforcing these package solids. This needs a provider
clarification before changing module material dimensions or coordinates.

The provider SVGs narrow this ambiguity. Their unfilled `baseplate footprint`
rectangles are exactly the catalogue package dimensions: for example, ED0 B5
rectangles are 138.330 x 32.000 mm and B6 rectangles are 159.996 x 32.000 mm.
Thus the current full-package x extent is not an implementation invention. The
combined SVG explicitly draws overlapping inward/outward baseplate rectangles;
ED0 row 1 B511 and B501 reproduce the same 71.0688 mm overlap found by the
audit. Their active-area rectangles overlap in projection as well.

Each disk is also delivered as four filtered views:
`front/back x inward/outward`. The legend defines FRONT as active area bridging
ridges and BACK as active area bridging furrows; this classification is distinct
from inward/outward facing. The filtered SVGs show that projected overlaps are
intended to be evaluated in separate facing categories. However, SVG contains
only x-y rectangles and no component thickness or z interval, so it does not
specify how the mechanically overlapping baseplates/adhesive are separated in
z. The gap reports likewise group row spacing by FRONT/BACK and do not provide
an inward/outward baseplate z-clearance rule. The provider must still clarify
the 3D stack relationship; changing the x footprint is no longer a supported
resolution.

## Step 7j: provisional outward-module z separation

Pending the provider's 3D stack clarification, outward-facing modules receive
a named, editable `SiEndcapTilingOutwardModule_z_offset = 0.230 mm`. The value
is exactly the modeled 0.080 mm external epoxy plus 0.150 mm carbon baseplate.
The complete outward module is translated farther away from its corrugation
surface along its source sensor-facing normal; inward modules and all source
CSV values remain unchanged.

In a surface-relative coordinate pointing into the channel, the inward epoxy
and baseplate remain at `-0.080..0.000` and `-0.230..-0.080 mm`. The translated
outward epoxy and baseplate occupy `-0.310..-0.230` and
`-0.460..-0.310 mm`. Thus the two backing stacks have no volume overlap and
touch only at `-0.230 mm`. This is explicitly provisional: it moves every
outward sensor reference 0.230 mm away from its supplied `z_sensor_mm` and
leaves the translated module away from the nominal corrugation contact.

Module parameters retain `source_corrugation_surface_z_mm` and
`source_sensor_reference_z_mm`, and additionally record
`applied_outward_z_offset_mm` and `placed_sensor_reference_z_mm` so the
temporary assumption is visible in exported geometry. Command:
`python3 logs/validate_svt_provisional_z_offset_step7j.py`; result PASS. The
recorded compile/export/overlap rerun is
`logs/check_svt_provisional_z_offset_step7j.sh`.

**Step 7j geometry result (2026-10-08): provisional module separation
passes.** The build, export, and both overlap commands completed. The Geant4
log contains 12 reports representing ten unique pairs, all classified as
beampipe/layer-envelope intersections. There are zero module/module,
module/support, support/support, epoxy, baseplate, sensor, electronics, or
fallback-support pairs. Thus the 0.230 mm provisional outward translation
removes the backing-stack intersections exactly as intended, while preserving
the Step 7h support clearance.

This result validates the temporary simulation geometry, not the underlying
mechanical interpretation. Keep the offset visibly provisional until the
provider confirms the intended facing-dependent z relationship. The remaining
ten beampipe/layer-envelope reports are independent of module placement and
should be handled as a separate disk-solid/opening task.

## Step 7k: beampipe/layer-envelope overlap audit

The twelve Step 7j Geant4 reports collapse to ten unique pairs, all between a
10 mm-thick disk-layer air mother and a beampipe volume. They affect HD2, HD3b,
HD4, ED3, and ED4. No physical module or corrugation component appears in
these reports, so this is a boolean layer-envelope problem rather than a new
module-placement or support-clearance problem.

The audit found that `DiskBoundary::lepton_opening` and
`DiskBoundary::hadron_opening` currently serve three different purposes:

1. rejecting delivered module rectangles that cross a source opening;
2. clipping corrugated carbon-support x intervals; and
3. subtracting holes from the full 10 mm-thick layer mother.

When tiling metadata is loaded, its opening primitives overwrite the XML
beampipe openings for all three purposes. That conflates a two-dimensional
source-layout exclusion with the physical clearance required by the installed
EPIC beampipe across the complete layer thickness.

There are two independent numerical mismatches:

| Disk | delivered metadata circle `(cx, r)` [mm] | physical XML opening `(cx, r)` [mm] | observation |
|---|---:|---:|---|
| HD2 | `(+3.136, 37.136)` | `(-1.440, 40.488)` | metadata is displaced toward positive x; overlap samples are on negative x |
| HD3b | `(+8.736, 42.736)` | `(-7.066, 46.111)` | opposite center direction; largest sampled penetration is 4.9619 mm |
| HD4 | `(+14.336, 48.336)` | `(-12.692, 51.734)` | opposite center direction |
| ED3 | centered circle `(0, 34.900)` | centered circle `(0, 38.500)` | sampled point has radius 35.207 mm, outside metadata but inside XML opening |
| ED4 | centered circle `(0, 34.900)` | centered circle `(0, 38.500)` | sampled points reach radius about 38.10 mm, outside metadata but inside XML opening |

The HD metadata signs are source-exact and also appear in the provider SVGs;
they must not be silently edited. They nevertheless do not follow the current
EPIC downstream beampipe direction. On the electron side, the metadata opening
is simply smaller than the pipe-aware layer opening. The existing XML comments
state that its openings are evaluated at the disk face nearest the IP for the
10 mm layer thickness and include `Beampipe_bakeout_buffer = 5 mm`. That is the
appropriate physical-envelope convention.

Recommended next change: retain separate **source placement openings** from
metadata for `module_inside_disk`, and retain separate **physical layer
openings** from XML for `build_disk_solid` and corrugation clipping. This keeps
the supplied module selection and CSVs exact while preventing the air mother
or carbon support from occupying the beampipe clearance. Do not add 5 mm to
the metadata or change its center signs. No source data or geometry was changed
in this audit step.

### Step 7k upstream-main correction (2026-10-08)

After refreshing `origin/main` to `25fb5c1d1` (2026-10-05), comparison showed
that the two-circle XML openings above are branch-local and must not be treated
as the latest upstream dimensions. Current upstream main still uses the legacy
trap-endcap geometry with a single centered scalar `rmin`; it has no
`beampipe_opening` elements or per-pipe opening centers. Its own XML describes
these radii as enlarged approximations for the beampipe offset that need to be
fixed later.

The evaluated upstream-main design `rmin` values for the affected disks are:

| Disk | upstream-main `rmin` [mm] | expression before evaluation |
|---|---:|---|
| HD2 | 38.460 | `32.86 + 5.0 + 0.6` |
| HD3 | 50.920 | `39.29 + 5.0 + 6.63` |
| HD4 | 63.420 | `45.72 + 5.0 + 12.70` |
| ED3 | 40.061 | `34.244 + 5.0 + 0.817` |
| ED4 | 46.353 | `38.043 + 5.0 + 3.310` |

The legacy layer mother subtracts a further 1.001 mm through its layer and ACTS
envelope allowances, so its actual centered solid boundary is approximately
37.459, 49.919, 62.419, 39.060, and 45.352 mm respectively. Those dimensions
are useful conservative cross-checks, but they cannot specify the exact
two-circle opening required by the new tiled geometry.

Therefore the architectural recommendation to separate source-placement and
physical-envelope openings remains valid, but the branch-local two-circle XML
numbers are not an authoritative physical input. Before implementing that
separation, derive the physical opening from the current beampipe geometry or
adopt and clearly label an upstream-main-style conservative centered opening.
The fetched remote reference changed only `.git` metadata; no worktree geometry
or source data was changed.
