# SVT disk layouts

Each scenario is self-contained so the geometry can select one package without
mixing its catalogue, module placements, or disk metadata with another.

- `all_6rsu/`: 2,160 six-RSU modules.
- `rsu_opt/`: 265 five-RSU and 1,899 six-RSU modules.

Each directory contains `catalog.csv`, ten `*_modules.csv` placement files,
and ten `*_metadata.txt` files. Disk IDs 0–9 run from ED4 at negative z to HD4
at positive z. Coordinates and lengths are in millimetres. The catalogue's
`active_*` dimensions are nominal tiling metrics; the detailed DD4hep
sensitive silicon excludes bias, backbone, and powerswitch strips.

The package validation and SHA-256 checksums are recorded in
`logs/validate_svt_handoff_step2.txt` and `logs/stage_svt_handoff_step3.txt`.
Geometry configuration files do not select these packages yet; the scenario
switch is a later implementation step.
