# SVT disk tiling deliveries

This directory holds immutable, versioned SVT endcap tiling deliveries. The
active delivery is `20261001_ES/`.

From Ernst S. Module tiling 20261001 -
https://conferences.lbl.gov/event/2600/

Each delivery directory owns its placement scenarios, metadata, and any shared
support table derived for that delivery. Geometry XML selects an explicit
delivery and scenario path. There is intentionally no `latest` alias: changing
the selected geometry must be a visible, reviewable XML change.

New deliveries belong in parallel directories named by delivery date and
provider, for example `YYYYMMDD_provider/`. Do not overwrite an older delivery
or mix catalogues, placements, metadata, or support tables between versions.

Within a delivery, each scenario directory contains one `catalog.csv`, ten
`*_modules.csv` placement files, and ten `*_metadata.txt` verification files.
Coordinates and lengths are in millimetres unless a file states otherwise.
