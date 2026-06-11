# RSU Drawing Summary

Source: `rsu.pdf`

The PDF is a single raster drawing, so the values below were read visually from the rendered page. All dimensions appear to be in millimeters.

## Clearly Readable Dimensions

### RSU / Sensor Strip

```text
RSU pitch / length along module long axis: 21.67
RSU active/sensor width:                  19.56
Six-RSU strip length:                    136.02
```

The six-RSU strip length is consistent with:

```text
136.02 = 4.50 + 6 * 21.67 + 1.50
```

Interpreted as:

```text
left margin before RSU chain   = 4.50
six RSUs                       = 6 * 21.67 = 130.02
right margin after RSU chain   = 1.50
total strip length             = 136.02
```

### Full Module / Package Envelope

```text
full package length: 152.02
full package width:   30.00
inner width marker:   29.00
```

The full package length is consistent with:

```text
152.02 = 11.00 + 136.02 + 5.00
```

Interpreted as:

```text
left external/package extension   = 11.00
central six-RSU strip             = 136.02
right external/package extension  = 5.00
full package length               = 152.02
```

### Width Direction / Carbon Fiber Support

The carbon fiber support has a main width of `29.00`, plus `3.00` mm associated with the bottom-edge cut features.

The bottom edge cutouts are:

```text
bottom-left cutout:   3.00 x 11.00
bottom-right cutout:  3.00 x 5.00
```

So the support is not a simple rectangle if these edge cuts are represented explicitly. A simplified envelope model could use the full package outline, while a more detailed support model should subtract or omit these two bottom-edge rectangles.

### RSU Width Within Support

The RSUs are narrower than the carbon fiber support. The drawing gives top and bottom margins of:

```text
top margin:     6.22
bottom margin:  6.22
```

Using the existing/read RSU width:

```text
RSU width = 19.56
```

the vertical span is consistent with:

```text
6.22 + 19.56 + 6.22 = 32.00
```

This supports the earlier `32 mm` corrugated/support envelope assumption for the RSU region, while the drawing also shows a `30.00 mm` package/support width marker in another view. This should be resolved before finalizing the module envelope convention in code.

## Other Visible Callouts

The following dimensions are visible in smaller end/detail views and should be treated as drawing-read values until cross-checked against the CAD/source model:

```text
edge clearance values:       0.50, repeated at ends
small end feature:           2.40
side/end view width:        30.00
side/end active width:      19.56
RSU top/bottom margins:      6.22 each
curvature radius:            R5.12
small rectangular feature:   4.00
```

Additional left-side FPC/readout-looking feature callouts:

```text
10.00
9.30
17.00
15.00
7.00
5.39
5.00
1.00
diameter 1.50
```

## FPC Geometry Considerations

The bridge and main FPC dimensions should be treated separately in the geometry model.

Bridge FPCs are good candidates for the next simple module-local approximation because their dimensions are tied to the local end/readout features:

```text
left bridge FPC:
  implemented rectangular approximation: 27.432 mm x 10.0 mm
  stack approximation:    80 um Kapton
                          + (0.69 * 15 um) effective aluminum
                          + (0.94 * 15 um) effective aluminum

right bridge FPC:
  implemented rectangular approximation: 19.7612 mm x 4.0 mm
  stack approximation:    80 um Kapton
                          + (0.42 * 15 um) effective aluminum
                          + no bottom metal layer in the current table
```

The left-bridge width was rounded to `10.0 mm` so it fits the `11.0 mm` LEC-side span with `0.5 mm` clearance to the LEC and `0.5 mm` clearance to the module edge. The right-bridge width was rounded to `4.0 mm` so it fits the `5.0 mm` REC-side span with the same two `0.5 mm` clearances. These are rectangular material-equivalent approximations, not exact outline models.

Main FPC should be deferred from the module-local RSU implementation. Page 2 gives useful reference values:

```text
connector-end width:          12.5476 mm
body width:                    9.1186 mm
reference total length:      514.83 mm
sensor-to-main-FPC distance:  17.385 mm theoretical, 17.392 mm practical
RSU overlap:                   1.0 mm
top/bottom fill factors:       0.45 / 0.68
```

Because the main FPC length varies row by row with the number of modules, it should probably be implemented later as row-level disk geometry, either derived from the production placement CSV or controlled by a dedicated FPC/corrugation-style CSV.

## Constants Already Matching the Drawing

The existing XML constants in `compact/tracking/silicon_disks_modules.xml` match the drawing well:

```xml
<constant name="SiEndcapRSU_width"  value="19.564*mm"/>
<constant name="SiEndcapRSU_length" value="21.666*mm"/>
```

These correspond to the visible drawing values:

```text
RSU width  ~= 19.56
RSU length ~= 21.67
```

## RSU Inactive-Region Assumptions

The current geometry implementation uses a simplified tile-level RSU approximation informed by `1RSU_top_half.png`. This is still an approximation: fine routing traces and detailed metal features are not modeled.

Each RSU is represented as:

```text
2 divisions along local x  x  3 tile columns per x-half  x  2 divisions along local y
= 12 active rectangular silicon tile regions
```

The inactive pieces are represented as passive silicon around or between those active regions:

```text
backbone width:     0.06
power-switch width: 0.02
bias width:         0.06
periphery width:    0.525
```

Interpretation:

- `backbone`: narrow inactive strip at the start of each local-x half-RSU. In the current model it runs across the full RSU width in local `y`.
- `power-switch`: narrow inactive strip after each tile column within a local-x half-RSU.
- `bias`: narrow inactive strip at the internal boundary between the two local-y halves. It follows the vertex-barrel convention where bias regions sit near the central boundary.
- `periphery`: inactive strip near the outer local-y edges of the RSU, representing the readout periphery plus pads/dicing lane. The current `0.525 mm` value corresponds to `0.200 + 0.325 mm` from the RSU architecture drawing.

With the current constants, each active rectangle is:

```text
tile_x = (RSU_length/2 - backbone_width - 3 * power-switch_width) / 3
       = (21.666/2 - 0.06 - 3 * 0.02) / 3
       = 3.571

active_y = RSU_width/2 - bias_width - periphery_width
         = 19.564/2 - 0.06 - 0.525
         = 9.197
```

This approximation captures the tracking-relevant active versus inactive silicon pattern while representing the drawing's 12-tile RSU structure. A 6-RSU module therefore has `6 * 12 = 72` sensitive silicon regions.

## Implementation Notes

For a detailed six-RSU corrugated/endcap module, the drawing suggests that the old simplified `130 mm x 30 mm` module footprint is only the approximate six-RSU active pitch region. A more complete package model should likely distinguish:

```text
RSU chain length:       6 * 21.666 = 129.996
sensor/readout strip:   136.02
full package length:    152.02
full package width:      30.00
```

Candidate constants for implementation:

```xml
<constant name="SiEndcapModule6RSU_package_length"      value="152.02*mm"/>
<constant name="SiEndcapModule_width_corrugated"        value="32.0*mm"/>
<constant name="SiEndcapModule6RSU_left_extension"      value="11.0*mm"/>
<constant name="SiEndcapModule6RSU_right_extension"     value="5.0*mm"/>
<constant name="SiEndcapModule6RSU_sensor_left_margin"  value="4.5*mm"/>
<constant name="SiEndcapModule6RSU_sensor_right_margin" value="1.5*mm"/>
```

Drawing-only values not currently represented as active XML constants:

```text
sensor/readout strip length: 136.02 mm
carbon fiber inner width:    29.0 mm
bottom cut width:             3.0 mm
left/right bottom cuts:      11.0 mm / 5.0 mm
RSU y margin:                 6.22 mm
```

Open validation items before final geometry implementation:

- Confirm whether the drawing's `30.00 mm` package width supersedes the earlier `32 mm` corrugated-width assumption.
- Confirm whether the RSU-region width convention should use `32.00 = 6.22 + 19.56 + 6.22`, while the carbon fiber support/package view uses `29.00 + 3.00`.
- Confirm whether the carbon fiber support bottom-edge cutouts should be represented explicitly or approximated by a rectangular envelope.
- Confirm which end features correspond to LEC, REC, FPC, or mechanical support.
- Confirm material choices for the non-sensitive package extensions.
- Confirm whether the small left-side detail dimensions should be represented in DD4hep or treated as simplified envelope material.
- Confirm bridge FPC orientation and anchoring relative to LEC/REC before implementing the rectangular approximations.
- Defer main FPC implementation until row-level length and placement conventions are defined.
