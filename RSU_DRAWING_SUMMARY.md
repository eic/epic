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

The current geometry implementation uses a simplified vertex-barrel-style RSU approximation. This is an implementation assumption, not a fully detailed extraction from `rsu.pdf`.

Each RSU is represented as:

```text
2 divisions along local x  x  2 divisions along local y
= 4 active rectangular silicon regions
```

The inactive pieces are represented as passive silicon around or between those active regions:

```text
backbone width:  0.09
bias width:      0.06
periphery width: 0.398
```

Interpretation:

- `backbone`: narrow inactive strip separating the two local-x halves of an RSU. In the current model it runs across the full RSU width in local `y`.
- `bias`: narrow inactive strip at the internal boundary between the two local-y halves. It follows the vertex-barrel convention where bias regions sit near the central boundary.
- `periphery`: inactive strip near the outer local-y edges of the RSU, representing non-sensitive edge/readout/periphery area.

With the current constants, each active rectangle is:

```text
active_x = RSU_length/2 - backbone_width
         = 21.666/2 - 0.09
         = 10.743

active_y = RSU_width/2 - bias_width - periphery_width
         = 19.564/2 - 0.06 - 0.398
         = 9.324
```

This approximation captures the tracking-relevant active versus inactive silicon pattern without modeling every physical tile, readout, or edge-detail feature.

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
<constant name="SiEndcapModule6RSU_strip_length"        value="136.02*mm"/>
<constant name="SiEndcapModule6RSU_package_length"      value="152.02*mm"/>
<constant name="SiEndcapModule_width_corrugated"        value="30.0*mm"/>
<constant name="SiEndcapModule_inner_width_corrugated"  value="29.0*mm"/>
<constant name="SiEndcapModule_support_bottom_cut_width" value="3.0*mm"/>
<constant name="SiEndcapModule6RSU_left_cut_length"     value="11.0*mm"/>
<constant name="SiEndcapModule6RSU_right_cut_length"    value="5.0*mm"/>
<constant name="SiEndcapRSU_y_margin"                   value="6.22*mm"/>
<constant name="SiEndcapModule6RSU_left_extension"      value="11.0*mm"/>
<constant name="SiEndcapModule6RSU_right_extension"     value="5.0*mm"/>
<constant name="SiEndcapModule6RSU_sensor_left_margin"  value="4.5*mm"/>
<constant name="SiEndcapModule6RSU_sensor_right_margin" value="1.5*mm"/>
```

Open validation items before final geometry implementation:

- Confirm whether the drawing's `30.00 mm` package width supersedes the earlier `32 mm` corrugated-width assumption.
- Confirm whether the RSU-region width convention should use `32.00 = 6.22 + 19.56 + 6.22`, while the carbon fiber support/package view uses `29.00 + 3.00`.
- Confirm whether the carbon fiber support bottom-edge cutouts should be represented explicitly or approximated by a rectangular envelope.
- Confirm which end features correspond to LEC, REC, FPC, or mechanical support.
- Confirm material choices for the non-sensitive package extensions.
- Confirm whether the small left-side detail dimensions should be represented in DD4hep or treated as simplified envelope material.
