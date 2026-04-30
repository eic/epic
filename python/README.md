# Python Geometry Plugins

DD4hep geometry plugins are normally written in C++ and registered via
`DECLARE_DETELEMENT`.  This directory provides infrastructure that lets you
write geometry plugins in **Python** instead, using DD4hep's cppyy bindings.

## How it works

The C++ plugin `epic_PythonDetector` (in `src/PythonDetector_geo.cpp`) acts as
a *trampoline*: it is declared to DD4hep as a normal C++ plugin but immediately
delegates to a Python function via ROOT's `TPython::Exec`.  C++ objects
(`Detector&`, the XML handle, and `SensitiveDetector`) are passed to Python by
address and reconstructed as cppyy proxies; the resulting `DetElement` is
transferred back to C++ via `ROOT::Internal::SwapWithObjAtAddr`.

## Compact XML usage

```xml
<detector id="1"
          name="MyPythonDet"
          type="epic_PythonDetector"
          module="my_detector_geo"
          function="create_detector"
          readout="MyDetHits"
          pythonpath="${DETECTOR_PATH}/python">
  <!-- any sub-elements your Python function reads -->
</detector>
```

| Attribute      | Required | Description |
|----------------|----------|-------------|
| `module`       | yes      | Python module name (`import module`) |
| `function`     | yes      | Callable inside that module |
| `pythonpath`   | no       | Directory prepended to `sys.path`; `${ENV_VAR}` is expanded |

## Python function signature

```python
def create_detector(description, xml_element, sens) -> dd4hep.DetElement:
    ...
```

The arguments arrive as cppyy proxies for:

| Argument      | C++ type                    |
|---------------|-----------------------------|
| `description` | `dd4hep::Detector&`         |
| `xml_element` | `dd4hep::xml::Handle_t`     |
| `sens`        | `dd4hep::SensitiveDetector` |

## Shared helper module

All the cppyy boilerplate is collected in `python/epic_geo_helpers.py`.
Start every plugin with a single import:

```python
from epic_geo_helpers import *
```

This provides:

| Name | Description |
|------|-------------|
| `TMath` | ROOT math constants |
| `_handle(el)` | Cast any XML element to `Handle_t` (required by `xml_coll_t`) |
| `_collValid(c)` | Test whether an `xml_coll_t` iterator is still valid |
| `_hasAttr(el, name)` | Test for an XML attribute (char* overload) |
| `_getAttrBool(el, name)` | Read a bool XML attribute |
| `Assembly`, `DetElement`, `PlacedVolume`, `Position`, `Transform3D`, `RotationY`, `Tube`, `Volume` | Common DD4hep geometry types |
| `_U`, `_toString` | XML macro equivalents (`Strng_t`, `_toString`) |
| `xml_det_t`, `xml_coll_t`, `xml_comp_t` | XML handle types |

## Writing a plugin — key cppyy patterns

### Minimal skeleton

```python
from epic_geo_helpers import *

def create_detector(description, xml_element, sens):
    x_det    = xml_det_t(xml_element)
    det_name = x_det.nameStr()
    sdet     = DetElement(det_name, x_det.id())
    assembly = Assembly(det_name)
    # ... build geometry ...
    mother_pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(0, 0, 0))
    mother_pv.addPhysVolID("system", x_det.id())
    sdet.setPlacement(mother_pv)
    return sdet
```

### Iterating over child XML elements

```python
# xml_coll_t constructor (const char* overload works fine)
i = xml_coll_t(_handle(parent_element), "child_tag")
while _collValid(i):                # NOT: while i: — always True!
    child = xml_comp_t(i.current())
    ...
    i.__preinc__()
```

### Reading attributes

```python
# Regular string attributes work via nameStr() / materialStr() / visStr() etc.
name = x_det.nameStr()

# Numeric/bool attributes: use _hasAttr/_getAttrBool from epic_geo_helpers
if _hasAttr(x_det, "reflect"):
    reflect = _getAttrBool(x_det, "reflect")

# Child element access via _U (Strng_t) works fine:
pos = xml_comp_t(x_det.child(_U("position")))
```

## Example plugin

`python/SimpleDisk_python_geo.py` is a complete Python port of
`src/SimpleDiskDetector_geo.cpp`.  It builds a layered endcap disk detector
with full layer/slice/PhysVolID support.

Test it with:

```bash
source <install_prefix>/bin/thisepic.sh
echo ".q" | geoDisplay -compact compact/tests/python_disk.xml -no-vis
```

Expected output:
```
Compact  INFO  ++ Converted subdetector:PythonDisk of type epic_PythonDetector
```

## Known limitations

* **No re-entrant calls**: the trampoline uses a single static `DetElement` to
  pass the result back to C++.  Nested Python detectors are not supported.
* **`xml_coll_t` operator bool**: cppyy does not map `operator bool()` to
  Python `__bool__` for this class.  Always use `_collValid()`.
* **Attribute access**: `hasAttr` / `attr<T>` take `const char16_t*`; Python
  `str` is incompatible.  Always use `_hasAttr()`, `_getAttrBool()`, and the
  underlying `epic_python::getAttrDbl` helper from `epic_geo_helpers`.
