# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2025 ePIC Collaboration
#
# Python port of SimpleDiskDetector_geo.cpp using DD4hep cppyy bindings.
# Serves as the plugin target for the C++ trampoline epic_PythonDetector.
#
# Compact XML usage:
#   <detector id="X_ID" name="X" type="epic_PythonDetector"
#             module="SimpleDisk_python_geo" function="create_detector"
#             reflect="false" readout="XHits"
#             pythonpath="${DETECTOR_PATH}/python">
#     <position x="0" y="0" z="0*cm"/>
#     <layer inner_z="..." inner_r="..." outer_r="...">
#       <slice material="Silicon" thickness="2*mm" sensitive="true"/>
#     </layer>
#   </detector>

from epic_geo_helpers import *  # noqa: F401,F403


def create_detector(description, xml_element, sens):
    """
    Create a simple disk detector (layered endcap).

    Parameters
    ----------
    description : dd4hep::Detector &
    xml_element : dd4hep::xml::Handle_t  (the <detector> node)
    """
    try:
        return _create_detector_impl(description, xml_element, sens)
    except Exception:
        import traceback
        traceback.print_exc()
        raise


def _create_detector_impl(description, xml_element, sens):
    x_det    = xml_det_t(xml_element)
    air      = description.air()
    det_name = x_det.nameStr()
    reflect  = _getAttrBool(x_det, "reflect") if _hasAttr(x_det, "reflect") else False
    sdet     = DetElement(det_name, x_det.id())
    assembly = Assembly(det_name)

    pos   = xml_comp_t(x_det.child(_U("position")))
    l_num = 0
    pv    = PlacedVolume()

    i = xml_coll_t(_handle(x_det), "layer")
    while _collValid(i):
        x_layer     = xml_comp_t(i.current())
        l_nam       = det_name + _toString(l_num, "_layer%d")
        zmin        = x_layer.inner_z()
        rmin        = x_layer.inner_r()
        rmax        = x_layer.outer_r()
        layer_width = 0.0
        s_num       = 0

        j = xml_coll_t(_handle(x_layer), "slice")
        while _collValid(j):
            x_slice = xml_comp_t(j.current())
            layer_width += x_slice.thickness()
            j.__preinc__()

        l_tub = Tube(rmin, rmax, layer_width / 2.0, 2 * TMath.Pi())
        l_vol = Volume(l_nam, l_tub, air)
        l_vol.setVisAttributes(description, x_layer.visStr())

        if not reflect:
            layer    = DetElement(sdet, l_nam + "_pos", l_num)
            layer_pv = assembly.placeVolume(l_vol, Position(0, 0, zmin + layer_width / 2.0))
            layer_pv.addPhysVolID("barrel", 3).addPhysVolID("layer", l_num)
            layer.setPlacement(layer_pv)
        else:
            layer    = DetElement(sdet, l_nam + "_neg", l_num)
            layer_pv = assembly.placeVolume(
                l_vol, Transform3D(RotationY(TMath.Pi()), Position(0, 0, -zmin - layer_width / 2.0))
            )
            layer_pv.addPhysVolID("barrel", 2).addPhysVolID("layer", l_num)
            layer.setPlacement(layer_pv)

        tot_thickness = -layer_width / 2.0
        j = xml_coll_t(_handle(x_layer), "slice")
        while _collValid(j):
            x_slice = xml_comp_t(j.current())
            thick   = x_slice.thickness()
            mat     = description.material(x_slice.materialStr())
            s_nam   = l_nam + _toString(s_num, "_slice%d")
            s_vol   = Volume(s_nam, Tube(rmin, rmax, thick / 2.0), mat)
            s_nam  += "_pos" if not reflect else "_neg"
            slice_de = DetElement(layer, s_nam, s_num)
            if x_slice.isSensitive():
                sens.setType("tracker")
                s_vol.setSensitiveDetector(sens)
            s_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr())
            pv = l_vol.placeVolume(s_vol, Position(0, 0, tot_thickness + thick / 2.0))
            pv.addPhysVolID("slice", s_num)
            slice_de.setPlacement(pv)
            tot_thickness += thick
            s_num += 1
            j.__preinc__()

        i.__preinc__()
        l_num += 1

    if _hasAttr(x_det, "combineHits"):
        sdet.setCombineHits(_getAttrBool(x_det, "combineHits"), sens)

    mother_pv = description.pickMotherVolume(sdet).placeVolume(
        assembly, Position(pos.x(), pos.y(), pos.z())
    )
    mother_pv.addPhysVolID("system", x_det.id())
    sdet.setPlacement(mother_pv)
    return sdet
