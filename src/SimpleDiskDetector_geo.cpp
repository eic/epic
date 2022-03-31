//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t SimpleDiskDetector_create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  xml_det_t      x_det    = e;
  Material       air      = description.air();
  string         det_name = x_det.nameStr();
  bool           reflect  = x_det.reflect();
  DetElement     sdet(det_name, x_det.id());
  Assembly       assembly(det_name);
  PlacedVolume   pv;
  int            l_num = 0;
  xml::Component pos   = x_det.position();

  Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();
  detWorldExt->addType("endcap", "detector");
  sdet.addExtension<Acts::ActsExtension>(detWorldExt);

  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++l_num) {
    xml_comp_t x_layer    = i;
    string     l_nam      = det_name + _toString(l_num, "_layer%d");
    double     zmin       = x_layer.inner_z();
    double     rmin       = x_layer.inner_r();
    double     rmax       = x_layer.outer_r();
    double     z          = zmin;
    double     layerWidth = 0.;
    int        s_num      = 0;

    for (xml_coll_t j(x_layer, _U(slice)); j; ++j) {
      double thickness = xml_comp_t(j).thickness();
      layerWidth += thickness;
    }
    Tube   l_tub(rmin, rmax, layerWidth/2.0, 2 * M_PI);
    Volume l_vol(l_nam, l_tub, air);
    l_vol.setVisAttributes(description, x_layer.visStr());
    DetElement layer;
    PlacedVolume layer_pv;
    if (!reflect) {
      layer = DetElement(sdet, l_nam + "_pos", l_num);
      layer_pv = assembly.placeVolume(l_vol, Position(0, 0, zmin + layerWidth / 2.));
      layer_pv.addPhysVolID("barrel", 3).addPhysVolID("layer", l_num);
      layer.setPlacement(layer_pv);
      Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
      layerExtension->addType("sensitive disk", "layer");
      //layerExtension->addType("axes", "definitions", "XzY");
      // need all four of these or else it is ignored.
      //layerExtension->addValue(0, "r_min", "envelope");
      //layerExtension->addValue(0, "r_max", "envelope");
      //layerExtension->addValue(0, "z_min", "envelope");
      //layerExtension->addValue(0, "z_max", "envelope");
      // layerExtension->addType("axes", "definitions", "XZY");

      layer.addExtension<Acts::ActsExtension>(layerExtension);
    } else {
      layer = DetElement(sdet, l_nam + "_neg", l_num);
      layer_pv = assembly.placeVolume(l_vol, Transform3D(RotationY(M_PI), Position(0, 0, -zmin - layerWidth / 2)));
      layer_pv.addPhysVolID("barrel", 2).addPhysVolID("layer", l_num);
      layer.setPlacement(layer_pv);
      // DetElement layerR = layer.clone(l_nam+"_neg");
      // sdet.add(layerR.setPlacement(pv));
      Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
      layerExtension->addType("sensitive disk", "layer");
      //layerExtension->addValue(0, "r_min", "envelope");
      //layerExtension->addValue(0, "r_max", "envelope");
      //layerExtension->addValue(0, "z_min", "envelope");
      //layerExtension->addValue(0, "z_max", "envelope");
      layer.addExtension<Acts::ActsExtension>(layerExtension);
    }

    double tot_thickness = -layerWidth / 2.0;
    for (xml_coll_t j(x_layer, _U(slice)); j; ++j, ++s_num) {
      xml_comp_t x_slice = j;
      double     thick   = x_slice.thickness();
      Material   mat     = description.material(x_slice.materialStr());
      string     s_nam   = l_nam + _toString(s_num, "_slice%d");
      Volume     s_vol(s_nam, Tube(rmin, rmax, thick/2.0), mat);
      if(!reflect){
        s_nam += "_pos";
      } else {
        s_nam += "_neg";
      }
      DetElement slice_de(layer, s_nam , s_num);
      if (x_slice.isSensitive()) {
        sens.setType("tracker");
        s_vol.setSensitiveDetector(sens);
        Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
        //sensorExtension->addType("sensor", "detector");
        slice_de.addExtension<Acts::ActsExtension>(sensorExtension);
      }
      s_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = l_vol.placeVolume(s_vol, Position(0, 0, tot_thickness + thick / 2));
      pv.addPhysVolID("slice", s_num);
      slice_de.setPlacement(pv);
      tot_thickness = tot_thickness + thick;
    }

  }
  if (x_det.hasAttr(_U(combineHits))) {
    sdet.setCombineHits(x_det.attr<bool>(_U(combineHits)), sens);
  }
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(pos.x(), pos.y(), pos.z()));
  pv.addPhysVolID("system", x_det.id()); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(ref_SolenoidEndcap, SimpleDiskDetector_create_detector)
DECLARE_DETELEMENT(athena_SolenoidEndcap, SimpleDiskDetector_create_detector)
