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

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = description.air();
  DetElement sdet        (det_name,x_det.id());

  // dimensions
  xml_comp_t dims = x_det.dimensions();
  auto rmin   = dims.rmin();
  auto rmax   = dims.rmax();
  auto length = dims.length();
  Assembly assembly(det_name + "_assembly");
  PlacedVolume pv;

  int i_layer = 0;
  for (xml_coll_t l_iter(x_det,_U(layer)); l_iter; ++l_iter, ++i_layer)  {
    xml_comp_t x_layer = l_iter;
    string l_name = getAttrOrDefault<string>(x_layer, _U(name), det_name + _toString(i_layer, "_layer%d"));
    double outer_z = x_layer.outer_z();
    double inner_r = x_layer.inner_r();
    double outer_r = inner_r;
    DetElement layer(sdet, _toString(i_layer, "layer%d"), x_layer.id());
    Tube   l_tub(inner_r, 2*inner_r, outer_z); // outer_r will be updated
    Volume l_vol(l_name, l_tub, air);

    int i_slice = 0;
    double l_thickness = 0.0;
    for (xml_coll_t s_iter(x_layer,_U(slice)); s_iter; ++s_iter, ++i_slice)  {
      // If slices are only given a thickness attribute, they are radially concentric slices
      // If slices are given an inner_z attribute, they are longitudinal slices with equal rmin
      xml_comp_t x_slice = s_iter;
      Material mat  = description.material(x_slice.materialStr());
      string s_name = getAttrOrDefault<string>(x_slice, _U(name), l_name + _toString(i_slice, "_slice%d"));
      double thickness = x_slice.thickness();
      if (thickness > l_thickness) l_thickness = thickness;
      double s_outer_z = dd4hep::getAttrOrDefault(x_slice, _Unicode(outer_z), outer_z);
      double s_inner_z = dd4hep::getAttrOrDefault(x_slice, _Unicode(inner_z), 0.0*cm);
      Tube   s_tub(inner_r, inner_r + thickness, (s_inner_z > 0? 0.5 * (s_outer_z - s_inner_z): s_outer_z));
      Volume s_vol(s_name, s_tub, mat);

      s_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      if (s_inner_z > 0) {
        // Place off-center volumes twice
        Position s_pos(0, 0, 0.5 * (s_outer_z + s_inner_z));
        pv = l_vol.placeVolume(s_vol, -s_pos);
        pv = l_vol.placeVolume(s_vol, +s_pos);
      } else {
        outer_r += l_thickness;
        pv = l_vol.placeVolume(s_vol);
      }
      // Slices have no extra id. Take the ID of the layer!
      pv.addPhysVolID("slice",i_slice);
    }
    l_tub.setDimensions(inner_r, outer_r, outer_z);
    l_vol.setVisAttributes(description, x_layer.visStr());

    pv = assembly.placeVolume(l_vol);
    pv.addPhysVolID("layer",i_layer);
    layer.setPlacement(pv);
  }

  // Get position and place volume
  xml::Component x_pos = x_det.position();
  Position pos(x_pos.x(), x_pos.y(), x_pos.z());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, pos);
  pv.addPhysVolID("system", sdet.id()).addPhysVolID("barrel",0);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(ecce_Solenoid,create_detector)
