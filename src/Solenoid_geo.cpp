// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong

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
// Unified solenoid detector constructor with barrel and endcap components
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e,
                             [[maybe_unused]] SensitiveDetector sens) {
  xml_det_t x_det = e;
  string det_name = x_det.nameStr();
  Material air    = description.air();
  DetElement sdet(det_name, x_det.id());

  // Read envelope dimensions
  xml_comp_t x_envelope = x_det.child(_U(envelope));
  double env_rmin       = x_envelope.rmin();
  double env_rmax       = x_envelope.rmax();
  double env_length     = x_envelope.length();

  // Create envelope tube
  Tube envelope_tube(env_rmin, env_rmax, env_length / 2.0);
  Volume envelope_vol(det_name + "_envelope", envelope_tube, air);
  envelope_vol.setVisAttributes(description, x_det.visStr());

  int component_num = 0;

  // Iterate over components
  for (xml_coll_t comp_iter(x_det, _Unicode(component)); comp_iter; ++comp_iter, ++component_num) {
    xml_comp_t x_comp = comp_iter;
    string comp_name  = x_comp.nameStr();
    string comp_type  = x_comp.attr<string>(_Unicode(type));

    DetElement comp_det(sdet, comp_name, component_num);
    Assembly comp_assembly(comp_name + "_assembly");

    if (comp_type == "barrel") {
      // ===== BARREL COMPONENT LOGIC =====
      int i_layer = 0;
      for (xml_coll_t l_iter(x_comp, _U(layer)); l_iter; ++l_iter, ++i_layer) {
        xml_comp_t x_layer = l_iter;
        string l_name =
            getAttrOrDefault<string>(x_layer, _U(name), comp_name + _toString(i_layer, "_layer%d"));
        double outer_z = x_layer.outer_z();
        double inner_r = x_layer.inner_r();
        double outer_r = inner_r;
        DetElement layer(comp_det, _toString(i_layer, "layer%d"), x_layer.id());
        Tube l_tub(inner_r, 2 * inner_r, outer_z);
        Volume l_vol(l_name, l_tub, air);

        int i_slice        = 0;
        double l_thickness = 0.0;
        for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter, ++i_slice) {
          xml_comp_t x_slice = s_iter;
          Material mat       = description.material(x_slice.materialStr());
          string s_name =
              getAttrOrDefault<string>(x_slice, _U(name), l_name + _toString(i_slice, "_slice%d"));
          double thickness = x_slice.thickness();
          if (thickness > l_thickness)
            l_thickness = thickness;
          double s_outer_z = dd4hep::getAttrOrDefault(x_slice, _Unicode(outer_z), outer_z);
          double s_inner_z = dd4hep::getAttrOrDefault(x_slice, _Unicode(inner_z), 0.0 * cm);
          Tube s_tub(inner_r, inner_r + thickness,
                     (s_inner_z > 0 ? 0.5 * (s_outer_z - s_inner_z) : s_outer_z));
          Volume s_vol(s_name, s_tub, mat);

          s_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(),
                              x_slice.visStr());
          if (s_inner_z > 0) {
            Position s_pos(0, 0, 0.5 * (s_outer_z + s_inner_z));
            PlacedVolume pv1 = l_vol.placeVolume(s_vol, -s_pos);
            PlacedVolume pv2 = l_vol.placeVolume(s_vol, +s_pos);
            pv1.addPhysVolID("slice", i_slice);
            pv2.addPhysVolID("slice", i_slice);
          } else {
            outer_r += l_thickness;
            PlacedVolume pv = l_vol.placeVolume(s_vol);
            pv.addPhysVolID("slice", i_slice);
          }
        }
        l_tub.setDimensions(inner_r, outer_r, outer_z);
        l_vol.setVisAttributes(description, x_layer.visStr());

        PlacedVolume layer_pv = comp_assembly.placeVolume(l_vol);
        layer_pv.addPhysVolID("layer", i_layer);
        layer.setPlacement(layer_pv);
      }

    } else if (comp_type == "disk") {
      // ===== DISK (ENDCAP) COMPONENT LOGIC =====
      bool reflect = x_comp.attr<bool>(_Unicode(reflect));
      int l_num    = 0;

      for (xml_coll_t l_iter(x_comp, _U(layer)); l_iter; ++l_iter, ++l_num) {
        xml_comp_t x_layer = l_iter;
        string l_nam       = comp_name + _toString(l_num, "_layer%d");
        double zmin        = x_layer.inner_z();
        double rmin        = x_layer.inner_r();
        double rmax        = x_layer.outer_r();
        double layerWidth  = 0.;

        for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter) {
          double thickness = xml_comp_t(s_iter).thickness();
          layerWidth += thickness;
        }

        Tube l_tub(rmin, rmax, layerWidth / 2.0, 2 * M_PI);
        Volume l_vol(l_nam, l_tub, air);
        l_vol.setVisAttributes(description, x_layer.visStr());

        DetElement layer;
        PlacedVolume layer_pv;
        if (!reflect) {
          layer    = DetElement(comp_det, l_nam + "_pos", l_num);
          layer_pv = comp_assembly.placeVolume(l_vol, Position(0, 0, zmin + layerWidth / 2.));
          layer_pv.addPhysVolID("layer", l_num);
          layer.setPlacement(layer_pv);
        } else {
          layer    = DetElement(comp_det, l_nam + "_neg", l_num);
          layer_pv = comp_assembly.placeVolume(
              l_vol, Transform3D(RotationY(M_PI), Position(0, 0, -zmin - layerWidth / 2)));
          layer_pv.addPhysVolID("layer", l_num);
          layer.setPlacement(layer_pv);
        }

        double tot_thickness = -layerWidth / 2.0;
        int s_num            = 0;
        for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter, ++s_num) {
          xml_comp_t x_slice = s_iter;
          double thick       = x_slice.thickness();
          Material mat       = description.material(x_slice.materialStr());
          string s_nam       = l_nam + _toString(s_num, "_slice%d");
          Volume s_vol(s_nam, Tube(rmin, rmax, thick / 2.0), mat);
          if (!reflect) {
            s_nam += "_pos";
          } else {
            s_nam += "_neg";
          }
          DetElement slice_de(layer, s_nam, s_num);
          s_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(),
                              x_slice.visStr());
          PlacedVolume pv = l_vol.placeVolume(s_vol, Position(0, 0, tot_thickness + thick / 2));
          pv.addPhysVolID("slice", s_num);
          slice_de.setPlacement(pv);
          tot_thickness = tot_thickness + thick;
        }
      }
    }

    // Place component assembly in envelope
    PlacedVolume comp_pv = envelope_vol.placeVolume(comp_assembly);
    comp_pv.addPhysVolID("component", component_num);
    comp_det.setPlacement(comp_pv);
  }

  // Get position and place envelope volume in world
  xml::Component x_pos = x_det.position();
  Position pos(x_pos.x(), x_pos.y(), x_pos.z());
  PlacedVolume pv = description.pickMotherVolume(sdet).placeVolume(envelope_vol, pos);
  pv.addPhysVolID("system", sdet.id());
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_Solenoid, create_detector)
