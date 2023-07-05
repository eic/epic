// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Maria Zurek, Whitney Armstrong

// Detector plugin to support a hybrid central barrel calorimeter
// The detector consists of interlayers of Pb/ScFi (segmentation in global r, phi) and W/Si (segmentation in local x, y)
// Assembly is used as the envelope so two different detectors can be interlayered with each other
//
//
// 06/19/2021: Implementation of the Sci Fiber geometry. M. Żurek
// 07/09/2021: Support interlayers between multiple detectors. C. Peng
// 07/23/2021: Add assemblies as mother volumes of fibers to reduce the number of daughter volumes. C. Peng, M. Żurek
//     Reference: TGeo performance issue with large number of daughter volumes
//     https://indico.cern.ch/event/967418/contributions/4075358/attachments/2128099/3583278/201009_shKo_dd4hep.pdf
// 07/24/2021: Changed support implementation to avoid too many uses of boolean geometries. DAWN view seems to have
//     issue dealing with it. C. Peng

#include "DD4hep/DetFactoryHelper.h"
#include "Math/Point2D.h"
#include "TGeoPolygon.h"
#include "XML/Layering.h"
#include <functional>

using namespace std;
using namespace dd4hep;

typedef ROOT::Math::XYPoint Point;

// barrel ecal layers contained in an assembly
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  Layering   layering(e);
  xml_det_t  x_det    = e;
  int        det_id   = x_det.id();
  string     det_name = x_det.nameStr();
  double     offset   = x_det.attr<double>(_Unicode(offset));
  xml_comp_t x_dim    = x_det.dimensions();
  int        nsides   = x_dim.numsides();
  double     inner_r  = x_dim.rmin();
  double     dphi     = (2 * M_PI / nsides);
  double     hphi     = dphi / 2;

  DetElement sdet(det_name, det_id);
  Volume     motherVol = desc.pickMotherVolume(sdet);

  Assembly     envelope(det_name);
  Transform3D  tr_global = Translation3D(0, 0, offset) * RotationZ(hphi);
  PlacedVolume env_phv   = motherVol.placeVolume(envelope, tr_global);
  sens.setType("calorimeter");

  env_phv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_phv);

  // build a single sector
  DetElement sector_det("sector0", det_id);
  Assembly   mod_vol("sector");

  // keep tracking of the total thickness
  double l_pos_z = inner_r;
  { // =====  buildBarrelStave(desc, sens, module_volume) =====
    // Parameters for computing the layer X dimension:
    double tan_hphi = std::tan(hphi);
    double l_dim_y  = x_dim.z() / 2.;

    // Loop over the sets of layer elements in the detector.
    int l_num = 1;
    for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
      xml_comp_t x_layer         = li;
      int        repeat          = x_layer.repeat();
      double     l_space_between = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_between), 0.);
      double     l_space_before  = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_before), 0.);
      bool       framing         = x_layer.hasChild(_Unicode(frame));
      l_pos_z += l_space_before;
      // Loop over number of repeats for this layer.
      for (int j = 0; j < repeat; j++) {
        // make an envelope for this layer
        string l_name         = Form("layer%d", l_num);
        double l_sthickness   = layering.layer(l_num - 1)->thickness(); // total thickness of slices
        double l_dim_x        = tan_hphi * l_pos_z;
        double l_thickness    = l_sthickness;
        auto   lfill_mat      = desc.air();
        auto   lenv_mat       = desc.air();
        double wt             = 0.;
        // update values if there will be a frame
        if (framing) {
          auto x_frame = x_layer.child(_Unicode(frame));
          l_thickness  = x_frame.attr<double>(_Unicode(height));
          wt           = x_frame.attr<double>(_Unicode(thickness));
          lenv_mat     = desc.material(x_frame.attr<std::string>(_Unicode(material)));
          lfill_mat    = desc.material(x_frame.attr<std::string>(_Unicode(fill)));
          // sanity check
          if (l_thickness < (l_sthickness + 2*wt)) {
            std::cerr << "EcalBarrelImaging_geo Error: frame available space " << l_thickness - 2.*wt
                      << "is less than total thickness of slices!" << l_sthickness
                      << std::endl;
          }
        }

        Position   l_pos(0, 0, l_pos_z + l_thickness / 2.); // Position of the layer.
        double     l_trd_x1 = l_dim_x;
        double     l_trd_x2 = l_dim_x + l_thickness * tan_hphi;
        double     l_trd_y1 = l_dim_y;
        double     l_trd_y2 = l_trd_y1;
        double     l_trd_z  = l_thickness / 2;
        Trapezoid  l_shape(l_trd_x1, l_trd_x2, l_trd_y1, l_trd_y2, l_trd_z);
        Volume     l_vol(l_name, l_shape, lenv_mat);
        DetElement layer(sector_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
        int    s_num   = 1;
        double s_pos_z = -(l_thickness / 2.) + wt;
        for (xml_coll_t si(x_layer, _U(slice)); si; ++si) {
          xml_comp_t x_slice  = si;
          string     s_name   = Form("slice%d", s_num);
          double     s_thick  = x_slice.thickness();
          double     s_trd_x1 = l_dim_x + (s_pos_z + l_thickness / 2) * tan_hphi;
          double     s_trd_x2 = l_dim_x + (s_pos_z + l_thickness / 2 + s_thick) * tan_hphi;
          double     s_trd_y1 = l_trd_y1;
          double     s_trd_y2 = s_trd_y1;
          double     s_trd_z  = s_thick / 2.;
          Trapezoid  s_shape(s_trd_x1 - wt, s_trd_x2 - wt, s_trd_y1, s_trd_y2, s_trd_z);
          Volume     s_vol(s_name, s_shape, desc.material(x_slice.materialStr()));
          DetElement slice(layer, s_name, det_id);

          if (x_slice.isSensitive()) {
            s_vol.setSensitiveDetector(sens);
          }
          s_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

          // Slice placement.
          PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick / 2));
          slice_phv.addPhysVolID("slice", s_num);
          slice.setPlacement(slice_phv);
          // Increment Z position of slice.
          s_pos_z += s_thick;
          ++s_num;
        }

        // filling gap
        if (framing && (l_thickness > (l_sthickness + 2.*wt))) {
          string     s_name   = Form("layer_filling%d", l_num);
          double     s_thick  = l_thickness - l_sthickness -2 *wt;
          double     s_trd_x1 = l_dim_x + (s_pos_z + l_thickness / 2) * tan_hphi;
          double     s_trd_x2 = l_dim_x + (s_pos_z + l_thickness / 2 + s_thick) * tan_hphi;
          double     s_trd_y1 = l_trd_y1;
          double     s_trd_y2 = s_trd_y1;
          double     s_trd_z  = s_thick / 2.;
          Trapezoid  s_shape(s_trd_x1, s_trd_x2, s_trd_y1, s_trd_y2, s_trd_z);
          Volume     s_vol(s_name, s_shape, lfill_mat);
          s_vol.setVisAttributes(desc.invisible());
          l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick / 2));
        }

        // Set region, limitset, and vis of layer.
        l_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

        PlacedVolume layer_phv = mod_vol.placeVolume(l_vol, l_pos);
        layer_phv.addPhysVolID("layer", l_num);
        layer.setPlacement(layer_phv);
        // Increment to next layer Z position. Do not add space_between for the last layer
        l_pos_z += l_thickness;
        if (j < repeat - 1) {
          l_pos_z += l_space_between;
        }
        ++l_num;
      }
    }
  }
  // Phi start for a sector.
  double phi = M_PI / nsides;
  // Create nsides sectors.
  for (int i = 0; i < nsides; i++, phi -= dphi) { // i is module number
    // Compute the sector position
    Transform3D  tr(RotationZYX(0, phi, M_PI * 0.5), Translation3D(0, 0, 0));
    PlacedVolume pv = envelope.placeVolume(mod_vol, tr);
    pv.addPhysVolID("module", i + 1);
    DetElement sd = (i == 0) ? sector_det : sector_det.clone(Form("sector%d", i));
    sd.setPlacement(pv);
    sdet.add(sd);
  }

  // Set envelope volume attributes.
  envelope.setAttributes(desc, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  return sdet;
}

DECLARE_DETELEMENT(epic_EcalBarrelImaging, create_detector)
