// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include <XML/Layering.h>
//////////////////////////////////////////////////
// Far Forward Ion Zero Degree Calorimeter - Hcal
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  string detName  = x_det.nameStr();
  int detID       = x_det.id();

  xml_dim_t dim = x_det.dimensions();
  double Width  = dim.x();
  // double     Length     = dim.z();

  xml_dim_t pos = x_det.position();
  double z      = pos.z();
  xml_dim_t rot = x_det.rotation();

  Material Vacuum = desc.material("Vacuum");

  double totWidth = Layering(x_det).totalThickness();

  Box envelope(Width / 2.0, Width / 2.0, totWidth / 2.0);
  Volume envelopeVol(detName + "_envelope", envelope, Vacuum);
  envelopeVol.setVisAttributes(desc.visAttributes(x_det.visStr()));
  PlacedVolume pv;

  int layer_num = 1;
  // Read layers
  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    int repeat         = x_layer.repeat();
    double layerWidth  = 0;

    for (xml_coll_t l(x_layer, _U(slice)); l; ++l)
      layerWidth += xml_comp_t(l).thickness();

    // Loop over repeat#
    for (int i = 0; i < repeat; i++) {
      double zlayer     = z;
      string layer_name = detName + _toString(layer_num, "_layer%d");
      Volume layer_vol(layer_name, Box(Width / 2.0, Width / 2.0, layerWidth / 2.0), Vacuum);

      int slice_num = 1;
      // Loop over slices
      for (xml_coll_t l(x_layer, _U(slice)); l; ++l) {
        xml_comp_t x_slice = l;
        double w           = x_slice.thickness();
        string slice_name  = layer_name + _toString(slice_num, "slice%d");
        Material slice_mat = desc.material(x_slice.materialStr());
        Volume slice_vol(slice_name, Box(Width / 2.0, Width / 2.0, w / 2.0), slice_mat);

        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }

        slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        pv = layer_vol.placeVolume(
            slice_vol, Transform3D(RotationZYX(0, 0, 0),
                                   Position(0.0, 0.0, z - zlayer - layerWidth / 2.0 + w / 2.0)));
        pv.addPhysVolID("slice", slice_num);
        z += w;
        ++slice_num;
      }

      string layer_vis =
          dd4hep::getAttrOrDefault<std::string>(x_layer, _Unicode(vis), "InvisibleWithDaughters");
      layer_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), layer_vis);
      pv = envelopeVol.placeVolume(
          layer_vol,
          Transform3D(RotationZYX(0, 0, 0),
                      Position(0, 0, zlayer - pos.z() - totWidth / 2.0 + layerWidth / 2.0)));
      pv.addPhysVolID("layer", layer_num);
      ++layer_num;
    }
  }

  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);
  Transform3D tr(RotationZYX(rot.z(), rot.y(), rot.x()),
                 Position(pos.x(), pos.y(), pos.z() + totWidth / 2.0));
  PlacedVolume phv = motherVol.placeVolume(envelopeVol, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(ZDC_Sampling, createDetector)
