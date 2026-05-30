// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong
// Copyright (C) 2026 Wouter Deconinck

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include <XML/Layering.h>
#include <XML/Utilities.h>
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

  xml_dim_t pos = x_det.position();
  xml_dim_t rot = x_det.rotation();

  double totWidth = Layering(x_det).totalThickness();

  Box envelope(Width / 2.0, Width / 2.0, totWidth / 2.0);
  Volume envelopeVol(detName + "_envelope", envelope, desc.material("Vacuum"));
  envelopeVol.setVisAttributes(desc.visAttributes(x_det.visStr()));
  PlacedVolume pv;

  // z of the front face of the current layer group, relative to envelope centre
  double z_env = -totWidth / 2.0;
  // 1-based layer counter used for physVolID (preserves existing cell IDs)
  int layer_num = 1;
  int lg        = 0; // layer-group counter used for unique volume naming

  // Read layer groups from XML
  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    int repeat         = x_layer.repeat();
    if (repeat <= 0)
      continue;
    ++lg;

    double layerWidth = 0;
    for (xml_coll_t l(x_layer, _U(slice)); l; ++l)
      layerWidth += xml_comp_t(l).thickness();

    // Build ONE layer template volume per layer group (not one per repetition).
    // All `repeat` copies share this single TGeoVolume, reducing unique
    // TGeoVolume count from  repeat * (1 + nSlices)  to  1 + nSlices.
    string layer_name = detName + _toString(lg, "_lg%d");
    Volume layer_vol(layer_name, Box(Width / 2.0, Width / 2.0, layerWidth / 2.0),
                     desc.material("Vacuum"));

    // Place slice volumes into the template once; their positions are in
    // layer-local coordinates, identical for every repetition.
    int slice_num  = 1;
    double z_slice = -layerWidth / 2.0; // front face of current slice (layer-local)
    for (xml_coll_t l(x_layer, _U(slice)); l; ++l, ++slice_num) {
      xml_comp_t x_slice = l;
      double w           = x_slice.thickness();
      string slice_name  = layer_name + _toString(slice_num, "_sl%d");
      Material slice_mat = desc.material(x_slice.materialStr());
      Volume slice_vol(slice_name, Box(Width / 2.0, Width / 2.0, w / 2.0), slice_mat);

      if (x_slice.isSensitive()) {
        sens.setType("calorimeter");
        slice_vol.setSensitiveDetector(sens);
      }
      slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = layer_vol.placeVolume(slice_vol, Position(0.0, 0.0, z_slice + w / 2.0));
      pv.addPhysVolID("slice", slice_num);
      z_slice += w;
    }

    string layer_vis =
        dd4hep::getAttrOrDefault<std::string>(x_layer, _Unicode(vis), "InvisibleWithDaughters");
    layer_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), layer_vis);

    // Place the shared template `repeat` times.  Using the same TGeoVolume
    // for all copies avoids creating `repeat` identical TGeoVolume objects.
    for (int i = 0; i < repeat; ++i) {
      double z_center = z_env + (i + 0.5) * layerWidth;
      pv              = envelopeVol.placeVolume(layer_vol, Position(0.0, 0.0, z_center));
      pv.addPhysVolID("layer", layer_num + i);
    }

    z_env += repeat * layerWidth;
    layer_num += repeat;
  }

  DetElement det(detName, detID);

  // apply any detector type flags set in XML
  dd4hep::xml::setDetectorTypeFlag(x_det, det);

  Volume motherVol = desc.pickMotherVolume(det);
  Transform3D tr(RotationZYX(rot.z(), rot.y(), rot.x()),
                 Position(pos.x(), pos.y(), pos.z() + totWidth / 2.0));
  PlacedVolume phv = motherVol.placeVolume(envelopeVol, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(ZDC_Sampling, createDetector)
