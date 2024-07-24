// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
//////////////////////////////////////////////////
// Far Forward Ion Zero Degree Calorimeter - Ecal
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  string detName  = x_det.nameStr();
  int detID       = x_det.id();

  xml_dim_t dim    = x_det.dimensions();
  double Width     = dim.x();
  double Thickness = dim.z();

  xml_dim_t pos = x_det.position();
  xml_dim_t rot = x_det.rotation();

  Material Vacuum = desc.material("Vacuum");

  xml_comp_t mod    = x_det.child(_Unicode(module));
  string modName    = mod.nameStr();
  Material mPbWO4   = desc.material(mod.materialStr());
  double mThickness = mod.attr<double>(_Unicode(thickness));
  double mWidth     = mod.attr<double>(_Unicode(width));
  double mGap       = mod.attr<double>(_Unicode(gap));
  int mNTowers      = mod.attr<double>(_Unicode(ntower));

  // Create Global Volume
  Box ffi_ZDC_GVol_Solid(Width * 0.5, Width * 0.5, Thickness * 0.5);
  Volume detVol("ffi_ZDC_GVol_Logic", ffi_ZDC_GVol_Solid, Vacuum);
  detVol.setVisAttributes(desc.visAttributes(x_det.visStr()));

  // Construct Tower
  // Single Module
  Box ffi_ZDC_ECAL_Solid_Tower(mWidth * 0.5, mWidth * 0.5, mThickness * 0.5);
  Volume modVol("ffi_ZDC_ECAL_Logic_Tower", ffi_ZDC_ECAL_Solid_Tower, mPbWO4);
  modVol.setVisAttributes(desc.visAttributes(mod.visStr()));
  sens.setType("calorimeter");
  modVol.setSensitiveDetector(sens);

  // Module Position
  double mod_x = 0.0 * mm;
  double mod_y = 0.0 * mm;
  double mod_z = -1.0 * Thickness / 2.0 + mThickness / 2.0 + 2.0 * mm;

  int k = -1;
  // Place Modules
  for (int j = 0; j < mNTowers; j++) {
    if (j == 0)
      mod_y = Width / 2.0 - mWidth / 2.0 - mGap;
    else
      mod_y -= (mWidth + mGap);

    if (abs(mod_y + mWidth / 2.0) > Width / 2.0)
      continue;

    mod_x = Width / 2.0 - (mWidth + mGap) * 0.5;

    for (int i = 0; i < mNTowers; i++) {
      if (i > 0)
        mod_x -= (mWidth + mGap);
      if (abs(mod_x + mWidth / 2.0) > Width / 2.0)
        continue;
      k++;
      string module_name  = detName + _toString(k, "_ECAL_Phys_%d");
      PlacedVolume pv_mod = detVol.placeVolume(modVol, Position(mod_x, mod_y, mod_z));
      pv_mod.addPhysVolID("module", k + 1);
    }
  }

  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);
  Transform3D tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);
  return det;
}
DECLARE_DETELEMENT(ZDC_ECAL, createDetector)
