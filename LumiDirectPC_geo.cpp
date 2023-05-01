// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Yasir Ali

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <tuple>
//////////////////////////////////////////////////
// Reference from ATHENA ScFiCalorimeter_geo.cpp
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

// main
static Ref_t create_detector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string     detName = detElem.nameStr();
  int             detID   = detElem.id();
  DetElement      det(detName, detID);
  sens.setType("calorimeter");
  auto      dim    = detElem.dimensions();
  auto      width  = dim.x();
  //auto      height = dim.y();
  auto      length = dim.z();
  
  xml_dim_t pos    = detElem.position();
  xml_dim_t rot    = detElem.rotation();

  // envelope
  Box    box(width , width , length );
  Volume vol(detName + "_vol", box, desc.material("Air"));
  vol.setVisAttributes(desc.visAttributes(detElem.visStr()));

  // build module
  xml_comp_t mod_x = detElem.child(_Unicode(module));
  auto       sx    = mod_x.attr<double>(_Unicode(sizex));
  auto       sy    = mod_x.attr<double>(_Unicode(sizey));
  auto       sz    = mod_x.attr<double>(_Unicode(sizez));

  Box    modShape(sx/2, sy/2, sz/2);
  auto   modMat = desc.material(mod_x.attr<std::string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);
  modVol.setVisAttributes(desc.visAttributes(mod_x.visStr()));
  // modVol.setSensitiveDetector(sens);

  if (mod_x.hasChild(_Unicode(fiber))) {
    auto   fiber_x  = mod_x.child(_Unicode(fiber));
    auto   fr       = fiber_x.attr<double>(_Unicode(radius));
    auto   fsx      = fiber_x.attr<double>(_Unicode(spacex));
    auto   fsy      = fiber_x.attr<double>(_Unicode(spacey));
    auto   foff     = dd4hep::getAttrOrDefault<double>(fiber_x, _Unicode(offset), 0.5 * mm);
    auto   fiberMat = desc.material(fiber_x.attr<std::string>(_Unicode(material)));
    Tube   fiberShape(0., fr, sz / 2. - 1. * mm);
    Volume fiberVol("fiber_vol", fiberShape, fiberMat);
    fiberVol.setSensitiveDetector(sens);

    // Fibers are placed in a honeycomb with the radius = sqrt(3)/2. * hexagon side length
    // So each fiber is fully contained in a regular hexagon, which are placed as
    // the parameters space x and space y are used to add additional spaces between the hexagons
    double fside  = 2. / std::sqrt(3.) * fr;
    double fdistx = 2. * fside + fsx;
    double fdisty = 2. * fr + fsy;

    // maximum numbers of the fibers, help narrow the loop range
    int nx = int(sx / (2. * fr)) + 1;
    int ny = int(sy / (2. * fr)) + 1;

    // place the fibers
    double y0      = (foff + fside);
    int    nfibers = 0;
    for (int iy = 0; iy < ny; ++iy) {
      double y = y0 + fdisty * iy;
      // about to touch the boundary
      if ((sy - y) < y0) {
        break;
      }
      double x0 = (iy % 2) ? (foff + fside) : (foff + fside + fdistx / 2.);
      for (int ix = 0; ix < nx; ++ix) {
        double x = x0 + fdistx * ix;
        // about to touch the boundary
        if ((sx - x) < x0) {
          break;
        }
        auto fiberPV = modVol.placeVolume(fiberVol, nfibers++, Position{x - sx / 2., y - sy / 2., 0});
        fiberPV.addPhysVolID("fiber_x", ix + 1).addPhysVolID("fiber_y", iy + 1);
      }
    }
    // if no fibers we make the module itself sensitive
  } else {
    modVol.setSensitiveDetector(sens);
  }

  // Module Position
  double mod_x_pos = 0.0;
  double mod_y_pos = 0.0;
  double mod_z_pos = 0.0 * mm;
  double mgap      = 0.000001 * mm;
  int    mNTowers  = floor(width / (sx + mgap));
  // std::cout << "mNTowers: " << mNTowers << std::endl;

  int k = 0;
  // Place Modules
  for (int j = 0; j < mNTowers; j++) {
    if (j == 0)
      mod_y_pos = width * 0.5 - (sy + mgap) * 0.5;
    else
      mod_y_pos -= (sy + mgap);

    for (int i = 0; i < mNTowers; i++) {
      if (i == 0)
        mod_x_pos = width * 0.5 - (sx + mgap) * 0.5;
      else
        mod_x_pos -= (sx + mgap);

      PlacedVolume pv_mod = vol.placeVolume(modVol, Position(mod_x_pos, mod_y_pos, mod_z_pos));
      pv_mod.addPhysVolID("module", k++);
      // std::cout << "j: " << j << "i: " << i << " Position: " << mod_x_pos << " " << mod_y_pos << " " << mod_z_pos <<
      // std::endl;
    }
  }

  // detector position and rotation
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume envPV = motherVol.placeVolume(vol, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  return det;
}

DECLARE_DETELEMENT(LumiDirectPCAL, create_detector)
