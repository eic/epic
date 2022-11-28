// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Shima Shimizu, Jihee Kim

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
// Far Forward Ion Zero Degree Calorimeter - Ecal
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
  auto      xwidth = dim.x();
  auto      ywidth = dim.y();
  auto      length = dim.z();
  xml_dim_t pos    = detElem.position();
  xml_dim_t rot    = detElem.rotation();

  // envelope
  Box    envShape(xwidth * 0.5, ywidth * 0.5, length * 0.5);
  Volume env(detName + "_envelope", envShape, desc.material("Air"));
  env.setVisAttributes(desc.visAttributes(detElem.visStr()));

  // build frame
  xml_comp_t fr         = detElem.child(_Unicode(support));
  auto       fx         = xwidth;
  auto       fy         = ywidth;
  auto       fz         = fr.attr<double>(_Unicode(sizez));
  auto       fthickness = fr.attr<double>(_Unicode(thickness));

  Box    frShape(fx / 2., fy / 2., fz / 2.);
  auto   frMat = desc.material(fr.attr<std::string>(_Unicode(material)));
  Volume frVol("frame_vol", frShape, frMat);
  frVol.setVisAttributes(desc.visAttributes(fr.visStr()));

  xml_comp_t mod_x = detElem.child(_Unicode(module));
  auto       nx    = mod_x.attr<int>(_Unicode(nx));
  auto       ny    = mod_x.attr<int>(_Unicode(ny));

  // crystal tower
  xml_comp_t twr    = mod_x.child(_Unicode(tower));
  double     tsx    = twr.attr<double>(_Unicode(cellx));
  double     tsy    = twr.attr<double>(_Unicode(celly));
  double     tsz    = twr.thickness();
  Material   t_mat  = desc.material(twr.materialStr());
  string     t_name = twr.nameStr();

  Box    t_Shape(tsx / 2., tsy / 2., tsz / 2.);
  Volume t_Vol("tower_vol", t_Shape, t_mat);
  t_Vol.setVisAttributes(desc.visAttributes(twr.visStr()));
  if (twr.isSensitive())
    t_Vol.setSensitiveDetector(sens);

  // readout socket
  xml_comp_t sct    = mod_x.child(_Unicode(socket));
  double     ssx    = sct.attr<double>(_Unicode(cellx));
  double     ssy    = sct.attr<double>(_Unicode(celly));
  double     ssz    = sct.thickness();
  Material   s_mat  = desc.material(sct.materialStr());
  string     s_name = sct.nameStr();

  Box    s_Shape(ssx / 2., ssy / 2., ssz / 2.);
  Volume s_Vol("socket_vol", s_Shape, s_mat);
  s_Vol.setVisAttributes(desc.visAttributes(sct.visStr()));

  PlacedVolume pv;
  double       x_pos_0          = -(nx * tsx + (nx - 1) * fthickness) / 2.;
  double       y_pos_0          = -(ny * tsy + (ny - 1) * fthickness) / 2.;
  double       twr_z_pos_in_fr  = -fz / 2. + tsz / 2.;
  double       sct_z_pos_in_env = -length / 2. + tsz + ssz / 2.;
  int          mod_i            = 0;
  for (int ix = 0; ix < nx; ix++) {
    double x_pos = x_pos_0 + ix * (tsx + fthickness) + tsx / 2.;

    for (int iy = 0; iy < ny; iy++) {
      double y_pos = y_pos_0 + iy * (tsy + fthickness) + tsy / 2.;

      mod_i++;

      Position twr_pos(x_pos, y_pos, twr_z_pos_in_fr);
      pv = frVol.placeVolume(t_Vol, twr_pos);
      pv.addPhysVolID(t_name, mod_i);

      Position sct_pos(x_pos, y_pos, sct_z_pos_in_env);
      pv = env.placeVolume(s_Vol, sct_pos);
    }
  }

  double   f_zpos = -length / 2. + fz / 2.;
  Position fr_pos(0, 0, f_zpos);
  pv = env.placeVolume(frVol, fr_pos);

  // detector position and rotation
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(rot.z(), -rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume envPV = motherVol.placeVolume(env, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  return det;
}

DECLARE_DETELEMENT(ZDC_Crystal, create_detector)
