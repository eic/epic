// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Jaroslav Adam

//DD4hep
#include "DD4hep/DetFactoryHelper.h"

using namespace dd4hep;

static Ref_t create_element(Detector& lccdd, xml_h e, SensitiveDetector /*sens*/)
{
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();

  Info("ConeBeamEl", "Creating the segment: %s", name.c_str());

  //center position along z, mm
  double zpos = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(zpos), 0);

  //outer radius (rmax) and wall inner radius (rmin) at lower z (0) and higher z (1), mm
  double rmin0 = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(rmin0), 0);
  double rmax0 = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(rmax0), 0);
  double rmin1 = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(rmin1), 0);
  double rmax1 = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(rmax1), 0);

  //full length along z, mm
  double dz = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(dz), 0);

  //outer shape with vacuum
  ConeSegment shape_outer(dz/2., 0, rmax0, 0, rmax1);

  //outer volume
  Volume vol_outer("vol_outer", shape_outer, lccdd.material("Vacuum"));
  vol_outer.setVisAttributes(lccdd.visAttributes("VisFwElInvisible"));

  //Aluminum wall in outer volume
  ConeSegment shape_wall(dz/2., rmin0, rmax0, rmin1, rmax1);
  Volume vol_wall("vol_wall", shape_wall, lccdd.material("Aluminum"));
  vol_wall.setVisAttributes(lccdd.visAttributes("GrayVis"));

  //wall placement in outer volume
  vol_outer.placeVolume(vol_wall, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, 0)));

  //placement in top
  DetElement det(x_det.nameStr(), x_det.id());

  Transform3D pos(RotationZYX(0, 0, 0), Position(0, 0, zpos));
  PlacedVolume pv = lccdd.pickMotherVolume(det).placeVolume(vol_outer, pos);
  det.setPlacement(pv);

  return det;
}//create_element

DECLARE_DETELEMENT(ConeBeamEl, create_element)
