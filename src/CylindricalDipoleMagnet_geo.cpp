// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Alex Jentsch, Wouter Deconinck, Whitney Armstrong

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "TMath.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace ROOT::Math;

static Ref_t build_magnet(Detector& dtor, xml_h e, SensitiveDetector /* sens */) {
  xml_det_t x_det     = e;
  int det_id          = x_det.id();
  string det_name     = x_det.nameStr();
  xml_dim_t pos       = x_det.child(_U(placement));
  double pos_x        = pos.x();
  double pos_y        = pos.y();
  double pos_z        = pos.z();
  double pos_theta    = pos.attr<double>(_U(theta));
  xml_dim_t dims      = x_det.dimensions();
  double dim_r        = dims.r();
  double dim_z        = dims.z();
  xml_dim_t apperture = x_det.child(_Unicode(apperture));
  double app_r        = apperture.r();
  Material iron       = dtor.material("Iron");

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name + "_assembly");

  const string module_name = "Quad_magnet";

  const string yoke_vis =
      dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "FFMagnetVis");

  sdet.setAttributes(dtor, assembly, x_det.regionStr(), x_det.limitsStr(), yoke_vis);

  // -- yoke
  Tube yoke_tube(app_r, dim_r, 0.5 * dim_z);
  Volume yoke_vol("yoke_vol", yoke_tube, iron);
  auto yoke_pv = assembly.placeVolume(yoke_vol);
  yoke_pv.addPhysVolID("element", 1);
  DetElement yoke_de(sdet, "yoke_de", 1);
  yoke_de.setPlacement(yoke_pv);
  yoke_de.setAttributes(dtor, yoke_vol, x_det.regionStr(), x_det.limitsStr(), yoke_vis);

  // -- finishing steps
  auto final_pos = Transform3D(Translation3D(pos_x, pos_y, pos_z) * RotationY(pos_theta));
  auto pv        = dtor.pickMotherVolume(sdet).placeVolume(assembly, final_pos);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(ip6_CylindricalDipoleMagnet, build_magnet)
