// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Dhevan Gangadharan

//==========================================================================
//
// Places a chain of far-backward beamline magnets
//
//==========================================================================

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

static Ref_t build_magnet(Detector& description, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t x_det     = e;
  string     det_name  = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Iron      = description.material("Iron");
  string     vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "FFMagnetVis");


  for( xml_coll_t magnet_coll(x_det, _Unicode(magnet)); magnet_coll; magnet_coll++ ) { // magnets

    xml_comp_t magnet( magnet_coll );

    string name = getAttrOrDefault<string>(magnet, _Unicode(name), "");
    double x = getAttrOrDefault<double>(magnet, _Unicode(x), 0);
    double y = getAttrOrDefault<double>(magnet, _Unicode(y), 0);
    double z = getAttrOrDefault<double>(magnet, _Unicode(z), 0);
    double theta = getAttrOrDefault<double>(magnet, _Unicode(theta), 0);
    double length = getAttrOrDefault<double>(magnet, _Unicode(length), 0);
    double rin = getAttrOrDefault<double>(magnet, _Unicode(rin), 0);
    double rout = getAttrOrDefault<double>(magnet, _Unicode(rout), 0);

    // -- yoke
    Tube   yoke_tube( rin, rout, 0.5 * length );
    Volume v_yoke( "v_yoke_" + name, yoke_tube, m_Iron );

    v_yoke.setVisAttributes(description.visAttributes( vis_name ) );

    assembly.placeVolume(v_yoke, Transform3D( RotationY(theta), Position(x, y, z)));
  }

  // Final placement
  auto pv_assembly =
    description.pickMotherVolume(sdet).placeVolume( assembly, Position(0.0, 0.0, 0.0));

  sdet.setPlacement(pv_assembly);

  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(CylindricalDipoleMagnetChain, build_magnet)
