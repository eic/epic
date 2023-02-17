// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Anna Kowalewska

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static dd4hep::Ref_t create_detector(dd4hep::Detector& description, xml_h e, dd4hep::SensitiveDetector sens)
{
  sens.setType("calorimeter");

  xml_det_t     x_det           = e;
  xml_comp_t    x_dim           = x_det.dimensions();
  xml_comp_t    x_pos           = x_det.position();
  xml_comp_t    x_rot           = x_det.rotation();
  //
  string        det_name        = x_det.nameStr();
  string        mat_name        = dd4hep::getAttrOrDefault<string>( x_det, _U(material), "PbWO4" );
  //
  double        sizeX           = x_dim.x();
  double        sizeY           = x_dim.y();
  double        sizeZ           = x_dim.z();
  double        posX            = x_pos.x();
  double        posY            = x_pos.y();
  double        posZ            = x_pos.z();
  double        rotX            = x_rot.x();
  double        rotY            = x_rot.y();
  double        rotZ            = x_rot.z();

  Box box( sizeX, sizeY, sizeZ );
  Volume vol( det_name + "_vol", box, desc.material( mat_name ) );
  //vol.setVisAttributes( x_det.visStr() );
  vol.setVisAttributes(desc.visAttributes(x_det.visStr()));
  vol.setSensitiveDetector(sens);

  Transform3D  pos( RotationZYX(rotX, rotY, rotZ), Position(posX, posY, posZ) );

  Volume motherVol = desc.pickMotherVolume( det );
  PlacedVolume phv = motherVol.placeVolume( vol, pos );
  phv.addPhysVolID("system", det_ID);
  det.setPlacement(phv);
  return det;
}

DECLARE_DETELEMENT(LumiDirect_PCAL, create_detector)
