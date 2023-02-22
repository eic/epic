// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Yasir Ali

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <tuple>

using namespace std;
using namespace dd4hep;
// Definition of function to build the modules
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /*sens*/)
{
  //sens.setType("calorimeter");
  xml_det_t     x_det           = e;
  //xml_comp_t    x_dim           = x_det.dimensions();
  xml_comp_t    x_pos           = x_det.position();
  xml_comp_t    x_rot           = x_det.rotation();

  //
  string        det_name        = x_det.nameStr();
  string        mat_name        = dd4hep::getAttrOrDefault<string>( x_det, _U(material), "Graphite" );
  //double        filterThickness = 5*cm;
  //double        filterDistanceZ = 10*cm;
  //
  double        sizeX           = 10*cm;
  double        sizeY           = 10*cm;
  double        sizeZ           = 5*cm;
  double        posX            = x_pos.x();
  double        posY            = x_pos.y();
  double        posZ            = x_pos.z();
  double        rotX            = x_rot.x();
  double        rotY            = x_rot.y();
  double        rotZ            = x_rot.z();

  Box slice( sizeX, sizeY, sizeZ );
  Volume vol( det_name + "_vol", slice, description.material( mat_name ) );
  vol.setVisAttributes( x_det.visStr() );

  Transform3D  filter_pos( RotationZYX(rotX, rotY, rotZ), Position(posX, posY, posZ ) );

  DetElement det(det_name, x_det.id());
  Volume motherVol = description.pickMotherVolume( det );
  PlacedVolume phv = motherVol.placeVolume( vol, filter_pos );

  det.setPlacement(phv);

  return det;
}

DECLARE_DETELEMENT(LumiDirectPCFilter, create_detector)
