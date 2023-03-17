// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Dhevan Gangadharan

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /*sens*/)
{

  xml_det_t     x_det           = e;
  xml_comp_t    x_dim           = x_det.dimensions();
  xml_comp_t    x_pos           = x_det.position();
  xml_comp_t    x_rot           = x_det.rotation();
  //
  string        det_name        = x_det.nameStr();
  string        mat_name        = dd4hep::getAttrOrDefault<string>( x_det, _U(material), "Aluminum" );
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
  Volume vol( det_name + "_vol", box, description.material( mat_name ) );
  vol.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Transform3D  pos( RotationZYX(rotX, rotY, rotZ), Position(posX, posY, posZ) );

  DetElement det(det_name, x_det.id());
  Volume motherVol = description.pickMotherVolume( det );
  PlacedVolume phv = motherVol.placeVolume( vol, pos );

  det.setPlacement(phv);

  return det;
}

DECLARE_DETELEMENT(LumiWindow, create_detector)
