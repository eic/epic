// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Dhevan Gangadharan

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
  xml_comp_t    x_dim           = x_det.child( _Unicode(dimensions) );
  xml_comp_t    x_pos           = x_det.child( _Unicode(position) );
  //
  string        det_name        = x_det.nameStr();
  string        mat_name        = getAttrOrDefault<string>( x_det, _Unicode(material), "Air" );
  //
  double        posZ1            = x_pos.attr<double>(_Unicode(z1));
  double        posZ2            = x_pos.attr<double>(_Unicode(z2));
  double        sizeX1           = x_dim.attr<double>(_Unicode(dx1));
  double        sizeX2           = x_dim.attr<double>(_Unicode(dx2));
  double        sizeY1           = x_dim.attr<double>(_Unicode(dy1));
  double        sizeY2           = x_dim.attr<double>(_Unicode(dy2));

  Trap trapezoid( fabs(posZ2-posZ1)/2., 0, 0, sizeY2, sizeX2, sizeX2, 0, sizeY1, sizeX1, sizeX1, 0);
  Volume vol( det_name + "_vol", trapezoid, description.material( mat_name ) );
  vol.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Transform3D  pos( RotationZYX(0, 0, 0), Position(0, 0, (posZ1+posZ2)/2.) );

  DetElement det(det_name, x_det.id());
  Volume motherVol = description.pickMotherVolume( det );
  PlacedVolume phv = motherVol.placeVolume( vol, pos );

  det.setPlacement(phv);

  return det;
}

DECLARE_DETELEMENT(LumiPhotonChamber, create_detector)
