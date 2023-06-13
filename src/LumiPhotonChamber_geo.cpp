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
  string        mat_pipeCap     = getAttrOrDefault<string>( x_det, _Unicode(pipeAndCapMaterial), "Aluminum" );
  string        mat_conv        = getAttrOrDefault<string>( x_det, _Unicode(convMaterial), "Aluminum" );
  string        mat_fill        = getAttrOrDefault<string>( x_det, _Unicode(fillMaterial), "Vacuum" );
  //
  double        posZ1            = x_pos.attr<double>(_Unicode(z1));
  double        posZ2            = x_pos.attr<double>(_Unicode(z2));
  double        posZconv         = x_pos.attr<double>(_Unicode(z_conv));
  double        rmin             = x_dim.attr<double>(_Unicode(rmin));
  double        pipe_DR          = x_dim.attr<double>(_Unicode(pipe_dr));
  double        cap1_DZ          = x_dim.attr<double>(_Unicode(cap1_dz));
  double        cap2_DZ          = x_dim.attr<double>(_Unicode(cap2_dz));
  double        conv_DZ          = x_dim.attr<double>(_Unicode(conv_dz));

  // Create main detector element to be returned at the end
  DetElement    det(det_name, x_det.id());
  
  // Mother volume
  Volume        motherVol = description.pickMotherVolume( det );

  // chamber assembly
  Assembly      assembly( det_name );
  assembly.setVisAttributes( description.invisible() );

  ////////////
  // G4 solids
  // main tube
  Tube tube( rmin, rmin + pipe_DR, (posZ1 - posZ2)/2.0, 0, 2*TMath::Pi() );
  // vacuum regions inside tube
  Tube vac1( 0, rmin, (posZ1 - posZconv - conv_DZ/2.0)/2.0, 0, 2*TMath::Pi() );
  Tube vac2( 0, rmin, (posZconv - posZ2 - conv_DZ/2.0)/2.0, 0, 2*TMath::Pi() );

  // end cap closest to IP
  Tube cap1( 0, rmin + pipe_DR, cap1_DZ/2.0, 0, 2*TMath::Pi() );
  // end cap farthest from IP
  Tube cap2( 0, rmin + pipe_DR, cap2_DZ/2.0, 0, 2*TMath::Pi() );
  // conversion foil
  Tube convFoil( 0, rmin, conv_DZ/2.0, 0, 2*TMath::Pi() );

  //////////
  // volumes
  Volume vol_vac1( det_name + "_vol_vac1", vac1, description.material( mat_fill ) );
  vol_vac1.setVisAttributes( description.invisible() );
  Volume vol_vac2( det_name + "_vol_vac2", vac2, description.material( mat_fill ) );
  vol_vac2.setVisAttributes( description.invisible() );

  Volume vol_tube( det_name + "_vol_tube", tube, description.material( mat_pipeCap ) );
  vol_tube.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Volume vol_cap1( det_name + "_vol_cap1", cap1, description.material( mat_pipeCap ) );
  vol_cap1.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Volume vol_cap2( det_name + "_vol_cap2", cap2, description.material( mat_pipeCap ) );
  vol_cap2.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Volume vol_conv( det_name + "_vol_conversionFoil", convFoil, description.material( mat_conv ) );
  vol_conv.setVisAttributes( description.visAttributes(x_det.visStr()) );

  // place each volume into assembly
  assembly.placeVolume(
          vol_tube, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, (posZ1 + posZ2)/2.)) );
  assembly.placeVolume(
          vol_cap1, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, posZ1 + cap1_DZ/2.)) );
  assembly.placeVolume(
          vol_cap2, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, posZ2 - cap2_DZ/2.)) );
  assembly.placeVolume(
          vol_conv, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, posZconv)) );
  assembly.placeVolume(
          vol_vac1, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, (posZ1 + posZconv)/2.)) );
  assembly.placeVolume(
          vol_vac2, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, (posZ2 + posZconv)/2.)) );

  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume phv = motherVol.placeVolume( assembly, Position(0.0, 0.0, 0.0) );

  det.setPlacement(phv);

  return det;
}

DECLARE_DETELEMENT(LumiPhotonChamber, create_detector)
