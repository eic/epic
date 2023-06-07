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
  string        matpipe_name    = getAttrOrDefault<string>( x_det, _Unicode(pipematerial), "Al" );
  string        matfill_name    = getAttrOrDefault<string>( x_det, _Unicode(fillmaterial), "Vacuum" );
  //
  double        posZ1            = x_pos.attr<double>(_Unicode(z1));
  double        posZ2            = x_pos.attr<double>(_Unicode(z2));
  double        rmin             = x_dim.attr<double>(_Unicode(rmin));
  double        pipe_DR          = x_dim.attr<double>(_Unicode(pipe_dr));
  double        cap1_DZ          = x_dim.attr<double>(_Unicode(cap1_dz));
  double        cap2_DZ          = x_dim.attr<double>(_Unicode(cap2_dz));

  // Create main detector element to be returned at the end
  DetElement    det(det_name, x_det.id());
  
  // Mother volume
  Volume        motherVol = description.pickMotherVolume( det );

  // chamber assembly
  Assembly      assembly( det_name );
  assembly.setVisAttributes( description.invisible() );

  //////////
  // G4 solids
  // photon chamber that encapsulates tube + caps
  Trap chamber( posZ1 - posZ2 + cap1_DZ + cap2_DZ, 2*(rmin + pipe_DR), 2*(rmin + pipe_DR), 2*(rmin + pipe_DR) );
  // main tube
  Tube tube( rmin, rmin + pipe_DR, (posZ1 - posZ2)/2.0, 0, 2*TMath::Pi() );
  // end cap closest to IP
  Tube cap1( 0, rmin + pipe_DR, cap1_DZ/2.0, 0, 2*TMath::Pi() );
  // end cap farthest from IP
  Tube cap2( 0, rmin + pipe_DR, cap2_DZ/2.0, 0, 2*TMath::Pi() );

  //////////
  // volumes
  Volume vol_chamber( det_name + "_vol_chamber", chamber, description.material( matfill_name ) );
  vol_chamber.setVisAttributes( description.invisible() );
  
  Volume vol_tube( det_name + "_vol_tube", tube, description.material( matpipe_name ) );
  vol_tube.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Volume vol_cap1( det_name + "_vol_cap1", cap1, description.material( matpipe_name ) );
  vol_cap1.setVisAttributes( description.visAttributes(x_det.visStr()) );

  Volume vol_cap2( det_name + "_vol_cap2", cap2, description.material( matpipe_name ) );
  vol_cap2.setVisAttributes( description.visAttributes(x_det.visStr()) );

  // place each volume into assembly
  assembly.placeVolume(
          vol_chamber, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, (posZ1 + posZ2)/2.)) );
  assembly.placeVolume(
          vol_tube, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, (posZ1 + posZ2)/2.)) );
  assembly.placeVolume(
          vol_cap1, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, posZ1 + cap1_DZ/2.)) );
  assembly.placeVolume(
          vol_cap2, Transform3D( RotationZYX(0, 0, 0), Position(0, 0, posZ2 - cap2_DZ/2.)) );
  
  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume phv = motherVol.placeVolume( assembly, Position(0.0, 0.0, 0.0) );

  det.setPlacement(phv);

  return det;
}

DECLARE_DETELEMENT(LumiPhotonChamber, create_detector)
