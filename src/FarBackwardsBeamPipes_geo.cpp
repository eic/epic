// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Dhevan Gangadharan

//==========================================================================
//
// Places far-backward beam pipes within and between the beamline magnets.
//
// Approximation used for beam pipes in between magnets.
// Right-angled ends with a small air gap to avoid G4 overlap errors
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Al     = description.material("Aluminum");
  Material   m_Vacuum = description.material("Vacuum");
  string     vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");
  double thickness = getAttrOrDefault<double>(x_det, _Unicode(wall_thickness), 0);

  vector<string> names;
  vector<double> xCenters;
  vector<double> zCenters;
  vector<double> lengths;
  vector<double> thetas;
  vector<double> rOuters1;
  vector<double> rOuters2;

  // Grab info for beamline magnets
  for( xml_coll_t pipe_coll(x_det, _Unicode(pipe)); pipe_coll; pipe_coll++ ) { // pipes

    xml_comp_t pipe( pipe_coll );

    names.push_back( getAttrOrDefault<string>(pipe, _Unicode(name), "") );
    int IsMagnetPipe = getAttrOrDefault<double>(pipe, _Unicode(magnetpipe), 0);

    if( IsMagnetPipe ) { // Beam pipes inside the magnets
      xCenters.push_back( getAttrOrDefault<double>(pipe, _Unicode(xcenter), 0) );
      zCenters.push_back( getAttrOrDefault<double>(pipe, _Unicode(zcenter), 0) );
      lengths.push_back( getAttrOrDefault<double>(pipe, _Unicode(length), 0) );
      thetas.push_back( getAttrOrDefault<double>(pipe, _Unicode(theta), 0) );
      rOuters1.push_back( getAttrOrDefault<double>(pipe, _Unicode(rout1), 0) );
      rOuters2.push_back( getAttrOrDefault<double>(pipe, _Unicode(rout2), 0) );
    }
    else{ // Left blank in xml and to be calculated below
      xCenters.push_back( 0 );
      zCenters.push_back( 0 );
      lengths.push_back( 0 );
      thetas.push_back( 0 );
      rOuters1.push_back( 0 );
      rOuters2.push_back( 0 );
    }
  }

  // Calculate parameters for connecting pipes in between magnets
  for( uint pipeN = 0; pipeN < names.size(); pipeN++ ) {

    if( ! lengths[pipeN] == 0 ) { continue; }
    if( int(pipeN-1) < 0 ) { continue; }

    double x = ( xCenters[pipeN-1] - lengths[pipeN-1]/2.*sin(thetas[pipeN-1]) + xCenters[pipeN+1] + lengths[pipeN+1]/2.*sin(thetas[pipeN+1]) ) / 2.;
    double z = ( zCenters[pipeN-1] - lengths[pipeN-1]/2.*cos(thetas[pipeN-1]) + zCenters[pipeN+1] + lengths[pipeN+1]/2.*cos(thetas[pipeN+1]) ) / 2.;
    double deltaX = (xCenters[pipeN-1] - lengths[pipeN-1]/2.*sin(thetas[pipeN-1])) - (xCenters[pipeN+1] + lengths[pipeN+1]/2.*sin(thetas[pipeN+1]));
    double deltaZ = (zCenters[pipeN-1] - lengths[pipeN-1]/2.*cos(thetas[pipeN-1])) - (zCenters[pipeN+1] + lengths[pipeN+1]/2.*cos(thetas[pipeN+1]));
    double l = sqrt( pow(deltaX, 2) + pow(deltaZ, 2) );

    // Small air gap between connecting and magnet beam pipes to avoid G4 overlap errors
    l -= 0.5;

    double theta = atan( deltaX / deltaZ );

    xCenters[pipeN] = x;
    zCenters[pipeN] = z;
    lengths[pipeN] = l;
    thetas[pipeN] = theta;
    rOuters1[pipeN] = rOuters2[pipeN-1];
    rOuters2[pipeN] = rOuters1[pipeN+1];
  }

  // Add all pipes to the assembly
  for( uint pipeN = 0; pipeN < xCenters.size(); pipeN++ ) {

    ConeSegment s_tube( lengths[pipeN] / 2.0, rOuters1[pipeN] - thickness, rOuters1[pipeN], rOuters2[pipeN] - thickness, rOuters2[pipeN] );
    ConeSegment s_vacuum( lengths[pipeN] / 2.0, 0, rOuters1[pipeN] - thickness, 0, rOuters2[pipeN] - thickness );

    Volume v_tube("v_tube_" + names[pipeN], s_tube, m_Al);
    Volume v_vacuum("v_vacuum_" + names[pipeN], s_vacuum, m_Vacuum);

    v_tube.setVisAttributes(description.visAttributes( vis_name ) );

    assembly.placeVolume(v_tube, Transform3D( RotationY(thetas[pipeN]), Position(xCenters[pipeN], 0, zCenters[pipeN])));
    assembly.placeVolume(v_vacuum, Transform3D( RotationY(thetas[pipeN]), Position(xCenters[pipeN], 0, zCenters[pipeN])));

  }

  // Final placement
  auto pv_assembly =
    description.pickMotherVolume(sdet).placeVolume( assembly, Position(0.0, 0.0, 0.0));

  sdet.setPlacement(pv_assembly);

  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(FarBackwardsBeamPipes, create_detector)
