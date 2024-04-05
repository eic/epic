// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Dhevan Gangadharan

//==========================================================================
//
// Places a chain of beam pipe segments within and between the beamline magnets.
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
  Material   m_SS     = description.material("StainlessSteel");
  Material   m_Cu     = description.material("Copper");
  Material   m_Vacuum = description.material("Vacuum");
  string     vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");
  double thickness = getAttrOrDefault<double>(x_det, _Unicode(wall_thickness), 0);
  double coating = getAttrOrDefault<double>(x_det, _Unicode(coating_thickness), 0);

  // beam pipe
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
    // Vectors momentarily filled with zeros for pipes in between magnets
    xCenters.push_back( getAttrOrDefault<double>(pipe, _Unicode(xcenter), 0) );
    zCenters.push_back( getAttrOrDefault<double>(pipe, _Unicode(zcenter), 0) );
    lengths.push_back( getAttrOrDefault<double>(pipe, _Unicode(length), 0) );
    thetas.push_back( getAttrOrDefault<double>(pipe, _Unicode(theta), 0) );
    rOuters1.push_back( getAttrOrDefault<double>(pipe, _Unicode(rout1), 0) );
    rOuters2.push_back( getAttrOrDefault<double>(pipe, _Unicode(rout2), 0) );
  }

  // Calculate parameters for connecting pipes in between magnets
  for( uint pipeN = 0; pipeN < names.size(); pipeN++ ) {

    if( lengths[pipeN] > 0 ) { continue; } // pipe parameters already set to nonzero values
    if( pipeN == 0 ) { continue; } // can't create pipe for an empty starting slot
    if( (pipeN+1) == names.size() ) { continue; } // can't create pipe for an empty end slot

    double x = ( xCenters[pipeN-1] - lengths[pipeN-1]/2.*sin(thetas[pipeN-1]) + xCenters[pipeN+1] + lengths[pipeN+1]/2.*sin(thetas[pipeN+1]) ) / 2.;
    double z = ( zCenters[pipeN-1] - lengths[pipeN-1]/2.*cos(thetas[pipeN-1]) + zCenters[pipeN+1] + lengths[pipeN+1]/2.*cos(thetas[pipeN+1]) ) / 2.;
    double deltaX = (xCenters[pipeN-1] - lengths[pipeN-1]/2.*sin(thetas[pipeN-1])) - (xCenters[pipeN+1] + lengths[pipeN+1]/2.*sin(thetas[pipeN+1]));
    double deltaZ = (zCenters[pipeN-1] - lengths[pipeN-1]/2.*cos(thetas[pipeN-1])) - (zCenters[pipeN+1] + lengths[pipeN+1]/2.*cos(thetas[pipeN+1]));
    double l = sqrt( pow(deltaX, 2) + pow(deltaZ, 2) );
    double theta = atan( deltaX / deltaZ );

    // Small air gap between connecting and magnet beam pipes to avoid G4 overlap errors
    if( (theta != thetas[pipeN-1]) || (theta != thetas[pipeN+1]) ) {
      l -= 0.5;
    }

    xCenters[pipeN] = x;
    zCenters[pipeN] = z;
    lengths[pipeN] = l;
    thetas[pipeN] = theta;
    rOuters1[pipeN] = rOuters2[pipeN-1];
    rOuters2[pipeN] = rOuters1[pipeN+1];
  }

  // Add all pipes to the assembly
  for( uint pipeN = 0; pipeN < xCenters.size(); pipeN++ ) {
    // beam pipe matter
    ConeSegment s_tube(  
                         lengths[pipeN] / 2.0,         // length 
                         rOuters2[pipeN] - thickness,  // r2 in
                         rOuters2[pipeN],              // r2 out
                         rOuters1[pipeN] - thickness,  // r1 in
                         rOuters1[pipeN]               // r1 out 
                      );
    // beam pipe coating
    ConeSegment s_coat(  
                         lengths[pipeN] / 2.0,                   // length 
                         rOuters2[pipeN] - thickness - coating,  // r2 in
                         rOuters2[pipeN] - thickness,            // r2 out
                         rOuters1[pipeN] - thickness - coating,  // r1 in
                         rOuters1[pipeN] - thickness             // r1 out 
                      );
    // beam pipe vacuum
    ConeSegment s_vacm(  
                         lengths[pipeN] / 2.0,                  // length 
                         0,                                     // r2 in
                         rOuters2[pipeN] - thickness - coating, // r2 out
                         0,                                     // r1 in
                         rOuters1[pipeN] - thickness - coating  // r1 out
                      );

    Volume v_tube("v_tube_"    + names[pipeN], s_tube, m_SS);
    Volume v_coat("v_coating_" + names[pipeN], s_coat, m_Cu);
    Volume v_vacm("v_vacuum_"  + names[pipeN], s_vacm, m_Vacuum);

    v_tube.setVisAttributes(description.visAttributes( vis_name ) );
    v_coat.setVisAttributes(description.visAttributes( vis_name + "Coating") );

    assembly.placeVolume(v_tube, Transform3D( RotationY(thetas[pipeN]), Position(xCenters[pipeN], 0, zCenters[pipeN])));
    assembly.placeVolume(v_coat, Transform3D( RotationY(thetas[pipeN]), Position(xCenters[pipeN], 0, zCenters[pipeN])));
    assembly.placeVolume(v_vacm, Transform3D( RotationY(thetas[pipeN]), Position(xCenters[pipeN], 0, zCenters[pipeN])));
  }

  // Final placement
  auto pv_assembly =
    description.pickMotherVolume(sdet).placeVolume( assembly, Position(0.0, 0.0, 0.0));

  sdet.setPlacement(pv_assembly);

  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(BeamPipeChain, create_detector)
