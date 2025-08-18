// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Dhevan Gangadharan, Simon Gardner

//==========================================================================
//
// Places a chain of beam pipe segments within and between the beamline magnets.
//
// Approximation used for beam pipes in between magnets.
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /* sens */) {

  using namespace ROOT::Math;
  xml_det_t x_det = e;
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly assembly(det_name + "_assembly");
  Material m_Al     = description.material("Aluminum");
  Material m_Vacuum = description.material("Vacuum");
  string vis_name   = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");
  double thickness  = getAttrOrDefault<double>(x_det, _Unicode(wall_thickness), 0);
  double bendRadius = 0.00001; //Small bend radius to allow the construction of toruses

  vector<string> names;
  vector<int> ids;
  vector<double> xCenters;
  vector<double> zCenters;
  vector<double> lengths;
  vector<double> thetas;
  vector<double> rOuters1;
  vector<double> rOuters2;
  vector<double> bendLengths;

  // Grab info for beamline magnets
  for (xml_coll_t pipe_coll(x_det, _Unicode(pipe)); pipe_coll; pipe_coll++) { // pipes

    xml_comp_t pipe(pipe_coll);

    names.push_back(getAttrOrDefault<string>(pipe, _Unicode(name), ""));
    ids.push_back(getAttrOrDefault<int>(pipe, _Unicode(id), 0));

    // Vectors momentarily filled with zeros for pipes in between magnets
    xCenters.push_back(getAttrOrDefault<double>(pipe, _Unicode(xcenter), 0));
    zCenters.push_back(getAttrOrDefault<double>(pipe, _Unicode(zcenter), 0));
    lengths.push_back(getAttrOrDefault<double>(pipe, _Unicode(length), 0));
    thetas.push_back(getAttrOrDefault<double>(pipe, _Unicode(theta), 0));
    rOuters1.push_back(getAttrOrDefault<double>(pipe, _Unicode(rout1), 0));
    rOuters2.push_back(getAttrOrDefault<double>(pipe, _Unicode(rout2), 0));
  }

  // Calculate parameters for connecting pipes in between magnets
  for (uint pipeN = 0; pipeN < names.size(); pipeN++) {

    if (lengths[pipeN] > 0) {
      continue;
    } // pipe parameters already set to nonzero values
    if (pipeN == 0) {
      continue;
    } // can't create pipe for an empty starting slot
    if ((pipeN + 1) == names.size()) {
      continue;
    } // can't create pipe for an empty end slot

    double previousPipeEndX =
        xCenters[pipeN - 1] - lengths[pipeN - 1] / 2. * sin(thetas[pipeN - 1]);
    double previousPipeEndZ =
        zCenters[pipeN - 1] - lengths[pipeN - 1] / 2. * cos(thetas[pipeN - 1]);
    double nextPipeStartX = xCenters[pipeN + 1] + lengths[pipeN + 1] / 2. * sin(thetas[pipeN + 1]);
    double nextPipeStartZ = zCenters[pipeN + 1] + lengths[pipeN + 1] / 2. * cos(thetas[pipeN + 1]);

    double x      = (previousPipeEndX + nextPipeStartX) / 2.;
    double z      = (previousPipeEndZ + nextPipeStartZ) / 2.;
    double deltaX = previousPipeEndX - nextPipeStartX;
    double deltaZ = previousPipeEndZ - nextPipeStartZ;
    double l      = sqrt(deltaX * deltaX + deltaZ * deltaZ);
    double theta  = atan2(deltaX, deltaZ);

    xCenters[pipeN] = x;
    zCenters[pipeN] = z;
    lengths[pipeN]  = l;
    thetas[pipeN]   = theta;
    rOuters1[pipeN] = rOuters2[pipeN - 1];
    rOuters2[pipeN] = rOuters1[pipeN + 1];
  }

  // If there is an bend in the pipe, calculate the length reduction of the pipe and joint length
  for (uint i = 1; i < thetas.size(); i++) {
    // Start at the join between the first two pipes ending at the join between the last two pipes N-1
    if (thetas[i - 1] == thetas[i]) {
      bendLengths.push_back(0);
    } else // Correct for tubes, not yet cones so imperfect
    {
      double bendAngle  = thetas[i] - thetas[i - 1];
      double bendLength = abs(rOuters1[i] * tan(bendAngle / 2));
      bendLengths.push_back(bendLength + bendRadius);
    }
  }

  // Add all pipes to the assembly
  for (uint pipeN = 0; pipeN < xCenters.size(); pipeN++) {

    double length  = lengths[pipeN];
    double theta   = thetas[pipeN];
    double xCenter = xCenters[pipeN];
    double zCenter = zCenters[pipeN];
    double rOuter1 = rOuters1[pipeN];
    double rOuter2 = rOuters2[pipeN];

    //Change straight pipe length to account for bend
    if (pipeN != 0) {
      length -= bendLengths[pipeN - 1];
      xCenter -= bendLengths[pipeN - 1] / 2 * sin(theta);
      zCenter -= bendLengths[pipeN - 1] / 2 * cos(theta);
    }
    if (pipeN != xCenters.size() - 1) {
      length -= bendLengths[pipeN];
      xCenter += bendLengths[pipeN] / 2 * sin(theta);
      zCenter += bendLengths[pipeN] / 2 * cos(theta);
    }

    ConeSegment s_tube(length / 2.0, rOuter2 - thickness, rOuter2, rOuter1 - thickness, rOuter1);
    ConeSegment s_vacuum(length / 2.0, 0, rOuter2 - thickness, 0, rOuter1 - thickness);

    Volume v_tube("v_tube_" + names[pipeN], s_tube, m_Al);
    Volume v_vacuum("v_vacuum_" + names[pipeN], s_vacuum, m_Vacuum);

    v_tube.setVisAttributes(description.visAttributes(vis_name));

    assembly.placeVolume(v_tube, Transform3D(RotationY(theta), Position(xCenter, 0, zCenter)));
    auto placed_vacuum = assembly.placeVolume(
        v_vacuum, Transform3D(RotationY(theta), Position(xCenter, 0, zCenter)));

    DetElement vacuum_element(sdet, names[pipeN] + "_vacuum", ids[pipeN]);
    vacuum_element.setPlacement(placed_vacuum);

    // // Add joining bend to next pipe
    if (pipeN != xCenters.size() - 1 && bendLengths[pipeN] != 0) {

      // double bendLength = bendLengths[pipeN];
      double bendAngle = thetas[pipeN + 1] - thetas[pipeN];

      // Create a torus to join the two pipes
      Torus s_bend_solid(rOuter2 + bendRadius, rOuter2 - thickness, rOuter2, 0.0, abs(bendAngle));

      //Create a vacuum torus to place inside the solid segment
      Torus s_bend_vac(rOuter2 + bendRadius, 0, rOuter2 - thickness, 0, abs(bendAngle));

      // //Create the volumes
      Volume v_bend_soild("v_bend_solid_" + names[pipeN], s_bend_solid, m_Al);
      Volume v_bend_vac("v_bend_vac_" + names[pipeN], s_bend_vac, m_Vacuum);

      // Calculate the position and rotation to place bend at the joint
      double bendCenterX = xCenter + rOuter2 * cos(theta) - (length / 2) * sin(theta);
      double bendCenterZ = zCenter - rOuter2 * sin(theta) - (length / 2) * cos(theta);
      double rotation    = pi - theta;
      if (bendAngle > 0) {
        bendCenterX = xCenter - rOuter2 * cos(theta) - (length / 2) * sin(theta);
        bendCenterZ = zCenter + rOuter2 * sin(theta) - (length / 2) * cos(theta);
        rotation    = -bendAngle - theta;
      }

      // Place the bend in the assembly
      assembly.placeVolume(v_bend_soild, Transform3D(RotationZYX(rotation, 0, pi / 2),
                                                     Position(bendCenterX, 0, bendCenterZ)));
      assembly.placeVolume(v_bend_vac, Transform3D(RotationZYX(rotation, 0, pi / 2),
                                                   Position(bendCenterX, 0, bendCenterZ)));

      // Set vis attributes
      v_bend_soild.setVisAttributes(description.visAttributes(vis_name));
    }
  }

  // Final placement
  auto pv_assembly =
      description.pickMotherVolume(sdet).placeVolume(assembly, Position(0.0, 0.0, 0.0));

  sdet.setPlacement(pv_assembly);

  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(BeamPipeChain, create_detector)
