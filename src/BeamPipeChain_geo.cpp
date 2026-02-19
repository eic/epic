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
#include <algorithm>

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
  double bendRadius = 0.00001 * mm; //Small bend radius to allow the construction of toruses

  vector<string> names;
  vector<int> ids;
  vector<double> xCenters;
  vector<double> zCenters;
  vector<double> lengths;
  vector<double> thetas;
  vector<double> rOuters1;
  vector<double> rOuters2;
  vector<double> bendLengths;

  // Optional boolean-combine groups, configured via XML:
  //   <combine name="..." start_id="..." end_id="..."/>
  // All pieces (straight and bend) whose pipe id falls in the range can be merged
  // into one unified pipe/vacuum boolean solid.
  vector<int> combineStartIds;
  vector<int> combineEndIds;
  vector<std::string> combineNames;

  for (xml_coll_t combine_coll(x_det, _Unicode(combine)); combine_coll; ++combine_coll) {
    xml_comp_t combine(combine_coll);
    int start_id = getAttrOrDefault<int>(combine, _Unicode(start_id), 0);
    int end_id   = getAttrOrDefault<int>(combine, _Unicode(end_id), 0);
    if (end_id < start_id) {
      std::swap(start_id, end_id);
    }
    std::string name =
        getAttrOrDefault<std::string>(combine, _Unicode(name), "combined_" + std::to_string(start_id) + "_" + std::to_string(end_id));

    combineStartIds.push_back(start_id);
    combineEndIds.push_back(end_id);
    combineNames.push_back(name);
  }

  // Return the combine-range index for a pipe id, or -1 if not in a combine group.
  auto find_range_index = [&](int pipe_id) {
    for (size_t range_n = 0; range_n < combineStartIds.size(); ++range_n) {
      if (pipe_id >= combineStartIds[range_n] && pipe_id <= combineEndIds[range_n]) {
        return static_cast<int>(range_n);
      }
    }
    return -1;
  };

  // Parallel vectors storing prebuilt solids and transforms for straight segments.
  // We keep these separate from placement so combined and non-combined paths can
  // share the same geometry source without duplication.
  vector<std::string> straightNames;
  vector<int> straightIds;
  vector<int> straightRangeIndices;
  vector<Solid> straightTubes;
  vector<Solid> straightVacuums;
  vector<Transform3D> straightTransforms;

  // Parallel vectors for bend-join torus pieces between neighboring straight pipes.
  vector<std::string> bendNames;
  vector<int> bendIds;
  vector<int> bendRangeIndices;
  vector<Solid> bendTubes;
  vector<Solid> bendVacuums;
  vector<Transform3D> bendTransforms;

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
    double bendAngle = thetas[i] - thetas[i - 1];
    if (std::abs(bendAngle) < 0.01 * mrad) {
      bendLengths.push_back(0);
    } else // Correct for tubes, not yet cones so imperfect
    {
      double bendLength = abs(rOuters1[i] * tan(bendAngle / 2));
      bendLengths.push_back(bendLength + bendRadius);
    }
  }

  // Build all straight and bend pieces first, then decide later whether each piece
  // is placed as an individual volume or absorbed into a combined boolean group.
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

    Transform3D straightTransform(RotationY(theta), Position(xCenter, 0, zCenter));
    straightNames.push_back(names[pipeN]);
    straightIds.push_back(ids[pipeN]);
    straightRangeIndices.push_back(find_range_index(ids[pipeN]));
    straightTubes.push_back(Solid(s_tube));
    straightVacuums.push_back(Solid(s_vacuum));
    straightTransforms.push_back(straightTransform);

    // // Add joining bend to next pipe
    if (pipeN != xCenters.size() - 1 && bendLengths[pipeN] != 0) {

      // double bendLength = bendLengths[pipeN];
      double bendAngle = thetas[pipeN + 1] - thetas[pipeN];

      // Create a torus to join the two pipes
      Torus s_bend_solid(rOuter2 + bendRadius, rOuter2 - thickness, rOuter2, 0.0, abs(bendAngle));

      //Create a vacuum torus to place inside the solid segment
      Torus s_bend_vac(rOuter2 + bendRadius, 0, rOuter2 - thickness, 0, abs(bendAngle));

      // Calculate the position and rotation to place bend at the joint
      double bendCenterX = xCenter + rOuter2 * cos(theta) - (length / 2) * sin(theta);
      double bendCenterZ = zCenter - rOuter2 * sin(theta) - (length / 2) * cos(theta);
      double rotation    = pi - theta;
      if (bendAngle > 0) {
        bendCenterX = xCenter - rOuter2 * cos(theta) - (length / 2) * sin(theta);
        bendCenterZ = zCenter + rOuter2 * sin(theta) - (length / 2) * cos(theta);
        rotation    = -bendAngle - theta;
      }

      Transform3D bendTransform(RotationZYX(rotation, 0, pi / 2),
                                Position(bendCenterX, 0, bendCenterZ));

      // A bend belongs to a combine group only if BOTH neighboring straight
      // segments are in the same range. This keeps range boundaries clean.
      int bend_range_index = -1;
      int this_range       = find_range_index(ids[pipeN]);
      int next_range       = find_range_index(ids[pipeN + 1]);
      if (this_range >= 0 && this_range == next_range) {
        bend_range_index = this_range;
      }

      bendNames.push_back(names[pipeN] + "_to_" + names[pipeN + 1]);
      bendIds.push_back(ids[pipeN]);
      bendRangeIndices.push_back(bend_range_index);
      bendTubes.push_back(Solid(s_bend_solid));
      bendVacuums.push_back(Solid(s_bend_vac));
      bendTransforms.push_back(bendTransform);
    }
  }

  // Place non-combined straight segments
  for (size_t piece_n = 0; piece_n < straightNames.size(); ++piece_n) {
    if (straightRangeIndices[piece_n] >= 0) {
      continue;
    }

    Volume v_tube("v_tube_" + straightNames[piece_n], straightTubes[piece_n], m_Al);
    Volume v_vacuum("v_vacuum_" + straightNames[piece_n], straightVacuums[piece_n], m_Vacuum);
    v_tube.setVisAttributes(description.visAttributes(vis_name));

    assembly.placeVolume(v_tube, straightTransforms[piece_n]);
    auto placed_vacuum = assembly.placeVolume(v_vacuum, straightTransforms[piece_n]);

    DetElement vacuum_element(sdet, straightNames[piece_n] + "_vacuum", straightIds[piece_n]);
    vacuum_element.setPlacement(placed_vacuum);
  }

  // Place non-combined bends
  for (size_t piece_n = 0; piece_n < bendNames.size(); ++piece_n) {
    if (bendRangeIndices[piece_n] >= 0) {
      continue;
    }

    Volume v_bend_soild("v_bend_solid_" + bendNames[piece_n], bendTubes[piece_n], m_Al);
    Volume v_bend_vac("v_bend_vac_" + bendNames[piece_n], bendVacuums[piece_n], m_Vacuum);
    v_bend_soild.setVisAttributes(description.visAttributes(vis_name));

    assembly.placeVolume(v_bend_soild, bendTransforms[piece_n]);
    assembly.placeVolume(v_bend_vac, bendTransforms[piece_n]);
  }

  // Build and place combined ranges as unified boolean solids.
  // Strategy:
  //  1) pick one piece as the local boolean reference frame (base_transform)
  //  2) union all matching tube pieces in that local frame
  //  3) union all matching vacuum pieces in that local frame
  //  4) subtract vacuum-union from tube-union for the final pipe shell
  for (size_t range_n = 0; range_n < combineStartIds.size(); ++range_n) {
    bool has_piece = false;
    bool base_from_straight = false;
    size_t base_index = 0;

    for (size_t piece_n = 0; piece_n < straightNames.size(); ++piece_n) {
      if (straightRangeIndices[piece_n] == static_cast<int>(range_n)) {
        has_piece = true;
        base_from_straight = true;
        base_index = piece_n;
        break;
      }
    }
    if (!has_piece) {
      for (size_t piece_n = 0; piece_n < bendNames.size(); ++piece_n) {
        if (bendRangeIndices[piece_n] == static_cast<int>(range_n)) {
          has_piece = true;
          base_from_straight = false;
          base_index = piece_n;
          break;
        }
      }
    }

    if (!has_piece) {
      printout(WARNING, "BeamPipeChain",
               "Combine range '%s' (%d-%d) did not match any segments.",
               combineNames[range_n].c_str(), combineStartIds[range_n], combineEndIds[range_n]);
      continue;
    }

    // Use first matching piece as boolean base frame.
    const Transform3D& base_transform =
        base_from_straight ? straightTransforms[base_index] : bendTransforms[base_index];
    Transform3D inv_base              = base_transform.Inverse();
    Solid range_tube = base_from_straight ? straightTubes[base_index] : bendTubes[base_index];
    Solid range_vacuum =
        base_from_straight ? straightVacuums[base_index] : bendVacuums[base_index];

    // Add all straight pieces in this range into the base frame.
    for (size_t piece_n = 0; piece_n < straightNames.size(); ++piece_n) {
      if (straightRangeIndices[piece_n] != static_cast<int>(range_n)) {
        continue;
      }
      if (base_from_straight && piece_n == base_index) {
        continue;
      }
      Transform3D rel = inv_base * straightTransforms[piece_n];
      range_tube      = UnionSolid(range_tube, straightTubes[piece_n], rel);
      range_vacuum    = UnionSolid(range_vacuum, straightVacuums[piece_n], rel);
    }
    // Add all bend pieces in this range into the same base frame.
    for (size_t piece_n = 0; piece_n < bendNames.size(); ++piece_n) {
      if (bendRangeIndices[piece_n] != static_cast<int>(range_n)) {
        continue;
      }
      if (!base_from_straight && piece_n == base_index) {
        continue;
      }
      Transform3D rel = inv_base * bendTransforms[piece_n];
      range_tube      = UnionSolid(range_tube, bendTubes[piece_n], rel);
      range_vacuum    = UnionSolid(range_vacuum, bendVacuums[piece_n], rel);
    }

    SubtractionSolid range_pipe(range_tube, range_vacuum);

    Volume v_combined_pipe("v_pipe_combined_" + combineNames[range_n], range_pipe, m_Al);
    Volume v_combined_vacuum("v_vacuum_combined_" + combineNames[range_n], range_vacuum,
                             m_Vacuum);
    v_combined_pipe.setVisAttributes(description.visAttributes(vis_name));

    assembly.placeVolume(v_combined_pipe, base_transform);
    auto placed_vacuum = assembly.placeVolume(v_combined_vacuum, base_transform);

    DetElement vacuum_element(sdet, combineNames[range_n] + "_vacuum", combineStartIds[range_n]);
    vacuum_element.setPlacement(placed_vacuum);
  }

  // Final placement
  auto pv_assembly =
      description.pickMotherVolume(sdet).placeVolume(assembly, Position(0.0, 0.0, 0.0));

  sdet.setPlacement(pv_assembly);

  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(BeamPipeChain, create_detector)
