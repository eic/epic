// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Dhevan Gangadharan, Simon Gardner

//==========================================================================
//
// BEAMPIPECHAIN DETECTOR GEOMETRY FACTORY
//
// This file constructs a chain of beam pipe segments connecting beamline magnets.
// The geometry is built from XML configuration with automatic bend joins, optional
// boolean combination groups, and optional split-branch extensions.
//
// BASIC PIPE CHAIN:
//   - Reads series of conical pipe segments (names, ids, positions, angles, radii)
//   - Creates torus transition pieces between consecutive pipes with different angles
//   - Automatically adjusts pipe lengths to account for torus volumes at junctions
//   - Each pipe is modeled as Aluminum shell with inner Vacuum
//   - Vacuums are tracked separately to enable boolean subtraction
//
// COMBINE GROUPS (boolean ranges):
//   - Optional ranges to resolve overlaps between individually placed pipe components
//   - Specified in XML as <combine start_id="A" end_id="B"/> (inclusive range by id)
//   - Split components precede main-chain components in overlap ownership
//   - Earlier components own overlapping outer envelopes and vacuum regions
//   - Each wall subtracts earlier outer envelopes and later vacuum primitives
//   - Each vacuum subtracts earlier vacuum primitives; no unions are constructed
//   - Conservative bounding-box rejection omits only provably disjoint cuts
//   - Combines also include adjacent bend pieces if BOTH neighbors are in same range
//   - Pipes/bends outside any combine group are placed individually
//
// SPLIT BRANCHES (junction extensions):
//   - Allows straight-pipe branches to spawn from downstream end of a chosen source pipe
//   - Specified in XML as: <split from_id="N" id="M" length="L" theta="T"
//                                   rout1="R1" rout2="R2" name="Name"/>
//   - Split branches created as simple conical segments extending from source pipe endpoint
//   - Uses source outer radii if rout1/rout2 not specified in XML
//   - Split pieces inherit combine-group membership from their id
//   - Split pieces participate in vacuum subtraction within combine groups
//
// VACUUM SUBTRACTION:
//   - All pipe solids (straight, bend, split) have associated vacuum solids
//   - Individual pipes: vacuum is placed separately alongside tube
//   - Combine groups: individual wall/vacuum components partition the same material
//     regions as the former union construction, retaining the terminal end cuts
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include "TGeoBBox.h"
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
  vector<int> splitFromIds;
  vector<int> splitIds;
  vector<string> splitNames;
  vector<double> splitLengths;
  vector<double> splitThetas;
  vector<double> splitROuters1;
  vector<double> splitROuters2;
  vector<double> splitJunctionOverlaps;

  // Optional boolean-combine groups, configured via XML:
  //   <combine name="..." start_id="..." end_id="..."/>
  // Pieces in a range are placed separately with overlap ownership resolved by subtraction.
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
    std::string name = getAttrOrDefault<std::string>(combine, _Unicode(name),
                                                     "combined_" + std::to_string(start_id) + "_" +
                                                         std::to_string(end_id));

    combineStartIds.push_back(start_id);
    combineEndIds.push_back(end_id);
    combineNames.push_back(name);
  }

  // Optional branch segments starting at the downstream end of an existing pipe id:
  //   <split from_id="..." id="..." name="..." length="..." theta="..." rout1="..." rout2="..." junction_overlap="..."/>
  for (xml_coll_t split_coll(x_det, _Unicode(split)); split_coll; ++split_coll) {
    xml_comp_t split(split_coll);

    int from_id   = getAttrOrDefault<int>(split, _Unicode(from_id), -1);
    int id        = getAttrOrDefault<int>(split, _Unicode(id), -1);
    string name   = getAttrOrDefault<string>(split, _Unicode(name), "split_" + std::to_string(id));
    double length = getAttrOrDefault<double>(split, _Unicode(length), 0.0);
    double theta  = getAttrOrDefault<double>(split, _Unicode(theta), 0.0);
    double rout1  = getAttrOrDefault<double>(split, _Unicode(rout1), 0.0);
    double rout2  = getAttrOrDefault<double>(split, _Unicode(rout2), 0.0);
    double junction_overlap = getAttrOrDefault<double>(split, _Unicode(junction_overlap), 0.0);

    splitFromIds.push_back(from_id);
    splitIds.push_back(id);
    splitNames.push_back(name);
    splitLengths.push_back(length);
    splitThetas.push_back(theta);
    splitROuters1.push_back(rout1);
    splitROuters2.push_back(rout2);
    splitJunctionOverlaps.push_back(junction_overlap);
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
  vector<Solid> straightOuters;
  vector<Solid> straightVacuums;
  vector<Transform3D> straightTransforms;
  vector<double> straightXCenters;
  vector<double> straightZCenters;
  vector<double> straightLengthsPlaced;
  vector<double> straightThetasPlaced;
  vector<double> straightOuterRadiiMax;
  vector<size_t> clipPieceIndices;
  vector<int> clipEndSigns;

  // Parallel vectors for bend-join torus pieces between neighboring straight pipes.
  vector<std::string> bendNames;
  vector<int> bendIds;
  vector<int> bendRangeIndices;
  vector<Solid> bendTubes;
  vector<Solid> bendOuters;
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
  // needs overlap subtraction within a combine group.
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
    ConeSegment s_outer(length / 2.0, 0, rOuter2, 0, rOuter1);

    Transform3D straightTransform(RotationY(theta), Position(xCenter, 0, zCenter));
    straightNames.push_back(names[pipeN]);
    straightIds.push_back(ids[pipeN]);
    straightRangeIndices.push_back(find_range_index(ids[pipeN]));
    straightTubes.push_back(Solid(s_tube));
    straightOuters.push_back(Solid(s_outer));
    straightVacuums.push_back(Solid(s_vacuum));
    straightTransforms.push_back(straightTransform);
    straightXCenters.push_back(xCenter);
    straightZCenters.push_back(zCenter);
    straightLengthsPlaced.push_back(length);
    straightThetasPlaced.push_back(theta);
    straightOuterRadiiMax.push_back(std::max(rOuter1, rOuter2));

    // // Add joining bend to next pipe
    if (pipeN != xCenters.size() - 1 && bendLengths[pipeN] != 0) {

      // double bendLength = bendLengths[pipeN];
      double bendAngle = thetas[pipeN + 1] - thetas[pipeN];

      // Create a torus to join the two pipes
      Torus s_bend_solid(rOuter2 + bendRadius, rOuter2 - thickness, rOuter2, 0.0, abs(bendAngle));

      //Create a vacuum torus to place inside the solid segment
      Torus s_bend_vac(rOuter2 + bendRadius, 0, rOuter2 - thickness, 0, abs(bendAngle));
      Torus s_bend_outer(rOuter2 + bendRadius, 0, rOuter2, 0, abs(bendAngle));

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
      bendOuters.push_back(Solid(s_bend_outer));
      bendVacuums.push_back(Solid(s_bend_vac));
      bendTransforms.push_back(bendTransform);
    }
  }

  // Build optional split branches from downstream ends of configured source ids.
  // Split branches are created as simple conical segments extending from the source pipe.
  for (size_t split_n = 0; split_n < splitIds.size(); ++split_n) {
    int source_index = -1;
    for (size_t piece_n = 0; piece_n < straightIds.size(); ++piece_n) {
      if (straightIds[piece_n] == splitFromIds[split_n]) {
        source_index = static_cast<int>(piece_n);
        break;
      }
    }

    if (source_index < 0) {
      printout(WARNING, "BeamPipeChain", "Split '%s' references unknown from_id=%d. Skipping.",
               splitNames[split_n].c_str(), splitFromIds[split_n]);
      continue;
    }

    double branch_length = splitLengths[split_n];
    if (branch_length <= 0) {
      printout(WARNING, "BeamPipeChain", "Split '%s' has non-positive length. Skipping.",
               splitNames[split_n].c_str());
      continue;
    }

    double source_theta    = straightThetasPlaced[source_index];
    double source_x_center = straightXCenters[source_index];
    double source_z_center = straightZCenters[source_index];
    double source_length   = straightLengthsPlaced[source_index];
    double branch_theta    = splitThetas[split_n];
    double branch_r1 = splitROuters1[split_n] > 0 ? splitROuters1[split_n] : rOuters1[source_index];
    double branch_r2 = splitROuters2[split_n] > 0 ? splitROuters2[split_n] : branch_r1;

    // Calculate position at downstream end of source pipe
    double source_end_x = source_x_center - 0.5 * source_length * sin(source_theta);
    double source_end_z = source_z_center - 0.5 * source_length * cos(source_theta);

    // Configurable junction overlap for connecting to source or external geometry.
    // Positive values create overlap (for internal connections to avoid gaps).
    // Zero or negative values create separation (for external geometry connections).
    double junction_overlap = splitJunctionOverlaps[split_n];
    if (junction_overlap > 0 && junction_overlap > 0.45 * branch_length) {
      junction_overlap = 0.45 * branch_length;
    }
    double branch_length_effective = branch_length + junction_overlap;
    double source_attach_x         = source_end_x + junction_overlap * sin(branch_theta);
    double source_attach_z         = source_end_z + junction_overlap * cos(branch_theta);

    // Position split branch center: halfway along its length from source endpoint
    double branch_center_x = source_attach_x - 0.5 * branch_length_effective * sin(branch_theta);
    double branch_center_z = source_attach_z - 0.5 * branch_length_effective * cos(branch_theta);

    // Create split tube and vacuum solids
    ConeSegment s_split_tube(branch_length_effective / 2.0, branch_r2 - thickness, branch_r2,
                             branch_r1 - thickness, branch_r1);
    ConeSegment s_split_vac(branch_length_effective / 2.0, 0, branch_r2 - thickness, 0,
                            branch_r1 - thickness);
    ConeSegment s_split_outer(branch_length_effective / 2.0, 0, branch_r2, 0, branch_r1);
    Transform3D split_transform(RotationY(branch_theta),
                                Position(branch_center_x, 0, branch_center_z));

    // Add split branch to straight piece vectors
    straightNames.push_back(splitNames[split_n]);
    straightIds.push_back(splitIds[split_n]);
    straightRangeIndices.push_back(find_range_index(splitIds[split_n]));
    straightTubes.push_back(Solid(s_split_tube));
    straightOuters.push_back(Solid(s_split_outer));
    straightVacuums.push_back(Solid(s_split_vac));
    straightTransforms.push_back(split_transform);
    straightXCenters.push_back(branch_center_x);
    straightZCenters.push_back(branch_center_z);
    straightLengthsPlaced.push_back(branch_length_effective);
    straightThetasPlaced.push_back(branch_theta);
    straightOuterRadiiMax.push_back(std::max(branch_r1, branch_r2));

    size_t branch_index = straightNames.size() - 1;
    int branch_range    = straightRangeIndices[branch_index];
    if (branch_range >= 0) {
      clipPieceIndices.push_back(branch_index);
      clipEndSigns.push_back(-1);
    }
  }

  // Clip the terminal downstream end of the original straight chain for each combine
  // range. This is the last non-split piece in that range.
  size_t base_straight_count = names.size();
  for (size_t range_n = 0; range_n < combineStartIds.size(); ++range_n) {
    int last_piece = -1;
    for (size_t piece_n = 0; piece_n < base_straight_count; ++piece_n) {
      if (straightRangeIndices[piece_n] == static_cast<int>(range_n)) {
        last_piece = static_cast<int>(piece_n);
      }
    }
    if (last_piece >= 0) {
      clipPieceIndices.push_back(static_cast<size_t>(last_piece));
      clipEndSigns.push_back(-1);
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

  // Partition each range without unions. Cutters always remain original primitives.
  // Split vacuums own junction overlaps, keeping the luminosity bore shallow.
  struct RangePiece {
    string name;
    int id;
    Solid tube;
    Solid outer;
    Solid vacuum;
    Transform3D transform;
  };
  for (size_t range_n = 0; range_n < combineStartIds.size(); ++range_n) {
    vector<RangePiece> pieces;
    auto append_straight = [&](size_t begin, size_t end) {
      for (size_t piece_n = begin; piece_n < end; ++piece_n) {
        if (straightRangeIndices[piece_n] == static_cast<int>(range_n)) {
          pieces.push_back({straightNames[piece_n], straightIds[piece_n], straightTubes[piece_n],
                            straightOuters[piece_n], straightVacuums[piece_n],
                            straightTransforms[piece_n]});
        }
      }
    };
    append_straight(base_straight_count, straightNames.size());
    append_straight(0, base_straight_count);
    for (size_t piece_n = 0; piece_n < bendNames.size(); ++piece_n) {
      if (bendRangeIndices[piece_n] == static_cast<int>(range_n)) {
        pieces.push_back({bendNames[piece_n], bendIds[piece_n], bendTubes[piece_n],
                          bendOuters[piece_n], bendVacuums[piece_n], bendTransforms[piece_n]});
      }
    }
    if (pieces.empty()) {
      printout(WARNING, "BeamPipeChain", "Combine range '%s' (%d-%d) did not match any segments.",
               combineNames[range_n].c_str(), combineStartIds[range_n], combineEndIds[range_n]);
      continue;
    }

    const Transform3D& base_transform = pieces.front().transform;
    Transform3D inv_base              = base_transform.Inverse();
    Assembly range_pipe("v_pipe_combined_" + combineNames[range_n]);
    Assembly range_vacuum("v_vacuum_combined_" + combineNames[range_n]);
    DetElement vacuum_element(sdet, combineNames[range_n] + "_vacuum", combineStartIds[range_n]);

    vector<Solid> cuts;
    vector<Transform3D> cut_transforms;
    vector<string> cut_names;
    for (size_t cut_n = 0; cut_n < clipPieceIndices.size(); ++cut_n) {
      size_t piece_n = clipPieceIndices[cut_n];
      if (straightRangeIndices[piece_n] != static_cast<int>(range_n)) {
        continue;
      }

      double theta  = straightThetasPlaced[piece_n];
      double halfL  = 0.5 * straightLengthsPlaced[piece_n];
      double rOuter = straightOuterRadiiMax[piece_n];
      double sign   = static_cast<double>(clipEndSigns[cut_n]);
      double axis_x = sign * sin(theta);
      double axis_z = sign * cos(theta);

      Position end_pos(straightXCenters[piece_n] + sign * halfL * sin(theta), 0.0,
                       straightZCenters[piece_n] + sign * halfL * cos(theta));
      Position cut_center(end_pos.x() + halfL * axis_x, 0.0, end_pos.z() + halfL * axis_z);

      Tube cut_tube(0.0, rOuter, halfL);
      cuts.push_back(Solid(cut_tube));
      cut_transforms.emplace_back(RotationY(theta), cut_center);
      cut_names.push_back(straightNames[piece_n]);
    }

    for (size_t piece_n = 0; piece_n < pieces.size(); ++piece_n) {
      const auto& piece     = pieces[piece_n];
      Transform3D inv_piece = piece.transform.Inverse();
      Solid tube            = piece.tube;
      Solid vacuum          = piece.vacuum;
      auto subtract = [](Solid left, Solid right, const Transform3D& transform) -> Solid {
        auto left_box = dynamic_cast<TGeoBBox*>(left.ptr());
        auto right_box = dynamic_cast<TGeoBBox*>(right.ptr());
        if (!left_box || !right_box) {
          except("BeamPipeChain", "Cannot bound a subtraction operand.");
        }
        double low[3] = {TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};
        double high[3] = {-TGeoShape::Big(), -TGeoShape::Big(), -TGeoShape::Big()};
        for (int ix : {-1, 1}) {
          for (int iy : {-1, 1}) {
            for (int iz : {-1, 1}) {
              XYZPoint corner(right_box->GetOrigin()[0] + ix * right_box->GetDX(),
                              right_box->GetOrigin()[1] + iy * right_box->GetDY(),
                              right_box->GetOrigin()[2] + iz * right_box->GetDZ());
              XYZPoint point = transform * corner;
              const double coordinates[3] = {point.x(), point.y(), point.z()};
              for (int axis = 0; axis < 3; ++axis) {
                low[axis] = std::min(low[axis], coordinates[axis]);
                high[axis] = std::max(high[axis], coordinates[axis]);
              }
            }
          }
        }
        const double half_lengths[3] = {left_box->GetDX(), left_box->GetDY(), left_box->GetDZ()};
        // Retain touching or numerically ambiguous cuts.
        const double margin = 1.e-7 * mm;
        for (int axis = 0; axis < 3; ++axis) {
          if (high[axis] < left_box->GetOrigin()[axis] - half_lengths[axis] - margin ||
              low[axis] > left_box->GetOrigin()[axis] + half_lengths[axis] + margin) {
            return left;
          }
        }
        return SubtractionSolid(left, right, transform);
      };
      for (size_t other_n = 0; other_n < pieces.size(); ++other_n) {
        if (other_n == piece_n) {
          continue;
        }
        const auto& other = pieces[other_n];
        Transform3D rel   = inv_piece * other.transform;
        if (other_n < piece_n) {
          tube   = subtract(tube, other.outer, rel);
          vacuum = subtract(vacuum, other.vacuum, rel);
        } else {
          tube = subtract(tube, other.vacuum, rel);
        }
      }
      // Distribute the former group-wide terminal cuts over every component.
      for (size_t cut_n = 0; cut_n < cuts.size(); ++cut_n) {
        // A finite primitive already ends at its own terminal clipping plane.
        if (piece.name == cut_names[cut_n]) {
          continue;
        }
        Transform3D rel_cut = inv_piece * cut_transforms[cut_n];
        tube                = subtract(tube, cuts[cut_n], rel_cut);
        vacuum              = subtract(vacuum, cuts[cut_n], rel_cut);
      }
      Volume v_tube("v_tube_" + piece.name, tube, m_Al);
      Volume v_vacuum("v_vacuum_" + piece.name, vacuum, m_Vacuum);
      v_tube.setVisAttributes(description.visAttributes(vis_name));
      Transform3D rel = inv_base * piece.transform;
      range_pipe.placeVolume(v_tube, rel);
      auto placed_vacuum = range_vacuum.placeVolume(v_vacuum, rel);
      DetElement piece_element(vacuum_element, piece.name + "_vacuum", piece.id);
      piece_element.setPlacement(placed_vacuum);
    }

    assembly.placeVolume(range_pipe, base_transform);
    auto placed_vacuum = assembly.placeVolume(range_vacuum, base_transform);

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
