// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2025 Simon Gardner

//==========================================================================
//
// Places thin sensitive slices of vacuum along a beam pipe
// Component bits distinguish fragments; system/pipe/end identify the full plane.
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>
#include <algorithm>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {

  using namespace ROOT::Math;
  xml_det_t x_det   = e;
  string det_name   = x_det.nameStr();
  int det_id        = x_det.id();
  Material m_Vacuum = description.material("Vacuum");
  string vis_name   = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");

  sens.setType("tracker");

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name + "_assembly");

  // Loop over each requested slice from the geometry description
  for (xml_coll_t slice_coll(x_det, _Unicode(slice)); slice_coll; slice_coll++) { // pipes

    string grandmotherName = slice_coll.attr<string>(_Unicode(grandmother));
    string motherName      = slice_coll.attr<string>(_Unicode(mother));
    bool detStart          = getAttrOrDefault<bool>(slice_coll, _Unicode(end), true);
    int pipe_id            = getAttrOrDefault<int>(slice_coll, _Unicode(pipe_id), 0);
    string slice_name      = slice_coll.attr<string>(_Unicode(name));
    DetElement mother      = description.detector(grandmotherName).child(motherName);
    vector<DetElement> slice_mothers;
    if (mother.volume().isAssembly()) {
      for (const auto& child : mother.children()) {
        slice_mothers.push_back(child.second);
      }
      if (slice_mothers.empty()) {
        except("BeamPipeTracking", "Vacuum assembly '%s' has no component DetElements.",
               motherName.c_str());
      }
    } else {
      slice_mothers.push_back(mother);
    }

    if (slice_mothers.size() > 256) {
      except("BeamPipeTracking", "Slice '%s' exceeds the 8-bit component field capacity.",
             slice_name.c_str());
    }
    for (size_t component_id = 0; component_id < slice_mothers.size(); ++component_id) {
      const auto& slice_mother = slice_mothers[component_id];
      // Get the mother volume
      Volume mother_vol     = slice_mother.volume();
      string component_name = slice_mothers.size() > 1 ? slice_mother.name() : motherName;

      double sensitive_thickness = 0.1 * mm;
      Solid sensitive_solid;
      Transform3D disk_transform; // Full placement transform

      // Check if we should use plane-based cross section (default: false, use cone)
      bool use_plane_xs = getAttrOrDefault<bool>(slice_coll, _Unicode(use_cross_section), false);
      if (mother.volume().isAssembly() && !use_plane_xs) {
        except("BeamPipeTracking", "Vacuum assembly '%s' requires use_cross_section=true.",
               motherName.c_str());
      }

      if (use_plane_xs) {

        // Get plane definition from XML (in world frame)
        double plane_x    = getAttrOrDefault<double>(slice_coll, _Unicode(plane_x), 0.0);
        double plane_y    = getAttrOrDefault<double>(slice_coll, _Unicode(plane_y), 0.0);
        double plane_z    = getAttrOrDefault<double>(slice_coll, _Unicode(plane_z), 0.0);
        double plane_yrot = getAttrOrDefault<double>(slice_coll, _Unicode(plane_yrot), 0.0);

        // Get the mother volume's transformation to world frame
        auto mother_matrix          = slice_mother.nominal().worldTransformation();
        const Double_t* translation = mother_matrix.GetTranslation();
        const Double_t* rotation    = mother_matrix.GetRotationMatrix();

        // Build Transform3D from TGeoMatrix
        Transform3D mother_to_world(rotation[0], rotation[1], rotation[2], translation[0],
                                    rotation[3], rotation[4], rotation[5], translation[1],
                                    rotation[6], rotation[7], rotation[8], translation[2]);

        // Transform the plane position and rotation from world frame to mother's local frame
        Position plane_pos_world(plane_x, plane_y, plane_z);
        RotationY plane_rot_world(plane_yrot);
        Transform3D plane_transform_world(plane_rot_world, plane_pos_world);

        // Apply inverse transformation to get plane in mother's local coordinates
        disk_transform = mother_to_world.Inverse() * plane_transform_world;

        // Get maximum dimension for cutting box
        double max_dim = getAttrOrDefault<double>(slice_coll, _Unicode(max_dimension), 1.0 * m);

        // Create a thin box perpendicular to the plane normal
        Box cutting_box(max_dim, max_dim, sensitive_thickness / 2);

        // Assembly components far from this plane need no sensitive daughter.
        if (mother.volume().isAssembly()) {
          auto bounds = dynamic_cast<TGeoBBox*>(mother_vol.solid().ptr());
          if (!bounds) {
            except("BeamPipeTracking", "Vacuum component '%s' has no bounding box.",
                   component_name.c_str());
          }
          Transform3D to_disk = disk_transform.Inverse();
          double min_x = max_dim, max_x = -max_dim;
          double min_y = max_dim, max_y = -max_dim;
          double min_z = sensitive_thickness / 2, max_z = -sensitive_thickness / 2;
          for (int ix : {-1, 1}) {
            for (int iy : {-1, 1}) {
              for (int iz : {-1, 1}) {
                XYZPoint corner(bounds->GetOrigin()[0] + ix * bounds->GetDX(),
                                bounds->GetOrigin()[1] + iy * bounds->GetDY(),
                                bounds->GetOrigin()[2] + iz * bounds->GetDZ());
                XYZPoint point = to_disk * corner;
                min_x          = std::min(min_x, point.x());
                max_x          = std::max(max_x, point.x());
                min_y          = std::min(min_y, point.y());
                max_y          = std::max(max_y, point.y());
                min_z          = std::min(min_z, point.z());
                max_z          = std::max(max_z, point.z());
              }
            }
          }
          if (min_x >= max_dim || max_x <= -max_dim || min_y >= max_dim || max_y <= -max_dim ||
              min_z >= sensitive_thickness / 2 || max_z <= -sensitive_thickness / 2) {
            continue;
          }
        }

        // Create intersection of cutting box with mother volume (reversed order to preserve rotation)
        // Apply the inverse transform to the second operand (mother) to express it in the box frame
        sensitive_solid =
            IntersectionSolid(cutting_box, mother_vol.solid(), disk_transform.Inverse());

      } else {
        // Use cone segment (default behavior before outline of xs is implemented in benchmark)
        ConeSegment mother_shape = mother_vol.solid();

        // Get the parameters of the mother volume
        double rOuter1 = mother_shape.rMax1();
        double rOuter2 = mother_shape.rMax2();
        double length  = 2 * mother_shape.dZ();

        //Calculate R of cone after sensitive layer
        double rEnd = rOuter2 - (rOuter2 - rOuter1) * sensitive_thickness / length;
        double zPos = length / 2.0 - sensitive_thickness / 2.0;
        if (detStart) {
          rEnd = rOuter1 - (rOuter1 - rOuter2) * sensitive_thickness / length;
          zPos = -length / 2.0 + sensitive_thickness / 2.0;
        }

        sensitive_solid = ConeSegment(sensitive_thickness / 2, 0.0, rOuter2, 0.0, rEnd);
        disk_transform  = Transform3D(Rotation3D(), Position(0.0, 0.0, zPos));
      }

      string disk_name =
          mother.volume().isAssembly() ? component_name + "_" + slice_name : motherName;
      Volume v_start_disk("v_start_disk_" + disk_name, sensitive_solid, m_Vacuum);
      v_start_disk.setSensitiveDetector(sens);

      auto disk_placement = mother_vol.placeVolume(v_start_disk, disk_transform);
      disk_placement.addPhysVolID("end", detStart);
      disk_placement.addPhysVolID("pipe", pipe_id);
      disk_placement.addPhysVolID("system", det_id);
      disk_placement.addPhysVolID("component", component_id);

      DetElement slice_element(slice_mother, slice_name, pipe_id);

      slice_element.setPlacement(disk_placement);
    }
  }

  auto pv_assembly = description.worldVolume().placeVolume(assembly);
  pv_assembly.addPhysVolID("system", det_id);
  sdet.setPlacement(pv_assembly);

  return sdet;
}

DECLARE_DETELEMENT(BeamPipeTracking, create_detector)
