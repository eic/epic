// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Simon Gardner

//==========================================================================
//
// Places a small sensitive disk of vacuum at the end of beam pipes
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

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

  // Grab info for beamline magnets
  for (xml_coll_t slice_coll(x_det, _Unicode(slice)); slice_coll; slice_coll++) { // pipes

    string grandmotherName = slice_coll.attr<string>(_Unicode(grandmother));
    string motherName      = slice_coll.attr<string>(_Unicode(mother));
    bool detStart          = getAttrOrDefault<bool>(slice_coll, _Unicode(end), true);
    int pipe_id            = getAttrOrDefault<int>(slice_coll, _Unicode(pipe_id), 0);
    string slice_name      = slice_coll.attr<string>(_Unicode(name));
    DetElement mother      = description.detector(grandmotherName).child(motherName);

    // Get the mother volume
    Volume mother_vol = mother.volume();

    double sensitive_thickness = 0.1 * mm;
    Solid sensitive_solid;
    Position disk_position;

    // Check if we should use plane-based cross section (default: false, use cone)
    bool use_plane_xs = getAttrOrDefault<bool>(slice_coll, _Unicode(use_cross_section), false);

    if (use_plane_xs) {
      // Use plane-based cross section
      
      // Get plane definition from XML (normal vector and distance from origin)
      double plane_nx = getAttrOrDefault<double>(slice_coll, _Unicode(plane_nx), 0.0);
      double plane_ny = getAttrOrDefault<double>(slice_coll, _Unicode(plane_ny), 0.0);
      double plane_nz = getAttrOrDefault<double>(slice_coll, _Unicode(plane_nz), 1.0);
      double plane_d  = getAttrOrDefault<double>(slice_coll, _Unicode(plane_d), 0.0);
      
      // Normalize the plane normal vector
      double norm = sqrt(plane_nx*plane_nx + plane_ny*plane_ny + plane_nz*plane_nz);
      plane_nx /= norm;
      plane_ny /= norm;
      plane_nz /= norm;
      
      // Get maximum dimension for cutting box
      double max_dim = getAttrOrDefault<double>(slice_coll, _Unicode(max_dimension), 1.0 * m);
      
      // Create a thin box perpendicular to the plane normal
      Box cutting_box(max_dim, max_dim, sensitive_thickness / 2);
      
      // Calculate rotation to align box with plane normal
      // Default box normal is along Z, rotate to align with plane normal
      RotationZYX rotation;
      if (fabs(plane_nz - 1.0) > 1e-6 || fabs(plane_nx) > 1e-6 || fabs(plane_ny) > 1e-6) {
        // Calculate rotation angles to align Z-axis with plane normal
        double theta = acos(plane_nz); // angle from Z axis
        double phi = atan2(plane_ny, plane_nx); // angle in XY plane
        rotation = RotationZYX(0.0, theta, phi);
      }
      
      // Position the cutting plane at distance plane_d along the normal
      Position plane_pos(plane_nx * plane_d, plane_ny * plane_d, plane_nz * plane_d);
      
      // Create intersection of mother volume with the cutting box
      sensitive_solid = IntersectionSolid(mother_vol.solid(), cutting_box, 
                                          Transform3D(rotation, plane_pos));
      disk_position = Position(0.0, 0.0, 0.0);

    } else {
      // Use cone segment (default behavior)
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
      disk_position = Position(0.0, 0.0, zPos);
    }

    Volume v_start_disk("v_start_disk_" + motherName, sensitive_solid, m_Vacuum);
    v_start_disk.setSensitiveDetector(sens);

    auto disk_placement = mother_vol.placeVolume(v_start_disk, disk_position);
    disk_placement.addPhysVolID("end", detStart);
    disk_placement.addPhysVolID("pipe", pipe_id);
    disk_placement.addPhysVolID("system", det_id);

    DetElement slice_element(sdet, slice_name, pipe_id);

    slice_element.setPlacement(disk_placement);
    description.declareParent(slice_name, mother);
  }

  auto pv_assembly = description.worldVolume().placeVolume(assembly, Position(0.0, 0.0, 0.0));
  pv_assembly.addPhysVolID("system", det_id);
  sdet.setPlacement(pv_assembly);

  return sdet;
}

DECLARE_DETELEMENT(BeamPipeTracking, create_detector)
