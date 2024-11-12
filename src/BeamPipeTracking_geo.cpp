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

    // Get mother volume shape as cone segment
    ConeSegment mother_shape = mother_vol.solid();

    // Get the parameters of the mother volume
    double rOuter1 = mother_shape.rMax1();
    double rOuter2 = mother_shape.rMax2();
    double length  = 2 * mother_shape.dZ();

    double sensitive_thickness = 0.1 * mm;

    //Calculate R or cone after sensitive layer

    double rEnd = rOuter2 - (rOuter2 - rOuter1) * sensitive_thickness / length;
    double zPos = length / 2.0 - sensitive_thickness / 2.0;
    if (detStart) {
      rEnd = rOuter1 - (rOuter1 - rOuter2) * sensitive_thickness / length;
      zPos = -length / 2.0 + sensitive_thickness / 2.0;
    }

    ConeSegment s_start_disk(sensitive_thickness / 2, 0.0, rOuter2, 0.0, rEnd);
    Volume v_start_disk("v_start_disk_" + motherName, s_start_disk, m_Vacuum);
    v_start_disk.setSensitiveDetector(sens);

    auto disk_placement = mother_vol.placeVolume(v_start_disk, Position(0.0, 0.0, zPos));
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
