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

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /* sens */) {

  using namespace ROOT::Math;
  xml_det_t x_det   = e;
  string det_name   = x_det.nameStr();
  int det_id        = x_det.id();
  Material m_Vacuum = description.material("Vacuum");
  string vis_name   = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");

  string grandmotherName = x_det.attr<string>(_Unicode(grandmother));
  string motherName      = x_det.attr<string>(_Unicode(mother));
  bool detStart          = getAttrOrDefault<bool>(x_det, _Unicode(end), true);
  DetElement mother      = description.detector(grandmotherName).child(motherName);

  DetElement sdet(det_name, det_id);

  // Get the mother volume
  Volume mother_vol = mother.volume();

  // Get mother volume shape as cone segment
  ConeSegment mother_shape = mother_vol.solid();

  // Get the parameters of the mother volume
  double rOuter1 = mother_shape.rMax1();
  double rOuter2 = mother_shape.rMax2();
  double length  = 2 * mother_shape.dZ();

  double sensitive_thickness = 100 * mm;

  //Calculate R or cone after sensitive layer
  double rEnd = rOuter2 - (rOuter2 - rOuter1) * sensitive_thickness / length;
  double zPos = length / 2.0 - sensitive_thickness / 2.0;
  if (detStart) {
    rEnd = rOuter1 - (rOuter1 - rOuter2) * sensitive_thickness / length;
    zPos = -length / 2.0 + sensitive_thickness / 2.0;
  }

  ConeSegment s_start_disk(sensitive_thickness / 2, 0.0, rOuter2, 0.0, rEnd);
  Volume v_start_disk("stop_disk_" + motherName, s_start_disk, m_Vacuum);

  v_start_disk.setLimitSet(description, "kill_limits");

  auto disk_placement = mother_vol.placeVolume(v_start_disk, Position(0.0, 0.0, zPos));

  sdet.setPlacement(disk_placement);
  description.declareParent(det_name, mother);

  return sdet;
}

DECLARE_DETELEMENT(BeamPipeStop, create_detector)
