// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Leszek Kosarzewski

//==========================================================================
//  Implementation of backward endcap flux return
//--------------------------------------------------------------------------
//  Author: Leszek Kosarzewski (OSU)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"
#include "XML/Layering.h"
/*using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
*/

static dd4hep::Ref_t create_detector(dd4hep::Detector& description, xml_h e, [[maybe_unused]] dd4hep::SensitiveDetector sens)
{
  xml_det_t x_det    = e;
  int       det_id   = x_det.id();
  std::string    det_name = x_det.nameStr();
  dd4hep::Material  air      = description.air();
  xml_comp_t x_pos  = x_det.position();


  dd4hep::Assembly     assembly(det_name);
  dd4hep::DetElement   sdet(det_name, det_id);
  dd4hep::PlacedVolume pv;

  double disksGap = 0.0;

  // Looping through all the different layer sections
  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    xml_comp_t x_layer         = xc;

    int     layer_id = x_layer.id();
    double     layer_rmax = x_layer.rmax();
    double     layer_rmin = x_layer.rmin();
    double     layer_thickness = x_layer.thickness();
    double     layer_zpos = x_layer.zpos();

    dd4hep::Material   l_mat   = description.material(x_layer.materialStr());


    dd4hep::DetElement disk_ele("disk_ele", layer_id);
    dd4hep::Volume disk(x_layer.nameStr(), dd4hep::Tube(layer_rmin, layer_rmax, layer_thickness / 2, 0.0, 2.0 * M_PI), air);
    disk.setVisAttributes(description.visAttributes(x_layer.visStr()));

    dd4hep::Volume halfdisk("halfdisk", dd4hep::Tube(layer_rmin, layer_rmax, layer_thickness / 2, M_PI / 2, M_PI * 3 / 2), l_mat);
    halfdisk.setVisAttributes(description.visAttributes(x_layer.visStr()));
    dd4hep::PlacedVolume s_phv1 = disk.placeVolume(halfdisk, dd4hep::Position(-disksGap / 2, 0, 0));
    s_phv1.addPhysVolID("halfdisk", 0);
    dd4hep::PlacedVolume s_phv2 = disk.placeVolume(halfdisk, dd4hep::Transform3D(dd4hep::RotationZYX(M_PI, 0, 0), dd4hep::Position(+disksGap / 2, 0, 0)));
    s_phv2.addPhysVolID("halfdisk", 1);


    pv = assembly.placeVolume(disk, dd4hep::Position(0, 0, -layer_thickness/2-layer_zpos));
    pv.addPhysVolID("layer", layer_id);
    disk_ele.setPlacement(pv);
  }

  // Get position and place volume
  dd4hep::Position       pos(x_pos.x(), x_pos.y(), x_pos.z());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, pos);
  pv.addPhysVolID("system", det_id).addPhysVolID("barrel", 0);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_EndcapFluxReturnN, create_detector)
