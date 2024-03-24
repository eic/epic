// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Dmitry Romanov

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////
// Central Barrel DIRC
//////////////////////////////////

using namespace std;
using namespace dd4hep;

static dd4hep::Trap MakeTrap(const std::string& pName, double pZ, double pY, double pX, double pLTX);

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t xml_det = e;

  // Detector element
  string     det_name = xml_det.nameStr();
  int        det_id   = xml_det.id();
  DetElement det(det_name, det_id);

  // Detector dimension, position, rotation
  xml_dim_t dirc_dim = xml_det.dimensions();
  xml_dim_t dirc_pos = xml_det.position();
  double    det_rmin = dirc_dim.rmin();
  double    det_rmax = dirc_dim.rmax();
  double    det_ravg = (det_rmin + det_rmax) / 2;

  // Detector type
  sens.setType("tracker");

  // Entire DIRC assembly
  Assembly det_volume("DIRC");
  det_volume.setVisAttributes(desc.visAttributes(xml_det.visStr()));
  Transform3D det_tr(RotationY(0), Position(0.0, 0.0, dirc_pos.z()));
  det.setPlacement(desc.pickMotherVolume(det).placeVolume(det_volume, det_tr).addPhysVolID("system", det_id));

  // Construct module
  xml_comp_t xml_module = xml_det.child(_U(module));

  Assembly dirc_module("DIRCModule");
  dirc_module.setVisAttributes(desc.visAttributes(xml_module.visStr()));

  // Bar
  xml_comp_t xml_bar    = xml_module.child(_Unicode(bar));
  double     bar_height = xml_bar.height();
  double     bar_width  = xml_bar.width();
  double     bar_length = xml_bar.length();
  Box        bar_box("bar_box", bar_height / 2, bar_width / 2, bar_length / 2);
  Volume     bar_vol("bar_vol", bar_box, desc.material(xml_bar.materialStr()));
  bar_vol.setVisAttributes(desc.visAttributes(xml_bar.visStr()));

  // Glue
  xml_comp_t xml_glue       = xml_module.child(_Unicode(glue));
  double     glue_thickness = xml_glue.thickness();
  Box        glue_box("glue_box", bar_height / 2, bar_width / 2, glue_thickness / 2);
  Volume     glue_vol("glue_vol", glue_box, desc.material(xml_glue.materialStr()));
  glue_vol.setVisAttributes(desc.visAttributes(xml_glue.visStr()));

  auto bar_repeat_y    = xml_bar.attr<int>(_Unicode(repeat_y));
  auto bar_repeat_z    = xml_bar.attr<int>(_Unicode(repeat_z));
  auto bar_gap         = xml_bar.gap();
  auto bar_assm_width  = (bar_width + bar_gap) * bar_repeat_y - bar_gap;
  auto bar_assm_length = (bar_length + glue_thickness) * bar_repeat_z;

  // Mirror construction
  xml_comp_t xml_mirror       = xml_module.child(_Unicode(mirror));
  auto       mirror_width     = xml_mirror.width();
  auto       mirror_height    = xml_mirror.height();
  auto       mirror_thickness = xml_mirror.thickness();
  Box        mirror_box("mirror_box", mirror_height / 2, mirror_width / 2, mirror_thickness / 2);
  Volume     mirror_vol("mirror_vol", mirror_box, desc.material(xml_mirror.materialStr()));
  mirror_vol.setVisAttributes(desc.visAttributes(xml_mirror.visStr()));

  // Mirror optical surface
  auto        surfMgr = desc.surfaceManager();
  auto        surf    = surfMgr.opticalSurface("DIRC_MirrorOpticalSurface");
  SkinSurface skin(desc, det, Form("dirc_mirror_optical_surface"), surf, mirror_vol);
  skin.isValid();

  // Envelope for bars + mirror
  Box Envelope_box("Envelope_box", (mirror_height + 1*mm)/2, 5*(bar_width + 0.15*mm), 2*(bar_length + glue_thickness) + 0.5*mirror_thickness);
  Volume Envelope_box_vol("Envelope_box_vol", Envelope_box, desc.material("AirOptical"));
  dirc_module.placeVolume(Envelope_box_vol, Position(0, 0, 0.5*mirror_thickness));

  for (int y_index = 0; y_index < bar_repeat_y; y_index++) {
    double y = 0.5 * bar_assm_width - 0.5 * bar_width - (bar_width + bar_gap) * y_index;
    for (int z_index = 0; z_index < bar_repeat_z; z_index++) {
      double z = 0.5 * bar_assm_length - 0.5 * mirror_thickness - 0.5 * bar_length - (bar_length + glue_thickness) * z_index;
      Envelope_box_vol.placeVolume(glue_vol, Position(0, y, z - 0.5 * (bar_length + glue_thickness)));
      Envelope_box_vol.placeVolume(bar_vol, Position(0, y, z)).addPhysVolID("section", z_index).addPhysVolID("bar", y_index);
    }
  }

  // Place mirror
  Envelope_box_vol.placeVolume(mirror_vol, Position(0, 0, 0.5 * bar_assm_length));

  // Prism variables
  xml_comp_t xml_prism        = xml_module.child(_Unicode(prism));
  double     prism_angle      = xml_prism.angle();
  double     prism_width      = xml_prism.width();
  double     prism_length     = xml_prism.length();
  double     prism_short_edge = getAttrOrDefault(xml_prism, _Unicode(short_edge), 50 * mm);
  double     prism_long_edge  = prism_short_edge + prism_length * tan(prism_angle);

  // Lens variables
  xml_comp_t xml_lens    = xml_module.child(_Unicode(lens));
  double     lens_shift  = getAttrOrDefault(xml_lens, _Unicode(shift), 0 * mm);
  double lens_width = getAttrOrDefault(xml_lens, _Unicode(width), 35 * mm);
  double lens_thickness = getAttrOrDefault(xml_lens, _Unicode(thickness), 12 * mm);
  double lens_r1        = getAttrOrDefault(xml_lens, _Unicode(r1), 62 * mm);
  double lens_r2        = getAttrOrDefault(xml_lens, _Unicode(r2), 36 * mm);

  // Lens construction
  // 3-layer spherical lens ---

  double lens_radius = sqrt(lens_width * lens_width / 4. + bar_height * bar_height / 4.);

  double lens_min_thickness = 2.0 * mm;

  double ztrans1 = -lens_thickness / 2. - sqrt(lens_r1 * lens_r1 - lens_radius * lens_radius) + lens_min_thickness;
  double ztrans2 = -lens_thickness / 2. - sqrt(lens_r2 * lens_r2 - lens_radius * lens_radius) + lens_min_thickness * 2;

  Box  lens_symm_box("lens_symm_box", 0.5 * prism_short_edge, 0.5 * lens_width, 0.5 * lens_thickness);
  Volume Envelope_lens_vol("Envelope_lens_vol", lens_symm_box, desc.material("AirOptical"));

  Tube lens_symm_tube(0, lens_radius, 0.5 * lens_thickness);

  Sphere lens_sphere1(0, lens_r1);
  Sphere lens_sphere2(0, lens_r2);

  IntersectionSolid lens_box("lens_box", lens_symm_box, lens_symm_box, Position(0, 0, -lens_min_thickness * 2));
  IntersectionSolid lens_tube("lens_tube", lens_symm_tube, lens_symm_box, Position(0, 0, lens_min_thickness * 2));
  UnionSolid        lens_box_tube("lens_box_tube", lens_box, lens_tube);

  IntersectionSolid lens_layer1_solid("lens_layer1_solid", lens_box_tube, lens_sphere1, Position(0, 0, -ztrans1));
  SubtractionSolid  lens_layer23_solid("lens_layer23_solid", lens_box_tube, lens_sphere1, Position(0, 0, -ztrans1));

  IntersectionSolid lens_layer2_solid("lens_layer2_solid", lens_layer23_solid, lens_sphere2, Position(0, 0, -ztrans2));
  SubtractionSolid  lens_layer3_solid("lens_layer3_solid", lens_layer23_solid, lens_sphere2, Position(0, 0, -ztrans2));

  Volume lens_layer1_vol("lens_layer1_vol", lens_layer1_solid,
                         desc.material(xml_lens.attr<std::string>(_Unicode(material1))));
  Volume lens_layer2_vol("lens_layer2_vol", lens_layer2_solid,
                         desc.material(xml_lens.attr<std::string>(_Unicode(material2))));
  Volume lens_layer3_vol("lens_layer3_vol", lens_layer3_solid,
                         desc.material(xml_lens.attr<std::string>(_Unicode(material3))));

  lens_layer1_vol.setVisAttributes(desc.visAttributes(xml_lens.attr<std::string>(_Unicode(vis1))));
  lens_layer2_vol.setVisAttributes(desc.visAttributes(xml_lens.attr<std::string>(_Unicode(vis2))));
  lens_layer3_vol.setVisAttributes(desc.visAttributes(xml_lens.attr<std::string>(_Unicode(vis3))));

  double   lens_position_x = lens_shift;
  double   lens_position_z = -0.5 * (bar_assm_length + lens_thickness);

  for(int y_index = 0; y_index < bar_repeat_y; y_index++)
    {
      double lens_position_y = y_index*lens_width - 0.5*(prism_width - lens_width);

      Position lens_position(lens_position_x, lens_position_y, lens_position_z);
      dirc_module.placeVolume(Envelope_lens_vol, lens_position);
    }

  Envelope_lens_vol.placeVolume(lens_layer1_vol);
  Envelope_lens_vol.placeVolume(lens_layer2_vol);
  Envelope_lens_vol.placeVolume(lens_layer3_vol);

  // Prism construction
  Trap   prism_trap = MakeTrap("prism_trap", prism_width, prism_length, prism_long_edge, prism_short_edge);
  Volume prism_vol("prism_vol", prism_trap, desc.material(xml_prism.materialStr()));
  prism_vol.setVisAttributes(desc.visAttributes(xml_prism.visStr()));

  double    prism_position_x = (prism_long_edge + prism_short_edge) / 4. - 0.5 * prism_short_edge + lens_shift;
  double    prism_position_z = -0.5 * (bar_assm_length + prism_length) - lens_thickness;
  RotationX prism_rotation(M_PI / 2.);
  Position  prism_position(prism_position_x, 0, prism_position_z);

  // Envelope for prism + mcp

  double Envelope_trap_width = prism_width + 1*mm;
  double Envelope_trap_length = prism_length + 1*mm; // mcp thickness is 1 mm
  double Envelope_trap_short_edge = prism_short_edge + 1*mm;
  double Envelope_trap_long_edge = Envelope_trap_short_edge + Envelope_trap_length * tan(prism_angle);

  Trap Envelope_trap = MakeTrap("Envelope_trap", Envelope_trap_width, Envelope_trap_length, Envelope_trap_long_edge, Envelope_trap_short_edge);
  Position Envelope_trap_position(prism_position_x, 0, prism_position_z - 0.5*mm);

  Volume Envelope_trap_vol("Envelope_trap_vol", Envelope_trap, desc.material("AirOptical"));
  dirc_module.placeVolume(Envelope_trap_vol, Transform3D(prism_rotation, Envelope_trap_position));

  Envelope_trap_vol.placeVolume(prism_vol, Position(0, 0.5*mm, 0));

  // MCP variables
  xml_comp_t xml_mcp       = xml_module.child(_Unicode(mcp));
  double     mcp_thickness = xml_mcp.thickness();
  double     mcp_height    = xml_mcp.height();
  double     mcp_width     = xml_mcp.width();

  // MCP construction
  Box    mcp_box("mcp_box", mcp_height / 2, mcp_width / 2, mcp_thickness / 2);
  Volume mcp_vol("mcp_vol", mcp_box, desc.material(xml_mcp.materialStr()));
  mcp_vol.setVisAttributes(desc.visAttributes(xml_mcp.visStr())).setSensitiveDetector(sens);

  double mcp_position_z = -0.5 * prism_length;
  Position mcp_position(prism_position_x, mcp_position_z, 0);
  RotationX mcp_rotation(-M_PI / 2.);

  Envelope_trap_vol.placeVolume(mcp_vol, Transform3D(mcp_rotation, mcp_position));

  // Place modules
  const int    module_repeat = xml_module.repeat();
  const double dphi          = 2. * M_PI / module_repeat;
  for (int i = 0; i < module_repeat; i++) {
    double phi = dphi * i;
    double x   = det_ravg * cos(phi);
    double y   = det_ravg * sin(phi);

    Transform3D tr(RotationZ(phi), Position(x, y, 0));
    det_volume.placeVolume(dirc_module, tr).addPhysVolID("module", i);
  }

  // Construct support
  xml_comp_t xml_support = xml_det.child(_U(support));
  Assembly dirc_support("DIRCSupport");
  dirc_support.setVisAttributes(desc.visAttributes(xml_support.visStr()));

  // Rail
  xml_comp_t xml_rail    = xml_support.child(_Unicode(rail));
  xml_dim_t  rail_pos = xml_rail.position();
  double     rail_height = xml_rail.height();
  double     rail_width2 = xml_rail.width();
  double     rail_distance_to_chord2 = rail_width2/2 / tan(dphi/2);
  double     rail_distance_to_chord1 = rail_distance_to_chord2 - rail_height;
  double     rail_width1 = 2*rail_distance_to_chord1 * tan(dphi/2);
  double     rail_length = xml_rail.length();
  Trap       rail_trap("rail_trap", rail_length / 2, 0, 0,
                       rail_height / 2, rail_width1 / 2, rail_width2 / 2, 0,
                       rail_height / 2, rail_width1 / 2, rail_width2 / 2, 0);
  Volume     rail_vol("rail_vol", rail_trap, desc.material(xml_rail.materialStr()));
  rail_vol.setVisAttributes(desc.visAttributes(xml_rail.visStr()));

  // Place rail
  Position rail_position(rail_pos.x(), rail_pos.y(), rail_pos.z());
  RotationZ rail_rotation(-M_PI / 2.);
  dirc_support.placeVolume(rail_vol, Transform3D(rail_rotation, rail_position));

  // Place support
  for (int i = 0; i < module_repeat; i++) {
    double phi = dphi * i + dphi / 2;
    double x   = det_ravg * cos(phi);
    double y   = det_ravg * sin(phi);

    Transform3D tr(RotationZ(phi), Position(x, y, 0));
    det_volume.placeVolume(dirc_support, tr);
  }


  return det;
}

static dd4hep::Trap MakeTrap(const std::string& pName, double pZ, double pY, double pX, double pLTX)
{
  // Fixed Trap constructor. This function is a workaround of this bug:
  // https://github.com/AIDASoft/DD4hep/issues/850
  // Should be used instead of dd4hep::Trap(pName, pZ, pY, pX, pLTX) constructor

  double fDz         = 0.5 * pZ;
  double fTthetaCphi = 0;
  double fTthetaSphi = 0;
  double fDy1        = 0.5 * pY;
  double fDx1        = 0.5 * pX;
  double fDx2        = 0.5 * pLTX;
  double fTalpha1    = atan(0.5 * (pLTX - pX) / pY);
  double fDy2        = fDy1;
  double fDx3        = fDx1;
  double fDx4        = fDx2;
  double fTalpha2    = fTalpha1;

  return Trap(pName, fDz, fTthetaCphi, fTthetaSphi, fDy1, fDx1, fDx2, fTalpha1, fDy2, fDx3, fDx4, fTalpha2);
}

DECLARE_DETELEMENT(epic_DIRC, createDetector)
