// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022-2025 Simon Gardner

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

//////////////////////////////////////////////////
// Far backwards vacuum drift volume
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

// Helper function to make the tagger tracker detectors
static void Make_Tagger(Detector& desc, xml_coll_t& mod, Assembly& env);

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector /* sens */) {

  xml_det_t x_det = e;
  string detName  = x_det.nameStr();
  int detID       = x_det.id();

  DetElement det(detName, detID);

  string vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BackwardsBox");

  // Dimensions of main beamline pipe
  xml::Component dim = x_det.child(_Unicode(dimensions));
  double WidthL      = dim.attr<double>(_Unicode(xL));
  double WidthR      = dim.attr<double>(_Unicode(xR));

  double Width     = (WidthL + WidthR) / 2;
  double Height    = dim.y();
  double Thickness = dim.z();

  // Materials
  Material Vacuum = desc.material("Vacuum");
  Material Steel  = desc.material("StainlessSteelSAE304");

  // Central focal point of the geometry
  xml::Component pos = x_det.child(_Unicode(focus));
  double off         = pos.z();

  // Beamline rotation
  xml_dim_t rot = x_det.rotation();

  // Beampipe thickness
  double wall = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(wall), 1 * mm);

  // Make bounding box to make IntersectionSolid with other components
  xml::Component BB = x_det.child(_Unicode(bounding));
  double BB_MinX    = BB.xmin();
  double BB_MinY    = BB.ymin();
  double BB_MinZ    = BB.zmin();
  double BB_MaxX    = BB.xmax();
  double BB_MaxY    = BB.ymax();
  double BB_MaxZ    = BB.zmax();

  double BB_X = abs(BB_MaxX - BB_MinX);
  double BB_Y = abs(BB_MaxY - BB_MinY);
  double BB_Z = abs(BB_MaxZ - BB_MinZ);

  Box Far_Backwards_Box(BB_X, BB_Y, BB_Z);

  // Entry box geometry description joining magnet, taggers and lumi
  xml::Component EB = x_det.child(_Unicode(exitdim));
  double ED_X       = EB.x();
  double ED_Y       = EB.y();
  double ED_Z       = off - EB.attr<double>(_Unicode(lumiZ));
  double Lumi_R     = EB.attr<double>(_Unicode(lumiR));

  // Maximum theta to exit the dipole from
  double exitTheta = EB.attr<double>(_Unicode(maxTheta));

  // Generic box for making intersection solid with
  double xbox = 10 * m;
  double ybox = 10 * m;
  double zbox = 50 * m;

  Box Cut_Box(xbox, ybox, zbox);

  // Central pipe box
  //Tube Extended_Beam_Box(Width,Width+wall,Thickness); // More realistic tube pipe
  Box Extended_Beam_Box(Width + wall, Height + wall, Thickness); // Simpler box pipe

  // Central vacuum box
  //Tube Extended_Vacuum_Box(0,Width,Thickness); // More realistic tube pipe
  Box Extended_Vacuum_Box(Width, Height, Thickness); // Simpler box pipe

  Solid Wall_Box   = Extended_Beam_Box;
  Solid Vacuum_Box = Extended_Vacuum_Box;

  Assembly DetAssembly("Tagger_vacuum_assembly");

  //-----------------------------------------------------------------
  // Add Tagger box containers and vacuum box extension for modules
  //-----------------------------------------------------------------
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod) {

    int moduleID      = dd4hep::getAttrOrDefault<int>(mod, _Unicode(id), 0);
    string moduleName = dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(name), "Tagger0");

    xml_dim_t mod_pos_global = mod.child(_U(position));
    xml_dim_t mod_rot_global = mod.child(_U(rotation));
    Position mod_pos(mod_pos_global.x(), mod_pos_global.y(), mod_pos_global.z());
    Position vac_pos(mod_pos_global.x() + wall / 2 * cos(mod_rot_global.theta()),
                     mod_pos_global.y(),
                     mod_pos_global.z() + wall / 2 * sin(mod_rot_global.theta()));
    RotationY mod_rot(mod_rot_global.theta() - rot.theta());

    // Size f the actual tagger box, replicated in BackwardsTagger
    xml_dim_t moddim = mod.child(_Unicode(dimensions));
    double vac_w     = moddim.x() / 2;
    double vac_h     = moddim.y() / 2;
    double vac_l     = moddim.z() / 2;

    // Width and height of box volume
    auto box_w = vac_w + wall;
    auto box_h = vac_h + wall;

    // auto l          = 100 * cm; // Length of the tagger box, arbitrary as long as > thickness of tagger

    Box TagWallBox(box_w, box_h, vac_l);
    Box TagVacBox(vac_w + wall / 2, vac_h, vac_l);

    Wall_Box   = UnionSolid(Wall_Box, TagWallBox, Transform3D(mod_rot, mod_pos));
    Vacuum_Box = UnionSolid(Vacuum_Box, TagVacBox, Transform3D(mod_rot, vac_pos));

    Assembly TaggerAssembly("Tagger_module_assembly");

    PlacedVolume pv_mod2 = DetAssembly.placeVolume(
        TaggerAssembly,
        Transform3D(mod_rot,
                    mod_pos + Position(-vac_l * sin(mod_rot_global.theta() - rot.theta()), 0,
                                       -vac_l * cos(mod_rot_global.theta() - rot.theta()))));
    DetElement moddet(det, moduleName, moduleID);
    pv_mod2.addPhysVolID("module", moduleID);
    moddet.setPlacement(pv_mod2);

    Make_Tagger(desc, mod, TaggerAssembly);
  }

  //-----------------------------------------------------------------
  // Cut off any vacuum right of the main beamline
  //-----------------------------------------------------------------

  Wall_Box   = IntersectionSolid(Wall_Box, Cut_Box, Position(-xbox + Width + wall, 0, 0));
  Vacuum_Box = IntersectionSolid(Vacuum_Box, Cut_Box, Position(-xbox + Width, 0, 0));

  //-----------------------------------------------------------------
  // Luminosity connecting box
  //-----------------------------------------------------------------
  bool addLumi = dd4hep::getAttrOrDefault<bool>(x_det, _Unicode(lumi), true);

  if (addLumi) {

    Box Entry_Beam_Box(ED_X + wall, ED_Y + wall, ED_Z);
    Box Entry_Vacuum_Box(ED_X, ED_Y, ED_Z - wall);
    Tube Lumi_Exit(0, Lumi_R, ED_Z);

    // Future angled exit window and more realistic tube shaped pipe.
    // double angle = -pi/4;
    // CutTube Entry_Beam_Box  (ED_X, ED_X + wall, ED_Z,        0,2*pi, sin(angle),0,cos(angle), 0,0,1);
    // CutTube Entry_Vacuum_Box(0,    ED_X,        ED_Z - wall, 0,2*pi, sin(angle),0,cos(angle), 0,0,1);
    // CutTube Lumi_Exit       (0,    Lumi_R,      ED_Z,        0,2*pi, sin(angle),0,cos(angle), 0,0,1);

    // Add entry boxes to main beamline volume
    Wall_Box   = UnionSolid(Wall_Box, Entry_Beam_Box, Transform3D(RotationY(-rot.theta())));
    Vacuum_Box = UnionSolid(Vacuum_Box, Entry_Vacuum_Box, Transform3D(RotationY(-rot.theta())));
    Vacuum_Box = UnionSolid(Vacuum_Box, Lumi_Exit, Transform3D(RotationY(-rot.theta())));
  }

  //-----------------------------------------------------------------
  // Restrict tagger boxes into region defined by exitTheta from the dipole magnet
  //-----------------------------------------------------------------
  double exitDist = BB_MinZ - off;
  double cutX     = (ED_X - exitDist * tan(-rot.theta())) * cos(rot.theta());
  double cutZ =
      (ED_X - exitDist * tan(-rot.theta())) * sin(rot.theta()) + exitDist * cos(rot.theta());
  double cutXwall = (ED_X - wall - exitDist * tan(-rot.theta())) * cos(rot.theta());
  double cutZwall =
      (ED_X - wall - exitDist * tan(-rot.theta())) * sin(rot.theta()) + exitDist * cos(rot.theta());

  Wall_Box = IntersectionSolid(Wall_Box, Cut_Box,
                               Transform3D(RotationY(exitTheta), Position(xbox - cutX, 0, cutZ)));
  Vacuum_Box =
      IntersectionSolid(Vacuum_Box, Cut_Box,
                        Transform3D(RotationY(exitTheta), Position(xbox - cutXwall, 0, cutZwall)));

  //-----------------------------------------------------------------
  // Cut solids so they are only in the far backwards box
  //-----------------------------------------------------------------
  RotationY rotate2(-rot.theta());
  Position position(0, 0, (exitDist - BB_Z) / cos(rot.theta()));

  IntersectionSolid Wall_Box_Sub(Wall_Box, Far_Backwards_Box, Transform3D(rotate2, position));
  IntersectionSolid Vacuum_Box_Sub(Vacuum_Box, Far_Backwards_Box, Transform3D(rotate2, position));
  SubtractionSolid Wall_Box_Out(Wall_Box_Sub, Vacuum_Box_Sub);

  Volume vacVol("TaggerStation_Vacuum", Vacuum_Box_Sub, Vacuum);
  vacVol.setVisAttributes(desc.visAttributes("BackwardsVac"));
  vacVol.placeVolume(DetAssembly);

  Volume wallVol("TaggerStation_Container", Wall_Box_Out, Steel);
  wallVol.setVisAttributes(desc.visAttributes(vis_name));

  Assembly backAssembly(detName + "_assembly");
  backAssembly.placeVolume(wallVol);
  backAssembly.placeVolume(vacVol);

  // placement in mother volume
  Transform3D tr(RotationY(rot.theta()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(backAssembly, tr);
  detPV.addPhysVolID("system", detID);

  det.setPlacement(detPV);

  return det;
}

static void Make_Tagger(Detector& desc, xml_coll_t& mod, Assembly& env) {

  xml_dim_t moddim = mod.child(_Unicode(dimensions));
  double tag_w     = moddim.x() / 2;
  double tag_h     = moddim.y() / 2;

  double window_thickness = 0;

  // Add window layer and air-vacuum boxes
  for (xml_coll_t lay(mod, _Unicode(windowLayer)); lay; ++lay) {

    string layerType = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "window");
    string layerVis =
        dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerShieldingVis");
    double layerRot = dd4hep::getAttrOrDefault<double>(lay, _Unicode(angle), 0);
    double layerThickness =
        dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 1 * mm);
    string layerMaterial = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(material), "Copper");
    
    window_thickness = layerThickness;

    Material WindowMaterial = desc.material(layerMaterial);

    RotationY rotate(layerRot);

    Box Window_Box(tag_w, tag_h, layerThickness / 2);
    Volume layVol("WindowVolume", Window_Box, WindowMaterial);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    env.placeVolume(layVol, Position(0, 0, layerThickness / 2));

    // Currently only one "window" layer implemented
    break;
  }

  // Add window layer and air-vacuum boxes
  for (xml_coll_t lay(mod, _Unicode(foilLayer)); lay; ++lay) {

    string layerType = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "foil");
    string layerVis =
        dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerShieldingVis");
    // double layerZ   = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerRot = dd4hep::getAttrOrDefault<double>(lay, _Unicode(angle), 45 * deg);
    double layerThickness =
        dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 100 * um);
    string layerMaterial = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(material), "Copper");

    Material FoilMaterial = desc.material(layerMaterial);

    RotationY rotate(layerRot);

    Box Foil_Box(tag_w / cos(layerRot) - 0.5 * layerThickness * tan(layerRot), tag_h,
                 layerThickness / 2);
    Volume layVol("FoilVolume", Foil_Box, FoilMaterial);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    env.placeVolume(layVol, Transform3D(rotate, Position(0, 0, window_thickness+tag_w * tan(layerRot))));

    // Currently only one "foil" layer implemented
    break;
  }

}

DECLARE_DETELEMENT(FarBackwardVacuum, create_detector)
