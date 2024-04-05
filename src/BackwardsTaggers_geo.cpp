// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Simon Gardner

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
#include "DD4hepDetectorHelper.h"


//////////////////////////////////////////////////
// Far backwards vacuum drift volume
// Low Q2 tagger trackers placed either in or out of vacuum
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

// Helper function to make the tagger tracker detectors
static void Make_Tagger(Detector& desc, xml_coll_t& mod, Assembly& env, DetElement modElement, SensitiveDetector& sens);


static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)
{

  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  DetElement det(detName, detID);
  dd4hep::xml::setDetectorTypeFlag(x_det, det);

  string vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BackwardsBox");

  // Dimensions of main beamline pipe
  xml::Component dim    = x_det.child(_Unicode(dimensions));
  double         WidthL = dim.attr<double>(_Unicode(xL));
  double         WidthR = dim.attr<double>(_Unicode(xR));

  double Width     = (WidthL + WidthR) / 2;
  double Height    = dim.y();
  double Thickness = dim.z();

  // Materials
  Material Vacuum  = desc.material("Vacuum");
  Material Steel   = desc.material("StainlessSteel");

  // Central focal point of the geometry
  xml::Component pos = x_det.child(_Unicode(focus));
  double         off = pos.z();

  // Beamline rotation
  xml_dim_t rot = x_det.rotation();

  // Beampipe thickness
  double wall = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(wall), 1 * mm);

  // Make bounding box to make IntersectionSolid with other components
  xml::Component BB      = x_det.child(_Unicode(bounding));
  double         BB_MinX = BB.xmin();
  double         BB_MinY = BB.ymin();
  double         BB_MinZ = BB.zmin();
  double         BB_MaxX = BB.xmax();
  double         BB_MaxY = BB.ymax();
  double         BB_MaxZ = BB.zmax();

  double BB_X = abs(BB_MaxX - BB_MinX);
  double BB_Y = abs(BB_MaxY - BB_MinY);
  double BB_Z = abs(BB_MaxZ - BB_MinZ);

  Box Far_Backwards_Box(BB_X, BB_Y, BB_Z);

  // Entry box geometry description joining magnet, taggers and lumi
  xml::Component EB     = x_det.child(_Unicode(exitdim));
  double         ED_X   = EB.x();
  double         ED_Y   = EB.y();
  double         ED_Z   = off - EB.attr<double>(_Unicode(lumiZ));
  double         Lumi_R = EB.attr<double>(_Unicode(lumiR));

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
  Assembly DetAssemblyAir("Tagger_air_assembly");
  int      nVacuum = 0;
  int      nAir    = 0;

  //-----------------------------------------------------------------
  // Add Tagger box containers and vacuum box extension for modules
  //-----------------------------------------------------------------
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod) {

    int    moduleID   = dd4hep::getAttrOrDefault<int>(mod, _Unicode(id), 0);
    string moduleName = dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(name), "Tagger0");

    // Offset from the electron beam
    double tagoff  = dd4hep::getAttrOrDefault<double>(mod, _Unicode(offset_min), 50.0 * mm);

    // Overlap left beyond theta setting
    double overlap = dd4hep::getAttrOrDefault<double>(mod, _Unicode(overlap), 0.0 * mm);

    // Theta coverage expected
    double thetamin = dd4hep::getAttrOrDefault<double>(mod, _Unicode(theta_min), 0.030 * rad) - rot.theta();
    double thetamax = dd4hep::getAttrOrDefault<double>(mod, _Unicode(theta_max), 0.030 * rad) - rot.theta();

    // Align box to max or minimum theta expected at the tagger from focal point
    bool max_align = dd4hep::getAttrOrDefault<bool>(mod, _Unicode(max_align), false);

    // Extend the beam vacuum around the taggers
    bool extend_vacuum = dd4hep::getAttrOrDefault<bool>(mod, _Unicode(extend_vacuum), true);

    // Size f the actual tagger box, replicated in BackwardsTagger
    xml_dim_t moddim  = mod.child(_Unicode(dimensions));
    double    w       = moddim.x();
    double    h       = moddim.y();
    double    tagboxL = moddim.z();

    // Width and height of vacuum volume
    auto vac_w = w;
    auto vac_h = h;

    // Width and height of box volume
    auto box_w = w + wall;
    auto box_h = h + wall;

    // Angle in relation to the main beam
    auto theta      = thetamin;


    auto offsetx    = -(box_w - wall) * (cos(theta));
    auto offsetz    = (box_w - wall) * (sin(theta));
    auto vacoffsetx = -vac_w * (cos(theta));
    auto vacoffsetz = vac_w * (sin(theta));
    auto l          = (tagoff) / (sin(theta)) + tagboxL;

    auto tagoffsetx = vacoffsetx - (l) * sin(theta);
    auto tagoffsetz = vacoffsetz - (l) * cos(theta);

    if (max_align) {
      theta      = thetamax;
      offsetx    = (overlap+box_w - wall) * (cos(theta));
      offsetz    = -(overlap+box_w - wall) * (sin(theta));
      vacoffsetx = (overlap+vac_w) * (cos(theta));
      vacoffsetz = -(overlap+vac_w) * (sin(theta));
      l          = (2 * offsetx + tagoff) / sin(theta);
      tagoffsetx = vacoffsetx - (l) * sin(theta);
      tagoffsetz = vacoffsetz - (l) * cos(theta);
    }

    Box TagWallBox(box_w, box_h, l + wall);
    Box TagVacBox(vac_w, vac_h, l);

    RotationY rotate(theta);

    Volume mother = DetAssemblyAir;

    if (extend_vacuum) {
      Wall_Box   = UnionSolid(Wall_Box, TagWallBox, Transform3D(rotate, Position(offsetx, 0, offsetz)));
      Vacuum_Box = UnionSolid(Vacuum_Box, TagVacBox, Transform3D(rotate, Position(vacoffsetx, 0, vacoffsetz)));
      mother     = DetAssembly;
      nVacuum++;
    } else {
      nAir++;
    }

    //Assembly TaggerAssembly("Tagger_module_assembly");
    Assembly TaggerAssembly(moduleName);

    PlacedVolume pv_mod2 = mother.placeVolume(
        TaggerAssembly,
        Transform3D(rotate, Position(tagoffsetx, 0,
                                     tagoffsetz))); // Very strange y is not centered and offset needs correcting for...
    DetElement moddet(det,moduleName, moduleID);
    pv_mod2.addPhysVolID("module", moduleID);
    moddet.setPlacement(pv_mod2);
    
    auto &moduleParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(moddet);
    moduleParams.set<string>("layer_pattern", "Tagger_tracker_layer\\d");

    Make_Tagger(desc, mod, TaggerAssembly, moddet, sens);

    // //Loop through moddet checking the names of its children
    // std::cout << "Children of " << moddet.name() << std::endl;
    // for (auto child : moddet.children()) {
    //   std::cout << child.first << std::endl;
    // }

  }

  //Loop through det checking the names of its children
  std::cout << "Children of " << det.name() << std::endl;
  for (auto child : det.children()) {
    std::cout << child.first << std::endl;
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
  double cutZ     = (ED_X - exitDist * tan(-rot.theta())) * sin(rot.theta()) + exitDist * cos(rot.theta());
  double cutXwall = (ED_X - wall - exitDist * tan(-rot.theta())) * cos(rot.theta());
  double cutZwall = (ED_X - wall - exitDist * tan(-rot.theta())) * sin(rot.theta()) + exitDist * cos(rot.theta());

  Wall_Box = IntersectionSolid(Wall_Box, Cut_Box, Transform3D(RotationY(exitTheta), Position(xbox - cutX, 0, cutZ)));
  Vacuum_Box =
      IntersectionSolid(Vacuum_Box, Cut_Box, Transform3D(RotationY(exitTheta), Position(xbox - cutXwall, 0, cutZwall)));

  //-----------------------------------------------------------------
  // Cut solids so they are only in the far backwards box
  //-----------------------------------------------------------------
  RotationY rotate2(-rot.theta());
  Position  position(0, 0, (exitDist - BB_Z) / cos(rot.theta()));

  IntersectionSolid Wall_Box_Sub(Wall_Box, Far_Backwards_Box, Transform3D(rotate2, position));
  IntersectionSolid Vacuum_Box_Sub(Vacuum_Box, Far_Backwards_Box, Transform3D(rotate2, position));
  SubtractionSolid  Wall_Box_Out(Wall_Box_Sub, Vacuum_Box_Sub);

  Volume vacVol("TaggerStation_Vacuum", Vacuum_Box_Sub, Vacuum);
  vacVol.setVisAttributes(desc.visAttributes("BackwardsVac"));
  if (nVacuum > 0)
    vacVol.placeVolume(DetAssembly);

  Volume wallVol("TaggerStation_Container", Wall_Box_Out, Steel);
  wallVol.setVisAttributes(desc.visAttributes(vis_name));

  Assembly backAssembly(detName + "_assembly");
  backAssembly.placeVolume(wallVol);
  backAssembly.placeVolume(vacVol);

  if (nAir > 0)
    backAssembly.placeVolume(DetAssemblyAir);

  // placement in mother volume
  Transform3D  tr(RotationY(rot.theta()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(backAssembly, tr);
  detPV.addPhysVolID("system", detID);

  det.setPlacement(detPV);


  return det;
}

static void Make_Tagger(Detector& desc, xml_coll_t& mod, Assembly& env, DetElement modElement, SensitiveDetector& sens)
{

  sens.setType("tracker");

  Material Air     = desc.material("Air");
  Material Silicon = desc.material("Silicon");

  xml_dim_t moddim  = mod.child(_Unicode(dimensions));
  double    tag_w   = moddim.x();
  double    tag_h   = moddim.y();
  double    tagboxL = moddim.z();

  Volume Tagger_Air;

  double airThickness    = 0;
  double vacuumThickness = tagboxL;

  // Add window layer and air-vacuum boxes
  for (xml_coll_t lay(mod, _Unicode(foilLayer)); lay; ++lay) {

    string layerType      = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "foil");
    string layerVis       = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerShieldingVis");
    double layerZ         = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerRot       = dd4hep::getAttrOrDefault<double>(lay, _Unicode(angle), 45*deg);
    double layerThickness = dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 100 * um);
    string layerMaterial  = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(material), "Copper");

    Material FoilMaterial = desc.material(layerMaterial);

    RotationY rotate(layerRot);

    airThickness    = tagboxL - layerZ;
    vacuumThickness = tagboxL - airThickness;

    Box Box_Air(tag_w, tag_h, airThickness / 2);
    Tagger_Air = Volume("AirVolume", Box_Air, Air);
    Tagger_Air.setVisAttributes(desc.visAttributes("BackwardsAir"));

    Box    Foil_Box(tag_w/cos(layerRot)-0.5*layerThickness*tan(layerRot), tag_h, layerThickness / 2);
    Volume layVol("FoilVolume", Foil_Box, FoilMaterial);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    env.placeVolume(layVol, Transform3D(rotate,Position(0, 0, tagboxL+tag_w*tan(layerRot))));

    // Currently only one "foil" layer implemented
    break;
  }

  // Add window layer and air-vacuum boxes
  for (xml_coll_t lay(mod, _Unicode(windowLayer)); lay; ++lay) {

    string layerType      = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "window");
    string layerVis       = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerShieldingVis");
    double layerZ         = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerRot       = dd4hep::getAttrOrDefault<double>(lay, _Unicode(angle), 0);
    double layerThickness = dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 1 * mm);
    string layerMaterial  = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(material), "Copper");

    Material WindowMaterial = desc.material(layerMaterial);

    RotationY rotate(layerRot);

    airThickness    = tagboxL - layerZ;
    vacuumThickness = tagboxL - airThickness;

    Box Box_Air(tag_w, tag_h, airThickness / 2);
    Tagger_Air = Volume("AirVolume", Box_Air, Air);
    Tagger_Air.setVisAttributes(desc.visAttributes("BackwardsAir"));

    Box    Window_Box(tag_w, tag_h, layerThickness / 2);
    Volume layVol("WindowVolume", Window_Box, WindowMaterial);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    Tagger_Air.placeVolume(layVol, Position(0, 0, airThickness / 2 - layerThickness / 2));

    env.placeVolume(Tagger_Air, Position(0, 0, tagboxL - airThickness / 2));

    // Currently only one "window" layer implemented
    break;
  }

  // Add Hodoscope layers
  int N_layers = 0;
  for (xml_coll_t lay(mod, _Unicode(trackLayer)); lay; ++lay, ++N_layers) {

    int    layerID        = dd4hep::getAttrOrDefault<int>(lay, _Unicode(id), 0);
    string layerType      = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "timepix");
    string layerVis       = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerLayerVis");
    double layerRot       = dd4hep::getAttrOrDefault<double>(lay, _Unicode(angle), 0);
    double layerZ         = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerThickness = dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 200 * um);

    double envelope_r_min = dd4hep::getAttrOrDefault<double>(lay, _Unicode(envelope_r_min),  -1.0*mm);
    double envelope_r_max = dd4hep::getAttrOrDefault<double>(lay, _Unicode(envelope_r_max),  1.0*mm);
    double envelope_z_min = dd4hep::getAttrOrDefault<double>(lay, _Unicode(envelope_z_min),  -1.0*mm);
    double envelope_z_max = dd4hep::getAttrOrDefault<double>(lay, _Unicode(envelope_z_max),  1.0*mm);

    Volume mother          = env;
    double MotherThickness = tagboxL;

    RotationY rotate(layerRot);

    if (layerZ > vacuumThickness) {
      mother = Tagger_Air;
      layerZ -= vacuumThickness;
      MotherThickness = airThickness / 2;
    }

    Box    Layer_Box(tag_w, tag_h, layerThickness / 2);

    std::string layerName = "Tagger_tracker_layer" + std::to_string(layerID);
    Volume layVol(layerName, Layer_Box, Silicon);
    layVol.setSensitiveDetector(sens);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    
    PlacedVolume pv_layer = mother.placeVolume(layVol, Transform3D(rotate, Position(0, 0, MotherThickness - layerZ + layerThickness / 2)));
    pv_layer.addPhysVolID("layer", layerID);

    DetElement laydet(modElement,pv_layer.volume().name(), layerID);

    // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
    Vector3D u(-1., 0., 0.);
    Vector3D v( 0.,-1., 0.);
    Vector3D n( 0., 0., 1.);

    // Add surface to layer for acts reconstruction 
    SurfaceType type(SurfaceType::Sensitive);

    layVol->GetShape()->ComputeBBox();
    auto &layerParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(laydet);
    layerParams.set<string>("axis_definitions", "XZY");
    layerParams.set<double>("envelope_r_min", envelope_r_min);
    layerParams.set<double>("envelope_r_max", envelope_r_max);
    layerParams.set<double>("envelope_z_min", envelope_z_min);
    layerParams.set<double>("envelope_z_max", envelope_z_max);

    for (xml_coll_t lmat(mod, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }


    VolPlane surf(layVol, type, 1.0, 2.0, u, v, n); //,o ) ;
    volSurfaceList(laydet)->push_back(surf);

    laydet.setPlacement(pv_layer);

  }
}

DECLARE_DETELEMENT(BackwardsTagger, create_detector)
