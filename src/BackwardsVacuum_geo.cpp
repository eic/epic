#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Vacuum drift volume in far backwards region
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector /*sens*/)
{

  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  string vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BackwardsBox");

  // Dimensions of main beamline pipe
  xml::Component dim    = x_det.child(_Unicode(dimensions));
  double         WidthL = dim.attr<double>(_Unicode(xL));
  double         WidthR = dim.attr<double>(_Unicode(xR));

  double Width     = (WidthL + WidthR) / 2;
  double Height    = dim.y();
  double Thickness = dim.z();

  // Materials
  Material Vacuum = desc.material("Vacuum");
  Material Steel  = desc.material("StainlessSteel");

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
  Box Extended_Beam_Box(Width + wall, Height + wall, Thickness);

  // Central vacuum box
  Box Extended_Vacuum_Box(Width, Height, Thickness);

  Solid Wall_Box   = Extended_Beam_Box;
  Solid Vacuum_Box = Extended_Vacuum_Box;

  Assembly DetAssembly(detName + "_assembly");
  Assembly DetAssemblyAir(detName + "_assembly_air");

  DetElement det(detName, detID);

  //-----------------------------------------------------------------
  // Add Tagger box containers and vacuum box extension for modules
  //-----------------------------------------------------------------
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod) {

    int    moduleID   = dd4hep::getAttrOrDefault<int>(mod, _Unicode(id), 0);
    string moduleName = dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(modname), "Tagger0");

    // Offset from the electron beam
    double tagoff = dd4hep::getAttrOrDefault<double>(mod, _Unicode(offset_min), 50.0 * mm);

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
    //    auto vac_w = w+wall/2; // Adds space for wall to be added in tagger geometry
    auto vac_h = h;

    // Width and height of box volume
    auto box_w = w + wall;
    auto box_h = h + wall;

    auto theta      = thetamin;
    auto offsetx    = -(box_w - wall) * (cos(theta));
    auto offsetz    = (box_w - wall) * (sin(theta));
    auto vacoffsetx = -vac_w * (cos(theta));
    auto vacoffsetz = vac_w * (sin(theta));
    auto l          = (tagoff) / (sin(theta));
    //    auto l          = (tagoff)/(sin(theta));

    auto tagoffsetx = vacoffsetx - (l + tagboxL) * sin(theta);
    auto tagoffsetz = vacoffsetz - (l + tagboxL) * cos(theta);
    //     auto tagoffsetx = vacoffsetx-(l+tagboxL/2)*sin(theta);
    //     auto tagoffsetz = vacoffsetz-(l+tagboxL/2)*cos(theta);

    if (max_align) {
      theta      = thetamax;
      offsetx    = (box_w - wall) * (cos(theta));
      offsetz    = -(box_w - wall) * (sin(theta));
      vacoffsetx = vac_w * (cos(theta));
      vacoffsetz = -vac_w * (sin(theta));
      l          = (2 * offsetx + tagoff) / sin(theta);
      tagoffsetx = vacoffsetx - (l + tagboxL) * sin(theta);
      tagoffsetz = vacoffsetz - (l + tagboxL) * cos(theta);
      //       tagoffsetx = -wall+vacoffsetx-(l+tagboxL/2)*sin(theta);
      //       tagoffsetz = vacoffsetz-(l+tagboxL/2)*cos(theta);
    }

    Box TagWallBox(box_w, box_h, l + tagboxL);
    Box TagVacBox(vac_w, vac_h, l + tagboxL);

    RotationY rotate(theta);

    Volume mother = DetAssemblyAir;

    if (extend_vacuum) {
      Wall_Box   = UnionSolid(Wall_Box, TagWallBox, Transform3D(rotate, Position(offsetx, 0, offsetz)));
      Vacuum_Box = UnionSolid(Vacuum_Box, TagVacBox, Transform3D(rotate, Position(vacoffsetx, 0, vacoffsetz)));
      mother     = DetAssembly;
    }

    Assembly TaggerAssembly("tagAssembly");

    PlacedVolume pv_mod2 = mother.placeVolume(
        TaggerAssembly,
        Transform3D(rotate, Position(tagoffsetx, 0,
                                     tagoffsetz))); // Very strange y is not centered and offset needs correcting for...
    DetElement moddet(moduleName, moduleID);
    moddet.setPlacement(pv_mod2);
    det.add(moddet);
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

    Box  Entry_Beam_Box(ED_X + wall, ED_Y + wall, ED_Z);
    Box  Entry_Vacuum_Box(ED_X, ED_Y, ED_Z - wall);
    Tube Lumi_Exit(0, Lumi_R, ED_Z);

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

  Volume vacVol("Vacuum_Box", Vacuum_Box_Sub, Vacuum);
  // vacVol.setVisAttributes(desc.visAttributes("RedGreenVis"));
  vacVol.placeVolume(DetAssembly);
  Volume wallVol("Tagger_Box", Wall_Box_Sub, Steel);
  wallVol.setVisAttributes(desc.visAttributes(vis_name));
  wallVol.placeVolume(vacVol);

  Assembly backAssembly("assembly");
  backAssembly.placeVolume(wallVol);
  backAssembly.placeVolume(DetAssemblyAir);

  // placement in mother volume
  Transform3D  tr(RotationY(rot.theta()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(backAssembly, tr);
  det.setPlacement(detPV);

  return det;
}

DECLARE_DETELEMENT(BackwardsVacuum, create_detector)
