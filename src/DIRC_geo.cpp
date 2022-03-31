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

// Fixed Trap constructor. This function is a workaround of this bug:
// https://github.com/AIDASoft/DD4hep/issues/850
// Should be used instead of dd4hep::Trap(pName, pZ, pY, pX, pLTX) constructor
dd4hep::Trap MakeTrap( const std::string& pName, double pZ, double pY, double pX, double pLTX );

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t xml_det   = e;
  string    det_name = xml_det.nameStr();
  int       det_id   = xml_det.id();

  // Main detector xml element
  xml_dim_t dirc_dim   = xml_det.dimensions();
  xml_dim_t dirc_pos   = xml_det.position();
  xml_dim_t dirc_rot   = xml_det.rotation();
  double    det_rin   = dirc_dim.rmin();
  double    det_rout  = dirc_dim.rmax();
  double    SizeZ = dirc_dim.length();

  // DEBUG
  // double mirror_r1 = x_det.attr<double>(_Unicode(r1));

  // DIRC box:
  xml_comp_t xml_box_module = xml_det.child(_U(module));

  Material Vacuum = desc.material("Vacuum");
  Material air = desc.material("AirOptical");
  Material quartz = desc.material("Quartz");
  Material epotek = desc.material("Epotek");
  Material nlak33a = desc.material("Nlak33a");
  auto& bar_material = quartz;
  auto mirror_material = desc.material("Aluminum");  // mirror material

  Tube     det_geo(det_rin, det_rout, SizeZ / 2., 0., 360.0 * deg);
  //Volume   det_volume("DIRC", det_geo, Vacuum);
  Assembly   det_volume("DIRC");
  det_volume.setVisAttributes(desc.visAttributes(xml_det.visStr()));

  DetElement   det(det_name, det_id);
  Volume       mother_vol = desc.pickMotherVolume(det);

  Transform3D  tr(RotationZYX(0, dirc_rot.theta(), 0.0), Position(0.0, 0.0, dirc_pos.z()));
  PlacedVolume det_plvol = mother_vol.placeVolume(det_volume, tr);

  det_plvol.addPhysVolID("system", det_id);
  det.setPlacement(det_plvol);

  // Parts Dimentions

  int fLensId = 6; // focusing system
                   // 0    no lens
                   // 1    spherical lens
                   // 3    3-layer spherical lens
                   // 6    3-layer cylindrical lens
                   // 10   ideal lens (thickness = 0, ideal focusing)
  
  int fGeomType = 0;  // Full DIRC - 0, 1 only one plate
  int fRunType = 0;   // 0, 10 - simulation, 1, 5 - lookup table, 2,3,4 - reconstruction


  double fPrizm[4];
  fPrizm[0] = 360 * mm;
  fPrizm[1] = 300 * mm;
  fPrizm[3] = 50 * mm;
  fPrizm[2] = fPrizm[3] + 300 * tan(32 * deg) * mm;
  double fBarsGap = 0.15 * mm;

  std::cout << "DIRC: fPrizm[2] " << fPrizm[2] << std::endl;

  double fdTilt = 80 * deg;
  double fPrizmT[6];
  fPrizmT[0] = 390 * mm;
  fPrizmT[1] = (400 - 290 * cos(fdTilt)) * mm; //
  fPrizmT[2] = 290 * sin(fdTilt) * mm;       // hight
  fPrizmT[3] = 50 * mm;                      // face
  fPrizmT[4] = 290 * mm;
  fPrizmT[5] = 290 * cos(fdTilt)* mm;

  double fMirror[3];
  fMirror[0] = 20 * mm;
  fMirror[1] = fPrizm[0];
  fMirror[2] = 1 * mm;
  //  fPrizm[0] = 170; fPrizm[1] = 300; fPrizm[2] = 50+300*tan(45*deg); fPrizm[3] = 50;

//  double fBar[3];
//  fBar[0] = 17 * mm;
//  fBar[1] = (fPrizm[0] - (fNBar - 1) * fBarsGap) / fNBar;
//  fBar[2] = 1050 * mm; // 4200; //4200

  double fMcpTotal[3];
  double fMcpActive[3];
  fMcpTotal[0] = fMcpTotal[1] = 53 + 4;
  fMcpTotal[2] = 1*mm;
  fMcpActive[0] = fMcpActive[1] = 53;
  fMcpActive[2] = 1*mm;

  double fLens[4];
  fLens[0] = fLens[1] = 40 * mm;
  fLens[2] = 10 * mm;
  double fRadius = (det_rin + det_rout)/2;

  double fBoxWidth = fPrizm[0];

  double fFd[3];
  fFd[0] = fBoxWidth;
  fFd[1] = fPrizm[2];
  fFd[2] = 1 * mm;

  fLens[0] = fPrizm[3];
  fLens[1] = fPrizm[0];
  fLens[2] = 12 * mm;

  // Getting box XML
  const int fNBoxes = xml_box_module.repeat();
  const double box_width = xml_box_module.width();
  const double box_height = xml_box_module.height();
  const double box_length = xml_box_module.length() + 550*mm;


  // The DIRC
  Assembly dirc_module("DIRCModule");

  //Volume lDirc("lDirc", gDirc, air);
  dirc_module.setVisAttributes(desc.visAttributes(xml_box_module.visStr()));


  // FD... whatever F and D is
  xml_comp_t xml_fd = xml_box_module.child(_Unicode(fd));
  Box gFd("gFd", xml_fd.height()/2, xml_fd.width()/2, xml_fd.thickness()/2);
  Volume lFd ("lFd", gFd, desc.material(xml_fd.materialStr()));
  lFd.setVisAttributes(desc.visAttributes(xml_fd.visStr()));
  //lFd.setSensitiveDetector(sens);


  // The Bar
  xml_comp_t xml_bar = xml_box_module.child(_Unicode(bar));
  double bar_height = xml_bar.height();
  double bar_width = xml_bar.width();
  double bar_length = xml_bar.length();
  Box gBar("gBar", bar_height/2, bar_width/2, bar_length/2);
  Volume lBar("lBar", gBar, desc.material(xml_bar.materialStr()));
  lBar.setVisAttributes(desc.visAttributes(xml_bar.visStr()));

  // Glue
  xml_comp_t xml_glue = xml_box_module.child(_Unicode(glue));
  double glue_thickness = xml_glue.thickness();  // 0.05 * mm;
  Box gGlue("gGlue", bar_height/2, bar_width/2, glue_thickness/2);
  Volume lGlue("lGlue", gGlue, desc.material(xml_glue.materialStr()));
  lGlue.setVisAttributes(desc.visAttributes(xml_glue.visStr()));

  sens.setType("tracker");
  lBar.setSensitiveDetector(sens);

  int bars_repeat_z = 4;    // TODO parametrize!
  double bar_assm_length = (bar_length + glue_thickness) * bars_repeat_z;
  int fNBar = xml_bar.repeat();
  double bar_gap = xml_bar.gap();
  for (int y_index = 0; y_index < fNBar; y_index++) {
    double shift_y = y_index * (bar_width + bar_gap) - 0.5 * box_width + 0.5 * bar_width;
    for (int z_index = 0; z_index < bars_repeat_z; z_index++) {
      double z = -0.5 * bar_assm_length + 0.5 * bar_length + (bar_length + glue_thickness) * z_index;
      auto placed_bar = dirc_module.placeVolume(lBar, Position(0, shift_y, z));
      dirc_module.placeVolume(lGlue, Position(0, shift_y, z + 0.5 * (bar_length + glue_thickness)));
      placed_bar.addPhysVolID("section", z_index);
      placed_bar.addPhysVolID("bar", y_index);
    }
  }


  // The Mirror
  xml_comp_t xml_mirror = xml_box_module.child(_Unicode(mirror));
  Box gMirror("gMirror", xml_mirror.height()/2, xml_mirror.width()/2, xml_mirror.thickness()/2);
  Volume lMirror("lMirror", gMirror, desc.material(xml_mirror.materialStr()));
  dirc_module.placeVolume(lMirror, Position(0, 0, -0.5 * (bar_assm_length - xml_mirror.thickness())));
  lMirror.setVisAttributes(desc.visAttributes(xml_mirror.visStr()));

  // The mirror optical surface
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  auto surf = surfMgr.opticalSurface("MirrorOpticalSurface");
  SkinSurface skin(desc, det, Form("dirc_mirror_optical_surface"), surf, lMirror);
  skin.isValid();

  // LENS
  // Lens volumes
  Volume lLens1;
  Volume lLens2;
  Volume lLens3;

  double lensMinThikness = 2.0 * mm;
  double layer12 = lensMinThikness * 2;

  // r1 = (r1==0)? 27.45: r1;
  // r2 = (r2==0)? 20.02: r2;

  double r1 = 33 * mm;
  double r2 = 24 * mm;
  double shight = 25 * mm;

  Position zTrans1(0, 0, -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - shight / 2. * shight / 2.) + lensMinThikness);
  Position zTrans2(0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - shight / 2. * shight / 2.) + layer12);

  Box gfbox("fbox", 0.5 * fLens[0], 0.5 * fLens[1], 0.5 * fLens[2]);
  Box gcbox("cbox", 0.5 * fLens[0], 0.5 * fLens[1] + 1*mm, 0.5 * fLens[2]);

//  Volume gfbox_volume("gfbox_volume", gfbox, bar_material);
//  lDirc.placeVolume(gfbox_volume, Position(0, 0, 0));
//
//  Volume gcbox_volume("gcbox_volume", gcbox, bar_material);
//  lDirc.placeVolume(gcbox_volume, Position(0, 0, 50));


  Position tTrans1(0.5 * (fLens[0] + shight), 0, -fLens[2] + layer12);
  Position tTrans0(-0.5 * (fLens[0] + shight), 0, -fLens[2] + layer12);
  SubtractionSolid tubox("tubox", gfbox, gcbox, tTrans1);
  SubtractionSolid gubox("gubox", tubox, gcbox, tTrans0);

//  Volume tubox_volume("tubox_volume", tubox, bar_material);
//  lDirc.placeVolume(tubox_volume, Position(0, 0, 100));
//
//  Volume gubox_volume("gubox_volume", gubox, bar_material);
//  lDirc.placeVolume(gubox_volume, Position(0, 0, 150));

  Tube gcylinder1("Cylinder1", 0, r1, 0.5 * fLens[1], 0 * deg, 360 * deg);
  Tube gcylinder2("Cylinder2", 0, r2, 0.5 * fLens[1] - 0.5*mm, 0 * deg, 360 * deg);
  Tube gcylinder1c("Cylinder1c", 0, r1, 0.5 * fLens[1] + 0.5*mm, 0 * deg, 360 * deg);
  Tube gcylinder2c("Cylinder2c", 0, r2, 0.5 * fLens[1] + 0.5*mm, 0 * deg, 360 * deg);
  RotationX xRot(-M_PI / 2.);

  IntersectionSolid gLens1("Lens1", gubox, gcylinder1, Transform3D(xRot, zTrans1));
  SubtractionSolid gLenst("temp", gubox, gcylinder1c, Transform3D(xRot, zTrans1));

//  Volume gLens1_volume("gLens1_volume", gLens1, bar_material);
//  lDirc.placeVolume(gLens1_volume, Position(0, 0, 200));
//
//  Volume gLenst_volume("gLenst_volume", gLenst, bar_material);
//  lDirc.placeVolume(gLenst_volume, Position(0, 0, 250));

  IntersectionSolid gLens2("Lens2", gLenst, gcylinder2, Transform3D(xRot, zTrans2));
  SubtractionSolid gLens3("Lens3", gLenst, gcylinder2c, Transform3D(xRot, zTrans2));

  lLens1 = Volume("lLens1", gLens1, bar_material);
  lLens2 = Volume("lLens2", gLens2, nlak33a);
  lLens3 = Volume("lLens3", gLens3, bar_material);

  lLens1.setVisAttributes(desc.visAttributes("DIRCLens1"));
  lLens2.setVisAttributes(desc.visAttributes("DIRCLens2"));
  lLens3.setVisAttributes(desc.visAttributes("DIRCLens3"));


  double shifth = 0.5 * (bar_assm_length + fLens[2]);
  // fmt::print("LENS HERE shifth={}\n", shifth);

  lLens1.setVisAttributes(desc.visAttributes("AnlTeal"));
  dirc_module.placeVolume(lLens1, Position(0, 0, shifth));
  dirc_module.placeVolume(lLens2, Position(0, 0, shifth));
  dirc_module.placeVolume(lLens3, Position(0, 0, shifth));
  
  // The Prizm
  Trap gPrizm = MakeTrap("gPrizm", fPrizm[0], fPrizm[1], fPrizm[2], fPrizm[3]);
  Volume lPrizm("lPrizm", gPrizm, bar_material);
  lPrizm.setVisAttributes(desc.visAttributes("DIRCPrism"));

  //G4RotationMatrix *fdRot = new G4RotationMatrix();
  //G4RotationMatrix *fdrot = new G4RotationMatrix();
  double evshiftz = 0.5 * bar_assm_length + fPrizm[1] + fMcpActive[2] / 2. + fLens[2];
  double evshiftx = -3*mm;

  double prism_shift_x = (fPrizm[2] + fPrizm[3]) / 4. - 0.5 * fPrizm[3] + 1.5*mm;
  double prism_shift_z = 0.5 * (bar_assm_length + fPrizm[1]) + fLens[2];

  Position  fPrismShift(prism_shift_x, 0, prism_shift_z);
  dirc_module.placeVolume(lPrizm, Transform3D(xRot, fPrismShift));
  dirc_module.placeVolume(lFd, Position(0.5 * fFd[1] - 0.5 * fPrizm[3] - evshiftx, 0, evshiftz));

  double dphi = 2 * M_PI / (double)fNBoxes;
  for (int i = 0; i < fNBoxes; i++) {
    double phi = dphi * i;
    double dx = -fRadius * cos(phi);
    double dy = -fRadius * sin(phi);

    //G4RotationMatrix *tRot = new G4RotationMatrix();

    Transform3D  tr(RotationZ(phi+M_PI), Position(dx, dy, 0));
    PlacedVolume box_placement = det_volume.placeVolume(dirc_module, tr);
    box_placement.addPhysVolID("module", i);


    // fmt::print("placing dircbox # {} -tphi={:.0f} dx={:.0f}, dy={:.0f}\n", i, phi/deg, dx/cm, dy/cm);

    //new G4PVPlacement(tRot, G4ThreeVector(dx, dy, 0), lDirc, "wDirc", lExpHall, false, i);
  }

    
  //////////////////
  // DIRC Bars
  //////////////////

  // double bar_radius   = 83.65 * cm;
  // double bar_length   = SizeZ;
  // double bar_width    = 42. * cm;
  // double bar_thicknes = 1.7 * cm;
  // int    bar_count    = 2 * M_PI * bar_radius / bar_width;
  // double bar_dphi     = 2 * 3.1415926 / bar_count;
  // Material bar_material = desc.material("Quartz");

  // Box    bar_geo(bar_thicknes / 2., bar_width / 2., bar_length / 2.);
  // Volume bar_volume("cb_DIRC_bars_Logix", bar_geo, bar_material);
  // bar_volume.setVisAttributes(desc.visAttributes(xml_det.visStr()));
  // sens.setType("tracker");
  // bar_volume.setSensitiveDetector(sens);

  // for (int ia = 0; ia < bar_count; ia++) {
  //   double phi = (ia * (bar_dphi));
  //   double x   = -bar_radius * cos(phi);
  //   double y   = -bar_radius * sin(phi);

  //   Transform3D  tr(RotationZ(bar_dphi * ia), Position(x, y, 0));
  //   PlacedVolume barPV = det_volume.placeVolume(bar_volume, tr);
  //   barPV.addPhysVolID("module", ia);
  // }
  return det;
}

DECLARE_DETELEMENT(cb_DIRC, createDetector)

dd4hep::Trap MakeTrap( const std::string& pName, double pZ, double pY, double pX, double pLTX )
{
    // Fixed Trap constructor. This function is a workaround of this bug:
    // https://github.com/AIDASoft/DD4hep/issues/850
    // Should be used instead of dd4hep::Trap(pName, pZ, pY, pX, pLTX) constructor

    double fDz  = 0.5*pZ;
    double fTthetaCphi = 0;
    double fTthetaSphi = 0;
    double fDy1 = 0.5*pY;
    double fDx1 = 0.5*pX;
    double fDx2 = 0.5*pLTX;
    double fTalpha1 = 0.5*(pLTX - pX)/pY;
    double fDy2 = fDy1;
    double fDx3 = fDx1;
    double fDx4 = fDx2;
    double fTalpha2 = fTalpha1;

    return Trap(pName, fDz,  fTthetaCphi,  fTthetaSphi,
                fDy1,  fDx1,  fDx2,  fTalpha1,
                fDy2,  fDx3,  fDx4,  fTalpha2);
}

