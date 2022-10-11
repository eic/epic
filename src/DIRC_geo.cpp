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
dd4hep::Trap MakeTrap(const std::string& pName, double pZ, double pY, double pX, double pLTX);

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t xml_det  = e;
  string    det_name = xml_det.nameStr();
  int       det_id   = xml_det.id();

  // Main detector xml element
  xml_dim_t dirc_dim = xml_det.dimensions();
  xml_dim_t dirc_pos = xml_det.position();
  xml_dim_t dirc_rot = xml_det.rotation();
  double    det_rin  = dirc_dim.rmin();
  double    det_rout = dirc_dim.rmax();

  Material quartz       = desc.material("Quartz");
  Material nlak33a      = desc.material("Nlak33a");
  auto&    bar_material = quartz;

  Assembly det_volume("DIRC");
  det_volume.setVisAttributes(desc.visAttributes(xml_det.visStr()));

  DetElement det(det_name, det_id);
  Volume     mother_vol = desc.pickMotherVolume(det);

  Transform3D  tr_global(RotationZYX(0, dirc_rot.theta(), 0.0), Position(0.0, 0.0, dirc_pos.z()));
  PlacedVolume det_plvol = mother_vol.placeVolume(det_volume, tr_global);

  det_plvol.addPhysVolID("system", det_id);
  det.setPlacement(det_plvol);

  // Parts dimensions

  // double fdTilt = 80 * deg;

  double fMcpTotal[3];
  double fMcpActive[3];
  fMcpTotal[0] = fMcpTotal[1] = 53 + 4;
  fMcpTotal[2]                = 1 * mm;
  fMcpActive[0] = fMcpActive[1] = 53;
  fMcpActive[2]                 = 1 * mm;

  double box_placement_radius = (det_rin + det_rout) / 2;

  // Getting box XML
  xml_comp_t   xml_box_module = xml_det.child(_U(module));
  const int    box_repeat     = xml_box_module.repeat();
  const double box_width      = xml_box_module.width();
  // FIXME unused box height and length
  // const double box_height = xml_box_module.height();
  // const double box_length = xml_box_module.length() + 550 * mm;

  // The DIRC
  Assembly dirc_module("DIRCModule");

  // Volume lDirc("lDirc", gDirc, air);
  dirc_module.setVisAttributes(desc.visAttributes(xml_box_module.visStr()));

  // FD... whatever F and D is
  xml_comp_t xml_fd = xml_box_module.child(_Unicode(fd));
  Box        gFd("gFd", xml_fd.height() / 2, xml_fd.width() / 2, xml_fd.thickness() / 2);
  Volume     lFd("lFd", gFd, desc.material(xml_fd.materialStr()));
  lFd.setVisAttributes(desc.visAttributes(xml_fd.visStr()));
  // lFd.setSensitiveDetector(sens);

  // Bar
  xml_comp_t xml_bar    = xml_box_module.child(_Unicode(bar));
  double     bar_height = xml_bar.height();
  double     bar_width  = xml_bar.width();
  double     bar_length = xml_bar.length();
  Box        gBar("gBar", bar_height / 2, bar_width / 2, bar_length / 2);
  Volume     lBar("lBar", gBar, desc.material(xml_bar.materialStr()));
  lBar.setVisAttributes(desc.visAttributes(xml_bar.visStr()));

  // Glue
  xml_comp_t xml_glue       = xml_box_module.child(_Unicode(glue));
  double     glue_thickness = xml_glue.thickness(); // 0.05 * mm;
  Box        gGlue("gGlue", bar_height / 2, bar_width / 2, glue_thickness / 2);
  Volume     lGlue("lGlue", gGlue, desc.material(xml_glue.materialStr()));
  lGlue.setVisAttributes(desc.visAttributes(xml_glue.visStr()));

  sens.setType("tracker");
  lBar.setSensitiveDetector(sens);

  int    bars_repeat_z   = 4; // TODO parametrize!
  double bar_assm_length = (bar_length + glue_thickness) * bars_repeat_z;
  int    fNBar           = xml_bar.repeat();
  double bar_gap         = xml_bar.gap();
  for (int y_index = 0; y_index < fNBar; y_index++) {
    double shift_y = y_index * (bar_width + bar_gap) - 0.5 * box_width + 0.5 * bar_width;
    for (int z_index = 0; z_index < bars_repeat_z; z_index++) {
      double z          = -0.5 * bar_assm_length + 0.5 * bar_length + (bar_length + glue_thickness) * z_index;
      auto   placed_bar = dirc_module.placeVolume(lBar, Position(0, shift_y, z));
      dirc_module.placeVolume(lGlue, Position(0, shift_y, z + 0.5 * (bar_length + glue_thickness)));
      placed_bar.addPhysVolID("section", z_index);
      placed_bar.addPhysVolID("bar", y_index);
    }
  }

  // Mirror construction
  xml_comp_t xml_mirror = xml_box_module.child(_Unicode(mirror));
  Box        gMirror("gMirror", xml_mirror.height() / 2, xml_mirror.width() / 2, xml_mirror.thickness() / 2);
  Volume     lMirror("lMirror", gMirror, desc.material(xml_mirror.materialStr()));
  dirc_module.placeVolume(lMirror, Position(0, 0, -0.5 * (bar_assm_length - xml_mirror.thickness())));
  lMirror.setVisAttributes(desc.visAttributes(xml_mirror.visStr()));

  // Mirror optical surface
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  auto                  surf    = surfMgr.opticalSurface("MirrorOpticalSurface");
  SkinSurface           skin(desc, det, Form("dirc_mirror_optical_surface"), surf, lMirror);
  skin.isValid();

  // Prism variables
  xml_comp_t xml_prism        = xml_box_module.child(_Unicode(prism));
  double     prism_angle      = xml_prism.angle();
  double     prism_width      = xml_prism.width();
  double     prism_length     = xml_prism.length();
  double     prism_short_edge = getAttrOrDefault(xml_prism, _Unicode(short_edge), 50 * mm);
  double     prism_long_edge  = prism_short_edge + prism_length * tan(prism_angle);

  // Lens variables
  xml_comp_t xml_lens    = xml_box_module.child(_Unicode(lens));
  double     lens_height = getAttrOrDefault(xml_lens, _Unicode(height), 50 * mm);
  // double lens_width = getAttrOrDefault(xml_lens, _Unicode(width), 25 * mm);
  double lens_thickness = getAttrOrDefault(xml_lens, _Unicode(thickness), 12 * mm);
  double lens_r1        = getAttrOrDefault(xml_lens, _Unicode(r1), 62 * mm);
  double lens_r2        = getAttrOrDefault(xml_lens, _Unicode(r2), 36 * mm);

  // Lens construction

  // lensMinThickness is distance from face to r1, and from r1 to r2 at lens top edge
  double lensMinThickness = 0.5 * mm;
  double layer01          = 1.0 * lensMinThickness;
  double layer12          = 2.0 * lensMinThickness;

  Position zTrans1(0, 0,
                   -lens_thickness / 2. - sqrt(lens_r1 * lens_r1 - lens_height / 2. * lens_height / 2.) + layer01);
  Position zTrans2(0, 0,
                   -lens_thickness / 2. - sqrt(lens_r2 * lens_r2 - lens_height / 2. * lens_height / 2.) + layer12);

  Box gfbox("fbox", 0.5 * prism_short_edge, 0.5 * prism_width, 0.5 * lens_thickness);
  Box gcbox("cbox", 0.5 * prism_short_edge, 0.5 * prism_width + 1 * mm, 0.5 * lens_thickness);

  //  Volume gfbox_volume("gfbox_volume", gfbox, bar_material);
  //  lDirc.placeVolume(gfbox_volume, Position(0, 0, 0));
  //
  //  Volume gcbox_volume("gcbox_volume", gcbox, bar_material);
  //  lDirc.placeVolume(gcbox_volume, Position(0, 0, 50));

  Position         tTrans1(0.5 * (prism_short_edge + lens_height), 0, -lens_thickness + layer12);
  Position         tTrans0(-0.5 * (prism_short_edge + lens_height), 0, -lens_thickness + layer12);
  SubtractionSolid tubox("tubox", gfbox, gcbox, tTrans1);
  SubtractionSolid gubox("gubox", tubox, gcbox, tTrans0);

  //  Volume tubox_volume("tubox_volume", tubox, bar_material);
  //  lDirc.placeVolume(tubox_volume, Position(0, 0, 100));
  //
  //  Volume gubox_volume("gubox_volume", gubox, bar_material);
  //  lDirc.placeVolume(gubox_volume, Position(0, 0, 150));

  Tube      gcylinder1("Cylinder1", 0, lens_r1, 0.5 * prism_width, 0 * deg, 360 * deg);
  Tube      gcylinder2("Cylinder2", 0, lens_r2, 0.5 * prism_width - 0.5 * mm, 0 * deg, 360 * deg);
  Tube      gcylinder1c("Cylinder1c", 0, lens_r1, 0.5 * prism_width + 0.5 * mm, 0 * deg, 360 * deg);
  Tube      gcylinder2c("Cylinder2c", 0, lens_r2, 0.5 * prism_width + 0.5 * mm, 0 * deg, 360 * deg);
  RotationX xRot(-M_PI / 2.);

  IntersectionSolid gLens1("Lens1", gubox, gcylinder1, Transform3D(xRot, zTrans1));
  SubtractionSolid  gLenst("temp", gubox, gcylinder1c, Transform3D(xRot, zTrans1));

  //  Volume gLens1_volume("gLens1_volume", gLens1, bar_material);
  //  lDirc.placeVolume(gLens1_volume, Position(0, 0, 200));
  //
  //  Volume gLenst_volume("gLenst_volume", gLenst, bar_material);
  //  lDirc.placeVolume(gLenst_volume, Position(0, 0, 250));

  IntersectionSolid gLens2("Lens2", gLenst, gcylinder2, Transform3D(xRot, zTrans2));
  SubtractionSolid  gLens3("Lens3", gLenst, gcylinder2c, Transform3D(xRot, zTrans2));

  Volume lLens1("lLens1", gLens1, bar_material);
  Volume lLens2("lLens2", gLens2, nlak33a);
  Volume lLens3("lLens3", gLens3, bar_material);

  lLens1.setVisAttributes(desc.visAttributes("DIRCLens1"));
  lLens2.setVisAttributes(desc.visAttributes("DIRCLens2"));
  lLens3.setVisAttributes(desc.visAttributes("DIRCLens3"));

  double shifth = 0.5 * (bar_assm_length + lens_thickness);

  lLens1.setVisAttributes(desc.visAttributes("AnlTeal"));
  dirc_module.placeVolume(lLens1, Position(0, 0, shifth));
  dirc_module.placeVolume(lLens2, Position(0, 0, shifth));
  dirc_module.placeVolume(lLens3, Position(0, 0, shifth));

  // Prism construction
  Trap   gPrism = MakeTrap("gPrism", prism_width, prism_length, prism_long_edge, prism_short_edge);
  Volume lPrism("lPrism", gPrism, bar_material);
  lPrism.setVisAttributes(desc.visAttributes("DIRCPrism"));

  double evshiftz = 0.5 * bar_assm_length + prism_length + fMcpActive[2] / 2. + lens_thickness;
  double evshiftx = -3 * mm;

  double prism_shift_x = (prism_long_edge + prism_short_edge) / 4. - 0.5 * prism_short_edge + 1.5 * mm;
  double prism_shift_z = 0.5 * (bar_assm_length + prism_length) + lens_thickness;

  Position fPrismShift(prism_shift_x, 0, prism_shift_z);
  dirc_module.placeVolume(lPrism, Transform3D(xRot, fPrismShift));
  dirc_module.placeVolume(lFd, Position(0.5 * prism_long_edge - 0.5 * prism_short_edge - evshiftx, 0, evshiftz));

  double dphi = 2 * M_PI / (double)box_repeat;
  for (int i = 0; i < box_repeat; i++) {
    double phi = dphi * i;
    double dx  = -box_placement_radius * cos(phi);
    double dy  = -box_placement_radius * sin(phi);

    // G4RotationMatrix *tRot = new G4RotationMatrix();

    Transform3D  tr(RotationZ(phi + M_PI), Position(dx, dy, 0));
    PlacedVolume box_placement = det_volume.placeVolume(dirc_module, tr);
    box_placement.addPhysVolID("module", i);

    // fmt::print("placing dircbox # {} -tphi={:.0f} dx={:.0f}, dy={:.0f}\n", i, phi/deg, dx/cm, dy/cm);

    // new G4PVPlacement(tRot, G4ThreeVector(dx, dy, 0), lDirc, "wDirc", lExpHall, false, i);
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

#ifdef EPIC_ECCE_LEGACY_COMPAT
DECLARE_DETELEMENT(ecce_DIRC, createDetector)
#endif
DECLARE_DETELEMENT(epic_DIRC, createDetector)

dd4hep::Trap MakeTrap(const std::string& pName, double pZ, double pY, double pX, double pLTX)
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
  double fTalpha1    = 0.5 * (pLTX - pX) / pY;
  double fDy2        = fDy1;
  double fDx3        = fDx1;
  double fDx4        = fDx2;
  double fTalpha2    = fTalpha1;

  return Trap(pName, fDz, fTthetaCphi, fTthetaSphi, fDy1, fDx1, fDx2, fTalpha1, fDy2, fDx3, fDx4, fTalpha2);
}
