//==========================================================================
//  Implementation of longitudinally separated forward calorimeter
//--------------------------------------------------------------------------
//  Author: Friederike Bock (ORNL)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <XML/Helper.h>
using namespace dd4hep;

//************************************************************************************************************
//************************** create module assembly  ******************************************************
//************************************************************************************************************
Volume createModule(Detector& desc, int modID_x, int modID_y, double length, SensitiveDetector sens,
                    int doSamplingFractionStudies, bool isLargeMod, bool useMoreLayers)
{
  double modBox_face_width      = 1.5*cm;
  double modBox_length_tot      = length - 2*modBox_face_width;
  double modBox_length          = length - 10 * cm;
  double modBox_steel_length    = 120 * cm;
  double modBox_tungsten_length = 10 * cm;
  double modBox_width           = 20 * cm;
  if (isLargeMod)
    modBox_width = 30 * cm;
  double      modBox_height = 10 * cm;
  std::string modShrtName   = "8M";
  if (isLargeMod)
    modShrtName = "12M";
  std::string  baseName = "GFHCAL_" + modShrtName + _toString(modID_x, "_modx_%d") + _toString(modID_y, "_mody_%d");
  PlacedVolume pvm;

  int visDetails = 0;
  visDetails     = 1;
  // visDetails = 2;

  // std::cout << "create 12M module" << std::endl;
  // std::cout << "length = " << length << std::endl;
  Box    modBox(modBox_width / 2., modBox_height / 2., length / 2.);
  Volume vol_modBox(baseName, modBox, desc.material("Air"));
  if (visDetails) {
    vol_modBox.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
  } else {
    vol_modBox.setVisAttributes(desc.visAttributes("AnlOrange"));
  }

  double miniframe_thickness = 1.0 * mm;
  double miniframe_height    = modBox_height - 2 * miniframe_thickness;

  double tyvek_thickness   = 0.15 * mm;
  double sciSeg_length_tot = 10 * cm;
  double sciSeg_width_tot  = 4 * mm + tyvek_thickness * 2;
  // double sciSeg_height_tot = 5 * cm;

  double sciSeg_length = sciSeg_length_tot - tyvek_thickness * 2;
  double sciSeg_width  = 4 * mm;
  double sciSeg_height = 5 * cm - tyvek_thickness - miniframe_thickness;

  double modBox_sidewall_thickness = 2 * mm;
  // double modBox_topwall_thickness = 0.5* mm;
  double absorber_thickness = 16.8 * mm;
  if(useMoreLayers) {
    absorber_thickness = 12.1 * mm;
  }

  double pcb_gap       = 1.19 * mm;
  double pcb_thickness = 1.18 * mm;
  double pcb_height    = modBox_height - 3 * miniframe_thickness;

  double fullLayer_width = absorber_thickness + pcb_gap + sciSeg_width_tot + miniframe_thickness * 2;
  // std::cout << "gfhcal fullLayer_width = " << fullLayer_width << std::endl;
  double miniframe_width = pcb_gap + sciSeg_width_tot + miniframe_thickness * 2;

  Box scintBox(sciSeg_width / 2., sciSeg_height / 2., sciSeg_length / 2.);
  // Volume vol_scintBox(baseName + _toString(ilx, "_ix_%d") + _toString(ily, "_iy_%d") + _toString(ilz, "_iz_%d") ,
  // scintBox, desc.material("Polystyrene"));
  Volume vol_scintBox(baseName + "_Scintillator" + modShrtName, scintBox, desc.material("Polystyrene"));
  if (visDetails == 2) {
    vol_scintBox.setVisAttributes(desc.visAttributes("GFHCALLayerScintVis"));
  } else {
    vol_scintBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelWallBox(modBox_sidewall_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_steelWallBox(baseName + "_FeWall" + modShrtName, steelWallBox, desc.material("Steel235"));
  if (visDetails) {
    vol_steelWallBox.setVisAttributes(desc.visAttributes("AnlOrange"));
  } else {
    vol_steelWallBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelBoxFace(modBox_width / 2., modBox_height / 2., modBox_face_width / 2.);
  Volume vol_steelBoxFace(baseName + "_FeAbsorber" + modShrtName, steelBoxFace, desc.material("Steel235"));
  if (visDetails) {
    vol_steelBoxFace.setVisAttributes(desc.visAttributes("AnlOrange"));
  } else {
    vol_steelBoxFace.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelAbsorberBox(absorber_thickness / 2., modBox_height / 2., modBox_steel_length / 2.);
  Volume vol_steelAbsorberBox(baseName + "_FeAbsorber" + modShrtName, steelAbsorberBox, desc.material("Steel235"));
  if (doSamplingFractionStudies == 3)
    vol_steelAbsorberBox.setMaterial(desc.material("Tungsten"));
  if (visDetails) {
    vol_steelAbsorberBox.setVisAttributes(desc.visAttributes("AnlLight_Gray"));
  } else {
    vol_steelAbsorberBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    tungstenAbsorberBox(absorber_thickness / 2., modBox_height / 2., modBox_tungsten_length / 2.);
  Volume vol_tungstenAbsorberBox(baseName + "_WAbsorber" + modShrtName, tungstenAbsorberBox, desc.material("Tungsten"));
  if (doSamplingFractionStudies == 2)
    vol_tungstenAbsorberBox.setMaterial(desc.material("Steel235"));
  if (visDetails) {
    if (isLargeMod) {
      vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("AnlBlue"));
    } else {
      vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("AnlViolet"));
    }

  } else {
    vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelMiniFrameBox(miniframe_thickness / 2., miniframe_height / 2., modBox_length_tot / 2.);
  Volume vol_steelMiniFrameBox(baseName + "_MiniFrame" + modShrtName, steelMiniFrameBox, desc.material("Steel235"));
  if (visDetails == 2) {
    vol_steelMiniFrameBox.setVisAttributes(desc.visAttributes("AnlDarkRed"));
  } else {
    vol_steelMiniFrameBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelMiniFrameBox_cover(miniframe_width / 2., miniframe_thickness / 2., modBox_length_tot / 2.);
  Volume vol_steelMiniFrameBox_cover(baseName + "_MiniFrameCover" + modShrtName, steelMiniFrameBox_cover,
                                     desc.material("Steel235"));
  if (visDetails == 2) {
    // vol_steelMiniFrameBox_cover.setVisAttributes(desc.visAttributes("AnlDarkRed"));
  } else {
    vol_steelMiniFrameBox_cover.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    pcbBox(pcb_thickness / 2., pcb_height / 2., modBox_length_tot / 2.);
  Volume vol_pcbBox(baseName + "_PCB" + modShrtName, pcbBox, desc.material("Fr4"));
  if (visDetails == 2) {
    vol_pcbBox.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  } else {
    vol_pcbBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    tyvekBox(tyvek_thickness / 2., miniframe_height / 2., modBox_length / 2.);
  Volume vol_tyvekBox(baseName + "_Tyvek" + modShrtName, tyvekBox, desc.material("Tyvek"));
  if (visDetails == 2) {
    vol_tyvekBox.setVisAttributes(desc.visAttributes("AnlGold"));
  } else {
    vol_tyvekBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  // std::cout << "fullLayer_width 8M = " << fullLayer_width << std::endl;
  int nLayers_x = (int)((modBox_width - 2 * modBox_sidewall_thickness) / fullLayer_width);
  int nLayers_z = (int)(modBox_length / sciSeg_length_tot);
  // int nLayers_y = 2;

  double addAbsorber_thickness = (modBox_width - 2 * modBox_sidewall_thickness) - nLayers_x * fullLayer_width;
  Box    steelAbsorberBoxAdd(addAbsorber_thickness / 2., modBox_height / 2., modBox_steel_length / 2.);
  Volume vol_steelAbsorberBoxAdd(baseName + "_FeAbsorberAdd" + modShrtName, steelAbsorberBoxAdd,
                                 desc.material("Steel235"));
  if (visDetails) {
    vol_steelAbsorberBoxAdd.setVisAttributes(desc.visAttributes("AnlBrown"));
  } else {
    vol_steelAbsorberBoxAdd.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }
  Box    tungstenAbsorberBoxAdd(addAbsorber_thickness / 2., modBox_height / 2., modBox_tungsten_length / 2.);
  Volume vol_tungstenAbsorberBoxAdd(baseName + "_WAbsorberAdd" + modShrtName, tungstenAbsorberBoxAdd,
                                    desc.material("Tungsten"));
  if (visDetails) {
    vol_tungstenAbsorberBoxAdd.setVisAttributes(desc.visAttributes("AnlTeal"));
  } else {
    vol_tungstenAbsorberBoxAdd.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }
  vol_modBox.placeVolume(
      vol_steelBoxFace,
      Transform3D(RotationZYX(0, 0, 0), Position(0, 0, -length/2 + modBox_face_width / 2.)));
  vol_modBox.placeVolume(
      vol_steelBoxFace,
      Transform3D(RotationZYX(0, 0, 0), Position(0, 0, length/2 - modBox_face_width / 2.)));

  vol_modBox.placeVolume(
      vol_steelWallBox,
      Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness / 2, 0, 0)));
  vol_modBox.placeVolume(
      vol_steelWallBox,
      Transform3D(RotationZYX(0, 0, 0), Position(modBox_width / 2 - modBox_sidewall_thickness / 2, 0, 0)));

  vol_modBox.placeVolume(
      vol_steelAbsorberBoxAdd,
      Transform3D(RotationZYX(0, 0, 0),
                  Position(modBox_width / 2 - modBox_sidewall_thickness - addAbsorber_thickness / 2, 0, modBox_face_width)));
  vol_modBox.placeVolume(vol_tungstenAbsorberBoxAdd,
                         Transform3D(RotationZYX(0, 0, 0),
                                     Position(modBox_width / 2 - modBox_sidewall_thickness - addAbsorber_thickness / 2,
                                              0, -length / 2 + modBox_face_width + modBox_tungsten_length / 2)));

  for (int ilx = 0; ilx < nLayers_x; ilx++) {

    if (doSamplingFractionStudies) {
      sens.setType("calorimeter");
      vol_steelAbsorberBox.setSensitiveDetector(sens);
    }
    pvm = vol_modBox.placeVolume(
        vol_steelAbsorberBox,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness +
                                                       absorber_thickness / 2 + ilx * fullLayer_width,
                                                   0, modBox_face_width)));
    if (doSamplingFractionStudies) {
      pvm.addPhysVolID("modulex", modID_x)
          .addPhysVolID("moduley", modID_y)
          .addPhysVolID("layerx", ilx)
          .addPhysVolID("passive", 1);
    }

    if (doSamplingFractionStudies) {
      sens.setType("calorimeter");
      vol_tungstenAbsorberBox.setSensitiveDetector(sens);
    }
    pvm = vol_modBox.placeVolume(
        vol_tungstenAbsorberBox,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness +
                                                       absorber_thickness / 2 + ilx * fullLayer_width,
                                                   0, -modBox_length_tot / 2 + modBox_tungsten_length / 2)));
    if (doSamplingFractionStudies) {
      pvm.addPhysVolID("modulex", modID_x)
          .addPhysVolID("moduley", modID_y)
          .addPhysVolID("layerx", ilx)
          .addPhysVolID("passive", 2);
    }

    vol_modBox.placeVolume(
        vol_pcbBox,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness + absorber_thickness +
                                                       miniframe_thickness + pcb_gap / 2 + ilx * fullLayer_width,
                                                   0, 0)));

    if (doSamplingFractionStudies) {
      sens.setType("calorimeter");
      vol_steelMiniFrameBox.setSensitiveDetector(sens);
    }
    pvm = vol_modBox.placeVolume(
        vol_steelMiniFrameBox,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness + absorber_thickness +
                                                       miniframe_thickness / 2 + ilx * fullLayer_width,
                                                   0, 0)));

    if (doSamplingFractionStudies) {
      pvm.addPhysVolID("modulex", modID_x)
          .addPhysVolID("moduley", modID_y)
          .addPhysVolID("layerx", ilx)
          .addPhysVolID("layery", 0)
          .addPhysVolID("passive", 3);
    }
    pvm = vol_modBox.placeVolume(
        vol_steelMiniFrameBox,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness + fullLayer_width -
                                                       miniframe_thickness / 2 + ilx * fullLayer_width,
                                                   0, 0)));

    if (doSamplingFractionStudies) {
      pvm.addPhysVolID("modulex", modID_x)
          .addPhysVolID("moduley", modID_y)
          .addPhysVolID("layerx", ilx)
          .addPhysVolID("layery", 1)
          .addPhysVolID("passive", 3);
    }

    if (doSamplingFractionStudies) {
      sens.setType("calorimeter");
      vol_steelMiniFrameBox_cover.setSensitiveDetector(sens);
    }
    pvm = vol_modBox.placeVolume(
        vol_steelMiniFrameBox_cover,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness + absorber_thickness +
                                                       miniframe_width / 2 + ilx * fullLayer_width,
                                                   miniframe_height/2 + miniframe_thickness / 2, 0)));

    if (doSamplingFractionStudies) {
      pvm.addPhysVolID("modulex", modID_x)
          .addPhysVolID("moduley", modID_y)
          .addPhysVolID("layerx", ilx)
          .addPhysVolID("layerz", 2)
          .addPhysVolID("passive", 3);
    }
    pvm = vol_modBox.placeVolume(
        vol_steelMiniFrameBox_cover,
        Transform3D(RotationZYX(0, 0, 0), Position(-modBox_width / 2 + modBox_sidewall_thickness + absorber_thickness +
                                                       miniframe_width / 2 + ilx * fullLayer_width,
                                                   -miniframe_height/2 - miniframe_thickness / 2, 0)));

    if (doSamplingFractionStudies) {
      pvm.addPhysVolID("modulex", modID_x)
          .addPhysVolID("moduley", modID_y)
          .addPhysVolID("layerx", ilx)
          .addPhysVolID("layerz", 3)
          .addPhysVolID("passive", 3);
    }
    vol_modBox.placeVolume(vol_tyvekBox, Transform3D(RotationZYX(0, 0, 0),
                                                     Position(-modBox_width / 2 + modBox_sidewall_thickness +
                                                                  absorber_thickness + miniframe_thickness + pcb_gap +
                                                                  tyvek_thickness / 2 + ilx * fullLayer_width,
                                                            0, -modBox_tungsten_length / 2 + modBox_face_width)));
    vol_modBox.placeVolume(vol_tyvekBox,
                           Transform3D(RotationZYX(0, 0, 0),
                                       Position(-modBox_width / 2 + modBox_sidewall_thickness +
                                                                  absorber_thickness + miniframe_thickness + pcb_gap +
                                                                  tyvek_thickness + sciSeg_width + tyvek_thickness/2 + ilx * fullLayer_width,
                                                0, -modBox_tungsten_length / 2 + modBox_face_width)));
    // for (int ily = 0; ily < nLayers_y; ily++) {
    for (int ilz = 0; ilz < nLayers_z; ilz++) {
      sens.setType("calorimeter");
      vol_scintBox.setSensitiveDetector(sens);

      pvm = vol_modBox.placeVolume(
          vol_scintBox,
          Transform3D(
              RotationZYX(0, 0, 0),
              Position(-modBox_width / 2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness +
                            pcb_gap + sciSeg_width_tot / 2 + ilx * fullLayer_width,
                        -miniframe_height / 2 + sciSeg_height / 2 ,
                        -modBox_length_tot / 2 + sciSeg_length_tot / 2 + ilz * sciSeg_length_tot)));
      pvm.addPhysVolID("modulex", modID_x);
      pvm.addPhysVolID("moduley", modID_y);
      pvm.addPhysVolID("layerx", ilx);
      pvm.addPhysVolID("layery", 0);
      pvm.addPhysVolID("layerz", ilz);
      pvm.addPhysVolID("passive", 0);

      sens.setType("calorimeter");
      vol_scintBox.setSensitiveDetector(sens);

      pvm = vol_modBox.placeVolume(
          vol_scintBox,
          Transform3D(
              RotationZYX(0, 0, 0),
              Position(-modBox_width / 2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness +
                            pcb_gap + sciSeg_width_tot / 2 + ilx * fullLayer_width,
                        miniframe_height / 2 - sciSeg_height / 2 ,
                        -modBox_length_tot / 2 + sciSeg_length_tot / 2 + ilz * sciSeg_length_tot)));
      pvm.addPhysVolID("modulex", modID_x);
      pvm.addPhysVolID("moduley", modID_y);
      pvm.addPhysVolID("layerx", ilx);
      pvm.addPhysVolID("layery", 1);
      pvm.addPhysVolID("layerz", ilz);
      pvm.addPhysVolID("passive", 0);
    }
    // }
  }

  return vol_modBox;
}

//********************************************************************************************
//*                                                                                          *
//*                              Create detector                                             *
//==============================  MAIN FUNCTION  =============================================
//*                                                                                          *
//********************************************************************************************
static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{

  // global detector variables
  xml_det_t   x_det    = handle;
  int         det_id   = x_det.id();
  std::string det_name = x_det.nameStr();

  int doSamplingFractionStudies = getAttrOrDefault(x_det, _Unicode(doSamplingFractionStudies), 0);
  std::cout << "doSamplingFractionStudies = " << doSamplingFractionStudies << std::endl;
  int useMoreLayers = getAttrOrDefault(x_det, _Unicode(useMoreLayers), 0);
  std::cout << "useMoreLayers = " << useMoreLayers << std::endl;

  xml_dim_t dim    = x_det.dimensions();
  double    length = dim.z(); // Size along z-axis
  xml_dim_t pos    = x_det.position();

  std::cout << "global GFHCAL position" << pos.x() << "\t" << pos.y() << "\t" << pos.z() << std::endl;

  DetElement sdet(det_name, det_id);

  Assembly     assembly(det_name);
  PlacedVolume phv;

  int moduleIDx = 0;
  int moduleIDy = 0;

  // phv = motherVol.placeVolume(eightMassembly, Position(0, 0, 0));

  // create 8M modules
  std::vector<double> xpos8M;
  std::vector<double> ypos8M;
  // std::vector<double> zpos8M;

  for (xml_coll_t i(handle, _Unicode(eightmodulepositions)); i; ++i) {
    xml_comp_t x_mtrx = i;

    std::string mtrx_name   = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string mtrx_values = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double>* aptr = NULL;

    if (mtrx_name == "xpos")
      aptr = &xpos8M;
    else if (mtrx_name == "ypos")
      aptr = &ypos8M;
    // else if (mtrx_name == "zpos")
    //   aptr = &zpos8M;
    else {
      printout(WARNING, "GFHCAL", "unknown <eightmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t      posC      = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  if (xpos8M.size() != ypos8M.size()) {
    std::cout << xpos8M.size() << "\t" << ypos8M.size() << std::endl;
    std::cout << "idiot you can't count" << std::endl;
  } else {
    for (int e = 0; e < (int)xpos8M.size(); e++) {
      if (e % 20 == 0)
        std::cout << "GFHCAL placing 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << xpos8M[e] << "\t"
                  << ypos8M[e] << std::endl;
      moduleIDx = ((xpos8M[e] + 270) / 10);
      moduleIDy = ((ypos8M[e] + 265) / 10);
      if (moduleIDx < 0 || moduleIDy < 0) {
        // std::cout << "GFHCAL placing 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << xpos8M[e] << "\t"
        //           << ypos8M[e] << "\t" << zpos8M[e] << std::endl;
        std::cout << "GFHCAL WRONG ID FOR 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << moduleIDx << "\t"
                  << moduleIDy << std::endl;
      }
      Volume eightMassembly = createModule(desc, moduleIDx, moduleIDy, length, sens, doSamplingFractionStudies, false, useMoreLayers);
      // Volume eightMassembly = createEightMModule(desc, moduleIDx, moduleIDy, length, sens,
      // doSamplingFractionStudies);

      auto tr8M = Transform3D(Position(pos.x() - xpos8M[e] * dd4hep::cm - 0.5 * 20 * cm,
                                       pos.y() - ypos8M[e] * dd4hep::cm, pos.z() + length / 2.));
      phv       = assembly.placeVolume(eightMassembly, tr8M);
      phv.addPhysVolID("modulex", moduleIDx);
      phv.addPhysVolID("moduley", moduleIDy);
      phv.addPhysVolID("passive", 0);
    }
  }

  // create 12M modules
  std::vector<double> xpos12M;
  std::vector<double> ypos12M;
  // std::vector<double> zpos12M;

  for (xml_coll_t i(handle, _Unicode(twelvemodulepositions)); i; ++i) {
    xml_comp_t x_mtrx = i;

    std::string mtrx_name   = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string mtrx_values = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double>* aptr = NULL;

    if (mtrx_name == "xpos")
      aptr = &xpos12M;
    else if (mtrx_name == "ypos")
      aptr = &ypos12M;
    // else if (mtrx_name == "zpos")
    //   aptr = &zpos12M;
    else {
      printout(WARNING, "GFHCAL", "unknown <TwelveModulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t      posC      = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  if (xpos12M.size() != ypos12M.size()) {
    std::cout << xpos12M.size() << "\t" << ypos12M.size() << std::endl;
    std::cout << "idiot you can't count" << std::endl;
  } else {
    for (int e = 0; e < (int)xpos12M.size(); e++) {
      if (e % 20 == 0)
        std::cout << "GFHCAL placing 12M module: " << e << "/" << (int)xpos12M.size() << "\t" << xpos12M[e] << "\t"
                  << ypos12M[e] << std::endl;

      moduleIDx = ((xpos12M[e] + 275) / 10);
      moduleIDy = ((ypos12M[e] + 265) / 10);
      if (moduleIDx < 0 || moduleIDy < 0) {
        // std::cout << "GFHCAL placing 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << xpos8M[e] << "\t"
        //           << ypos8M[e] << "\t" << zpos8M[e] << std::endl;
        std::cout << "GFHCAL WRONG ID FOR 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << moduleIDx << "\t"
                  << moduleIDy << std::endl;
      }
      Volume eightMassembly = createModule(desc, moduleIDx, moduleIDy, length, sens, doSamplingFractionStudies, true, useMoreLayers);
      // Volume eightMassembly = createTwelveMModule(desc, moduleIDx, moduleIDy, length, sens,
      // doSamplingFractionStudies);

      auto tr12M = Transform3D(Position(pos.x() - xpos12M[e] * dd4hep::cm - 10 * cm, pos.y() - ypos12M[e] * dd4hep::cm,
                                        pos.z() + length / 2.));
      phv        = assembly.placeVolume(eightMassembly, tr12M);
      phv.addPhysVolID("modulex", moduleIDx);
      phv.addPhysVolID("moduley", moduleIDy);
      phv.addPhysVolID("passive", 0);
    }
  }

  Volume motherVol = desc.pickMotherVolume(sdet);
  // std::cout << motherVol.name() << std::endl;

  Transform3D  tr    = Translation3D(0., 0., 0.) * RotationZYX(0., 0., 0.);
  PlacedVolume envPV = motherVol.placeVolume(assembly, tr);
  envPV.addPhysVolID("system", det_id);
  sdet.setPlacement(envPV);

  return sdet;
}
DECLARE_DETELEMENT(epic_GFHCAL, createDetector)
