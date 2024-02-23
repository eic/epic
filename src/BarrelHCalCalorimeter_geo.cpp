// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 John Lajoie

//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"

#include "TVector3.h"
#include "TGDMLParse.h"
#include "FileLoaderHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{

  // printout(WARNING, "BarrelHCalCalorimeter", "called create_detector ");

  xml_det_t x_det    = e;
  int       det_id   = x_det.id();
  string    det_name = x_det.nameStr();
  Material  air      = description.air();

  DetElement sdet(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(sdet);

  // Create envelope to hold HCAL barrel

  double rmin1   = x_det.rmin1();
  double rmin2   = x_det.rmin2();
  double rmax    = x_det.rmax();
  double length1 = x_det.z1();
  double length2 = x_det.z2();

  // 5mm buffer on inner radius to acount for comb angle/tilt
  // (allows the combs to be added in later)
  std::vector<double> rmins = {rmin2-0.5*cm, rmin2-0.5*cm, rmin1-1.5*cm, rmin1-1.5*cm, rmin2-0.5*cm, rmin2-0.5*cm};
  // 1cm buffer on outer radius to acount for comb angle/tilt
  // (allows the combs to be added in later)
  std::vector<double> rmaxs = {rmax+1.0*cm, rmax+1.0*cm, rmax+1.0*cm, rmax+1.0*cm, rmax+1.0*cm, rmax+1.0*cm};
  // leave room for dogbones at the ends (not part of det table)
  std::vector<double> zs    = {-length2 / 2. - 8.2*cm, -length1 / 2., -length1 / 2., length1 / 2., length1 / 2., length2 / 2. + 8.2*cm};

  // printout(WARNING, "BarrelHCalCalorimeter", "%f %f %f %f %f", rmin1, rmin2, rmax, length1/2., length2/2.);

  Polycone ptube(0.0, 2.0 * M_PI, rmins, rmaxs, zs);
  Volume envelope(det_name, ptube, air);

  PlacedVolume env_phv = motherVol.placeVolume(envelope);
  env_phv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_phv);

  xml_comp_t det_define = x_det.child("define");

  // Pick up the constants

  double tilePlaneRotate  = 0.0;
  double tile_tolerance = 0.2; // Tile tolerance in mm to avoid overlaps

  // Sector steel tessellated shape gdml file info
  xml_comp_t x_det_sec_gdmlfile = x_det.child("sec_gdmlfile");
  std::string sec_gdml_file = getAttrOrDefault<std::string>(x_det_sec_gdmlfile, _Unicode(file), " ");;
  std::string sec_gdml_material = getAttrOrDefault<std::string>(x_det_sec_gdmlfile, _Unicode(material), " ");
  std::string sec_gdml_url = getAttrOrDefault<std::string>(x_det_sec_gdmlfile, _Unicode(url), " ");
  std::string sec_gdml_cache = getAttrOrDefault<std::string>(x_det_sec_gdmlfile, _Unicode(cache), " ");

  xml_comp_t x_det_csec_gdmlfile = x_det.child("csec_gdmlfile");
  std::string csec_gdml_file = getAttrOrDefault<std::string>(x_det_csec_gdmlfile, _Unicode(file), " ");;
  std::string csec_gdml_material = getAttrOrDefault<std::string>(x_det_csec_gdmlfile, _Unicode(material), " ");
  std::string csec_gdml_url = getAttrOrDefault<std::string>(x_det_csec_gdmlfile, _Unicode(url), " ");
  std::string csec_gdml_cache = getAttrOrDefault<std::string>(x_det_csec_gdmlfile, _Unicode(cache), " ");

  xml_comp_t x_det_er_gdmlfile = x_det.child("er_gdmlfile");
  std::string er_gdml_file = getAttrOrDefault<std::string>(x_det_er_gdmlfile, _Unicode(file), " ");;
  std::string er_gdml_material = getAttrOrDefault<std::string>(x_det_er_gdmlfile, _Unicode(material), " ");
  std::string er_gdml_url = getAttrOrDefault<std::string>(x_det_er_gdmlfile, _Unicode(url), " ");
  std::string er_gdml_cache = getAttrOrDefault<std::string>(x_det_er_gdmlfile, _Unicode(cache), " ");

  // Loop over the defines section and pick up the tile offsets, ref location and angles

  for (xml_coll_t i(det_define, _Unicode(constant)); i; ++i) {
    xml_comp_t x_const = i;

    std::string const_name  = getAttrOrDefault<std::string>(x_const, _Unicode(name), " ");

    if (const_name == "tilePlaneRotate"){
      std::string const_value = getAttrOrDefault<std::string>(x_const, _Unicode(value), " ");
      tilePlaneRotate = atof(const_value.c_str());
    }
    else
      printout(WARNING, "BarrelHCalCalorimeter", "unrecognized <constant> data!");
  }

  std::vector<double> xposOuter;
  std::vector<double> yposOuter;

  std::vector<double> xposTileS;
  std::vector<double> yposTileS;
  std::vector<double> zposTileS;

  std::vector<double> xposTileN;
  std::vector<double> yposTileN;
  std::vector<double> zposTileN;

  std::vector<double> xposChimneyTileS;
  std::vector<double> yposChimneyTileS;
  std::vector<double> zposChimneyTileS;

  for (xml_coll_t i(det_define, _Unicode(matrix)); i; ++i) {
    xml_comp_t x_mtrx = i;

    std::string mtrx_name   = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string mtrx_values = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double>* aptr = NULL;

    if (mtrx_name == "xposOuter")
      aptr = &xposOuter;
    else if (mtrx_name == "yposOuter")
      aptr = &yposOuter;
    else if (mtrx_name == "xposTileS")
      aptr = &xposTileS;
    else if (mtrx_name == "yposTileS")
      aptr = &yposTileS;
    else if (mtrx_name == "zposTileS")
      aptr = &zposTileS;
    else if (mtrx_name == "xposTileN")
      aptr = &xposTileN;
    else if (mtrx_name == "yposTileN")
      aptr = &yposTileN;
    else if (mtrx_name == "zposTileN")
      aptr = &zposTileN;
    else if (mtrx_name == "xposChimneyTileS")
      aptr = &xposChimneyTileS;
    else if (mtrx_name == "yposChimneyTileS")
      aptr = &yposChimneyTileS;
    else if (mtrx_name == "zposChimneyTileS")
      aptr = &zposChimneyTileS;
    else {
      printout(WARNING, "BarrelHCalCalorimeter", "unknown <matrix> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t      pos       = 0;
    std::string token;
    while ((pos = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, pos);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, pos + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  // Read in the barrel structure GDML file
  // three structures - normal sector, chimney sector, and end rings
  Assembly BarrelHCAL("BarrelHCAL");
  TGDMLParse parser;

  // sector
  EnsureFileFromURLExists(sec_gdml_url, sec_gdml_file, sec_gdml_cache);
  if (!fs::exists(fs::path(sec_gdml_file))) {
    printout(ERROR, "BarrelHCalCalorimeter_geo", "file " + sec_gdml_file + " does not exist");
    printout(ERROR, "BarrelHCalCalorimeter_geo", "use a FileLoader plugin before the field element");
    std::_Exit(EXIT_FAILURE);
  }

  Volume barrel_sector_vol = parser.GDMLReadFile(sec_gdml_file.c_str());
  if(!barrel_sector_vol.isValid()){
    printout(WARNING, "BarrelHCalCalorimeter", "%s", sec_gdml_file.c_str());
    printout(WARNING, "BarrelHCalCalorimeter", "barrel_sector_vol invalid, GDML parser failed!");
    std::_Exit(EXIT_FAILURE);
  }
  barrel_sector_vol.import();
  barrel_sector_vol.setVisAttributes(description, x_det.visStr());
  TessellatedSolid barrel_sector_solid = barrel_sector_vol.solid();
  barrel_sector_solid->CloseShape(true, true, true); // tesselated solid not closed by import!
  Material sector_material = description.material(sec_gdml_material.c_str());
  barrel_sector_vol.setMaterial(sector_material);

  // chimney sector
  EnsureFileFromURLExists(csec_gdml_url, csec_gdml_file, csec_gdml_cache);
  if (!fs::exists(fs::path(csec_gdml_file))) {
    printout(ERROR, "BarrelHCalCalorimeter_geo", "file " + csec_gdml_file + " does not exist");
    printout(ERROR, "BarrelHCalCalorimeter_geo", "use a FileLoader plugin before the field element");
    std::_Exit(EXIT_FAILURE);
  }

  Volume barrel_csector_vol = parser.GDMLReadFile(csec_gdml_file.c_str());
  if(!barrel_csector_vol.isValid()){
    printout(WARNING, "BarrelHCalCalorimeter", "%s", csec_gdml_file.c_str());
    printout(WARNING, "BarrelHCalCalorimeter", "barrel_csector_vol invalid, GDML parser failed!");
    std::_Exit(EXIT_FAILURE);
  }
  barrel_csector_vol.import();
  barrel_csector_vol.setVisAttributes(description, x_det.visStr());
  TessellatedSolid barrel_csector_solid = barrel_csector_vol.solid();
  barrel_csector_solid->CloseShape(true, true, true); // tesselated solid not closed by import!
  Material csector_material = description.material(csec_gdml_material.c_str());
  barrel_csector_vol.setMaterial(csector_material);

  // end ring
  EnsureFileFromURLExists(er_gdml_url, er_gdml_file, er_gdml_cache);
  if (!fs::exists(fs::path(er_gdml_file))) {
    printout(ERROR, "BarrelHCalCalorimeter_geo", "file " + er_gdml_file + " does not exist");
    printout(ERROR, "BarrelHCalCalorimeter_geo", "use a FileLoader plugin before the field element");
    std::_Exit(EXIT_FAILURE);
  }

  Volume barrel_er_vol = parser.GDMLReadFile(er_gdml_file.c_str());
  if(!barrel_er_vol.isValid()){
    printout(WARNING, "BarrelHCalCalorimeter", "%s", er_gdml_file.c_str());
    printout(WARNING, "BarrelHCalCalorimeter", "barrel_er_vol invalid, GDML parser failed!");
    std::_Exit(EXIT_FAILURE);
  }
  barrel_er_vol.import();
  barrel_er_vol.setVisAttributes(description, x_det.visStr());
  TessellatedSolid barrel_er_solid = barrel_er_vol.solid();
  barrel_er_solid->CloseShape(true, true, true); // tesselated solid not closed by import!
  Material er_material = description.material(er_gdml_material.c_str());
  barrel_er_vol.setMaterial(er_material);

  // Place steel in envelope

  double sec_rot_angle = 360.0/32.0;

  for(int k=0; k<29; k++){
    BarrelHCAL.placeVolume(barrel_sector_vol, k, Transform3D(RotationZ(-k*sec_rot_angle*dd4hep::deg)*RotationY(180.0* dd4hep::deg),Translation3D(0, 0, 0)));
  }
  BarrelHCAL.placeVolume(barrel_csector_vol, 0, Transform3D(RotationZ(sec_rot_angle*dd4hep::deg)*RotationY(180.0* dd4hep::deg),Translation3D(0, 0, 0)));
  BarrelHCAL.placeVolume(barrel_csector_vol, 1, Transform3D(RotationY(180.0* dd4hep::deg),Translation3D(0, 0, 0)));
  BarrelHCAL.placeVolume(barrel_csector_vol, 2, Transform3D(RotationZ(-sec_rot_angle*dd4hep::deg)*RotationY(180.0* dd4hep::deg),Translation3D(0, 0, 0)));
  BarrelHCAL.placeVolume(barrel_er_vol, 0, Transform3D(RotationY(180.0* dd4hep::deg),Translation3D(0, 0, 0)));
  BarrelHCAL.placeVolume(barrel_er_vol, 1, Transform3D(RotationY(0.0* dd4hep::deg),Translation3D(0, 0, 0)));

  // Loop over the tile solids, create them and add them to the detector volume

  Volume Tile[12];
  Volume ChimneyTile[4];

  for(int j=1; j<17; j++){

    std::string gdmlname;
    std::string solid_name;

    if(j<13){

      // standard tiles

      gdmlname = _toString(j,"tile%d_gdmlfile");
      solid_name = _toString(j,"OuterHCalTile%02d");

    }
    else{

      // chimney tiles

      gdmlname = _toString(j-4,"ctile%d_gdmlfile");
      solid_name = _toString(j-4,"OuterHCalChimneyTile%02d");

    }

    // tile shape gdml file info
    xml_comp_t x_det_tgdmlfile = x_det.child(gdmlname);

    std::string tgdml_file = getAttrOrDefault<std::string>(x_det_tgdmlfile, _Unicode(file), " ");;
    std::string tgdml_material = getAttrOrDefault<std::string>(x_det_tgdmlfile, _Unicode(material), " ");
    std::string tgdml_url = getAttrOrDefault<std::string>(x_det_tgdmlfile, _Unicode(url), " ");
    std::string tgdml_cache = getAttrOrDefault<std::string>(x_det_tgdmlfile, _Unicode(cache), " ");

    EnsureFileFromURLExists(tgdml_url, tgdml_file, tgdml_cache);
    if (!fs::exists(fs::path(tgdml_file))) {
      printout(ERROR, "BarrelHCalCalorimeter_geo", "file " + tgdml_file + " does not exist");
      printout(ERROR, "BarrelHCalCalorimeter_geo", "use a FileLoader plugin before the field element");
      std::_Exit(EXIT_FAILURE);
    }

    Volume solidVolume = parser.GDMLReadFile(tgdml_file.c_str());
    if(!solidVolume.isValid()){
      printout(WARNING, "BarrelHCalCalorimeter_geo", "%s", tgdml_file.c_str());
      printout(WARNING, "BarrelHCalCalorimeter_geo", "solidVolume invalid, GDML parser failed!");
      std::_Exit(EXIT_FAILURE);
    }
    solidVolume.import();
    solidVolume.setVisAttributes(description, x_det.visStr());
    TessellatedSolid volume_solid = solidVolume.solid();
    volume_solid->CloseShape(true, true, true); // tesselated solid not closed by import!
    Material tile_material = description.material(tgdml_material.c_str());
    solidVolume.setMaterial(tile_material);

    solidVolume.setSensitiveDetector(sens);

    // For tiles we build an assembly to get the full array of tiles
    // Offsets and rotation are to properly orient the tiles in the assembly.

    if (solid_name.size() > 0) {

      std::string type = solid_name.substr(0, solid_name.size() - 2);

      if (type == "OuterHCalTile" || type == "OuterHCalChimneyTile") {

        std::string stnum = solid_name.substr(solid_name.size() - 2, solid_name.size());
        int         tnum  = atoi(stnum.c_str()) - 1;

        // Tile numbers are indexed by the center (eta=0) out, we want them starting zero at one end.

        if (type == "OuterHCalTile") {

          Tile[11-tnum] = solidVolume;

        }
        else if ((tnum > 7) && (type == "OuterHCalChimneyTile")) {

          ChimneyTile[11-tnum] = solidVolume;

        }

      } else
        printout(WARNING, "BarrelHCalCalorimeter", "invalid solid_name, not a tile type?");

    } else
      printout(WARNING, "BarrelHCalCalorimeter", "solid_name.size() invalid! ");
  }

  // Place the tiles into the calorimeter volume

  double increment_angle = (360.0/320.0)*dd4hep::deg;
  double increment_offset = -10.0*increment_angle;

  DetElement tile_det("tile0", det_id);
  sens.setType("calorimeter");

  for (int i_eta = 0; i_eta < 12; i_eta++) {  // eta ring

    for (int i_phi = 0; i_phi < 320; i_phi++) { // phi index

      if (i_eta > 3) {

        // ordinary sector tiles

        PlacedVolume phv1 = BarrelHCAL.placeVolume(
                                                   Tile[i_eta], i_phi + (12 + i_eta) * 320,
                                                   RotationZ(i_phi * increment_angle + increment_offset) *
                                                   Transform3D(RotationY(90.0 * dd4hep::deg),
                                                               Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0 * dd4hep::mm)) *
                                                   RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                   Transform3D(RotationY(180.0 * dd4hep::deg),
                                                               Translation3D((xposTileN[i_eta] - ((i_eta-11) + 1) * tile_tolerance) * dd4hep::mm,
                                                                             yposTileN[i_eta] * dd4hep::mm,
                                                                             zposTileN[i_eta] * dd4hep::mm)) *
                                                   Translation3D(-xposTileS[i_eta] * dd4hep::mm,-yposTileS[i_eta] * dd4hep::mm, -zposTileS[i_eta] * dd4hep::mm));

        phv1.addPhysVolID("tile", i_phi + (12 + i_eta) * 320);
        DetElement sd1 = tile_det.clone(_toString(i_phi + (12 + i_eta) * 320, "tile%d"));
        sd1.setPlacement(phv1);
        sdet.add(sd1);

        PlacedVolume phv0 = BarrelHCAL.placeVolume(
                                                   Tile[i_eta], i_phi + i_eta * 320,
                                                   RotationZ(i_phi * increment_angle + increment_offset) *
                                                   Transform3D(RotationY(90.0 * dd4hep::deg),
                                                               Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0 * dd4hep::mm)) *
                                                   RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                   Translation3D((((i_eta-11) + 1) * tile_tolerance) * dd4hep::mm, 0.0 * dd4hep::mm, 0.0 * dd4hep::mm));
        phv0.addPhysVolID("tile", i_phi + i_eta * 320);
        DetElement sd0 = tile_det.clone(_toString(i_phi + i_eta * 320, "tile%d"));
        sd0.setPlacement(phv0);
        sdet.add(sd0);

      } else {

        // first three sectors are chimney sectors

        if(i_phi>29){

          // ordinary sector tiles

          PlacedVolume phv1 = BarrelHCAL.placeVolume(
                                                     Tile[i_eta], i_phi + i_eta * 320,
                                                     RotationZ(i_phi * increment_angle + increment_offset) *
                                                     Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                 Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0 * dd4hep::mm)) *
                                                     RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                     Transform3D(RotationY(180.0 * dd4hep::deg),
                                                                 Translation3D((xposTileN[i_eta] - ((i_eta-11) + 1) * tile_tolerance) * dd4hep::mm,
                                                                               yposTileN[i_eta] * dd4hep::mm,
                                                                               zposTileN[i_eta] * dd4hep::mm)) *
                                                     Translation3D(-xposTileS[i_eta] * dd4hep::mm,-yposTileS[i_eta] * dd4hep::mm, -zposTileS[i_eta] * dd4hep::mm));

          phv1.addPhysVolID("tile", i_phi + i_eta * 320);
          DetElement sd1 = tile_det.clone(_toString(i_phi + i_eta * 320, "tile%d"));
          sd1.setPlacement(phv1);
          sdet.add(sd1);

          PlacedVolume phv0 = BarrelHCAL.placeVolume(
                                                     Tile[i_eta], i_phi + (12 + i_eta) * 320,
                                                     RotationZ(i_phi * increment_angle + increment_offset) *
                                                     Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                 Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0 * dd4hep::mm)) *
                                                     RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                     Translation3D((((i_eta-11) + 1) * tile_tolerance) * dd4hep::mm, 0.0 * dd4hep::mm, 0.0 * dd4hep::mm));
          phv0.addPhysVolID("tile", i_phi + (12 + i_eta) * 320);
          DetElement sd0 = tile_det.clone(_toString(i_phi + (12 + i_eta) * 320, "tile%d"));
          sd0.setPlacement(phv0);
          sdet.add(sd0);

        }
        else{

          // chimney sector tile

          PlacedVolume phv1 = BarrelHCAL.placeVolume(
                                                     ChimneyTile[i_eta], i_phi + i_eta * 320,
                                                     RotationZ(i_phi * increment_angle + increment_offset) *
                                                     Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                 Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0 * dd4hep::mm)) *
                                                     RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                     Translation3D((((i_eta-11) + 1) * tile_tolerance) * dd4hep::mm, 0.0 * dd4hep::mm, 0.0 * dd4hep::mm));

          phv1.addPhysVolID("tile", i_phi + i_eta * 320);
          DetElement sd1 = tile_det.clone(_toString(i_phi + i_eta * 320, "tile%d"));
          sd1.setPlacement(phv1);
          sdet.add(sd1);

          PlacedVolume phv0 = BarrelHCAL.placeVolume(
                                                     Tile[i_eta], i_phi + (12 + i_eta) * 320,
                                                     RotationZ(i_phi * increment_angle + increment_offset) *
                                                     Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                 Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0 * dd4hep::mm)) *
                                                     RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                     Translation3D((((i_eta-11) + 1) * tile_tolerance) * dd4hep::mm, 0.0 * dd4hep::mm, 0.0 * dd4hep::mm));
          phv0.addPhysVolID("tile", i_phi + (12 + i_eta) * 320);
          DetElement sd0 = tile_det.clone(_toString(i_phi + (12 + i_eta) * 320, "tile%d"));
          sd0.setPlacement(phv0);
          sdet.add(sd0);

        }
      }
    }
  }

  // Place the detector into the envelope

  envelope.placeVolume(BarrelHCAL, 0, Transform3D(RotationZ(0.0),Translation3D(0, 0, 0)));

  std::string env_vis = getAttrOrDefault<std::string>(x_det, _Unicode(env_vis), "HcalBarrelEnvelopeVis");
  envelope.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), env_vis);
  return sdet;
}

DECLARE_DETELEMENT(epic_HcalBarrelGDML, create_detector)
