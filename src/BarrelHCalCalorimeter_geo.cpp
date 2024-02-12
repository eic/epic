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
  // Tube ptube(rmin1,rmax,length2);
  Volume envelope(det_name, ptube, air);

  PlacedVolume env_phv = motherVol.placeVolume(envelope);
  env_phv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_phv);

  // Storage for sectors and tile assemblies
  Assembly ChimneyTower[4];
  Assembly Tower[24];

  xml_comp_t det_define = x_det.child("define");

  // Pick up the constants

  double ctilePlaneRotate = 0.0;
  double tilePlaneRotate  = 0.0;

  double tile_tolerance = 0.2; // Tile tolerance in mm to avoid overlaps

  // Sector steel tessellated shape gdml file info
  xml_comp_t x_det_gdmlfile = x_det.child("gdmlfile");

  std::string gdml_file = getAttrOrDefault<std::string>(x_det_gdmlfile, _Unicode(file), " ");;
  std::string gdml_material = getAttrOrDefault<std::string>(x_det_gdmlfile, _Unicode(material), " ");
  std::string gdml_url = getAttrOrDefault<std::string>(x_det_gdmlfile, _Unicode(url), " ");
  std::string gdml_cache = getAttrOrDefault<std::string>(x_det_gdmlfile, _Unicode(cache), " ");

  // Loop over the defines section and pick up the tile offsets, ref location and angles

  for (xml_coll_t i(det_define, _Unicode(constant)); i; ++i) {
    xml_comp_t x_const = i;

    std::string const_name  = getAttrOrDefault<std::string>(x_const, _Unicode(name), " ");

    if (const_name == "ctilePlaneRotate"){
      std::string const_value = getAttrOrDefault<std::string>(x_const, _Unicode(value), " ");
      ctilePlaneRotate = atof(const_value.c_str());
    }
    else if (const_name == "tilePlaneRotate"){
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

  double increment_angle = (360.0/320.0)*dd4hep::deg;

  // Read in the barrel structure GDML file
  Assembly BarrelHCAL("BarrelHCAL");

  EnsureFileFromURLExists(gdml_url, gdml_file, gdml_cache);
  if (!fs::exists(fs::path(gdml_file))) {
    printout(ERROR, "BarrelHCalCalorimeter_geo", "file " + gdml_file + " does not exist");
    printout(ERROR, "BarrelHCalCalorimeter_geo", "use a FileLoader plugin before the field element");
    std::_Exit(EXIT_FAILURE);
  }

  TGDMLParse parser;
  Volume barrel_steel_vol = parser.GDMLReadFile(gdml_file.c_str());
  if(!barrel_steel_vol.isValid()){
    printout(WARNING, "BarrelHCalCalorimeter", "%s", gdml_file.c_str());
    printout(WARNING, "BarrelHCalCalorimeter", "barrel_steel_vol invalid, GDML parser failed!");
    std::_Exit(EXIT_FAILURE);
  }
  barrel_steel_vol.import();
  barrel_steel_vol.setVisAttributes(description, x_det.visStr());
  TessellatedSolid barrel_steel_solid = barrel_steel_vol.solid();
  barrel_steel_solid->CloseShape(true, true, true); // tesselated solid not closed by import!
  Material sector_material = description.material(gdml_material.c_str());
  barrel_steel_vol.setMaterial(sector_material);

  // Place steel in envelope
  BarrelHCAL.placeVolume(barrel_steel_vol, 0, Transform3D(RotationY(180.0* dd4hep::deg),Translation3D(0, 0, 0)));

  // Loop over the tile solids, create them and add them to the detector volume

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

    // For tiles we build an assembly to get the full array of tiles
    // Offsets and rotation are to properly orient the tiles in the assembly.

    if (solid_name.size() > 0) {

      std::string type = solid_name.substr(0, solid_name.size() - 2);

      if (type == "OuterHCalTile" || type == "OuterHCalChimneyTile") {

        std::string stnum = solid_name.substr(solid_name.size() - 2, solid_name.size());
        int         tnum  = atoi(stnum.c_str()) - 1;

        Assembly TempTower1(_toString(11 - tnum, "Tower%i"));
        Assembly TempTower2(_toString(12 + tnum, "Tower%i"));

        solidVolume.setSensitiveDetector(sens);

        DetElement tile_det("tile0", det_id);

        if (type == "OuterHCalTile") {

          Tower[11 - tnum] = TempTower1;
          Tower[12 + tnum] = TempTower2;

          for (int i = 0; i < 5; i++) {

            if (tnum < 8) {

              PlacedVolume phv0 = Tower[11 - tnum].placeVolume(
                                                               solidVolume, i,
                                                               RotationZ(i * increment_angle) *
                                                               Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                           Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                                                               RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                               Transform3D(RotationY(0.0),
                                                                           Translation3D((xposTileS[tnum] + (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                                                         yposTileS[tnum] * dd4hep::mm, zposTileS[tnum] * dd4hep::mm)) *
                                                               Translation3D(-xposTileS[tnum] * dd4hep::mm, -yposTileS[tnum] * dd4hep::mm, -zposTileS[tnum] * dd4hep::mm));

              phv0.addPhysVolID("tile", i);
              DetElement sd0 = tile_det.clone(_toString(i + (11 - tnum) * 10, "tile%d"));
              sd0.setPlacement(phv0);
              sdet.add(sd0);

              PlacedVolume phv1 = Tower[12 + tnum].placeVolume(
                                                               solidVolume, i + 5,
                                                               RotationZ(i * increment_angle) *
                                                               Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                           Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                                                               RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                               Transform3D(RotationY(180.0 * dd4hep::deg),
                                                                           Translation3D((xposTileN[tnum] - (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                                                         yposTileN[tnum] * dd4hep::mm, zposTileN[tnum] * dd4hep::mm)) *
                                                               Translation3D(-xposTileS[tnum] * dd4hep::mm, -yposTileS[tnum] * dd4hep::mm, -zposTileS[tnum] * dd4hep::mm));

              phv1.addPhysVolID("tile", i);
              DetElement sd1 = tile_det.clone(_toString(i + 5 + (12 + tnum) * 10, "tile%d"));
              sd1.setPlacement(phv1);
              sdet.add(sd1);

            } else {

              PlacedVolume phv0 = Tower[11 - tnum].placeVolume(
                                                               solidVolume, i,
                                                               RotationZ(i * increment_angle) *
                                                               Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                           Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                                                               RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                               Transform3D(RotationY(180.0 * dd4hep::deg),
                                                                           Translation3D((xposTileS[tnum] + (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                                                         yposTileS[tnum] * dd4hep::mm, zposTileS[tnum] * dd4hep::mm)) *
                                                               Translation3D(xposTileS[tnum] * dd4hep::mm, -yposTileS[tnum] * dd4hep::mm, -zposTileS[tnum] * dd4hep::mm));

              phv0.addPhysVolID("tile", i);
              DetElement sd0 = tile_det.clone(_toString(i + (11 - tnum) * 10, "tile%d"));
              sd0.setPlacement(phv0);
              sdet.add(sd0);

              PlacedVolume phv1 = Tower[12 + tnum].placeVolume(
                                                               solidVolume, i + 5,
                                                               RotationZ(i * increment_angle) *
                                                               Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                           Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                                                               RotationX(-tilePlaneRotate * dd4hep::deg) *
                                                               Transform3D(RotationY(0.0),
                                                                           Translation3D((xposTileN[tnum] - (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                                                         yposTileN[tnum] * dd4hep::mm, zposTileN[tnum] * dd4hep::mm)) *
                                                               Translation3D(xposTileS[tnum] * dd4hep::mm, -yposTileS[tnum] * dd4hep::mm, -zposTileS[tnum] * dd4hep::mm));
              phv1.addPhysVolID("tile", i);
              DetElement sd1 = tile_det.clone(_toString(i + 5 + (12 + tnum) * 10, "tile%d"));
              sd1.setPlacement(phv1);
              sdet.add(sd1);
            }
          }
        }

        if ((tnum > 7) && (type == "OuterHCalChimneyTile")) {

          Assembly TempChimneyTower1(_toString(11 - tnum, "ChimneyTower%i"));
          ChimneyTower[11 - tnum] = TempChimneyTower1;

          for (int i = 0; i < 5; i++) {

            PlacedVolume phv = ChimneyTower[11 - tnum].placeVolume(
                                                                   solidVolume, i,
                                                                   RotationZ(i * increment_angle) *
                                                                   Transform3D(RotationY(90.0 * dd4hep::deg),
                                                                               Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                                                                   RotationX(-ctilePlaneRotate * dd4hep::deg) *
                                                                   Transform3D(RotationY(0.0),
                                                                               Translation3D((xposChimneyTileS[tnum - 8] + (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                                                             yposChimneyTileS[tnum - 8] * dd4hep::mm,
                                                                                             zposChimneyTileS[tnum - 8] * dd4hep::mm)) *
                                                                   Translation3D(-xposChimneyTileS[tnum-8] * dd4hep::mm, -yposChimneyTileS[tnum-8] * dd4hep::mm,
                                                                                 -zposChimneyTileS[tnum-8] * dd4hep::mm));
            phv.addPhysVolID("tile", i);
            DetElement sd = tile_det.clone(_toString(i + (11 - tnum) * 10 + 480, "tile%d"));
            sd.setPlacement(phv);
            sdet.add(sd);
          }
        }

      } else
        printout(WARNING, "BarrelHCalCalorimeter", "invalid solid_name, not a tile type?");

    } else
      printout(WARNING, "BarrelHCalCalorimeter", "solid_name.size() invalid! ");
  }


  // Place the sector tile assemblies into the sectors

  sens.setType("calorimeter");

  DetElement tower_det("tower0", det_id);

  for (int j = 0; j < 32; j++) {  // sector

    // chimney sectors
    if( j < 3) {

      // special chimney sector towers
      for (int i = 0; i < 4; i++) {

        PlacedVolume tower_phv0 =
          BarrelHCAL.placeVolume(ChimneyTower[i], i + j*48, Transform3D(RotationZ((10 * (j-1) * (360.0 / 320.0) - 0.042) * dd4hep::deg), Translation3D(0.0, 0.0, 0.0)));
        tower_phv0.addPhysVolID("tower", i + j*48);
        DetElement tt0 = tower_det.clone(_toString(i + j*48, "tower%d"));
        tt0.setPlacement(tower_phv0);
        sdet.add(tt0);

        PlacedVolume tower_phv1 =
          BarrelHCAL.placeVolume(ChimneyTower[i], i + 24 + j*48, Transform3D(RotationZ(((10 * (j-1) + 5) * (360.0 / 320.0) - 0.042) * dd4hep::deg), Translation3D(0.0, 0.0, 0.0)));
        tower_phv1.addPhysVolID("tower", i + 24 + j*48);
        DetElement tt1 = tower_det.clone(_toString(i + 24 + j*48, "tower%d"));
        tt1.setPlacement(tower_phv1);
        sdet.add(tt1);

      }

      // ordinary towers in chimney sectors
      for (int i = 4; i < 24; i++) {

        PlacedVolume tower_phv0 =
          BarrelHCAL.placeVolume(Tower[i], i + j*48, Transform3D(RotationZ((10 * (j-1) * (360.0 / 320.0) - 0.042) * dd4hep::deg), Translation3D(0.0, 0.0, 0.0)));
        tower_phv0.addPhysVolID("tower", i + j*48);
        DetElement tt0 = tower_det.clone(_toString(i + j*48, "tower%d"));
        tt0.setPlacement(tower_phv0);
        sdet.add(tt0);

        PlacedVolume tower_phv1 =
          BarrelHCAL.placeVolume(Tower[i], i + 24 + j*48, Transform3D(RotationZ(((10 * (j-1) + 5) * (360.0 / 320.0) - 0.042) * dd4hep::deg), Translation3D(0.0, 0.0, 0.0)));
        tower_phv1.addPhysVolID("tower", i + 24 + j*48);
        DetElement tt1 = tower_det.clone(_toString(i + 24 + j*48, "tower%d"));
        tt1.setPlacement(tower_phv1);
        sdet.add(tt1);

      }

    }
    // ordinary sectors
    else {

      for (int i = 0; i < 24; i++) {

        PlacedVolume tower_phv0 =
          BarrelHCAL.placeVolume(Tower[i], i + j*48, Transform3D(RotationZ((10 * (j-1) * (360.0 / 320.0) - 0.0525) * dd4hep::deg), Translation3D(0.0, 0.0, 0.0)));
        tower_phv0.addPhysVolID("tower", i + j*48);
        DetElement tt0 = tower_det.clone(_toString(i + j*48, "tower%d"));
        tt0.setPlacement(tower_phv0);
        sdet.add(tt0);

        PlacedVolume tower_phv1 =
          BarrelHCAL.placeVolume(Tower[i], i + 24 + j*48, Transform3D(RotationZ(((10 * (j-1) + 5) * (360.0 / 320.0) - 0.0525) * dd4hep::deg), Translation3D(0.0, 0.0, 0.0)));
        tower_phv1.addPhysVolID("tower", i + 24 + j*48);
        DetElement tt1 = tower_det.clone(_toString(i + 24 + j*48, "tower%d"));
        tt1.setPlacement(tower_phv1);
        sdet.add(tt1);

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
