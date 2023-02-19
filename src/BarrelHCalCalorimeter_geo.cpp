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

  // get the solids section for this detector
  xml_comp_t x_solids = x_det.child("solids");

  DetElement sdet(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(sdet);

  // Create envelope to hold HCAL barrel

  double rmin1   = x_det.rmin1();
  double rmin2   = x_det.rmin2();
  double rmax    = x_det.rmax();
  double length1 = x_det.z1();
  double length2 = x_det.z2();

  std::vector<double> rmins = {rmin2, rmin2, rmin1, rmin1, rmin2, rmin2};
  std::vector<double> rmaxs = {rmax, rmax, rmax, rmax, rmax, rmax};
  std::vector<double> zs    = {-length2 / 2., -length1 / 2., -length1 / 2., length1 / 2., length1 / 2., length2 / 2.};

  // printout(WARNING, "BarrelHCalCalorimeter", "%f %f %f %f %f", rmin1, rmin2, rmax, length1/2., length2/2.);

  Polycone ptube(0.0, 2.0 * M_PI, rmins, rmaxs, zs);
  // Tube ptube(rmin1,rmax,length2);
  Volume envelope(det_name, ptube, air);

  PlacedVolume env_phv = motherVol.placeVolume(envelope);
  env_phv.addPhysVolID("system", det_id);
  //env_phv.addPhysVolID("barrel", 0);
  sdet.setPlacement(env_phv);

  // Storage for sectors and tile assemblies
  Assembly ChimneySector("ChimneySector");
  Assembly Sector("Sector");
  Assembly ChimneyTower[4];
  Assembly Tower[24];

  xml_comp_t det_define = x_det.child("define");

  // Pick up the constants

  double ctilePlaneRotate = 0.0;
  double tilePlaneRotate  = 0.0;
  double csectorRotate    = 0.0;
  double sectorRotate     = 0.0;

  double tile_tolerance = 0.2; // Tile tolerance in mm to avoid overlaps

  // Tile rotation starting points to align with sector plates
  double ctileRotateStart  = 5.4420 * (360.0 / 320.0) * dd4hep::deg;
  double octileRotateStart = 5.4420 * (360.0 / 320.0) * dd4hep::deg;
  // double tileRotateStart  = 20.38875*(360.0/320.0)*dd4hep::deg + ctileRotateStart;
  double tileRotateStart = 20.38675 * (360.0 / 320.0) * dd4hep::deg + ctileRotateStart;

  for (xml_coll_t i(det_define, _Unicode(constant)); i; ++i) {
    xml_comp_t x_const = i;

    std::string const_name  = getAttrOrDefault<std::string>(x_const, _Unicode(name), " ");
    std::string const_value = getAttrOrDefault<std::string>(x_const, _Unicode(value), " ");

    if (const_name == "ctilePlaneRotate")
      ctilePlaneRotate = atof(const_value.c_str());
    else if (const_name == "tilePlaneRotate")
      tilePlaneRotate = atof(const_value.c_str());
    else if (const_name == "csectorRotate")
      csectorRotate = atof(const_value.c_str());
    else if (const_name == "sectorRotate")
      sectorRotate = atof(const_value.c_str());
    else
      printout(WARNING, "BarrelHCalCalorimeter", "unrecognized <constant> data!");
  }

  // Loop over the defines section and pick up the tile offsets

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

  std::vector<double> plates_x;
  std::vector<double> plates_y;
  std::vector<double> plates_z;

  std::vector<double> tweak_tiles;
  std::vector<double> tweak_chimney_tiles;
  std::vector<double> tweak_sectors;

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
    else if (mtrx_name == "plates_x")
      aptr = &plates_x;
    else if (mtrx_name == "plates_y")
      aptr = &plates_y;
    else if (mtrx_name == "plates_z")
      aptr = &plates_z;
    else if (mtrx_name == "tweak_tiles")
      aptr = &tweak_tiles;
    else if (mtrx_name == "tweak_chimney_tiles")
      aptr = &tweak_chimney_tiles;
    else if (mtrx_name == "tweak_sectors")
      aptr = &tweak_sectors;
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

  // Loop over the solids, create them and add them to the detector volume

  for (xml_coll_t k(x_solids, _Unicode(solid)); k; ++k) {

    xml_comp_t x_solid = k;

    // get the sector solid definitions
    xml_comp_t define      = x_solid.child("define");
    xml_comp_t tessellated = x_solid.child("tessellated");

    std::string solid_name     = getAttrOrDefault<std::string>(x_solid, _Unicode(name), " ");
    std::string solidMatString = getAttrOrDefault<std::string>(x_solid, _Unicode(material), " ");
    Material    solid_material = description.material(solidMatString);

    double offset_x = atof(getAttrOrDefault<std::string>(x_solid, _Unicode(x), "0").c_str()) * dd4hep::mm;
    double offset_y = atof(getAttrOrDefault<std::string>(x_solid, _Unicode(y), "0").c_str()) * dd4hep::mm;
    double offset_z = atof(getAttrOrDefault<std::string>(x_solid, _Unicode(z), "0").c_str()) * dd4hep::mm;

    // Get the vertices
    std::vector<Tessellated::Vertex_t> vertices;
    for (xml_coll_t j(define, _Unicode(position)); j; ++j) {
      xml_comp_t pos = j;

      // create the vertex point

      double xp = atof(getAttrOrDefault<std::string>(pos, _Unicode(x), "0").c_str()) * dd4hep::mm - offset_x;
      double yp = atof(getAttrOrDefault<std::string>(pos, _Unicode(y), "0").c_str()) * dd4hep::mm - offset_y;
      double zp = atof(getAttrOrDefault<std::string>(pos, _Unicode(z), "0").c_str()) * dd4hep::mm - offset_z;

      // for the sector plates  we perform a rotation around y - the chimney cutout should be in the
      // electron arm

      if ((solid_name == "HCAL_Chimney_Sector_Half_Plate") || (solid_name == "HCAL_Chimney_Sector_Plate") ||
          (solid_name == "HCAL_Sector_Half_Plate") || (solid_name == "HCAL_Sector_Plate")) {
        xp = -xp;
        zp = -zp;
      }

      Tessellated::Vertex_t thisPoint(xp, yp, zp);

      vertices.push_back(thisPoint);
    }

    TessellatedSolid solid(solid_name.c_str(), vertices);

    for (xml_coll_t i(tessellated, _Unicode(triangular)); i; ++i) {
      xml_comp_t triang = i;

      int vtx1 = -1;
      int vtx2 = -1;
      int vtx3 = -1;

      std::string facetName1 = getAttrOrDefault<std::string>(triang, _Unicode(vertex1), " ");
      std::string facetName2 = getAttrOrDefault<std::string>(triang, _Unicode(vertex2), " ");
      std::string facetName3 = getAttrOrDefault<std::string>(triang, _Unicode(vertex3), " ");

      // Search the define collection to match things up
      int idx = 0;
      for (xml_coll_t j(define, _Unicode(position)); j; ++j) {
        xml_comp_t  pos     = j;
        std::string posName = getAttrOrDefault<std::string>(pos, _Unicode(name), " ");

        if (posName == facetName1)
          vtx1 = idx;
        if (posName == facetName2)
          vtx2 = idx;
        if (posName == facetName3)
          vtx3 = idx;

        if ((vtx1 >= 0) && (vtx2 >= 0) && (vtx3 >= 0))
          break;

        idx++;
      }

      // Add the facet to the solid

      if ((vtx1 >= 0) && (vtx2 >= 0) && (vtx3 >= 0) && (vtx1 != vtx2) && (vtx1 != vtx3) && (vtx2 != vtx3)) {

        solid->AddFacet(vtx1, vtx2, vtx3);

      } else
        printout(WARNING, "BarrelHCalCalorimeter", "bad facet! %d %d %d", vtx1, vtx2, vtx3);
    }

    // Complete the shape
    solid->CloseShape(true, true, true);

    Volume solidVolume(solid_name, solid, solid_material);
    solidVolume.setVisAttributes(description, x_det.visStr());

    // printout(WARNING, "BarrelHCalCalorimeter", "tesselated solid name %s", solid_name.c_str());

    if (solid_name == "HCAL_Chimney_Sector_Half_Plate") {

      ChimneySector.placeVolume(
          solidVolume, 0,
          RotationZ(csectorRotate * dd4hep::deg) *
              Transform3D(RotationZ(0.0),
                          Translation3D(plates_x[0] * dd4hep::mm, plates_y[0] * dd4hep::mm, plates_z[0] * dd4hep::mm)));

      ChimneySector.placeVolume(
          solidVolume, 1,
          RotationZ((9.60 * 2 * M_PI / 320) + csectorRotate * dd4hep::deg) *
              Transform3D(RotationZ(0.0),
                          Translation3D(plates_x[0] * dd4hep::mm, plates_y[0] * dd4hep::mm, plates_z[0] * dd4hep::mm)));

    } else if (solid_name == "HCAL_Chimney_Sector_Plate") {

      for (int i = 0; i < 9; i++)
        ChimneySector.placeVolume(
            solidVolume, i,
            RotationZ((i * 2 * M_PI / 320) + csectorRotate * dd4hep::deg) *
                Transform3D(RotationZ(0.0), Translation3D(plates_x[1] * dd4hep::mm, plates_y[1] * dd4hep::mm,
                                                          plates_z[1] * dd4hep::mm)));

    } else if (solid_name == "HCAL_Sector_Half_Plate") {

      Sector.placeVolume(
          solidVolume, 0,
          RotationZ(sectorRotate * dd4hep::deg) *
              Transform3D(RotationZ(0.4 * dd4hep::deg),
                          Translation3D(plates_x[2] * dd4hep::mm, plates_y[2] * dd4hep::mm, plates_z[2] * dd4hep::mm)));

      Sector.placeVolume(
          solidVolume, 1,
          RotationZ((9.60 * 2 * M_PI / 320) + sectorRotate * dd4hep::deg) *
              Transform3D(RotationZ(0.4 * dd4hep::deg),
                          Translation3D(plates_x[2] * dd4hep::mm, plates_y[2] * dd4hep::mm, plates_z[2] * dd4hep::mm)));

    } else if (solid_name == "HCAL_Sector_Plate") {

      for (int i = 0; i < 9; i++)
        Sector.placeVolume(solidVolume, i,
                           RotationZ((i * 2 * M_PI / 320) + sectorRotate * dd4hep::deg) *
                               Transform3D(RotationZ(0.4 * dd4hep::deg),
                                           Translation3D(plates_x[3] * dd4hep::mm, plates_y[3] * dd4hep::mm,
                                                         plates_z[3] * dd4hep::mm)));

    } else {

      // If it's not sectors then it's a tile - for these we build an assembly to get the full array of tiles
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
                    RotationZ(octileRotateStart + i * (360.0 / 320.0) * dd4hep::deg) *
                        Transform3D(RotationY(90.0 * dd4hep::deg),
                                    Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                        RotationX(-tilePlaneRotate * dd4hep::deg) *
                        Transform3D(RotationY(0.0),
                                    Translation3D((xposTileS[tnum] + (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                  yposTileS[tnum] * dd4hep::mm, zposTileS[tnum] * dd4hep::mm)));

                phv0.addPhysVolID("tile", i);
                DetElement sd0 = tile_det.clone(_toString(i + (11 - tnum) * 10, "tile%d"));
                sd0.setPlacement(phv0);
                sdet.add(sd0);

                PlacedVolume phv1 = Tower[12 + tnum].placeVolume(
                    solidVolume, i + 5,
                    RotationZ(octileRotateStart + i * (360.0 / 320.0) * dd4hep::deg) *
                        Transform3D(RotationY(90.0 * dd4hep::deg),
                                    Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                        RotationX(-tilePlaneRotate * dd4hep::deg) *
                        Transform3D(RotationY(180.0 * dd4hep::deg),
                                    Translation3D((xposTileN[tnum] - (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                  yposTileN[tnum] * dd4hep::mm, zposTileN[tnum] * dd4hep::mm)));
                phv1.addPhysVolID("tile", i);
                DetElement sd1 = tile_det.clone(_toString(i + 5 + (12 + tnum) * 10, "tile%d"));
                sd1.setPlacement(phv1);
                sdet.add(sd1);

              } else {

                PlacedVolume phv0 = Tower[11 - tnum].placeVolume(
                    solidVolume, i,
                    RotationZ(octileRotateStart + i * (360.0 / 320.0) * dd4hep::deg) *
                        Transform3D(RotationY(90.0 * dd4hep::deg),
                                    Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                        RotationX(-tilePlaneRotate * dd4hep::deg) *
                        Transform3D(RotationY(180.0 * dd4hep::deg),
                                    Translation3D((xposTileS[tnum] + (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                  yposTileS[tnum] * dd4hep::mm, zposTileS[tnum] * dd4hep::mm)));

                phv0.addPhysVolID("tile", i);
                DetElement sd0 = tile_det.clone(_toString(i + (11 - tnum) * 10, "tile%d"));
                sd0.setPlacement(phv0);
                sdet.add(sd0);

                PlacedVolume phv1 = Tower[12 + tnum].placeVolume(
                    solidVolume, i + 5,
                    RotationZ(octileRotateStart + i * (360.0 / 320.0) * dd4hep::deg) *
                        Transform3D(RotationY(90.0 * dd4hep::deg),
                                    Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                        RotationX(-tilePlaneRotate * dd4hep::deg) *
                        Transform3D(RotationY(0.0),
                                    Translation3D((xposTileN[tnum] - (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                  yposTileN[tnum] * dd4hep::mm, zposTileN[tnum] * dd4hep::mm)));
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
                  RotationZ(ctileRotateStart + i * (360.0 / 320.0) * dd4hep::deg) *
                      Transform3D(RotationY(90.0 * dd4hep::deg),
                                  Translation3D(xposOuter[0] * dd4hep::mm, yposOuter[0] * dd4hep::mm, 0.0)) *
                      RotationX(-ctilePlaneRotate * dd4hep::deg) *
                      Transform3D(RotationY(0.0),
                                  Translation3D((xposChimneyTileS[tnum - 8] + (tnum + 1) * tile_tolerance) * dd4hep::mm,
                                                yposChimneyTileS[tnum - 8] * dd4hep::mm,
                                                zposChimneyTileS[tnum - 8] * dd4hep::mm)));
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
  }

  // Place the sector tile assemblies into the sectors

  sens.setType("calorimeter");

  DetElement tower_det("tower0", det_id);

  // special chimney sector towers
  for (int i = 0; i < 4; i++) {

    PlacedVolume tower_phv0 = ChimneySector.placeVolume(
        ChimneyTower[i], i, Transform3D(RotationZ(tweak_chimney_tiles[i]), Translation3D(0.0, 0.0, 0.0)));
    tower_phv0.addPhysVolID("tower", i); // lower phi
    DetElement sd0 = tower_det.clone(_toString(i, "tower%d"));
    sd0.setPlacement(tower_phv0);
    sdet.add(sd0);

    PlacedVolume tower_phv1 =
        ChimneySector.placeVolume(ChimneyTower[i], i + 24,
                                  Transform3D(RotationZ(5 * (360.0 / 320.0) * dd4hep::deg + tweak_chimney_tiles[i]),
                                              Translation3D(0.0, 0.0, 0.0)));
    tower_phv1.addPhysVolID("tower", i+24); // upper phi
    DetElement sd1 = tower_det.clone(_toString(i + 24, "tower%d"));
    sd1.setPlacement(tower_phv1);
    sdet.add(sd1);
  }

  // ordinary towers in chimney sectors
  for (int i = 4; i < 24; i++) {

    PlacedVolume tower_phv0 = ChimneySector.placeVolume(
        Tower[i], i, Transform3D(RotationZ(tweak_chimney_tiles[i]), Translation3D(0.0, 0.0, 0.0)));
    tower_phv0.addPhysVolID("tower", i); // lower phi
    DetElement sd0 = tower_det.clone(_toString(i, "tower%d"));
    sd0.setPlacement(tower_phv0);
    sdet.add(sd0);

    PlacedVolume tower_phv1 =
        ChimneySector.placeVolume(Tower[i], i + 24,
                                  Transform3D(RotationZ(5 * (360.0 / 320.0) * dd4hep::deg + tweak_chimney_tiles[i]),
                                              Translation3D(0.0, 0.0, 0.0)));
    tower_phv1.addPhysVolID("tower", i+24); // upper phi
    DetElement sd1 = tower_det.clone(_toString(i + 24, "tower%d"));
    sd1.setPlacement(tower_phv1);
    sdet.add(sd1);
  }

  // ordinary sectors
  for (int i = 0; i < 24; i++) {

    PlacedVolume tower_phv0 = Sector.placeVolume(
        Tower[i], i + 24,
        Transform3D(RotationZ(tileRotateStart - octileRotateStart + tweak_tiles[i]), Translation3D(0.0, 0.0, 0.0)));
    tower_phv0.addPhysVolID("tower", i); // lower phi
    DetElement sd0 = tower_det.clone(_toString(i + 48, "tower%d"));
    sd0.setPlacement(tower_phv0);
    sdet.add(sd0);

    PlacedVolume tower_phv1 = Sector.placeVolume(
        Tower[i], i + 48,
        Transform3D(RotationZ(tileRotateStart - octileRotateStart + 5 * (360.0 / 320.0) * dd4hep::deg + tweak_tiles[i]),
                    Translation3D(0.0, 0.0, 0.0)));
    tower_phv1.addPhysVolID("tower", i+24); // upper phi
    DetElement sd1 = tower_det.clone(_toString(i + 72, "tower%d"));
    sd1.setPlacement(tower_phv1);
    sdet.add(sd1);
  }

  // Place the sectors into the envelope

  DetElement sector_det("sector0", det_id);

  // Chimney sectors
  for (int i = -1; i < 2; i++) {
    PlacedVolume sect_phv = envelope.placeVolume(
        ChimneySector, i + 1,
        Transform3D(RotationZ(((i - 1) * 2 * M_PI / 32) + tweak_sectors[i + 1]), Translation3D(0, 0, 0)));
    sect_phv.addPhysVolID("system", det_id);
    //sect_phv.addPhysVolID("barrel", 0);
    sect_phv.addPhysVolID("sector", i + 1);
    DetElement sd = sector_det.clone(_toString(i + 1, "sector%d"));
    sd.setPlacement(sect_phv);
    sdet.add(sd);
  }

  // Normal sectors

  for (int i = 3; i < 32; i++) {
    PlacedVolume sect_phv =
        envelope.placeVolume(Sector, i,
                             Transform3D(RotationZ((-2.075 * M_PI / 32) + (i - 3) * (2 * M_PI / 32) + tweak_sectors[i]),
                                         Translation3D(0, 0, 0)));
    sect_phv.addPhysVolID("system", det_id);
    sect_phv.addPhysVolID("sector", i);
    DetElement sd = sector_det.clone(_toString(i, "sector%d"));
    sd.setPlacement(sect_phv);
    sdet.add(sd);
  }

  std::string env_vis = getAttrOrDefault<std::string>(x_det, _Unicode(env_vis), "HcalBarrelEnvelopeVis");
  envelope.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), env_vis);
  return sdet;
}

DECLARE_DETELEMENT(epic_HcalBarrelGDML, create_detector)
