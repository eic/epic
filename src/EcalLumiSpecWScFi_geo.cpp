// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Aranya Giri, Simon Gardner

/* Date : 12/12/2023
W Scifi (EM Calorimeter) Pair Spectrometer */

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <TVector3.h>

using namespace std;
using namespace dd4hep;

// Definition of function to build the modules
static tuple<Volume, Position> build_specScFiCAL_module(const Detector& description,
                                                        const xml::Component& mod_x,
                                                        SensitiveDetector& sens);

// Driver Function
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  sens.setType("calorimeter");

  xml_det_t x_det  = e;
  xml_comp_t x_mod = x_det.child(_Unicode(module));
  string det_name  = x_det.nameStr();
  int det_ID       = x_det.id();

  Material Air = description.material("Air");

  // Create main detector element to be returned at the end
  DetElement det(det_name, det_ID);

  // Mother volume
  Volume motherVol = description.pickMotherVolume(det);

  // Detector assembly
  Assembly assembly(det_name);
  assembly.setVisAttributes(description.visAttributes(x_det.attr<std::string>(_Unicode(vis))));

  double detSizeXY  = getAttrOrDefault(x_det, _Unicode(sizeXY), 180 * mm);
  double detSizeZ   = getAttrOrDefault(x_det, _Unicode(sizeZ), 180 * mm);
  int nmod_perlayer = getAttrOrDefault(x_det, _Unicode(nmod_perlayer), 3);
  int nlayer        = getAttrOrDefault(x_det, _Unicode(nlayer), 20);

  // Global detector position and resolution
  xml_comp_t pos = x_det.position();
  xml_comp_t rot = x_det.rotation();

  //function for module definition
  auto [modVol, modSize] = build_specScFiCAL_module(description, x_mod, sens);

  //layer definition before rotation
  xml_comp_t x_layer = x_det.child(_Unicode(layer));

  double laySizeX       = getAttrOrDefault(x_layer, _Unicode(sizeX), 1.0 * cm);
  double laySizeY       = getAttrOrDefault(x_layer, _Unicode(sizeY), 1.0 * cm);
  double laySizeZ       = getAttrOrDefault(x_layer, _Unicode(sizeZ), 1.0 * cm);
  double layerCoatSizeX = getAttrOrDefault(x_layer, _Unicode(coatSizeX), 1.0 * cm);
  double layerCoatSizeY = getAttrOrDefault(x_layer, _Unicode(coatSizeY), 1.0 * cm);

  Material layMat = description.material(x_layer.attr<std::string>(_Unicode(material)));

  Box layerBox(laySizeX / 2.0, laySizeY / 2.0, laySizeZ / 2.0);
  Volume layerVol("layer", layerBox, layMat);
  layerVol.setVisAttributes(description.visAttributes(x_det.attr<std::string>(_Unicode(vis))));

  // Position of first module, layer from center of CAL
  double mod_pos0   = -(detSizeXY / 2.0) + (modSize.x() / 2.0) + layerCoatSizeX;
  double layer_pos0 = -(detSizeZ / 2.0) + (modSize.y() / 2.0) + layerCoatSizeY;

  //Fill uncoated layer with modules
  for (int mod_id = 0; mod_id < nmod_perlayer; mod_id++) {

    //Build // to z-axis, then rotate
    double mod_pos_z = 0.0 * cm;
    double mod_pos_y = 0.0 * cm;
    double mod_pos_x = mod_id * modSize.x() + mod_pos0;

    PlacedVolume modPV = layerVol.placeVolume(modVol, Position(mod_pos_x, mod_pos_y, mod_pos_z));
    modPV.addPhysVolID("module", mod_id);
  } //imod-loop close

  //sector definition
  Box sectorBox(detSizeXY / 2, detSizeXY / 2, detSizeZ / 2);
  Volume sectorVol(det_name + "_sector", sectorBox, Air);
  sectorVol.setVisAttributes(description.visAttributes(x_mod.attr<std::string>(_Unicode(vis))));

  //Fill sector with layers (coated)
  for (int layer_id = 0; layer_id < nlayer; layer_id++) {

    double lay_pos_z = -layer_id * (modSize.y() + 2.0 * layerCoatSizeY) - layer_pos0;
    double lay_pos_y = 0.0 * cm;
    double lay_pos_x = 0.0 * cm;
    int orientation  = layer_id % 2 == 0;

    //rotation
    RotationZYX lay_rot = RotationZYX(0, 0, -90.0 * degree);
    if (orientation)
      lay_rot *= RotationY(-90.0 * degree);

    PlacedVolume layPV = sectorVol.placeVolume(
        layerVol, Transform3D(lay_rot, Position(lay_pos_x, lay_pos_y, lay_pos_z)));
    layPV.addPhysVolID("layer", layer_id);
  } //layer_id-loop close

  // loop over sectors(top, bottom)
  for (xml_coll_t si(x_det, _Unicode(sector)); si; si++) {

    xml_comp_t x_sector(si);
    int sector_id      = x_sector.id();
    xml_comp_t sec_pos = x_sector.position();
    xml_comp_t sec_rot = x_sector.rotation();

    PlacedVolume secPV = assembly.placeVolume(
        sectorVol, Transform3D(RotationZYX(sec_rot.z(), sec_rot.y(), sec_rot.x()),
                               Position(sec_pos.x(), sec_pos.y(), sec_pos.z())));
    secPV.addPhysVolID("sector", sector_id);
  } // sectors

  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume detPV =
      motherVol.placeVolume(assembly, Transform3D(RotationZYX(rot.z(), rot.y(), rot.x()),
                                                  Position(pos.x(), pos.y(), pos.z())));
  detPV.addPhysVolID("system", det_ID);
  // Connect to system ID
  det.setPlacement(detPV);

  return det;
} //Driver class close

static tuple<Volume, Position> build_specScFiCAL_module(const Detector& description,
                                                        const xml::Component& mod_x,
                                                        SensitiveDetector& sens) {

  //--------------------Module Setup---------------------------------------------------------------------
  double sx = mod_x.attr<double>(_Unicode(sizex));
  double sy = mod_x.attr<double>(_Unicode(sizey));
  double sz = mod_x.attr<double>(_Unicode(sizez));

  Position modSize(sx, sy, sz);

  Box modShape(modSize.x() / 2.0, modSize.y() / 2.0, modSize.z() / 2.0);
  auto modMat = description.material(mod_x.attr<std::string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);

  modVol.setVisAttributes(description.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));

  //--------------------------Block of fibers in a module------------------------------------------------------

  //block fibers
  auto fiber_block = mod_x.child(_Unicode(block));

  auto fb_sx      = fiber_block.attr<double>(_Unicode(sizeX));
  auto fb_sy      = fiber_block.attr<double>(_Unicode(sizeY));
  auto fb_sz      = fiber_block.attr<double>(_Unicode(sizeZ));
  auto fb_SpaceXY = fiber_block.attr<double>(_Unicode(SpaceXY));

  Position fbSize(fb_sx, fb_sy, fb_sz);

  //fibers
  auto fiber_tube = mod_x.child(_Unicode(fiber));

  auto fr  = fiber_tube.attr<double>(_Unicode(radius));
  auto fsx = fiber_tube.attr<double>(_Unicode(spacex));
  auto fsy = fiber_tube.attr<double>(_Unicode(spacey));

  //fiber block description and placement in module
  Box fbShape(fbSize.x() / 2.0, fbSize.y() / 2.0, fbSize.z() / 2.0);
  Volume fbVol("fiberblock_volume", fbShape, modMat);
  fbVol.setVisAttributes(
      description.visAttributes(mod_x.attr<std::string>(_Unicode(vis)))); //same as module

  int num_fbX = int(modSize.x() / (fbSize.x() + 2.0 * fb_SpaceXY));
  int num_fbY = int(modSize.y() / (fbSize.y() + 2.0 * fb_SpaceXY));

  double fb_xpos0 = -(modSize.x() / 2.0) + (fbSize.x() / 2.0) + fb_SpaceXY;
  double fb_ypos0 = -(modSize.y() / 2.0) + (fbSize.y() / 2.0) + fb_SpaceXY;
  int nblock      = 0;

  for (int iy = 0; iy < num_fbY; iy++) {
    for (int ix = 0; ix < num_fbX; ix++) {
      double fb_pos_x = fb_xpos0 + ix * (fbSize.x() + 2.0 * fb_SpaceXY); //mm
      double fb_pos_y = fb_ypos0 + iy * (fbSize.y() + 2.0 * fb_SpaceXY); //mm
      double fb_pos_z = 0 * mm;

      auto fbPV = modVol.placeVolume(fbVol, nblock, Position{fb_pos_x, fb_pos_y, fb_pos_z});
      fbPV.addPhysVolID("block", nblock++);
    }
  }

  //fiber placement and description in blocks
  auto fiberMat = description.material(fiber_tube.attr<std::string>(_Unicode(material)));
  Tube fiberShape(0., fr, fbSize.z() / 2.0);
  Volume fiberVol("fiber_vol", fiberShape, fiberMat);
  fiberVol.setVisAttributes(description.visAttributes(fiber_tube.attr<std::string>(_Unicode(vis))));
  fiberVol.setSensitiveDetector(sens);

  int num_fX = int(fbSize.x() / (2 * fr + 2.0 * fsx));
  int num_fY = int(fbSize.y() / (2 * fr + 2.0 * fsy));

  double fiber_xpos0 = -(fbSize.x() / 2.0) + fr + fsx;
  double fiber_ypos0 = -(fbSize.y() / 2.0) + fr + fsy;
  int nfibers        = 0;

  //Fiber Holder
  auto fiberholder_x = mod_x.child(_Unicode(fiberholder));
  double fh_dz       = 0.25 * mm; //thickness of fiber holder

  double fh_outerbox_y = 2.0 * fr + 2.0 * fsy;
  double fh_outerbox_x = 2.0 * fr + 2.0 * fsx;
  Box fh_outerbox(fh_outerbox_x / 2.0, fh_outerbox_y / 2.0, fh_dz / 2.0);

  double fh_innerbox_y = 2.0 * fr;
  double fh_innerbox_x = 2.0 * fr;
  Box fh_innerbox(fh_innerbox_x / 2.0, fh_innerbox_y / 2.0, fh_dz / 2.0);

  SubtractionSolid fiberholder_solid(fh_outerbox, fh_innerbox, Position(0.0, 0.0, 0.0));
  auto fiberholderMat = description.material(fiberholder_x.attr<std::string>(_Unicode(material)));
  Volume fiberholderVol("fiberholder_vol", fiberholder_solid, fiberholderMat);
  fiberholderVol.setVisAttributes(
      description.visAttributes(fiberholder_x.attr<std::string>(_Unicode(vis))));

  int nfh = 0;

  //placement of fibers and fiberholder
  for (int iy = 0; iy < num_fY; iy++) {
    for (int ix = 0; ix < num_fX; ix++) {

      double fiber_pos_x = fiber_xpos0 + ix * (2.0 * fr + 2.0 * fsx); //mm
      double fiber_pos_y = fiber_ypos0 + iy * (2.0 * fr + 2.0 * fsy); //mm
      double fiber_pos_z = 0 * mm;                                    //mm

      //placement of fiber
      auto fiberPV =
          fbVol.placeVolume(fiberVol, nfibers++, Position{fiber_pos_x, fiber_pos_y, fiber_pos_z});
      fiberPV.addPhysVolID("fiber_x", ix + 1).addPhysVolID("fiber_y", iy + 1);

      //placement of fiber holder 6.6*cm apart c-to-c
      int num_holders  = 4; // which means 4 regions
      double fh_pos_z0 = -1 * (fbSize.z() / 2.0) + (fh_dz / 2.0);

      for (int iz = 0; iz < num_holders; iz++) {
        double fh_pos_z = fh_pos_z0 + iz * ((fbSize.z() - fh_dz) / (num_holders - 1));
        fbVol.placeVolume(fiberholderVol, nfh++, Position{fiber_pos_x, fiber_pos_y, fh_pos_z});
      } //iz close

    } //ix close
  } //iy close

  return make_tuple(modVol, modSize);

} // build_specScifiCAL_module function close

DECLARE_DETELEMENT(EcalLumiSpecWScFi, create_detector)
