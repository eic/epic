// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Dmitry Romanov, Pu-Kai Wang
//==========================================================================
//  A general implementation for homogeneous calorimeter
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 06/09/2021
//==========================================================================
//==========================================================================
//  Date: 03/10/2022
//  Add the new geometry and the supporting structue
//  Adapted the single module with additional wrraper and supporting structure
//--------------------------------------------------------------------------
//  Date: 20/03/2023
//  Reorganize and optimize the scripts
//  Adapted the inner supporting structure for improving low Q2 measurements
//--------------------------------------------------------------------------
//  Author: WANG Pu-Kai, ZHU Yuwei (IJClab)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "GeometryHelpers.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <tuple>
#include <vector>

using namespace dd4hep;

/** \addtogroup calorimeters Calorimeters
 */

/** \addtogroup Homogeneous Calorimeter
 * \brief Type: **HomogeneousCalorimeter**.
 * \author C. Peng
 * \ingroup calorimeters
 *
 * @{
 */

// headers
static std::tuple<int, std::pair<int, int>> add_12surface_disk(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                               SensitiveDetector& sens, int id);

// helper function to get x, y, z if defined in a xml component
template <class XmlComp>
Position get_xml_xyz(XmlComp& comp, dd4hep::xml::Strng_t name)
{
  Position pos(0., 0., 0.);
  if (comp.hasChild(name)) {
    auto child = comp.child(name);
    pos.SetX(dd4hep::getAttrOrDefault<double>(child, _Unicode(x), 0.));
    pos.SetY(dd4hep::getAttrOrDefault<double>(child, _Unicode(y), 0.));
    pos.SetZ(dd4hep::getAttrOrDefault<double>(child, _Unicode(z), 0.));
  }
  return pos;
}

// main
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string     detName = detElem.nameStr();
  int             detID   = detElem.id();
  DetElement      det(detName, detID);
  sens.setType("calorimeter");

  // assembly
  Assembly assembly(detName);

  // module placement
  xml::Component     plm = detElem.child(_Unicode(placements));


  std::map<int, std::pair<int, int>> sectorModuleRowsColumns;
  auto addRowColumnNumbers = [&sectorModuleRowsColumns](int sector, std::pair<int, int> rowcolumn) {
                         auto it = sectorModuleRowsColumns.find(sector);
                         if (it != sectorModuleRowsColumns.end()) {
                           it->second = rowcolumn;
                         } else {
                           sectorModuleRowsColumns[sector] = rowcolumn;
                         }
                       };


  int sector_id = 1;
  for (xml::Collection_t disk_12surface(plm, _Unicode(disk_12surface)); disk_12surface; ++disk_12surface) {
    auto [sector, rowcolumn] = add_12surface_disk(desc, assembly, disk_12surface, sens, sector_id++);
    addRowColumnNumbers(sector, rowcolumn);
  }

  for (auto [sector, rowcolumn] : sectorModuleRowsColumns) {
    desc.add(Constant(Form((detName + "_NModules_Sector%d").c_str(), sector), std::to_string((rowcolumn.first)), std::to_string((rowcolumn.second)) ));
  }


  // detector position and rotation
  auto         pos       = get_xml_xyz(detElem, _Unicode(position));
  auto         rot       = get_xml_xyz(detElem, _Unicode(rotation));
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr        = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
  PlacedVolume envPV     = motherVol.placeVolume(assembly, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  return det;
}

// helper function to build module with or w/o wrapper
std::tuple<Volume, Position> build_module(Detector& desc, xml::Collection_t& plm, SensitiveDetector& sens)
{
  auto   mod = plm.child(_Unicode(module));
  auto   mx  = mod.attr<double>(_Unicode(modulex));
  auto   my  = mod.attr<double>(_Unicode(moduley));
  auto   mz  = mod.attr<double>(_Unicode(modulez));
  auto   mdz = mod.attr<double>(_Unicode(moduleshift));
  Box    modshape(mx / 2., my / 2., mz / 2.);
  auto   modMat = desc.material(mod.attr<std::string>(_Unicode(gmaterial)));
  Volume modVol("module_vol", modshape, modMat);
  modVol.setVisAttributes(desc.visAttributes(mod.attr<std::string>(_Unicode(vis))));

  auto   cry  = plm.child(_Unicode(crystal));
  auto   cryx = cry.attr<double>(_Unicode(sizex));
  auto   cryy = cry.attr<double>(_Unicode(sizey));
  auto   cryz = cry.attr<double>(_Unicode(sizez));
  Box    crystalshape(cryx / 2., cryy / 2., cryz / 2.);
  auto   crystalMat = desc.material(cry.attr<std::string>(_Unicode(material)));
  Volume crystalVol("crystal_vol", crystalshape, crystalMat);
  modVol.placeVolume(crystalVol, Position(0., 0., -mdz / 2.));
  crystalVol.setVisAttributes(desc.visAttributes(cry.attr<std::string>(_Unicode(cryvis))));
  crystalVol.setSensitiveDetector(sens);

  if (!plm.hasChild(_Unicode(wrapper))){  // no wrapper
    printout(DEBUG, "HomogeneousCalorimeter", "without wrapper");
    return std::make_tuple(modVol, Position{mx, my, mz});
  }
  else{  // build wrapper
    auto wrp              = plm.child(_Unicode(wrapper)); // Read all the contents in the wrapper block
    auto wrapcfthickness = wrp.attr<double>(_Unicode(carbonfiber_thickness));
    auto wrapcflength    = wrp.attr<double>(_Unicode(carbonfiber_length));
    auto wrapVMthickness = wrp.attr<double>(_Unicode(VM2000_thickness));
    auto carbonMat        = desc.material(wrp.attr<std::string>(_Unicode(material_carbon)));
    auto wrpMat           = desc.material(wrp.attr<std::string>(_Unicode(material_wrap)));
    auto gapMat           = desc.material(wrp.attr<std::string>(_Unicode(material_gap)));

    if (wrapcfthickness < 1e-12 * mm)
      return std::make_tuple(modVol, Position{mx, my, mz});

    Box carbonShape(mx / 2., my / 2., wrapcflength / 2.);
    Box carbonShape_sub((mx - 2. * wrapcfthickness) / 2., (my - 2. * wrapcfthickness) / 2., wrapcflength / 2.);
    SubtractionSolid carbon_subtract(carbonShape, carbonShape_sub, Position(0., 0., 0.));

    Box gapShape(mx / 2., my / 2., (cryz - 2. * wrapcflength) / 2.);
    Box gapShape_sub((mx - 2. * wrapcfthickness) / 2., (my - 2. * wrapcfthickness) / 2., (cryz - 2. * wrapcflength) / 2.);
    SubtractionSolid gap_subtract(gapShape, gapShape_sub, Position(0., 0., 0.));

    Box wrpVM2000((mx - 2. * wrapcfthickness) / 2., (my - 2. * wrapcfthickness) / 2., mz / 2.);
    Box wrpVM2000_sub((mx - 2. * wrapcfthickness - 2. * wrapVMthickness) / 2., (my - 2. * wrapcfthickness - 2. * wrapVMthickness) / 2., cryz / 2.);
    SubtractionSolid wrpVM2000_subtract(wrpVM2000, wrpVM2000_sub, Position(0., 0., -mdz / 2.));

    Volume carbonVol("carbon_vol", carbon_subtract, carbonMat);
    Volume gapVol("gap_vol", gap_subtract, gapMat);
    Volume wrpVol("wrapper_vol", wrpVM2000_subtract, wrpMat);

    modVol.placeVolume(carbonVol, Position(0., 0., (cryz - wrapcflength - wrapVMthickness) / 2.)); // put the wrap in the both ends of crystal
    modVol.placeVolume(carbonVol, Position(0., 0., (wrapcflength - cryz - wrapVMthickness) / 2.));
    modVol.placeVolume(gapVol, Position(0., 0., -wrapVMthickness / 2.)); // put the gap between two carbon fiber
    modVol.placeVolume(wrpVol, Position(0., 0., 0.));


    carbonVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_carbon))));
    gapVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_gap))));
    wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_wrap))));

    printout(DEBUG, "HomogeneousCalorimeter", "with wrapper");

    return std::make_tuple(modVol, Position{mx, my, mz});
  }
}

// place 12 surface disk of modules
static std::tuple<int, std::pair<int, int>> add_12surface_disk(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                                               SensitiveDetector& sens, int sid)
{
  auto [modVol, modSize]             = build_module(desc, plm, sens);
  int         sector_id              = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  double      rmin                   = plm.attr<double>(_Unicode(rmin));
  double      rmax                   = plm.attr<double>(_Unicode(rmax));
  double      r12min                 = plm.attr<double>(_Unicode(r12min));
  double      r12max                 = plm.attr<double>(_Unicode(r12max));
  double      structure_frame_length = plm.attr<double>(_Unicode(outerringlength));
  double      calo_module_length     = plm.attr<double>(_Unicode(modulelength));
  double      Prot                   = plm.attr<double>(_Unicode(protate));
  double      Nrot                   = plm.attr<double>(_Unicode(nrotate));
  double      Oring_shift            = plm.attr<double>(_Unicode(outerringshift));
  double      Innera                 = plm.attr<double>(_Unicode(inneradiusa));
  double      Innerb                 = plm.attr<double>(_Unicode(inneradiusb));
  double      phimin                 = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double      phimax                 = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);

  std::vector<double> pt_innerframe_x;  //The points information for inner supporting frame
  std::vector<double> pt_innerframe_y;



  //=========================================================
  // The modules' positions followd the fillRectangles function
  //=========================================================
  float half_modx = modSize.x() * 0.5, half_mody = modSize.y() * 0.5;
  auto points = epic::geo::fillRectangles({half_modx, half_mody}, modSize.x(), modSize.y(), rmin, rmax, phimin, phimax);


  //=========================================================
  // Read the positions information from xml file
  //=========================================================
  xml_coll_t pts_extrudedpolygon(plm, _Unicode(points_extrudedpolygon));
  for (xml_coll_t position_i(pts_extrudedpolygon, _U(position)); position_i; ++position_i){
    xml_comp_t position_comp = position_i;
    pt_innerframe_x.push_back((position_comp.x()));
    pt_innerframe_y.push_back((position_comp.y()));
  }

  xml_coll_t positions_addmodules(plm, _Unicode(addmodulespos));
  for (xml_coll_t position_i(positions_addmodules, _U(position)); position_i; ++position_i){
    xml_comp_t position_comp = position_i;
    auto add_point = epic::geo::Point((position_comp.x()), (position_comp.y()));
    points.push_back(add_point);
  }


  //=========================================================
  // optional envelope volume and the supporting frame
  //=========================================================

  // Material for the structure and mother space
  //
  Material outer_ring_material = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "StainlessSteel"));
  Material inner_ring_material = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Copper"));

  //==============================
  // Outer supporting frame
  //==============================

  PolyhedraRegular solid_ring12(12, r12min, r12max, structure_frame_length);
  Volume           ring12_vol("ring12", solid_ring12, outer_ring_material);
  Transform3D      tr_global_Oring = RotationZYX(Prot, 0., 0.) * Translation3D(0., 0., Oring_shift);
  ring12_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));

  //=============================
  // Inner supporting frame
  //=============================

  // Version3: solid with elliptical inside
  //
  std::vector<double> sec_z  = {-calo_module_length / 2., calo_module_length / 2.};
  std::vector<double> sec_x  = {0., 0.};
  std::vector<double> sec_y  = {0., 0.};
  std::vector<double> zscale = {1., 1.};

  ExtrudedPolygon  inner_support_main(pt_innerframe_x, pt_innerframe_y, sec_z, sec_x, sec_y, zscale);
  EllipticalTube   subtract_a(Innera, Innerb, calo_module_length / 2.);
  SubtractionSolid inner_support_substracta(inner_support_main, subtract_a, Position(0., 0., 0.));
  Volume           inner_support_vol("inner_support_vol", inner_support_substracta, inner_ring_material);
  inner_support_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));
  Transform3D tr_global_Iring_elli = RotationZYX(Nrot, 0., 0.) * Translation3D(0., 0., 0.);

  //=============================
  // The mother volume of modules
  //=============================
  bool             has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  PolyhedraRegular solid_world(12, 0., r12min, calo_module_length);
  EllipticalTube   solid_sub(Innera, Innerb, calo_module_length / 2.);
  Transform3D      subtract_pos = RotationZYX(Nrot, 0., 0.) * Translation3D(0., 0., 0.);
  SubtractionSolid calo_subtract(solid_world, solid_sub, subtract_pos);
  Volume           env_vol(std::string(env.name()) + "_envelope", calo_subtract, outer_ring_material);
  Transform3D      tr_global = RotationZYX(Prot, 0., 0.) * Translation3D(0., 0., 0.);
  env_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_steel_gap))));

  // Place frames and mother volume of modules into the world volume
  //
  if (has_envelope) {
    env.placeVolume(env_vol, tr_global);                          // Place the mother volume for all modules
    env.placeVolume(ring12_vol, tr_global_Oring);                 // Place the outer supporting frame
    env_vol.placeVolume(inner_support_vol, tr_global_Iring_elli); // Place the version3 inner supporting frame
  }




  //=====================================================================
  // Placing The Modules
  // Since the inner and outer porfile is not the circle shape,
  // which means the modules placements can't be simply done by
  // fillRectangles functions. I hardcode the additional positions
  // at backward_PbWO4.xml to fill all the gap between circular
  // placement and supporting structures.
  //=====================================================================

  // retrieve the max and min position of modules
  auto [minl_ptsx, maxl_ptsx] = std::minmax_element(points.begin(), points.end(), [](auto a, auto b){ return a.x() < b.x(); });
  auto [minl_ptsy, maxl_ptsy] = std::minmax_element(points.begin(), points.end(), [](auto a, auto b){ return a.y() < b.y(); });

  // Place the modules, use row and column ID for each modules
  // row and column ID start from the top left corner
  //
  int row = 0, column = 0;
  int N_row =  std::round((maxl_ptsy[0].y() - minl_ptsy[0].y()) / modSize.y());
  int N_column =  std::round((maxl_ptsx[0].x() - minl_ptsx[0].x()) / modSize.x());

  auto rowcolumn = std::make_pair(N_row, N_column);

  for (auto& p : points) {
    column = std::round((p.x() - minl_ptsx[0].x()) / modSize.x());
    row = std::round((maxl_ptsy[0].y() - p.y()) / modSize.y());

    Transform3D tr_local = RotationZYX(Nrot, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("row", row).addPhysVolID("column", column);
  }

  return {sector_id, rowcolumn};
}

//@}
DECLARE_DETELEMENT(epic_HomogeneousCalorimeter, create_detector)
