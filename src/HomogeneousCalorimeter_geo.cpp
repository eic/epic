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

  auto   roc  = plm.child(_Unicode(readout));
  auto   PCBx = roc.attr<double>(_Unicode(PCB_sizex));
  auto   PCBy = roc.attr<double>(_Unicode(PCB_sizex));
  auto   PCBz = roc.attr<double>(_Unicode(PCB_thickness));
  Box    PCBshape(PCBx / 2., PCBy / 2., PCBz / 2.);
  auto   PCBMat = desc.material(roc.attr<std::string>(_Unicode(material)));
  Volume PCBVol("PCB_vol", PCBshape, PCBMat);
  
  auto   sensorx = roc.attr<double>(_Unicode(Sensor_sizex));
  auto   sensory = roc.attr<double>(_Unicode(Sensor_sizey));
  auto   sensorz = roc.attr<double>(_Unicode(Sensor_thickness));
  auto   sensorNx = roc.attr<int>(_Unicode(Nsensor_X));
  auto   sensorNy = roc.attr<int>(_Unicode(Nsensor_Y));
  Box    sensorshape(sensorx / 2., sensory / 2., sensorz / 2.);
  auto   sensorMat = desc.material(roc.attr<std::string>(_Unicode(material)));
  Volume sensorVol("sensor_vol", sensorshape, sensorMat);
  
  modVol.placeVolume(crystalVol, Position(0., 0., (PCBz + sensorz - mdz) / 2.));
  crystalVol.setVisAttributes(desc.visAttributes(cry.attr<std::string>(_Unicode(cryvis))));
  crystalVol.setSensitiveDetector(sens);

  modVol.placeVolume(PCBVol, Position(0., 0., (PCBz - mz) / 2.));

  auto intervalx = (PCBx - sensorNx * sensorx) / (sensorNx - 1);
  auto intervaly = (PCBy - sensorNy * sensory) / (sensorNy - 1);
  auto x0 = (sensorx - PCBx) / 2.;
  auto y0 = (sensory - PCBy) / 2.;
  for(int i = 0 ; i < sensorNx ; i++)
    for(int j = 0 ; j < sensorNy ; j++)
      modVol.placeVolume(sensorVol, Position(x0 + (sensorx + intervalx) * i, y0 + (sensory + intervaly) * j, PCBz + (sensorz - mz) / 2.));
    

  
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

    modVol.placeVolume(carbonVol, Position(0., 0., sensorz + PCBz + (cryz - wrapcflength - wrapVMthickness) / 2.)); // put the wrap in the both ends of crystal
    modVol.placeVolume(carbonVol, Position(0., 0., sensorz + PCBz + (wrapcflength - cryz - wrapVMthickness) / 2.));
    modVol.placeVolume(gapVol, Position(0., 0., sensorz + PCBz - (wrapVMthickness / 2.))); // put the gap between two carbon fiber
    modVol.placeVolume(wrpVol, Position(0., 0., sensorz + PCBz));


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
  double half_modx = modSize.x() * 0.5, half_mody = modSize.y() * 0.5;


  //=========================================================
  // Read the positions information from xml file
  //=========================================================
  xml_coll_t pts_extrudedpolygon(plm, _Unicode(points_extrudedpolygon));
  for (xml_coll_t position_i(pts_extrudedpolygon, _U(position)); position_i; ++position_i){
    xml_comp_t position_comp = position_i;
    pt_innerframe_x.push_back((position_comp.x()));
    pt_innerframe_y.push_back((position_comp.y()));
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
  //=====================================================================

  auto test_points = epic::geo::fillRectangles({half_modx, half_mody}, modSize.x(), modSize.y(), 0., (rmax/std::cos(Prot)), phimin, phimax);

  std::pair<double, double> c1 (0., 0.);
  auto polyVertex = epic::geo::getPolygonVertices(c1, (rmax/std::cos(Prot)), M_PI/12., 12);
  std::vector<epic::geo::Point> Out_vertices, In_vertices;
  int index = 0;
  for( auto p : polyVertex ){
    epic::geo::Point a = {p.first, p.second};
    Out_vertices.push_back(a);
    // std::cout << "{" << Out_vertices[index].x() << ", " << Out_vertices[index].y() << "}, " << std::endl;
    index++;
  }
  // std::cout << std::endl << std::endl;

  for (xml_coll_t position_i(pts_extrudedpolygon, _U(position)); position_i; ++position_i){
    xml_comp_t position_comp = position_i;
    epic::geo::Point inpt = {position_comp.x(), position_comp.y()};
    In_vertices.push_back(inpt);
  }

  double minX = 0., maxX = 0., minY = 0., maxY = 0.;   
  for (auto &square : test_points) {
    epic::geo::Point box[4] = {{square.x() + half_modx, square.y() + half_mody}, {square.x() - half_modx, square.y() + half_mody}, {square.x() - half_modx, square.y() - half_mody}, {square.x() + half_modx, square.y() - half_mody}};
    if (epic::geo::isBoxInsidePolygon(box, Out_vertices)) {      
      if( square.x() < minX ) minX = square.x();
      if( square.y() < minY ) minY = square.x();
      if( square.x() > maxX ) maxX = square.x();
      if( square.y() > maxY ) maxY = square.x();     
    }          
  }

  int total_count = 0;
  int row = 0, column = 0;
  int N_row = std::round((maxY - minY) / modSize.y());
  int N_column = std::round((maxX - minX) / modSize.x());
  auto rowcolumn = std::make_pair(N_row, N_column);

  for (auto &square : test_points) {
    epic::geo::Point box[4] = {{square.x() + half_modx, square.y() + half_mody}, {square.x() - half_modx, square.y() + half_mody}, {square.x() - half_modx, square.y() - half_mody}, {square.x() + half_modx, square.y() - half_mody}};
    
    if (epic::geo::isBoxInsidePolygon(box, Out_vertices)) {      
      if(!epic::geo::isBoxInsidePolygon(box, In_vertices)) {
        // std::cout << "Square inside the inner frame: " << square.x() << ", " << square.y() << std::endl;
        column = std::round((square.x() - minX) / modSize.x());
        row = std::round((maxY - square.y()) / modSize.y());
        Transform3D tr_local = RotationZYX(Nrot, 0.0, 0.0) * Translation3D(square.x(), square.y(), 0.0);
        auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
        modPV.addPhysVolID("sector", sector_id).addPhysVolID("row", row).addPhysVolID("column", column);
        total_count++;
      }
    }
    
  }

  std::cout << "New Min: " << minX << ", " << minY << " || Max: " << maxX << ", " << maxY << std::endl;
  std::cout << "Total count new: " << total_count << std::endl;

  return {sector_id, rowcolumn};
}

//@}
DECLARE_DETELEMENT(epic_HomogeneousCalorimeter, create_detector)
