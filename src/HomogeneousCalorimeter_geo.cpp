// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Dmitry Romanov, Pu-Kai Wang
//==========================================================================
//  A general implementation for homogeneous calorimeter
//  it supports three types of placements
//  1. Individual module placement with module dimensions and positions
//  2. Array placement with module dimensions and numbers of row and columns
//  3. Disk placement with module dimensions and (Rmin, Rmax), and (Phimin, Phimax)
//  4. Lines placement with module dimensions and (mirrorx, mirrory)
//     (NOTE: anchor point is the 0th block of the line instead of line center)
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 06/09/2021
//==========================================================================
//==========================================================================
//  Add the new geometry and the supporting structue
//  Adapted the single module with additional wrraper and supporting structure
//--------------------------------------------------------------------------
//  Author: WANG Pu-Kai, ZHU Yuwei (IJClab)
//  Date: 03/10/2022
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
 *
 * \code
 *   <detector id="1" name="HyCal" type="HomogeneousCalorimeter" readout="EcalHits">
 *     <position x="0" y="0" z="0"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <array nrow="34" ncol="34"
 *              envelope="false" sector="1">
 *         <position x="0" y="0" z="-9.73*cm"/>
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="18*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *         <removal row="16" col="16"/>
 *         <removal row="16" col="17"/>
 *         <removal row="17" col="16"/>
 *         <removal row="17" col="17"/>
 *       </array>
 *       <array nrow="6" ncol="24" sector="2">
 *         <position x="-17*(2.05+0.015)*cm+12*(3.8+0.015)*cm" y="17*(2.05+0.015)*cm+3*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *       <array nrow="24" ncol="6" sector="3">
 *         <position x="17*(2.05+0.015)*cm+3*(3.8+0.015)*cm" y="17*(2.05+0.015)*cm-12*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *       <array nrow="6" ncol="24" sector="4">
 *         <position x="17*(2.05+0.015)*cm-12*(3.8+0.015)*cm" y="-17*(2.05+0.015)*cm-3*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *       <array nrow="24" ncol="6" sector="3">
 *         <position x="-17*(2.05+0.015)*cm-3*(3.8+0.015)*cm" y="-17*(2.05+0.015)*cm+12*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *     </placements>
 *   </detector>
 *
 *   <detector id="2" name="SomeBlocks" type="HomogeneousCalorimeter" readout="EcalHits">
 *     <position x="0" y="0" z="30*cm"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <individuals sector="1"/>
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="20*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *         <placement x="1*cm" y="1*cm" z="0" id="1"/>
 *         <placement x="-1*cm" y="1*cm" z="0" id="2"/>
 *         <placement x="1*cm" y="-1*cm" z="0" id="3"/>
 *         <placement x="-1*cm" y="-1*cm" z="0" id="4"/>
 *       </individuals>
 *     </placements>
 *   </detector>
 *
 *   <detector id="2" name="DiskShapeCalorimeter" type="HomogeneousCalorimeter" readout="EcalHits">
 *     <position x="0" y="0" z="-30*cm"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <disk rmin="25*cm" rmax="125*cm" length="20.5*cm" phimin="0" phimax="360*degree"
 *             envelope="false" sector="1">
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="20*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </disk>
 *     </placements>
 *   </detector>
 *
 *   <detector id="3" name="SomeLines" type="HomogeneousCalorimeter" readout="EcalHits">
 *     <position x="0" y="0" z="60*cm"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <lines sector="1" mirrorx="true" mirrory="true"/>
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="20*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *         <line x="10.25*mm" y="10.25*mm" axis="x" begin="8" nmods="16"/>
 *         <line x="10.25*mm" y="30.75*mm" axis="y" begin="8" nmods="16"/>
 *         <line x="10.25*mm" y="51.25*mm" axis="z" begin="8" nmods="16"/>
 *       </individuals>
 *     </placements>
 *   </detector>
 * \endcode
 *
 * @{
 */

// headers
static std::tuple<int, int> add_individuals(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                            SensitiveDetector& sens, int id);
static std::tuple<int, int> add_array(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                      int id);
static std::tuple<int, int> add_disk(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                     int id);
static std::tuple<int, int> add_12surface_disk(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                               SensitiveDetector& sens, int id);
static std::tuple<int, int> add_lines(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                      int id);

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
  std::map<int, int> sectorModuleNumbers;
  auto               addModuleNumbers = [&sectorModuleNumbers](int sector, int nmod) {
    auto it = sectorModuleNumbers.find(sector);
    if (it != sectorModuleNumbers.end()) {
      it->second += nmod;
    } else {
      sectorModuleNumbers[sector] = nmod;
    }
  };
  int sector_id = 1;
  for (xml::Collection_t mod(plm, _Unicode(individuals)); mod; ++mod) {
    auto [sector, nmod] = add_individuals(desc, assembly, mod, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
  for (xml::Collection_t arr(plm, _Unicode(array)); arr; ++arr) {
    auto [sector, nmod] = add_array(desc, assembly, arr, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
  for (xml::Collection_t disk(plm, _Unicode(disk)); disk; ++disk) {
    auto [sector, nmod] = add_disk(desc, assembly, disk, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
  for (xml::Collection_t disk_12surface(plm, _Unicode(disk_12surface)); disk_12surface; ++disk_12surface) {
    auto [sector, nmod] = add_12surface_disk(desc, assembly, disk_12surface, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }
  for (xml::Collection_t lines(plm, _Unicode(lines)); lines; ++lines) {
    auto [sector, nmod] = add_lines(desc, assembly, lines, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }

  for (auto [sector, nmods] : sectorModuleNumbers) {
    desc.add(Constant(Form((detName + "_NModules_Sector%d").c_str(), sector), std::to_string(nmods)));
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
  auto   sx  = mod.attr<double>(_Unicode(gx));
  auto   sy  = mod.attr<double>(_Unicode(gy));
  auto   sz  = mod.attr<double>(_Unicode(gz));
  auto   sdz = mod.attr<double>(_Unicode(gdz));
  Box    modshape(sx / 2., sy / 2., (sz + sdz) / 2.);
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
  modVol.placeVolume(crystalVol, Position(0., 0., -sdz / 2.));
  crystalVol.setVisAttributes(desc.visAttributes(cry.attr<std::string>(_Unicode(cryvis))));
  crystalVol.setSensitiveDetector(sens);

  if (!plm.hasChild(_Unicode(wrapper))) // no wrapper
  {
    printout(DEBUG, "HomogeneousCalorimeter", "without wrapper");

    return std::make_tuple(modVol, Position{sx, sy, sz});

  } else                                                  // build wrapper
  {
    auto wrp              = plm.child(_Unicode(wrapper)); // Read all the content in the wrapper block
    auto carbon_thickness = wrp.attr<double>(_Unicode(carbon_thickness));
    auto wrap_thickness   = wrp.attr<double>(_Unicode(wrap_thickness));
    auto length_wrapper   = wrp.attr<double>(_Unicode(length));
    auto carbonMat        = desc.material(wrp.attr<std::string>(_Unicode(material_carbon)));
    auto wrpMat           = desc.material(wrp.attr<std::string>(_Unicode(material_wrap)));
    auto gapMat           = desc.material(wrp.attr<std::string>(_Unicode(material_gap)));

    if (carbon_thickness < 1e-12 * mm)
      return std::make_tuple(modVol, Position{sx, sy, sz});

    Box carbonShape(sx / 2., sy / 2., length_wrapper / 2.);
    Box carbonShape_sub((sx - 2. * carbon_thickness) / 2., (sy - 2. * carbon_thickness) / 2., length_wrapper / 2.);
    SubtractionSolid carbon_subtract(carbonShape, carbonShape_sub, Position(0., 0., 0.));

    Box wrpShape((sx - 2. * carbon_thickness) / 2., (sy - 2. * carbon_thickness) / 2., (sz + wrap_thickness) / 2.);
    Box wrpShape_sub((sx - 2. * carbon_thickness - 2. * wrap_thickness) / 2.,
                     (sy - 2. * carbon_thickness - 2. * wrap_thickness) / 2., sz / 2.);
    SubtractionSolid wrp_subtract(wrpShape, wrpShape_sub, Position(0., 0., -wrap_thickness));

    Box              gapShape(sx / 2., sy / 2., (sz - 2. * length_wrapper) / 2.);
    Box              gapShape_sub((sx - 2. * carbon_thickness) / 2., (sy - 2. * carbon_thickness) / 2.,
                                  (sz - 2. * length_wrapper) / 2.);
    SubtractionSolid gap_subtract(gapShape, gapShape_sub, Position(0., 0., 0.));

    Volume carbonVol("carbon_vol", carbon_subtract, carbonMat);
    Volume wrpVol("wrapper_vol", wrp_subtract, wrpMat);
    Volume gapVol("gap_vol", gap_subtract, gapMat);

    modVol.placeVolume(
        carbonVol,
        Position(0., 0., (sz - length_wrapper - wrap_thickness) / 2.)); // put the wrap in the both ends of crystal
    modVol.placeVolume(carbonVol, Position(0., 0., (length_wrapper - sz - wrap_thickness) / 2.));
    modVol.placeVolume(wrpVol, Position(0., 0., wrap_thickness / 2.));
    modVol.placeVolume(gapVol, Position(0., 0., -wrap_thickness / 2.)); // put the gap in the middle of crystal

    carbonVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_carbon))));
    wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_wrap))));
    gapVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_gap))));

    printout(DEBUG, "HomogeneousCalorimeter", "with wrapper");

    return std::make_tuple(modVol, Position{sx, sy, sz});
  }
}

// place modules, id must be provided
static std::tuple<int, int> add_individuals(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                            SensitiveDetector& sens, int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int sector_id          = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int nmodules           = 0;
  for (xml::Collection_t pl(plm, _Unicode(placement)); pl; ++pl) {
    Position    pos(dd4hep::getAttrOrDefault<double>(pl, _Unicode(x), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(y), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(z), 0.));
    Position    rot(dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotx), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(roty), 0.),
                    dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotz), 0.));
    auto        mid   = pl.attr<int>(_Unicode(id));
    Transform3D tr    = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
    auto        modPV = env.placeVolume(modVol, tr);
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", mid);
    nmodules++;
  }

  return {sector_id, nmodules};
}

// place array of modules
static std::tuple<int, int> add_array(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                      int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int sector_id          = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int id_begin           = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  int nrow               = plm.attr<int>(_Unicode(nrow));
  int ncol               = plm.attr<int>(_Unicode(ncol));

  // compute array position
  double begx = -modSize.x() * ncol / 2. + modSize.x() / 2.;
  double begy = modSize.y() * nrow / 2. - modSize.y() / 2.;

  std::vector<std::pair<int, int>> removals;
  // get the removal list
  for (xml::Collection_t rm(plm, _Unicode(removal)); rm; ++rm) {
    removals.push_back({rm.attr<int>(_Unicode(row)), rm.attr<int>(_Unicode(col))});
  }

  // placement inside mother
  auto pos = get_xml_xyz(plm, _Unicode(position));
  auto rot = get_xml_xyz(plm, _Unicode(rotation));

  // optional envelope volume
  bool        has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  Material    material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Air"));
  Box         solid(ncol * modSize.x() / 2.0, nrow * modSize.y() / 2.0, modSize.z() / 2.0);
  Volume      env_vol(std::string(env.name()) + "_envelope", solid, material);
  Transform3D tr_global = RotationZYX(rot.z(), rot.y(), rot.x()) * Translation3D(pos.x(), pos.y(), pos.z());
  if (has_envelope) {
    env.placeVolume(env_vol, tr_global);
  }

  // local placement of modules
  int nmodules = 0;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (std::find(removals.begin(), removals.end(), std::pair<int, int>(i, j)) != removals.end()) {
        continue;
      }
      double      px       = begx + modSize.x() * j;
      double      py       = begy - modSize.y() * i;
      Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(px, py, 0.0);
      auto        modPV =
          (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", i * ncol + j + id_begin);
      nmodules++;
    }
  }
  return {sector_id, nmodules};
}

// place disk of modules
static std::tuple<int, int> add_disk(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                     int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int    sector_id       = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int    id_begin        = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  double rmin            = plm.attr<double>(_Unicode(rmin));
  double rmax            = plm.attr<double>(_Unicode(rmax));
  double phimin          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double phimax          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);

  // placement inside mother
  auto pos = get_xml_xyz(plm, _Unicode(position));
  auto rot = get_xml_xyz(plm, _Unicode(rotation));

  // optional envelope volume
  bool        has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  Material    material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Air"));
  Tube        solid(rmin, rmax, modSize.z() / 2.0, phimin, phimax);
  Volume      env_vol(std::string(env.name()) + "_envelope", solid, material);
  Transform3D tr_global = RotationZYX(rot.z(), rot.y(), rot.x()) * Translation3D(pos.x(), pos.y(), pos.z());
  if (has_envelope)
    env.placeVolume(env_vol, tr_global);

  // local placement of modules
  int  mid    = 0;
  auto points = epic::geo::fillRectangles({0., 0.}, modSize.x(), modSize.y(), rmin, rmax, phimin, phimax);
  for (auto& p : points) {
    Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
  }

  return {sector_id, mid};
}

// place 12 surface disk of modules
static std::tuple<int, int> add_12surface_disk(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                               SensitiveDetector& sens, int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int    sector_id       = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int    id_begin        = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  double rmin            = plm.attr<double>(_Unicode(rmin));
  double rmax            = plm.attr<double>(_Unicode(rmax));
  double r12min          = plm.attr<double>(_Unicode(r12min));
  double r12max          = plm.attr<double>(_Unicode(r12max));
  double structure_frame_length = plm.attr<double>(_Unicode(SFlength));
  double calo_module_length = plm.attr<double>(_Unicode(CMlength));
  double phimin          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double phimax          = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);





  //=========================================================
  // optional envelope volume and the supporting frame
  //=========================================================

  // Material for the structure and mother space
  //
  Material outer_ring_material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "StainlessSteel"));
  Material inner_ring_material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Copper"));
  // Material hole_material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Vacuum"));




  //==============================
  // Outer supporting frame
  //==============================

  PolyhedraRegular solid_ring12(12, r12min, r12max, structure_frame_length);
  Volume ring12_vol("ring12", solid_ring12, outer_ring_material);
  Transform3D tr_global_Oring = RotationZYX(15.*degree, 0., 0.) * Translation3D(0., 0., -20.*cm);
  ring12_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));




  //=============================
  // Inner supporting frame
  //=============================

  // // Version1: circle shape
  // //
  // Tube        Ssolid_ring12(8.5*cm, rmin, calo_module_length/2., phimin, phimax);
  // Volume      Sring12_vol("Sring12", Ssolid_ring12, inner_ring_material);
  // Sring12_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));
  // Transform3D tr_global_Iring = RotationZYX(0., 0., 0.) * Translation3D(0., 0., 0.);


  // // Version2: elliptical shape
  // //
  // EllipticalTube   inner_elliptical_support_a(9.*cm, 7.5*cm, calo_module_length/2.);
  // EllipticalTube   inner_elliptical_support_b(8.5*cm, 7.*cm, calo_module_length/2.);
  // SubtractionSolid inner_elliptical_support_substract(inner_elliptical_support_a, inner_elliptical_support_b, Position(0., 0., 0.));
  // Volume           inner_elliptical_vol("inner_elliptical", inner_elliptical_support_substract, inner_ring_material);
  // inner_elliptical_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));
  // Transform3D tr_global_Iring_elli = RotationZYX(-15.*degree, 0., 0.) * Translation3D(2.05*cm, 0., 0.);


  // Version3: solid with elliptical inside
  //
  Box inner_support_main(8.2*cm, 6.15*cm, calo_module_length/2.);  // Original size
  // Box inner_support_main(8.2*cm, 5.125*cm, calo_module_length/2.);  // Adapted size
  Box subtract_corner(1.025*cm, 1.025*cm, calo_module_length/2.);
  EllipticalTube   subtract_a(7.5*cm, 5.5*cm, calo_module_length/2.);
  SubtractionSolid inner_support_substractb1(inner_support_main, subtract_corner, Position(7.175*cm, 5.125*cm, 0.));
  SubtractionSolid inner_support_substractb2(inner_support_substractb1, subtract_corner, Position(7.175*cm, -5.125*cm, 0.));
  SubtractionSolid inner_support_substractb3(inner_support_substractb2, subtract_corner, Position(-7.175*cm, 5.125*cm, 0.));
  SubtractionSolid inner_support_substractb4(inner_support_substractb3, subtract_corner, Position(-7.175*cm, -5.125*cm, 0.));
  SubtractionSolid inner_support_substracta(inner_support_substractb4, subtract_a, Position(0., 0., 0.));
  Volume           inner_support_vol("inner_support_vol", inner_support_substracta, inner_ring_material);
  inner_support_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));
  Transform3D tr_global_Iring_elli = RotationZYX(-15.*degree, 0., 0.) * Translation3D(0., 0., 0.);


  // // The vacuum inside the inner structure
  // //
  // EllipticalTube   inner_elliptical_vacuum(7.5*cm, 4.5*cm, calo_module_length/2.);
  // Volume           inner_elliptical_vacuum_vol("inner_elliptical_vacuum_vol", inner_elliptical_vacuum, hole_material);
  // inner_elliptical_vacuum_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));
  // Transform3D tr_global_Iring_elli_vacuum = RotationZYX(-15.*degree, 0., 0.) * Translation3D(1.025*cm, 0., 0.);


  //=============================
  // The mother volume of modules
  //=============================
  bool has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  PolyhedraRegular solid_world(12, 0., r12min, calo_module_length);
  EllipticalTube  solid_sub(7.5*cm, 5.5*cm, calo_module_length/2.);
  Transform3D subtract_pos = RotationZYX(-15.*degree, 0., 0.) * Translation3D(0., 0., 0.);
  SubtractionSolid calo_subtract(solid_world, solid_sub, subtract_pos);
  Volume      env_vol(std::string(env.name()) + "_envelope", calo_subtract, outer_ring_material);
  Transform3D tr_global = RotationZYX(15.*degree, 0., 0.) * Translation3D(0., 0., 0.);
  env_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_steel_gap))));



  // // =============================
  // // Supporting frame for the cabling
  // // =============================
  // // Having overlapped with tracker barrelendcapsupporting structure, comment here until the size and position determined
  // //
  // PolyhedraRegular cabling_support12(12, r12max, 70.*cm, 2.54*cm);
  // Volume cabling_support12_V("cabling_support12_V", cabling_support12, inner_ring_material);
  // Transform3D tr_global_csfront = RotationZYX(15.*degree, 0., 0.) * Translation3D(0., 0., 9.*cm);
  // Transform3D tr_global_csback = RotationZYX(15.*degree, 0., 0.) * Translation3D(0., 0., -49.*cm);

  // Box hole(13.5/2.*cm, 2.54/2.*cm, 2.54/2.*cm);
  // Volume hole_V("hole_V", hole, hole_material);
  // Transform3D hole_pos = RotationZYX(0., 0., 0.) * Translation3D(0., 0., 0.);

  // Place frames and mother volume of modules into the world volume
  //
  if (has_envelope)
    {
      env.placeVolume(env_vol, tr_global);                          // Place the mother volume for all modules
      env.placeVolume(ring12_vol, tr_global_Oring);                 // Place the outer supporting frame
      env_vol.placeVolume(inner_support_vol, tr_global_Iring_elli);  // Place the version3 inner supporting frame


      // env.placeVolume(cabling_support12_V, tr_global_csfront);  // Place the front cabling frame(Only the holes are visible)
      // env.placeVolume(cabling_support12_V, tr_global_csback);   // Place the back cabling frame(Only the holes are visible)

      // for(int i = 0 ; i < 12 ; i++)
      //   {
      //     hole_pos = RotationZYX((15. + i * 30.)*degree, 0., 0.) * Translation3D(82.5*mm, 675.*mm, 0.);
      //     cabling_support12_V.placeVolume(hole_V, hole_pos);
      //     hole_pos = RotationZYX((15. + i * 30.)*degree, 0., 0.) * Translation3D(-82.5*mm, 675.*mm, 0.);
      //     cabling_support12_V.placeVolume(hole_V, hole_pos);
      //   }
    }







  //=====================================================================
  // Placing The Modules
  // Since the inner and outer porfile is not the circle shape,
  // which means the modules placements can't be simply done by
  // fillRectangles functions. I hardcode the additional positions
  // to fill all the gap between circular placement and supporting
  // structures.
  //=====================================================================

  // Add the modules followd the fillRectangles function
  //
  int mid = 0, total_id = 0;
  float half_modx = modSize.x() * 0.5, half_mody = modSize.y() * 0.5;
  auto points = epic::geo::fillRectangles({half_modx, half_mody}, modSize.x(), modSize.y(), rmin, rmax, phimin, phimax);
  for (auto& p : points)
    {
      Transform3D tr_local = RotationZYX(-15.*degree, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
      auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", total_id);
      total_id = id_begin + mid++;
    }

  // Add the modules manually in the gap [Outer]
  //
  const int add_N_mod_outer = 56;
  float addX_outer[add_N_mod_outer] = {9.225, 11.275, 13.325, 15.375, 17.425, 19.475, 44.075, 46.125, 60.475, 60.475, 62.525, 62.525, 62.525, 62.525,  -9.225, -11.275, -13.325, -15.375, -17.425, -19.475, -44.075, -46.125, -60.475, -60.475, -62.525, -62.525, -62.525, -62.525,  -9.225, -11.275, -13.325, -15.375, -17.425, -19.475, -44.075, -46.125, -60.475, -60.475, -62.525, -62.525, -62.525, -62.525,  9.225, 11.275, 13.325, 15.375, 17.425, 19.475, 44.075, 46.125, 60.475, 60.475, 62.525, 62.525, 62.525, 62.525};
  float addY_outer[add_N_mod_outer] = {62.525, 62.525, 62.525, 62.525, 60.475, 60.475, 46.125, 44.075, 19.475, 17.425, 15.375, 13.325, 11.275, 9.225,  62.525, 62.525, 62.525, 62.525, 60.475, 60.475, 46.125, 44.075, 19.475, 17.425, 15.375, 13.325, 11.275, 9.225,  -62.525, -62.525, -62.525, -62.525, -60.475, -60.475, -46.125, -44.075, -19.475, -17.425, -15.375, -13.325, -11.275, -9.225, -62.525, -62.525, -62.525, -62.525, -60.475, -60.475, -46.125, -44.075, -19.475, -17.425, -15.375, -13.325, -11.275, -9.225};

  for (int im = 0; im < add_N_mod_outer; im++)
    {
      Transform3D add_local = RotationZYX(-15.*degree, 0.0, 0.0) * Translation3D(addX_outer[im] * cm, addY_outer[im] * cm, 0.0);
      auto modPV = (has_envelope ? env_vol.placeVolume(modVol, add_local) : env.placeVolume(modVol, tr_global * add_local));
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", total_id);
      total_id++;
    }

  // Add the modules manually in the gap [Inner]
  //
  const int add_N_mod_inner = 36;
  double addX_inner[add_N_mod_inner] = {1.025,1.025,3.075,3.075,5.125,7.175,7.175,9.225,9.225, -1.025,-1.025,-3.075,-3.075,-5.125,-7.175,-7.175,-9.225,-9.225, -1.025,-1.025,-3.075,-3.075,-5.125,-7.175,-7.175,-9.225,-9.225, 1.025,1.025,3.075,3.075,5.125,7.175,7.175,9.225,9.225};
  double addY_inner[add_N_mod_inner] = {7.175,9.225,7.175,9.225,7.175,7.175,5.125,1.025,3.075, 7.175,9.225,7.175,9.225,7.175,7.175,5.125,1.025,3.075, -7.175,-9.225,-7.175,-9.225,-7.175,-7.175,-5.125,-1.025,-3.075, -7.175,-9.225,-7.175,-9.225,-7.175,-7.175,-5.125,-1.025,-3.075};

  for (int im = 0; im < add_N_mod_inner; im++)
    {
      Transform3D add_local = RotationZYX(-15.*degree, 0.0, 0.0) * Translation3D(addX_inner[im] * cm, addY_inner[im] * cm, 0.0);
      auto modPV = (has_envelope ? env_vol.placeVolume(modVol, add_local) : env.placeVolume(modVol, tr_global * add_local));
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", total_id);
      total_id++;
    }



  return {sector_id, mid};
}

// place lines of modules (anchor point is the 0th module of this line)
static std::tuple<int, int> add_lines(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens,
                                      int sid)
{
  auto [modVol, modSize] = build_module(desc, plm, sens);
  int  sector_id         = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int  id_begin          = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  bool mirrorx           = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(mirrorx), false);
  bool mirrory           = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(mirrory), false);
  bool mirrorz           = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(mirrorz), false);

  // line placement
  int mid = 0;
  for (xml::Collection_t pl(plm, _Unicode(line)); pl; ++pl) {
    Position pos(dd4hep::getAttrOrDefault<double>(pl, _Unicode(x), 0.),
                 dd4hep::getAttrOrDefault<double>(pl, _Unicode(y), 0.),
                 dd4hep::getAttrOrDefault<double>(pl, _Unicode(z), 0.));
    Position rot(dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotx), 0.),
                 dd4hep::getAttrOrDefault<double>(pl, _Unicode(roty), 0.),
                 dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotz), 0.));
    // determine axis
    std::string axis = dd4hep::getAttrOrDefault<std::string>(pl, _Unicode(axis), "x");
    std::transform(axis.begin(), axis.end(), axis.begin(), [](char c) { return std::tolower(c); });
    if ((axis != "x") && (axis != "y") && (axis != "z")) {
      std::cerr << "HomogeneousCalorimeter Error: cannot determine axis of line " << axis
                << ", abort placement of this line." << std::endl;
      continue;
    }

    int begin = dd4hep::getAttrOrDefault<int>(pl, _Unicode(begin), 0);
    int nmods = pl.attr<int>(_Unicode(nmods));

    std::vector<Position> trans;
    for (int i = 0; i < nmods; ++i) {
      Position tran{(axis == "x") ? pos.x() + (begin + i) * modSize.x() : pos.x(),
                    (axis == "y") ? pos.y() + (begin + i) * modSize.y() : pos.y(),
                    (axis == "z") ? pos.z() + (begin + i) * modSize.z() : pos.z()};
      trans.push_back(tran);
    }

    // mirror placement
    if (mirrorx) {
      size_t curr_size = trans.size();
      for (size_t i = 0; i < curr_size; ++i) {
        trans.push_back(Position{-trans[i].x(), trans[i].y(), trans[i].z()});
      }
    }
    if (mirrory) {
      size_t curr_size = trans.size();
      for (size_t i = 0; i < curr_size; ++i) {
        trans.push_back(Position{trans[i].x(), -trans[i].y(), trans[i].z()});
      }
    }
    if (mirrorz) {
      size_t curr_size = trans.size();
      for (size_t i = 0; i < curr_size; ++i) {
        trans.push_back(Position{trans[i].x(), trans[i].y(), -trans[i].z()});
      }
    }

    // place volume
    for (auto& p : trans) {
      Transform3D tr    = RotationZYX(rot.z(), rot.y(), rot.x()) * Translation3D(p.x(), p.y(), p.z());
      auto        modPV = env.placeVolume(modVol, tr);
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
    }
  }
  return {sector_id, mid};
}
//@}
DECLARE_DETELEMENT(epic_HomogeneousCalorimeter, create_detector)
