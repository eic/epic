// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Dmitry Romanov, Pu-Kai Wang
//==========================================================================
//  A general implementation for homogeneous calorimeter
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
 * @{
 */

// headers
static std::tuple<int, int> add_12surface_disk(Detector& desc, Assembly& env, xml::Collection_t& plm,
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
  for (xml::Collection_t disk_12surface(plm, _Unicode(disk_12surface)); disk_12surface; ++disk_12surface) {
    auto [sector, nmod] = add_12surface_disk(desc, assembly, disk_12surface, sens, sector_id++);
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

// place 12 surface disk of modules
static std::tuple<int, int> add_12surface_disk(Detector& desc, Assembly& env, xml::Collection_t& plm,
                                               SensitiveDetector& sens, int sid)
{
  auto [modVol, modSize]             = build_module(desc, plm, sens);
  int         sector_id              = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  // int         id_begin               = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  double      rmin                   = plm.attr<double>(_Unicode(rmin));
  double      rmax                   = plm.attr<double>(_Unicode(rmax));
  double      r12min                 = plm.attr<double>(_Unicode(r12min));
  double      r12max                 = plm.attr<double>(_Unicode(r12max));
  double      structure_frame_length = plm.attr<double>(_Unicode(SFlength));
  double      calo_module_length     = plm.attr<double>(_Unicode(CMlength));
  double      NEEMC_Prot             = plm.attr<double>(_Unicode(NEEMC_PR));
  double      NEEMC_Nrot             = plm.attr<double>(_Unicode(NEEMC_NR));
  double      NEEMC_OR_shift         = plm.attr<double>(_Unicode(NEEMC_OR_RS));
  double      NEEMC_IR_a             = plm.attr<double>(_Unicode(Inner_a));
  double      NEEMC_IR_b             = plm.attr<double>(_Unicode(Inner_b));
  double      phimin                 = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double      phimax                 = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);
  std::string polygonX               = plm.attr<std::string>(_Unicode(ptsX_extrudedpolygon));
  std::string polygonY               = plm.attr<std::string>(_Unicode(ptsY_extrudedpolygon));
  std::string iposx                  = plm.attr<std::string>(_Unicode(inner_outer_add_posx));
  std::string iposy                  = plm.attr<std::string>(_Unicode(inner_outer_add_posy));

  //=========================================================
  // optional envelope volume and the supporting frame
  //=========================================================

  // Material for the structure and mother space
  //
  Material outer_ring_material = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "StainlessSteel"));
  Material inner_ring_material = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Copper"));
  // Material hole_material     = desc.material(getAttrOrDefault<std::string>(plm, _U(material), "Vacuum"));

  //==============================
  // Outer supporting frame
  //==============================

  PolyhedraRegular solid_ring12(12, r12min, r12max, structure_frame_length);
  Volume           ring12_vol("ring12", solid_ring12, outer_ring_material);
  Transform3D      tr_global_Oring = RotationZYX(NEEMC_Prot, 0., 0.) * Translation3D(0., 0., NEEMC_OR_shift);
  ring12_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));

  //=============================
  // Inner supporting frame
  //=============================

  // Version3: solid with elliptical inside
  //
  std::vector<double> pt_x;
  std::vector<double> pt_y;
  std::string         delimiter = " ";
  size_t              pos       = 0;
  std::string         token;
  while ((pos = polygonX.find(delimiter)) != std::string::npos) {
    token = polygonX.substr(0, pos);
    pt_x.push_back(atof(token.c_str()));
    polygonX.erase(0, pos + delimiter.length());
  }
  pt_x.push_back(atof(polygonX.c_str()));

  pos = 0;
  while ((pos = polygonY.find(delimiter)) != std::string::npos) {
    token = polygonY.substr(0, pos);
    pt_y.push_back(atof(token.c_str()));
    polygonY.erase(0, pos + delimiter.length());
  }
  pt_y.push_back(atof(polygonY.c_str()));
  std::vector<double> sec_z  = {-calo_module_length / 2., calo_module_length / 2.};
  std::vector<double> sec_x  = {0., 0.};
  std::vector<double> sec_y  = {0., 0.};
  std::vector<double> zscale = {1., 1.};

  ExtrudedPolygon  inner_support_main(pt_x, pt_y, sec_z, sec_x, sec_y, zscale);
  EllipticalTube   subtract_a(NEEMC_IR_a, NEEMC_IR_b, calo_module_length / 2.);
  SubtractionSolid inner_support_substracta(inner_support_main, subtract_a, Position(0., 0., 0.));
  Volume           inner_support_vol("inner_support_vol", inner_support_substracta, inner_ring_material);
  inner_support_vol.setVisAttributes(desc.visAttributes(plm.attr<std::string>(_Unicode(vis_struc))));
  Transform3D tr_global_Iring_elli = RotationZYX(NEEMC_Nrot, 0., 0.) * Translation3D(0., 0., 0.);

  //=============================
  // The mother volume of modules
  //=============================
  bool             has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
  PolyhedraRegular solid_world(12, 0., r12min, calo_module_length);
  EllipticalTube   solid_sub(NEEMC_IR_a, NEEMC_IR_b, calo_module_length / 2.);
  Transform3D      subtract_pos = RotationZYX(NEEMC_Nrot, 0., 0.) * Translation3D(0., 0., 0.);
  SubtractionSolid calo_subtract(solid_world, solid_sub, subtract_pos);
  Volume           env_vol(std::string(env.name()) + "_envelope", calo_subtract, outer_ring_material);
  Transform3D      tr_global = RotationZYX(NEEMC_Prot, 0., 0.) * Translation3D(0., 0., 0.);
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

  // Add the modules followd the fillRectangles function
  //
  int   mid = 0;
  // int total_id = 0;
  float half_modx = modSize.x() * 0.5, half_mody = modSize.y() * 0.5;
  auto points = epic::geo::fillRectangles({half_modx, half_mody}, modSize.x(), modSize.y(), rmin, rmax, phimin, phimax);
  float min_ptsx = 0., max_ptsx = 0., min_ptsy = 0., max_ptsy = 0.; // Find the min/max position X and Y

  // Segment the position list[string] and save them as vector
  //
  std::vector<double> inner_outer_posx;
  std::vector<double> inner_outer_posy;

  pos = 0;
  while ((pos = iposx.find(delimiter)) != std::string::npos) {
    token = iposx.substr(0, pos);
    inner_outer_posx.push_back(atof(token.c_str()));
    iposx.erase(0, pos + delimiter.length());
  }
  inner_outer_posx.push_back(atof(iposx.c_str()));

  pos = 0;
  while ((pos = iposy.find(delimiter)) != std::string::npos) {
    token = iposy.substr(0, pos);
    inner_outer_posy.push_back(atof(token.c_str()));
    iposy.erase(0, pos + delimiter.length());
  }
  inner_outer_posy.push_back(atof(iposy.c_str()));

  for (auto& p : points) {
    if (p.x() < min_ptsx)
      min_ptsx = p.x();
    if (p.x() > max_ptsx)
      max_ptsx = p.x();
    if (p.y() < min_ptsy)
      min_ptsy = p.y();
    if (p.y() > max_ptsy)
      max_ptsy = p.y();
  }


  // Place the modules, use row and column ID for each modules
  // row and column ID start from the top left corner
  //
  int row = 0, column = 0;
  for (auto& p : points) {
    column = (p.x() - min_ptsx) / modSize.x();
    row = (max_ptsy - p.y()) / modSize.y();
    Transform3D tr_local = RotationZYX(NEEMC_Nrot, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("row", row).addPhysVolID("column", column);
  }

  int im = 0;
  for (auto&& value : inner_outer_posx) {
    column = (value - min_ptsx) / modSize.x();
    row = (max_ptsy - inner_outer_posy[im]) / modSize.y();
    Transform3D add_local = RotationZYX(NEEMC_Nrot, 0.0, 0.0) * Translation3D(value * cm, inner_outer_posy[im] * cm, 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, add_local) : env.placeVolume(modVol, tr_global * add_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("row", row).addPhysVolID("column", column);
    im++;
  }

  return {sector_id, mid};
}

//@}
DECLARE_DETELEMENT(epic_HomogeneousCalorimeter, create_detector)
