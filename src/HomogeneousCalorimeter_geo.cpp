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

#include "DD4hep/DetFactoryHelper.h"
#include "GeometryHelpers.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <tuple>

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
  auto   sx  = mod.attr<double>(_Unicode(sizex));
  auto   sy  = mod.attr<double>(_Unicode(sizey));
  auto   sz  = mod.attr<double>(_Unicode(sizez));
  Box    modShape(sx / 2., sy / 2., sz / 2.);
  auto   modMat = desc.material(mod.attr<std::string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);
  modVol.setSensitiveDetector(sens);
  modVol.setVisAttributes(desc.visAttributes(mod.attr<std::string>(_Unicode(vis))));

  // no wrapper
  if (!plm.hasChild(_Unicode(wrapper))) {
    return std::make_tuple(modVol, Position{sx, sy, sz});
    // build wrapper
  } else {
    auto wrp       = plm.child(_Unicode(wrapper));
    auto thickness = wrp.attr<double>(_Unicode(thickness));
    if (thickness < 1e-12 * mm) {
      return std::make_tuple(modVol, Position{sx, sy, sz});
    }
    auto   wrpMat = desc.material(wrp.attr<std::string>(_Unicode(material)));
    Box    wrpShape((sx + thickness) / 2., (sy + thickness) / 2., sz / 2.);
    Volume wrpVol("wrapper_vol", wrpShape, wrpMat);
    wrpVol.placeVolume(modVol, Position(0., 0., 0.));
    wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis))));
    return std::make_tuple(wrpVol, Position{sx + thickness, sy + thickness, sz});
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
  if (has_envelope) {
    env.placeVolume(env_vol, tr_global);
  }

  // local placement of modules
  int  mid    = 0;
  auto points = ecce::geo::fillRectangles({0., 0.}, modSize.x(), modSize.y(), rmin, rmax, phimin, phimax);
  for (auto& p : points) {
    Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
    auto modPV = (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
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
#ifdef EPIC_ECCE_LEGACY_COMPAT
DECLARE_DETELEMENT(ecce_HomogeneousCalorimeter, create_detector)
#endif
DECLARE_DETELEMENT(epic_HomogeneousCalorimeter, create_detector)
