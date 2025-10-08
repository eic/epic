// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
 *
 * \ingroup trackers
 *
 * @{
 */
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Printout.h>
#include <DD4hep/Shapes.h>
#include <XML/Layering.h>
#include <XML/Utilities.h>

#include <cassert>

using namespace std;
using namespace dd4hep;

namespace {
std::pair<Volume, Transform3D> build_shape(const Detector& descr, const xml_det_t& x_det,
                                           const xml_comp_t& x_support, const xml_comp_t& x_child,
                                           const double offset = 0) {
  // Get Initial rotation/translation info
  xml_dim_t x_pos(x_child.child(_U(position), false));
  xml_dim_t x_rot(x_child.child(_U(rotation), false));
  Position pos3D{0, 0, 0};
  Rotation3D rot3D;

  if (x_rot) {
    rot3D = RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0));
  }
  if (x_pos) {
    pos3D = Position(x_pos.x(0), x_pos.y(0), x_pos.z(0));
  }

  // handle different known shapes and create solids
  Solid solid;
  const std::string type = x_support.attr<std::string>(_U(type));
  if (type == "Tube" || type == "Cylinder") {
    const double thickness = getAttrOrDefault(x_child, _U(thickness), x_support.thickness());
    const double length    = getAttrOrDefault(x_child, _U(length), x_support.length());
    const double rmin      = getAttrOrDefault(x_child, _U(rmin), x_support.rmin()) + offset;
    const double phimin    = getAttrOrDefault(
        x_child, _Unicode(phimin), getAttrOrDefault(x_support, _Unicode(phimin), 0.0 * deg));
    const double phimax = getAttrOrDefault(
        x_child, _Unicode(phimax), getAttrOrDefault(x_support, _Unicode(phimax), 360.0 * deg));
    solid = Tube(rmin, rmin + thickness, length / 2, phimin, phimax);
  } else if (type == "Box") {
    const double box_x = getAttrOrDefault(x_child, _U(x), x_support.x());
    const double box_y = getAttrOrDefault(x_child, _U(y), x_support.y());
    const double box_z = getAttrOrDefault(x_child, _U(z), x_support.z());
    solid              = Box(box_x / 2.0, box_y / 2.0, box_z / 2.0);
  }
  // A disk is a cylinder, constructed differently
  else if (type == "Disk") {
    const double thickness = getAttrOrDefault(x_child, _U(thickness), x_support.thickness());
    const double rmin      = getAttrOrDefault(x_child, _U(rmin), x_support.rmin());
    const double rmax      = getAttrOrDefault(x_child, _U(rmax), x_support.rmax());
    const double phimin    = getAttrOrDefault(
        x_child, _Unicode(phimin), getAttrOrDefault(x_support, _Unicode(phimin), 0.0 * deg));
    const double phimax = getAttrOrDefault(
        x_child, _Unicode(phimax), getAttrOrDefault(x_support, _Unicode(phimax), 360.0 * deg));
    pos3D = pos3D + Position(0, 0, -x_support.thickness() / 2 + thickness / 2 + offset);
    solid = Tube(rmin, rmax, thickness / 2, phimin, phimax);
  } else if (type == "Cone") {
    const double base_rmin1 = getAttrOrDefault(x_child, _U(rmin1), x_support.rmin1());
    const double base_rmin2 = getAttrOrDefault(x_child, _U(rmin2), x_support.rmin2());
    const double length     = getAttrOrDefault(x_child, _U(length), x_support.length());
    // Account for the fact that the distance between base_rmin1 and rmax2 is the projection
    // of the thickness on the transverse direction
    const double thickness = getAttrOrDefault(x_child, _U(thickness), x_support.thickness());
    const double transverse_thickness =
        thickness / cos(atan2(fabs(base_rmin2 - base_rmin1), length));
    // also account that the same is true for the offset
    const double transverse_offset = offset / cos(atan2(fabs(base_rmin2 - base_rmin1), length));
    const double rmin1             = base_rmin1 + transverse_offset;
    const double rmin2             = base_rmin2 + transverse_offset;
    const double rmax1             = rmin1 + transverse_thickness;
    const double rmax2             = rmin2 + transverse_thickness;
    // Allow for optional hard rmin/rmax cutoffs
    const double rmin = getAttrOrDefault(
        x_child, _U(rmin), getAttrOrDefault(x_support, _Unicode(rmin), min(rmin1, rmin2)));
    const double rmax = getAttrOrDefault(
        x_child, _U(rmax), getAttrOrDefault(x_support, _Unicode(rmax), max(rmax1, rmax2)));
    if (rmin > min(rmax1, rmax2)) {
      printout(ERROR, x_det.nameStr(),
               "%s: rmin (%f mm) must be smaller than the smallest rmax (%f %f mm)",
               x_support.nameStr().c_str(), rmin / mm, rmax1 / mm, rmax2 / mm);
      std::exit(1);
    }
    if (rmax < max(base_rmin1, base_rmin2)) {
      printout(ERROR, x_det.nameStr(),
               "%s: rmax (%f mm) must be larger than the largest rmin (%f %f mm)",
               x_support.nameStr().c_str(), rmax / mm, base_rmin1 / mm, base_rmin2 / mm);
      std::exit(1);
    }
    const double zmin  = -length / 2 + length * (rmin - rmin1) / (rmin2 - rmin1);
    const double zmax  = -length / 2 + length * (rmax - rmax1) / (rmax2 - rmax1);
    const auto rmin_at = [&](const double z) {
      return rmin1 + (z + length / 2) * (rmin2 - rmin1) / length;
    };
    const auto rmax_at = [&](const double z) {
      return rmax1 + (z + length / 2) * (rmax2 - rmax1) / length;
    };
    // Allow for optional phimin/phimax
    const double phimin = getAttrOrDefault<double>(
        x_child, _Unicode(phimin), getAttrOrDefault(x_support, _Unicode(phimin), 0.0 * deg));
    const double phimax = getAttrOrDefault<double>(
        x_child, _Unicode(phimax), getAttrOrDefault(x_support, _Unicode(phimax), 360.0 * deg));
    const double deltaphi = phimax - phimin;
    const double epsilon{TGeoShape::Tolerance()};
    if (fabs(zmin) >= length / 2 - epsilon && fabs(zmax) >= length / 2 - epsilon) {
      if (fabs(phimax - phimin - 360 * deg) < epsilon) {
        solid = Cone(length / 2, rmin1, rmax1, rmin2, rmax2);
      } else {
        solid = ConeSegment(length / 2, rmin1, rmax1, rmin2, rmax2, phimin, phimax);
      }
    } else {
      std::vector<double> v_rmin{max(rmin1, rmin), max(rmin2, rmin)},
          v_rmax{min(rmax1, rmax), min(rmax2, rmax)}, v_z{-length / 2, +length / 2};
      for (const auto& z :
           (zmin < zmax ? std::vector<double>{zmin, zmax} : std::vector<double>{zmax, zmin})) {
        if (-length / 2 + epsilon < z && z < -epsilon + length / 2) {
          v_rmin.insert(std::prev(v_rmin.end()), std::max(rmin, rmin_at(z)));
          v_rmax.insert(std::prev(v_rmax.end()), std::min(rmax, rmax_at(z)));
          v_z.insert(std::prev(v_z.end()), z);
        }
      }
      solid = Polycone(phimin, deltaphi, v_rmin, v_rmax, v_z);
    }
  } else if (type == "Disk") {
    const double thickness = getAttrOrDefault(x_child, _U(thickness), x_support.thickness());
    const double rmin      = getAttrOrDefault(x_child, _U(rmin), x_support.rmin());
    const double rmax      = getAttrOrDefault(x_child, _U(rmax), x_support.rmax());
    const double phimin    = getAttrOrDefault(
        x_child, _Unicode(phimin), getAttrOrDefault(x_support, _Unicode(phimin), 0.0 * deg));
    const double phimax = getAttrOrDefault(
        x_child, _Unicode(phimax), getAttrOrDefault(x_support, _Unicode(phimax), 360.0 * deg));
    pos3D = pos3D + Position(0, 0, -x_support.thickness() / 2 + thickness / 2 + offset);
    solid = Tube(rmin, rmax, thickness / 2, phimin, phimax);
  } else {
    printout(ERROR, x_det.nameStr(), "Unknown support type: %s", type.c_str());
    std::exit(1);
  }
  // Materials
  Material mat = descr.material(getAttrOrDefault<std::string>(x_child, _U(material), "Air"));
  // Create our volume
  Volume vol{getAttrOrDefault<std::string>(x_child, _U(name), "support_vol"), solid, mat};

  // Create full transformation
  Transform3D tr(rot3D, pos3D);

  // visualization?
  if (x_child.hasAttr(_U(vis))) {
    vol.setVisAttributes(descr.visAttributes(x_child.visStr()));
  }
  return {vol, tr};
}

// Function to create a subtraction of two shapes
std::pair<Volume, Transform3D> build_subtraction(const Detector& descr, const xml_det_t& x_det,
                                                 const xml_comp_t& x_subtraction) {
  // Get the two shapes to subtract
  xml_comp_t x_shape1 = x_subtraction.child(_Unicode(shape1));
  xml_comp_t x_shape2 = x_subtraction.child(_Unicode(shape2));

  // Build the primary shape (without transformation)
  auto [vol1, tr1_unused] = build_shape(descr, x_det, x_shape1, x_shape1);
  Solid solid1            = vol1.solid();

  // Build the shape to subtract (without transformation)
  auto [vol2, tr2_unused] = build_shape(descr, x_det, x_shape2, x_shape2);
  Solid solid2            = vol2.solid();

  // Get relative transformation for shape2 (relative to shape1's center)
  xml_dim_t x_pos2(x_shape2.child(_U(position), false));
  xml_dim_t x_rot2(x_shape2.child(_U(rotation), false));
  Position pos2{0, 0, 0};
  Rotation3D rot2;

  if (x_rot2) {
    rot2 = RotationZYX(x_rot2.z(0), x_rot2.y(0), x_rot2.x(0));
  }
  if (x_pos2) {
    pos2 = Position(x_pos2.x(0), x_pos2.y(0), x_pos2.z(0));
  }

  // Create the subtraction solid
  Transform3D tr_relative(rot2, pos2);
  SubtractionSolid subtracted_solid(solid1, solid2, tr_relative);

  // Create volume with the subtracted solid
  Material mat = descr.material(getAttrOrDefault<std::string>(x_subtraction, _U(material), "Air"));
  std::string vol_name = getAttrOrDefault<std::string>(x_subtraction, _U(name), "subtraction_vol");
  Volume vol{vol_name, subtracted_solid, mat};

  // Get overall position/rotation for the subtracted volume
  xml_dim_t x_pos(x_subtraction.child(_U(position), false));
  xml_dim_t x_rot(x_subtraction.child(_U(rotation), false));
  Position pos3D{0, 0, 0};
  Rotation3D rot3D;

  if (x_rot) {
    rot3D = RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0));
  }
  if (x_pos) {
    pos3D = Position(x_pos.x(0), x_pos.y(0), x_pos.z(0));
  }

  Transform3D tr(rot3D, pos3D);

  // Set visualization if specified
  if (x_subtraction.hasAttr(_U(vis))) {
    vol.setVisAttributes(descr.visAttributes(x_subtraction.visStr()));
  }

  // Debug output
  printout(DEBUG, "SupportServiceMaterial", "Created subtraction volume: %s", vol_name.c_str());

  return {vol, tr};
}

std::pair<Volume, Transform3D> build_shape(const Detector& descr, const xml_det_t& x_det,
                                           const xml_comp_t& x_support, const double offset = 0) {
  return build_shape(descr, x_det, x_support, x_support, offset);
}
} // namespace

/** Generic tracker support implementation, can consist of arbitrary shapes
 *
 * @author Sylvester Joosten
 */
static Ref_t create_SupportServiceMaterial(Detector& description, xml_h e,
                                           [[maybe_unused]] SensitiveDetector sens) {
  const xml_det_t x_det = e;
  const int det_id      = x_det.id();
  const string det_name = x_det.nameStr();

  // global z-offset for the entire support assembly
  const double offset = getAttrOrDefault(x_det, _U(offset), 0.);

  DetElement det(det_name, det_id);
  Assembly assembly(det_name + "_assembly");

  // Loop over the supports
  for (xml_coll_t su{x_det, _U(support)}; su; ++su) {
    xml_comp_t x_sup   = su;
    const double rot_z = getAttrOrDefault(x_sup, _U(phi0), 0.);
    RotationZ rot(rot_z);
    Transform3D tr_rot(
        rot, Position(0, 0, 0)); // additional rotation of the module after position offset

    auto [vol, tr]           = build_shape(description, x_det, x_sup);
    [[maybe_unused]] auto pv = assembly.placeVolume(vol, tr_rot * tr);
    // Loop over support components, if any
    double cumulative_thickness = 0;
    for (xml_coll_t com{x_sup, _U(component)}; com; ++com) {
      xml_comp_t x_com = com;
      auto [cvol, ctr] = build_shape(description, x_det, x_sup, x_com, cumulative_thickness);
      vol.placeVolume(cvol, ctr);
      cumulative_thickness += x_com.thickness();
    }
  }
  // Loop over any subtraction volumes
  for (xml_coll_t sub{x_det, _Unicode(subtraction)}; sub; ++sub) {
    xml_comp_t x_sub = sub;
    auto [vol, tr]   = build_subtraction(description, x_det, x_sub);
    assembly.placeVolume(vol, tr);
  }

  // final placement
  Volume motherVol = description.pickMotherVolume(det);
  Position pos(0, 0, offset);
  PlacedVolume pv = motherVol.placeVolume(assembly, pos);
  pv.addPhysVolID("system", det.id());
  det.setPlacement(pv);

  return det;
}

// clang-format off
DECLARE_DETELEMENT(epic_SupportServiceMaterial, create_SupportServiceMaterial)
