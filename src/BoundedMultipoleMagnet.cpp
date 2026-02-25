// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Wouter Deconinck

/** \addtogroup Fields
 *  \brief Multipole magnet field with an axis-aligned bounding-box (AABB) pre-filter.
 *
 *  Drop-in replacement for DD4hep's built-in `MultipoleMagnet` field type.
 *  The upstream factory parses `<position>`, `<rotation>`, `<shape>`, and
 *  `<coefficient>` from the compact XML unchanged.  Before delegating to the
 *  full tube-containment check and multipole calculation, a cheap rectangular
 *  AABB test in global coordinates is applied so that the expensive per-step
 *  work is skipped for the vast majority of tracks that are nowhere near the
 *  far-forward magnets.
 *
 *  Usage – replace `type="MultipoleMagnet"` with `type="BoundedMultipoleMagnet"`:
 *  \code
 *  <field name="B0PF_Magnet" type="BoundedMultipoleMagnet">
 *    <position x="..." y="0" z="..."/>
 *    <rotation x="0" y="..." z="0"/>
 *    <shape type="Tube" rmin="0" rmax="..." dz="..."/>
 *    <coefficient coefficient="..." skew="..."/>
 *  </field>
 *  \endcode
 */

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/FieldTypes.h>
#include <DD4hep/Plugins.h>
#include <DD4hep/Printout.h>
#include <DD4hep/Shapes.h>

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace dd4hep;

// ---------------------------------------------------------------------------

/** Multipole field that skips evaluation outside a world-axis-aligned bounding box.
 *
 *  Inherits all multipole state from dd4hep::MultipoleField.  The six AABB
 *  limits are set by the factory function and never change afterwards.
 */
class BoundedMultipoleField : public MultipoleField {
public:
  double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};

  BoundedMultipoleField() : MultipoleField() {}

  virtual void fieldComponents(const double* pos, double* field) override {
    if (pos[0] < xmin || pos[0] > xmax || pos[1] < ymin || pos[1] > ymax || pos[2] < zmin ||
        pos[2] > zmax)
      return;
    MultipoleField::fieldComponents(pos, field);
  }
};

// ---------------------------------------------------------------------------

static Ref_t create_bounded_multipole_magnet(Detector& det, xml::Handle_t handle) {

  // 1. Let the upstream MultipoleMagnet factory handle all XML parsing
  //    (position, rotation, shape, coefficients → transform/inverse/rotation/volume/B_z).
  NamedObject* obj = PluginService::Create<NamedObject*>("MultipoleMagnet", &det, &handle);
  auto* mf         = dynamic_cast<MultipoleField*>(obj);
  if (!mf)
    throw std::runtime_error(
        "BoundedMultipoleMagnet: upstream MultipoleMagnet factory returned null or wrong type");

  // 2. Derive a tight AABB in world (global) coordinates from the parsed state.
  //
  //    The transform translation is in DD4hep/Geant4 units (mm).
  //    TGeoShape accessors return ROOT units (cm); multiply by dd4hep::cm (= 10) to convert to mm.
  //
  //    For a cylinder with half-length dz, radius r, and world-frame axis d:
  //      AABB half-extent_i = dz*|d_i| + r*sqrt(1 - d_i^2)
  //
  //    For a box with half-extents (dx, dy, dz_b) and rotation matrix R:
  //      AABB half-extent_i = |R[i,0]|*dx + |R[i,1]|*dy + |R[i,2]|*dz_b

  double xmin, xmax, ymin, ymax, zmin, zmax;

  if (mf->volume.isValid()) {
    Position center;
    mf->transform.GetTranslation(center); // world-frame center, mm

    double hx, hy, hz;
    std::string shape_class = mf->volume.ptr()->ClassName();

    if (shape_class == "TGeoTubeSeg" || shape_class == "TGeoTube") {
      // Tube: AABB via cylinder extent formula
      ROOT::Math::XYZVector axis =
          mf->rotation * ROOT::Math::XYZVector(0, 0, 1); // world-frame tube axis (unit vector)

      // Shape dimensions: ROOT cm -> DD4hep mm
      Tube tube_shape(mf->volume);
      double r  = tube_shape.rMax() * dd4hep::cm;
      double dz = tube_shape.dZ() * dd4hep::cm;

      auto half_extent = [&](double d_i) {
        return dz * std::abs(d_i) + r * std::sqrt(std::max(0.0, 1.0 - d_i * d_i));
      };
      hx = half_extent(axis.X());
      hy = half_extent(axis.Y());
      hz = half_extent(axis.Z());

    } else if (shape_class == "TGeoBBox") {
      // Box: AABB via rotation-matrix abs-sum formula
      Box box_shape(mf->volume);
      double dx   = box_shape.x() * dd4hep::cm; // half-length x, mm
      double dy   = box_shape.y() * dd4hep::cm; // half-length y, mm
      double dz_b = box_shape.z() * dd4hep::cm; // half-length z, mm

      double xx, xy, xz, yx, yy, yz, zx, zy, zz;
      mf->rotation.GetComponents(xx, xy, xz, yx, yy, yz, zx, zy, zz);
      hx = std::abs(xx) * dx + std::abs(xy) * dy + std::abs(xz) * dz_b;
      hy = std::abs(yx) * dx + std::abs(yy) * dy + std::abs(yz) * dz_b;
      hz = std::abs(zx) * dx + std::abs(zy) * dy + std::abs(zz) * dz_b;

    } else {
      throw std::runtime_error("BoundedMultipoleMagnet: unsupported shape type '" + shape_class +
                               "' for field '" + xml_comp_t(handle).nameStr() +
                               "'; only Tube and Box are supported");
    }

    xmin = center.X() - hx;
    xmax = center.X() + hx;
    ymin = center.Y() - hy;
    ymax = center.Y() + hy;
    zmin = center.Z() - hz;
    zmax = center.Z() + hz;

    printout(DEBUG, "BoundedMultipoleMagnet",
             "%s  AABB  x[%.1f,%.1f]  y[%.1f,%.1f]  z[%.1f,%.1f] mm",
             xml_comp_t(handle).nameStr().c_str(), xmin, xmax, ymin, ymax, zmin, zmax);
  } else {
    // No bounding shape: disable the pre-filter (field valid everywhere, same as upstream)
    printout(WARNING, "BoundedMultipoleMagnet",
             "%s has no bounding shape; AABB pre-filter disabled",
             xml_comp_t(handle).nameStr().c_str());
    constexpr double inf = std::numeric_limits<double>::infinity();
    xmin = ymin = zmin = -inf;
    xmax = ymax = zmax = +inf;
  }

  // 3. Move parsed data into a BoundedMultipoleField.
  //    All MultipoleField public members are copied; the private `flag` and
  //    `translation` members are left at their zero-initialised defaults and
  //    will be lazily populated on the first fieldComponents() call.
  auto* bounded        = new BoundedMultipoleField();
  bounded->field_type  = mf->field_type;
  bounded->coefficents = mf->coefficents;
  bounded->skews       = mf->skews;
  bounded->volume      = mf->volume;
  bounded->transform   = mf->transform;
  bounded->inverse     = mf->inverse;
  bounded->rotation    = mf->rotation;
  bounded->B_z         = mf->B_z;
  delete mf; // temporary; not registered with the Detector

  bounded->xmin = xmin;
  bounded->xmax = xmax;
  bounded->ymin = ymin;
  bounded->ymax = ymax;
  bounded->zmin = zmin;
  bounded->zmax = zmax;

  CartesianField field;
  field.assign(bounded, xml_comp_t(handle).nameStr(), "BoundedMultipoleMagnet");
  return field;
}

DECLARE_XMLELEMENT(BoundedMultipoleMagnet, create_bounded_multipole_magnet)
