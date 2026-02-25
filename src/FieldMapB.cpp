// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck, Dhevan Gangadharan, Aranya Giri

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/FieldTypes.h>
#include <DD4hep/Printout.h>
#include <DD4hep/DD4hepUnits.h>
#include <XML/Utilities.h>
#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/array.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/morton.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/parameter_pack.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <variant>
namespace fs = std::filesystem;

#include "FileLoaderHelper.h"

// Strided (default, matches on-disk file format)
using fieldBrBz_strided_t =
    covfie::field<covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::size2, covfie::backend::array<covfie::vector::float2>>>>>;
using fieldBxByBz_strided_t =
    covfie::field<covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::size3, covfie::backend::array<covfie::vector::float3>>>>>;

// Morton (space-filling curve layout for improved cache locality)
using fieldBrBz_morton_t =
    covfie::field<covfie::backend::affine<covfie::backend::linear<covfie::backend::morton<
        covfie::vector::size2, covfie::backend::array<covfie::vector::float2>>>>>;
using fieldBxByBz_morton_t =
    covfie::field<covfie::backend::affine<covfie::backend::linear<covfie::backend::morton<
        covfie::vector::size3, covfie::backend::array<covfie::vector::float3>>>>>;

// Variant types (monostate = not yet loaded)
using fieldBrBz_t   = std::variant<std::monostate, fieldBrBz_strided_t, fieldBrBz_morton_t>;
using fieldBxByBz_t = std::variant<std::monostate, fieldBxByBz_strided_t, fieldBxByBz_morton_t>;

using namespace dd4hep;

//  Two coordinate options:
//  coord_type="BrBz"
//    expected fields contained within <dimensions>:
//          <R step="" min="" max=""/>
//          <Z step="" min="" max=""/>
//  coord_type="BxByBz"
//    expected fields contained within <dimensions>:
//          <X step="" min="" max=""/>
//          <Y step="" min="" max=""/>
//          <Z step="" min="" max=""/>

// implementation of the field map
class FieldMapB : public dd4hep::CartesianField::Object {

  enum FieldCoord { BrBz, BxByBz };
  enum class BackendType { strided, morton };

public:
  FieldMapB(const std::string& field_type_str   = "magnetic",
            const std::string& coord_type_str   = "BrBz",
            const std::string& backend_type_str = "strided");

  void LoadMap(const std::string& map_file, float scale);

  void SetCoordTranslation(const Transform3D& tr) {
    coordTranslate     = tr;
    coordTranslate_inv = tr.Inverse();
  }
  void SetFieldRotation(const Transform3D& tr) { fieldRot = tr; }

  virtual void fieldComponents(const double* pos, double* field);

private:
  FieldCoord fieldCoord;                                   // field coordinate type
  BackendType backendType{BackendType::strided};           // memory layout backend
  std::optional<Transform3D> coordTranslate{std::nullopt}; // coord translation
  std::optional<Transform3D> coordTranslate_inv{std::nullopt};
  std::optional<Transform3D> fieldRot{std::nullopt}; // field rotation

  fieldBrBz_t fieldBrBz;
  fieldBxByBz_t fieldBxByBz;

  // Coordinate bounds extracted from the field map after loading.
  // For BrBz: bmin[0]=r_min, bmin[1]=z_min, bmax[0]=r_max, bmax[1]=z_max
  // For BxByBz: bmin[0..2]=x/y/z_min, bmax[0..2]=x/y/z_max
  std::array<float, 3> bmin{}, bmax{};
  float fieldScale{1.0f}; // scale factor applied to field values in fieldComponents
};

// constructor
FieldMapB::FieldMapB(const std::string& field_type_str, const std::string& coord_type_str,
                     const std::string& backend_type_str) {
  std::string ftype = field_type_str;
  for (auto& c : ftype) {
    c = tolower(c);
  }

  // set type
  if (ftype == "magnetic") {
    field_type = CartesianField::MAGNETIC;
  } else if (ftype == "electric") {
    field_type = CartesianField::ELECTRIC;
  } else {
    field_type = CartesianField::UNKNOWN;
    printout(ERROR, "FieldMap", "Unknown field type " + ftype);
  }

  if (coord_type_str.compare("BrBz") == 0) { // BrBz
    fieldCoord = FieldCoord::BrBz;
  } else { // BxByBz
    fieldCoord = FieldCoord::BxByBz;
  }

  if (backend_type_str == "morton") {
    backendType = BackendType::morton;
  } else {
    backendType = BackendType::strided;
  }
}

// load data
void FieldMapB::LoadMap(const std::string& map_file, float scale) {
  std::ifstream input(map_file);
  if (!input) {
    printout(ERROR, "FieldMapB", "FieldMapB Error: file " + map_file + " cannot be read.");
    std::_Exit(EXIT_FAILURE);
  }

  switch (fieldCoord) {
  case FieldCoord::BrBz: {
    fieldBrBz_strided_t strided(input);
    switch (backendType) {
    case BackendType::morton:
      fieldBrBz = fieldBrBz_morton_t(strided);
      break;
    default:
      fieldBrBz = std::move(strided);
      break;
    }
    break;
  }
  case FieldCoord::BxByBz: {
    fieldBxByBz_strided_t strided(input);
    switch (backendType) {
    case BackendType::morton:
      fieldBxByBz = fieldBxByBz_morton_t(strided);
      break;
    default:
      fieldBxByBz = std::move(strided);
      break;
    }
    break;
  }
  }

  // Store scale factor for application in fieldComponents.
  fieldScale = scale;

  // Extract coordinate bounds from the affine matrix and array sizes so that
  // fieldComponents can return 0 for out-of-bounds queries.
  auto extract_bounds = [&](auto& f) {
    using T = std::decay_t<decltype(f)>;
    if constexpr (!std::is_same_v<T, std::monostate>) {
      const auto& mat         = f.backend().get_configuration();
      const auto& sizes       = f.backend().get_backend().get_backend().get_configuration();
      constexpr std::size_t N = std::decay_t<decltype(sizes)>::dimensions;
      for (std::size_t i = 0; i < N; ++i) {
        bmin[i] = -mat(i, N) / mat(i, i);
        bmax[i] = bmin[i] + static_cast<float>(sizes[i] - 1) / mat(i, i);
      }
    }
  };
  if (fieldCoord == FieldCoord::BrBz)
    std::visit(extract_bounds, fieldBrBz);
  else
    std::visit(extract_bounds, fieldBxByBz);

  // Print extent and central field value via fieldComponents to verify units.
  if (fieldCoord == FieldCoord::BrBz) {
    const float z_c = 0.5f * (bmin[1] + bmax[1]);
    // Query at r=0 (on-axis), mid-z, in world coords (undo translation).
    auto origin   = coordTranslate.has_value()
                        ? coordTranslate.value() * ROOT::Math::XYZPoint(0, 0, z_c)
                        : ROOT::Math::XYZPoint(0, 0, z_c);
    double pos[3] = {origin.x(), origin.y(), origin.z()};
    double B[3]   = {0, 0, 0};
    fieldComponents(pos, B);
    // Convert back from DD4hep units to Tesla for display.
    printout(INFO, "FieldMapB",
             "Loaded BrBz map: r=[%.1f,%.1f] cm, z=[%.1f,%.1f] cm; "
             "B(r=0,z=%.1f) = (Bz=%.4f) T",
             bmin[0], bmax[0], bmin[1], bmax[1], z_c, B[2] / dd4hep::tesla);
  } else {
    const float x_c = 0.5f * (bmin[0] + bmax[0]);
    const float y_c = 0.5f * (bmin[1] + bmax[1]);
    const float z_c = 0.5f * (bmin[2] + bmax[2]);
    auto origin     = coordTranslate.has_value()
                          ? coordTranslate.value() * ROOT::Math::XYZPoint(x_c, y_c, z_c)
                          : ROOT::Math::XYZPoint(x_c, y_c, z_c);
    double pos[3]   = {origin.x(), origin.y(), origin.z()};
    double B[3]     = {0, 0, 0};
    fieldComponents(pos, B);
    printout(INFO, "FieldMapB",
             "Loaded BxByBz map: x=[%.1f,%.1f] cm, y=[%.1f,%.1f] cm, z=[%.1f,%.1f] cm; "
             "B(x=%.1f,y=%.1f,z=%.1f) = (Bx=%.4f, By=%.4f, Bz=%.4f) T",
             bmin[0], bmax[0], bmin[1], bmax[1], bmin[2], bmax[2], x_c, y_c, z_c,
             B[0] / dd4hep::tesla, B[1] / dd4hep::tesla, B[2] / dd4hep::tesla);
  }
}

// get field components
void FieldMapB::fieldComponents(const double* pos, double* field) {
  // coordinate conversion
  auto p = coordTranslate_inv.has_value()
               ? coordTranslate_inv.value() * ROOT::Math::XYZPoint(pos[0], pos[1], pos[2])
               : ROOT::Math::XYZPoint(pos[0], pos[1], pos[2]);

  if (fieldCoord == FieldCoord::BrBz) {
    // coordinates conversion
    const float r   = sqrt(p.x() * p.x() + p.y() * p.y());
    const float z   = p.z();
    const float phi = atan2(p.y(), p.x());

    // Return 0 field outside map extent
    if (r < bmin[0] || r > bmax[0] || z < bmin[1] || z > bmax[1])
      return;

    std::visit(
        [&](auto& f) {
          using T = std::decay_t<decltype(f)>;
          if constexpr (!std::is_same_v<T, std::monostate>) {
            typename T::view_t view(f);
            typename T::output_t b = view.at(r, z);
            const float Br = b[0], Bz = b[1];
            auto B = fieldRot.has_value()
                         ? fieldRot.value() * ROOT::Math::XYZPoint(Br * cos(phi), Br * sin(phi), Bz)
                         : ROOT::Math::XYZPoint(Br * cos(phi), Br * sin(phi), Bz);
            field[0] += B.x() * dd4hep::tesla * fieldScale;
            field[1] += B.y() * dd4hep::tesla * fieldScale;
            field[2] += B.z() * dd4hep::tesla * fieldScale;
          }
        },
        fieldBrBz);

  } else { // BxByBz

    // Return 0 field outside map extent
    if (p.x() < bmin[0] || p.x() > bmax[0] || p.y() < bmin[1] || p.y() > bmax[1] ||
        p.z() < bmin[2] || p.z() > bmax[2])
      return;

    std::visit(
        [&](auto& f) {
          using T = std::decay_t<decltype(f)>;
          if constexpr (!std::is_same_v<T, std::monostate>) {
            typename T::view_t view(f);
            typename T::output_t b = view.at(p.x(), p.y(), p.z());
            auto B                 = fieldRot.has_value()
                                         ? fieldRot.value() * ROOT::Math::XYZPoint(b[0], b[1], b[2])
                                         : ROOT::Math::XYZPoint(b[0], b[1], b[2]);
            field[0] += B.x() * dd4hep::tesla * fieldScale;
            field[1] += B.y() * dd4hep::tesla * fieldScale;
            field[2] += B.z() * dd4hep::tesla * fieldScale;
          }
        },
        fieldBxByBz);
  }
}

// assign the field map to CartesianField
static Ref_t create_field_map_b(Detector& /*lcdd*/, xml::Handle_t handle) {
  xml_comp_t x_par(handle);

  if (!x_par.hasAttr(_Unicode(field_map))) {
    throw std::runtime_error(
        "FieldMapB Error: must have an xml attribute \"field_map\" for the field map.");
  }

  CartesianField field;
  std::string field_type = x_par.attr<std::string>(_Unicode(field_type));

  std::string coord_type = x_par.attr<std::string>(_Unicode(coord_type));

  if (coord_type.compare("BrBz") != 0 && coord_type.compare("BxByBz") != 0) {
    printout(ERROR, "FieldMapB", "Coordinate type: " + coord_type + ", is not BrBz nor BxByBz");
    std::_Exit(EXIT_FAILURE);
  }

  std::string field_map_file  = x_par.attr<std::string>(_Unicode(field_map));
  std::string field_map_url   = x_par.attr<std::string>(_Unicode(url));
  std::string field_map_cache = getAttrOrDefault<std::string>(x_par, _Unicode(cache), "");

  EnsureFileFromURLExists(field_map_url, field_map_file, field_map_cache);

  float field_map_scale = x_par.attr<float>(_Unicode(scale));

  if (!fs::exists(fs::path(field_map_file))) {
    printout(ERROR, "FieldMapB", "file " + field_map_file + " does not exist");
    printout(ERROR, "FieldMapB", "use a FileLoader plugin before the field element");
    std::_Exit(EXIT_FAILURE);
  }

  auto map = new FieldMapB(field_type, coord_type,
                           getAttrOrDefault<std::string>(x_par, _Unicode(backend_type), "strided"));

  // Optional field rotation and coordinate translation
  if (x_par.hasChild(_Unicode(dimensions))) {
    xml_comp_t x_dim = x_par.dimensions();

    if (x_dim.hasChild(_Unicode(rotationField))) {
      static float deg2r = ROOT::Math::Pi() / 180.;
      xml_comp_t rot_dim = x_dim.child(_Unicode(rotationField));
      RotationZYX rot(rot_dim.z() * deg2r, rot_dim.y() * deg2r, rot_dim.x() * deg2r);
      map->SetFieldRotation(Transform3D(rot));
    }

    if (x_dim.hasChild(_Unicode(translationCoord))) {
      xml_comp_t trans_dim = x_dim.child(_Unicode(translationCoord));
      Translation3D trans(trans_dim.x(), trans_dim.y(), trans_dim.z());
      map->SetCoordTranslation(Transform3D(trans));
    }
  }

  map->LoadMap(field_map_file, field_map_scale);
  field.assign(map, x_par.nameStr(), "FieldMapB");

  return field;
}

DECLARE_XMLELEMENT(epic_FieldMapB, create_field_map_b)
