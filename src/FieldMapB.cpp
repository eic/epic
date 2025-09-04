// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck, Dhevan Gangadharan, Aranya Giri

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/FieldTypes.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>
#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/array.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
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
namespace fs = std::filesystem;

#include "FileLoaderHelper.h"

using fieldBrBz_t =
    covfie::field<covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::size2, covfie::backend::array<covfie::vector::float2>>>>>;

using fieldBxByBz_t =
    covfie::field<covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::size3, covfie::backend::array<covfie::vector::float3>>>>>;

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

public:
  FieldMapB(const std::string& field_type_str = "magnetic",
            const std::string& coord_type_str = "BrBz");

  void LoadMap(const std::string& map_file, float scale);

  void SetCoordTranslation(const Transform3D& tr) {
    coordTranslate     = tr;
    coordTranslate_inv = tr.Inverse();
  }
  void SetFieldRotation(const Transform3D& tr) {
    fieldRot     = tr;
    fieldRot_inv = tr.Inverse();
  }

  virtual void fieldComponents(const double* pos, double* field);

private:
  FieldCoord fieldCoord;                                   // field coordinate type
  std::optional<Transform3D> coordTranslate{std::nullopt}; // coord translation
  std::optional<Transform3D> coordTranslate_inv{std::nullopt};
  std::optional<Transform3D> fieldRot{std::nullopt}; // field rotation
  std::optional<Transform3D> fieldRot_inv{std::nullopt};
  std::vector<float> steps, mins, maxs; // B map cell info

  std::shared_ptr<fieldBrBz_t> fieldBrBz;
  std::shared_ptr<fieldBxByBz_t> fieldBxByBz;
};

// constructor
FieldMapB::FieldMapB(const std::string& field_type_str, const std::string& coord_type_str) {
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
}

// load data
void FieldMapB::LoadMap(const std::string& map_file, float) {
  std::string line;
  std::ifstream input(map_file);
  if (!input) {
    printout(ERROR, "FieldMapB", "FieldMapB Error: file " + map_file + " cannot be read.");
  }

  switch (fieldCoord) {
  case FieldCoord::BrBz:
    fieldBrBz = std::make_shared<fieldBrBz_t>(input);
    break;
  case FieldCoord::BxByBz:
    fieldBxByBz = std::make_shared<fieldBxByBz_t>(input);
    break;
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
    const float r = sqrt(p.x() * p.x() + p.y() * p.y());
    const float z = p.z();

    fieldBrBz_t::view_t view(*fieldBrBz);
    fieldBrBz_t::output_t b = view.at(r, z);
    float Br                = b[0];
    float Bz                = b[1];

    // convert Br Bz to Bx By Bz and rotate field
    const float phi = atan2(p.y(), p.x());
    auto B          = fieldRot.has_value()
                          ? fieldRot.value() * ROOT::Math::XYZPoint(Br * cos(phi), Br * sin(phi), Bz)
                          : ROOT::Math::XYZPoint(Br * cos(phi), Br * sin(phi), Bz);

    field[0] += B.x();
    field[1] += B.y();
    field[2] += B.z();

  } else { // BxByBz

    fieldBxByBz_t::view_t view(*fieldBxByBz);
    fieldBxByBz_t::output_t b = view.at(p.x(), p.y(), p.z());

    // field rotation done in LoadMap()
    field[0] += b[0];
    field[1] += b[1];
    field[2] += b[2];
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

  // dimensions
  xml_comp_t x_dim = x_par.dimensions();

  // vector of dimension parameters: step, min, max
  std::vector<xml_comp_t> dimensions;

  if (coord_type.compare("BrBz") == 0) {
    dimensions.push_back(x_dim.child(_Unicode(R)));
    dimensions.push_back(x_dim.child(_Unicode(Z)));
  } else if (coord_type.compare("BxByBz") == 0) {
    dimensions.push_back(x_dim.child(_Unicode(X)));
    dimensions.push_back(x_dim.child(_Unicode(Y)));
    dimensions.push_back(x_dim.child(_Unicode(Z)));
  } else {
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

  auto map = new FieldMapB(field_type, coord_type);

  // translation
  if (x_dim.hasChild(_Unicode(rotationField))) {
    static float deg2r = ROOT::Math::Pi() / 180.;
    xml_comp_t rot_dim = x_dim.child(_Unicode(rotationField));
    RotationZYX rot(rot_dim.z() * deg2r, rot_dim.y() * deg2r, rot_dim.x() * deg2r);
    map->SetFieldRotation(Transform3D(rot));
  }
  // rotation
  if (x_dim.hasChild(_Unicode(translationCoord))) {
    xml_comp_t trans_dim = x_dim.child(_Unicode(translationCoord));
    Translation3D trans(trans_dim.x(), trans_dim.y(), trans_dim.z());
    map->SetCoordTranslation(Transform3D(trans));
  }

  map->LoadMap(field_map_file, field_map_scale);
  field.assign(map, x_par.nameStr(), "FieldMapB");

  return field;
}

DECLARE_XMLELEMENT(epic_FieldMapB, create_field_map_b)
