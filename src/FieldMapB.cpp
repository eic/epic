// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck, Dhevan Gangadharan, Aranya Giri

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/FieldTypes.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
namespace fs = std::filesystem;

#include "FileLoaderHelper.h"

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

// In DD4hep 1.26, the name for the `field_type` enum changed (from `type`).
#if DD4HEP_VERSION_GE(1, 26)
#define DD4HEP_FIELD_TYPE field_type
#else
#define DD4HEP_FIELD_TYPE type
#endif

// This allows us to specify the name of the variable by hand, when patching
// the previous versions, by setting `DD4HEP_FIELD_TYPE_OVERRIDE`.
// TODO remove after DD4hep 1.26 release
#define XSTR(x) STR(x)
#define STR(x) #x
#ifdef DD4HEP_FIELD_TYPE_OVERRIDE
#undef DD4HEP_FIELD_TYPE
#define DD4HEP_FIELD_TYPE DD4HEP_FIELD_TYPE_OVERRIDE
#pragma message("DD4HEP_FIELD_TYPE overridden as " XSTR(DD4HEP_FIELD_TYPE))
#endif
#pragma message("DD4HEP_FIELD_TYPE is " XSTR(DD4HEP_FIELD_TYPE))

// implementation of the field map
class FieldMapB : public dd4hep::CartesianField::Object {
public:
  FieldMapB(const std::string& field_type_str = "magnetic", const std::string& coord_type_str = "BrBz");
  void Configure(std::vector<xml_comp_t> dimensions);
  void LoadMap(const std::string& map_file, double scale);
  void GetIndices(std::vector<double> coord, std::vector<int> *indices, std::vector<double> *deltas);
  void SetTransform(const Transform3D& tr)
  {
    trans     = tr;
    trans_inv = tr.Inverse();
  }

  virtual void fieldComponents(const double* pos, double* field);

private:
  std::string                                                   coord_type;
  Transform3D                                                   trans, trans_inv;
  std::vector<double>                                           steps, mins, maxs;
  std::vector<std::vector<std::vector<double>>>                 Bvals_RZ; // {R}, {Z}, {Br,Bz}
  std::vector<std::vector<std::vector<std::vector<double>>>>    Bvals_XYZ; // {X}, {Y}, {Z}, {Bx,By,Bz}
};

// constructor
FieldMapB::FieldMapB(const std::string& field_type_str, const std::string& coord_type_str)
{
  std::string ftype = field_type_str;
  for (auto& c : ftype) {
    c = tolower(c);
  }

  // set type
  if (ftype == "magnetic") {
    DD4HEP_FIELD_TYPE = CartesianField::MAGNETIC;
  } else if (ftype == "electric") {
    DD4HEP_FIELD_TYPE = CartesianField::ELECTRIC;
  } else {
    DD4HEP_FIELD_TYPE = CartesianField::UNKNOWN;
    std::cout << "FieldMapB Warning: Unknown field type " << field_type_str << "!" << std::endl;
  }

  coord_type = coord_type_str;
}

// fill field vector
void FieldMapB::Configure(std::vector<xml_comp_t> dimensions)
{
  // Fill vectors with step size, min, max for each dimension
  for( auto el : dimensions ) {
    steps.push_back( el.step() ); 
    mins.push_back( getAttrOrDefault<double>(el, _Unicode(min), 0) );
    maxs.push_back( getAttrOrDefault<double>(el, _Unicode(max), 0) );
  }

  if( coord_type.compare("BrBz") == 0 ) {
    int nr = int((maxs[0] - mins[0]) / steps[0]) + 2;
    int nz = int((maxs[1] - mins[1]) / steps[1]) + 2;
    
    Bvals_RZ.resize(nr);
    for (auto& B2 : Bvals_RZ) {
      B2.resize(nz);
      for (auto& B : B2) {
        B.resize(2, 0.); // Br, Bz
      }
    }
  }

  if( coord_type.compare("BxByBz") == 0 ) {
    int nx = int((maxs[0] - mins[0]) / steps[0]) + 2;
    int ny = int((maxs[1] - mins[1]) / steps[1]) + 2;
    int nz = int((maxs[2] - mins[2]) / steps[2]) + 2;
   
    Bvals_XYZ.resize(nx);
    for (auto& B3 : Bvals_XYZ) {
      B3.resize(ny);
      for (auto& B2 : B3) {
        B2.resize(nz);
        for (auto& B : B2) {
          B.resize(3, 0.); // Bx, By, Bz
        }
      }
    }
  }
}

// get cell indices corresponding to point of interest 
void FieldMapB::GetIndices(std::vector<double> coord, std::vector<int> *indices, std::vector<double> *deltas)
{
  for( uint i = 0; i < coord.size(); i++ ) {
     // Range check
     if( coord[i] < mins[i] || coord[i] > maxs[i] ) {
       indices->clear();
       indices->resize( coord.size(), -1 );
       return;
     }

     double temp_idx;
     deltas->push_back( std::modf( (coord[i] - mins[i]) / steps[i], &temp_idx ) );
     indices->push_back( static_cast<int>(temp_idx) );
  }
}

// load data
void FieldMapB::LoadMap(const std::string& map_file, double scale)
{
  std::string   line;
  std::ifstream input(map_file);
  if (!input) {
    std::cout << "FieldMapB Error: file \"" << map_file << "\" cannot be read." << std::endl;
  }

  std::vector<double> coord = {};
  std::vector<double> Bcomp = {};
  std::vector<int> indices = {};
  std::vector<double> deltas = {};

  while (std::getline(input, line).good()) {
    std::istringstream iss(line);

    coord.clear();
    Bcomp.clear();
    indices.clear();
    deltas.clear();

    if( coord_type.compare("BrBz") == 0 ) {
      coord.resize(2);
      Bcomp.resize(2);
      iss >> coord[0] >> coord[1] >> Bcomp[0] >> Bcomp[1];
    }
    if( coord_type.compare("BxByBz") == 0 ) {
      coord.resize(3);
      Bcomp.resize(3);
      iss >> coord[0] >> coord[1] >> coord[2] >> Bcomp[0] >> Bcomp[1] >> Bcomp[2];
    }

    GetIndices(coord, &indices, &deltas);

    if (std::count( indices.begin(), indices.end(), -1)) {
      std::cout << "FieldMapB Warning: coordinates out of range, skipped it."
        << std::endl;
    } else {
      if( coord_type.compare("BrBz") == 0 ) {
        Bvals_RZ[ indices[0] ][ indices[1] ] = {Bcomp[0] * scale, Bcomp[1] * scale};
        //ROOT::Math::XYZPoint p(coord[0], 0, coord[1]);
        //std::cout << p << " -> " << trans*p << std::endl;
        //std::cout << indices[0] << ", " << indices[1] << ", " << Bcomp[0] << ", " << Bcomp[1] << std::endl;
      }
      if( coord_type.compare("BxByBz") == 0 ) {
        Bvals_XYZ[ indices[0] ][ indices[1] ][ indices[2] ] = {Bcomp[0] * scale, Bcomp[1] * scale, Bcomp[2] * scale};
      }
    }
  }
}

// get field components
void FieldMapB::fieldComponents(const double* pos, double* field)
{
  // coordinate conversion
  auto p = trans_inv * ROOT::Math::XYZPoint(pos[0], pos[1], pos[2]);
  
  std::vector<int>  indices = {};
  std::vector<double> deltas = {};
  ROOT::Math::XYZPoint B = ROOT::Math::XYZPoint(0,0,0);

  if( coord_type.compare("BrBz") == 0 ) {
    // coordinates conversion
    const double r   = sqrt(p.x() * p.x() + p.y() * p.y());
    const double z   = p.z();
    const double phi = atan2(p.y(), p.x());

    GetIndices({r, z}, &indices, &deltas);

    // out of the range
    if (std::count( indices.begin(), indices.end(), -1)) {
      return;
    }

    // p1    p3
    //    p
    // p0    p2
    auto& p0 = Bvals_RZ[indices[0]][indices[1]];
    auto& p1 = Bvals_RZ[indices[0]][indices[1] + 1];
    auto& p2 = Bvals_RZ[indices[0] + 1][indices[1]];
    auto& p3 = Bvals_RZ[indices[0] + 1][indices[1] + 1];

    // Bilinear interpolation
    double Br = p0[0] * (1 - deltas[0]) * (1 - deltas[1]) 
      + p1[0] * (1 - deltas[0]) * deltas[1] 
      + p2[0] * deltas[0] * (1 - deltas[1])
      + p3[0] * deltas[0] * deltas[1];

    double Bz = p0[1] * (1 - deltas[0]) * (1 - deltas[1]) 
      + p1[1] * (1 - deltas[0]) * deltas[1] 
      + p2[1] * deltas[0] * (1 - deltas[1]) 
      + p3[1] * deltas[0] * deltas[1];

    // convert Br Bz to Bx By Bz
    B = ROOT::Math::XYZPoint(Br * sin(phi), Br * cos(phi), Bz);
    B = trans * B;
  }

  if( coord_type.compare("BxByBz") == 0 ) {

    GetIndices({p.x(), p.y(), p.z()}, &indices, &deltas);

    // out of the range
    if (std::count( indices.begin(), indices.end(), -1)) {
      return;
    }

    for(int comp = 0; comp < 3; comp++) { // field component loop
      // Trilinear interpolation
      // First along X, along 4 lines
      double b00 = Bvals_XYZ[indices[0]][indices[1]][indices[2]][comp]         * (1 - deltas[0])         
        + Bvals_XYZ[indices[0] + 1][indices[1]][indices[2]][comp]              * deltas[0];
      double b01 = Bvals_XYZ[indices[0]][indices[1]][indices[2] + 1][comp]     * (1 - deltas[0])     
        + Bvals_XYZ[indices[0] + 1][indices[1]][indices[2] + 1][comp]          * deltas[0];
      double b10 = Bvals_XYZ[indices[0]][indices[1] + 1][indices[2]][comp]     * (1 - deltas[0])     
        + Bvals_XYZ[indices[0] + 1][indices[1] + 1][indices[2]][comp]          * deltas[0];
      double b11 = Bvals_XYZ[indices[0]][indices[1] + 1][indices[2] + 1][comp] * (1 - deltas[0]) 
        + Bvals_XYZ[indices[0] + 1][indices[1] + 1][indices[2] + 1][comp]      * deltas[0];
      //// Next along Y, along 2 lines
      double b0 = b00 * (1 - deltas[1]) + b10 * deltas[1];
      double b1 = b01 * (1 - deltas[1]) + b11 * deltas[1];
      //// Finally along Z
      double b = b0 * (1 - deltas[2]) + b1 * deltas[2];
      if(comp == 0) { B.SetX( b ); }
      if(comp == 1) { B.SetY( b ); }
      if(comp == 2) { B.SetZ( b ); }
    }
  }

  field[0] += B.x() * tesla;
  field[1] += B.y() * tesla;
  field[2] += B.z() * tesla;

  return;  
}

// assign the field map to CartesianField
static Ref_t create_field_map_b(Detector& /*lcdd*/, xml::Handle_t handle)
{
  xml_comp_t x_par(handle);

  if (!x_par.hasAttr(_Unicode(field_map))) {
    throw std::runtime_error("FieldMapB Error: must have an xml attribute \"field_map\" for the field map.");
  }

  CartesianField field;
  std::string    field_type = x_par.attr<std::string>(_Unicode(field_type));
  
  // coordinate type: "BrBz" or "BxByBz"
  std::string    coord_type = x_par.attr<std::string>(_Unicode(coord_type));
  
  // dimensions
  xml_comp_t x_dim = x_par.dimensions();

  // vector of dimension parameters: step, min, max
  std::vector<xml_comp_t> dimensions;

  if( coord_type.compare("BrBz") == 0 ) {
     dimensions.push_back( x_dim.child(_Unicode(R)) );
     dimensions.push_back( x_dim.child(_Unicode(Z)) );
  }
  else if( coord_type.compare("BxByBz") == 0 ) {
     dimensions.push_back( x_dim.child(_Unicode(X)) );
     dimensions.push_back( x_dim.child(_Unicode(Y)) );
     dimensions.push_back( x_dim.child(_Unicode(Z)) );
  }
  else {
    printout(ERROR, "FieldMapB", "Coordinate type: " + coord_type + ", is not BrBz nor BxByBz");
    std::_Exit(EXIT_FAILURE);
  }
  
  std::string field_map_file  = x_par.attr<std::string>(_Unicode(field_map));
  std::string field_map_url   = x_par.attr<std::string>(_Unicode(url));
  std::string field_map_cache = getAttrOrDefault<std::string>(x_par, _Unicode(cache), "");

  EnsureFileFromURLExists(field_map_url, field_map_file, field_map_cache);

  double field_map_scale = x_par.attr<double>(_Unicode(scale));

  if (!fs::exists(fs::path(field_map_file))) {
    printout(ERROR, "FieldMapB", "file " + field_map_file + " does not exist");
    printout(ERROR, "FieldMapB", "use a FileLoader plugin before the field element");
    std::_Exit(EXIT_FAILURE);
  }

  auto map = new FieldMapB(field_type, coord_type);
  map->Configure( dimensions );

  // translation, rotation
  static double deg2r = ROOT::Math::Pi() / 180.;
  RotationZYX   rot(0., 0., 0.);
  if (x_dim.hasChild(_Unicode(rotation))) {
    xml_comp_t rot_dim = x_dim.child(_Unicode(rotation));
    rot                = RotationZYX(rot_dim.z() * deg2r, rot_dim.y() * deg2r, rot_dim.x() * deg2r);
  }

  Translation3D trans(0., 0., 0.);
  if (x_dim.hasChild(_Unicode(translation))) {
    xml_comp_t trans_dim = x_dim.child(_Unicode(translation));
    trans                = Translation3D(trans_dim.x(), trans_dim.y(), trans_dim.z());
  }
  map->SetTransform(trans * rot);

  map->LoadMap(field_map_file, field_map_scale);
  field.assign(map, x_par.nameStr(), "FieldMapB");

  return field;
}

DECLARE_XMLELEMENT(epic_FieldMapB, create_field_map_b)
