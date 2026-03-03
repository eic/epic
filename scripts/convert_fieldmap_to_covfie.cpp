// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 the ePIC collaboration

// Convert a text field map to covfie binary format for use with FieldMapB.
//
// The output field is stored with coordinates in DD4hep length units (cm, since
// dd4hep::cm == 1.0) so that FieldMapB can call view.at() directly with the
// values it already has. Text field map coordinates are assumed to be in cm.
// Field values are stored in Tesla; FieldMapB multiplies by dd4hep::tesla.
//
// Usage examples:
//   BrBz (solenoid):
//     convert_fieldmap_to_covfie --coord-type BrBz
//       --r-min 0 --r-max 998 --r-step 2
//       --z-min -800 --z-max 798 --z-step 2
//       -i field.txt -o field.covfie
//
//   BxByBz (dipole):
//     convert_fieldmap_to_covfie --coord-type BxByBz
//       --x-min -7.5 --x-max 7.5 --x-step 0.5
//       --y-min -34  --y-max 34  --y-step 2
//       --z-min -80  --z-max 80  --z-step 2
//       -i field.txt -o field.covfie

#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/array.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/parameter_pack.hpp>

#include <boost/program_options.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace po = boost::program_options;

// Field types matching FieldMapB.cpp (linear ≡ nearest_neighbour on disk)
using fieldBrBz_t = covfie::field<
    covfie::backend::affine<covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::size2, covfie::backend::array<covfie::vector::float2>>>>>;

using fieldBxByBz_t = covfie::field<
    covfie::backend::affine<covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::size3, covfie::backend::array<covfie::vector::float3>>>>>;

// Coordinates in the text field maps are in cm.
// DD4hep's internal length unit is cm (dd4hep::cm == 1.0), so the affine
// transform must map cm → grid index to match what FieldMapB::fieldComponents
// receives in pos[].
// Field values in the text maps are in Tesla. FieldMapB::fieldComponents
// converts to DD4hep internal field units (dd4hep::tesla = 1e-13) before
// accumulating into field[]. The converter therefore stores raw Tesla values.

void convert_BrBz(const std::string& input_file, const std::string& output_file, float r_min,
                  float r_max, float r_step, float z_min, float z_max, float z_step, float scale) {
  const std::size_t nr = static_cast<std::size_t>(std::lround((r_max - r_min) / r_step)) + 1;
  const std::size_t nz = static_cast<std::size_t>(std::lround((z_max - z_min) / z_step)) + 1;
  std::cout << "Grid: " << nr << " x " << nz << " (BrBz)" << std::endl;

  // Affine: map cm coordinates → grid indices
  covfie::algebra::affine<2> translation = covfie::algebra::affine<2>::translation(-r_min, -z_min);
  covfie::algebra::affine<2> scaling     = covfie::algebra::affine<2>::scaling(
      static_cast<float>(nr - 1) / (r_max - r_min), static_cast<float>(nz - 1) / (z_max - z_min));

  fieldBrBz_t field(covfie::make_parameter_pack(
      fieldBrBz_t::backend_t::configuration_t(scaling * translation),
      fieldBrBz_t::backend_t::backend_t::configuration_t{},
      fieldBrBz_t::backend_t::backend_t::backend_t::configuration_t{nr, nz}));
  fieldBrBz_t::view_t fv(field);

  std::ifstream f(input_file);
  if (!f.good()) {
    std::cerr << "Failed to open input file: " << input_file << std::endl;
    std::exit(1);
  }

  std::size_t count = 0;
  float r, z, Br, Bz;
  while (f >> r >> z >> Br >> Bz) {
    fieldBrBz_t::view_t::output_t& p = fv.at(r, z); // r, z already in cm
    p[0]                             = Br * scale;
    p[1]                             = Bz * scale;
    ++count;
  }
  std::cout << "Read " << count << " field values." << std::endl;

  std::ofstream fs(output_file, std::ofstream::binary);
  if (!fs.good()) {
    std::cerr << "Failed to open output file: " << output_file << std::endl;
    std::exit(1);
  }
  field.dump(fs);
  std::cout << "Wrote " << output_file << std::endl;
}

void convert_BxByBz(const std::string& input_file, const std::string& output_file, float x_min,
                    float x_max, float x_step, float y_min, float y_max, float y_step, float z_min,
                    float z_max, float z_step, float scale) {
  const std::size_t nx = static_cast<std::size_t>(std::lround((x_max - x_min) / x_step)) + 1;
  const std::size_t ny = static_cast<std::size_t>(std::lround((y_max - y_min) / y_step)) + 1;
  const std::size_t nz = static_cast<std::size_t>(std::lround((z_max - z_min) / z_step)) + 1;
  std::cout << "Grid: " << nx << " x " << ny << " x " << nz << " (BxByBz)" << std::endl;

  // Affine: map cm coordinates → grid indices
  covfie::algebra::affine<3> translation =
      covfie::algebra::affine<3>::translation(-x_min, -y_min, -z_min);
  covfie::algebra::affine<3> scaling = covfie::algebra::affine<3>::scaling(
      static_cast<float>(nx - 1) / (x_max - x_min), static_cast<float>(ny - 1) / (y_max - y_min),
      static_cast<float>(nz - 1) / (z_max - z_min));

  fieldBxByBz_t field(covfie::make_parameter_pack(
      fieldBxByBz_t::backend_t::configuration_t(scaling * translation),
      fieldBxByBz_t::backend_t::backend_t::configuration_t{},
      fieldBxByBz_t::backend_t::backend_t::backend_t::configuration_t{nx, ny, nz}));
  fieldBxByBz_t::view_t fv(field);

  std::ifstream f(input_file);
  if (!f.good()) {
    std::cerr << "Failed to open input file: " << input_file << std::endl;
    std::exit(1);
  }

  std::size_t count = 0;
  float x, y, z, Bx, By, Bz;
  while (f >> x >> y >> z >> Bx >> By >> Bz) {
    fieldBxByBz_t::view_t::output_t& p = fv.at(x, y, z); // x, y, z already in cm
    p[0]                               = Bx * scale;
    p[1]                               = By * scale;
    p[2]                               = Bz * scale;
    ++count;
  }
  std::cout << "Read " << count << " field values." << std::endl;

  std::ofstream fs(output_file, std::ofstream::binary);
  if (!fs.good()) {
    std::cerr << "Failed to open output file: " << output_file << std::endl;
    std::exit(1);
  }
  field.dump(fs);
  std::cout << "Wrote " << output_file << std::endl;
}

int main(int argc, char** argv) {
  po::options_description opts("convert_fieldmap_to_covfie options");
  // clang-format off
  opts.add_options()
    ("help,h", "produce help message")
    ("input,i",       po::value<std::string>()->required(), "input text field map file")
    ("output,o",      po::value<std::string>()->required(), "output covfie binary file")
    ("coord-type,t",  po::value<std::string>()->required(), "coordinate type: BrBz or BxByBz")
    ("scale,s",       po::value<float>()->default_value(1.0f), "scale factor applied to field values")
    ("r-min",  po::value<float>(), "R minimum [cm] (BrBz)")
    ("r-max",  po::value<float>(), "R maximum [cm] (BrBz)")
    ("r-step", po::value<float>(), "R step size [cm] (BrBz)")
    ("x-min",  po::value<float>(), "X minimum [cm] (BxByBz)")
    ("x-max",  po::value<float>(), "X maximum [cm] (BxByBz)")
    ("x-step", po::value<float>(), "X step size [cm] (BxByBz)")
    ("y-min",  po::value<float>(), "Y minimum [cm] (BxByBz)")
    ("y-max",  po::value<float>(), "Y maximum [cm] (BxByBz)")
    ("y-step", po::value<float>(), "Y step size [cm] (BxByBz)")
    ("z-min",  po::value<float>(), "Z minimum [cm]")
    ("z-max",  po::value<float>(), "Z maximum [cm]")
    ("z-step", po::value<float>(), "Z step size [cm]");
  // clang-format on

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opts), vm);
    if (vm.count("help")) {
      std::cout << opts << std::endl;
      return 0;
    }
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "Error: " << e.what() << "\n" << opts << std::endl;
    return 1;
  }

  const auto input      = vm["input"].as<std::string>();
  const auto output     = vm["output"].as<std::string>();
  const auto coord_type = vm["coord-type"].as<std::string>();
  const float scale     = vm["scale"].as<float>();

  if (coord_type == "BrBz") {
    for (const char* p : {"r-min", "r-max", "r-step", "z-min", "z-max", "z-step"}) {
      if (!vm.count(p)) {
        std::cerr << "Error: --" << p << " is required for BrBz" << std::endl;
        return 1;
      }
    }
    convert_BrBz(input, output, vm["r-min"].as<float>(), vm["r-max"].as<float>(),
                 vm["r-step"].as<float>(), vm["z-min"].as<float>(), vm["z-max"].as<float>(),
                 vm["z-step"].as<float>(), scale);
  } else if (coord_type == "BxByBz") {
    for (const char* p :
         {"x-min", "x-max", "x-step", "y-min", "y-max", "y-step", "z-min", "z-max", "z-step"}) {
      if (!vm.count(p)) {
        std::cerr << "Error: --" << p << " is required for BxByBz" << std::endl;
        return 1;
      }
    }
    convert_BxByBz(input, output, vm["x-min"].as<float>(), vm["x-max"].as<float>(),
                   vm["x-step"].as<float>(), vm["y-min"].as<float>(), vm["y-max"].as<float>(),
                   vm["y-step"].as<float>(), vm["z-min"].as<float>(), vm["z-max"].as<float>(),
                   vm["z-step"].as<float>(), scale);
  } else {
    std::cerr << "Error: --coord-type must be BrBz or BxByBz" << std::endl;
    return 1;
  }

  return 0;
}
