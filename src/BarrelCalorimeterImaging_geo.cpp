// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Maria Zurek, Whitney Armstrong

// Detector plugin to support a hybrid central barrel calorimeter
// The detector consists of interlayers of Pb/ScFi (segmentation in global r, phi) and W/Si (segmentation in local x, y)
// Assembly is used as the envelope so two different detectors can be interlayered with each other
//
//
// 06/19/2021: Implementation of the Sci Fiber geometry. M. Żurek
// 07/09/2021: Support interlayers between multiple detectors. C. Peng
// 07/23/2021: Add assemblies as mother volumes of fibers to reduce the number of daughter volumes. C. Peng, M. Żurek
//     Reference: TGeo performance issue with large number of daughter volumes
//     https://indico.cern.ch/event/967418/contributions/4075358/attachments/2128099/3583278/201009_shKo_dd4hep.pdf
// 07/24/2021: Changed support implementation to avoid too many uses of boolean geometries. DAWN view seems to have
//     issue dealing with it. C. Peng

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/Surface.h"
#include "Math/Point2D.h"
#include "TGeoPolygon.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

#include "DD4hepDetectorHelper.h"

using namespace dd4hep;

typedef ROOT::Math::XYPoint Point;

// geometry helpers
static void buildSupport(Detector& desc, Volume& mother, xml_comp_t x_support,
                         const std::tuple<double, double, double, double>& dimensions);

// barrel ecal layers contained in an assembly
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens) {
  xml_det_t x_detector      = e;
  int detector_id           = x_detector.id();
  std::string detector_name = x_detector.nameStr();
  double offset             = x_detector.attr<double>(_Unicode(offset));
  xml_comp_t x_dimensions   = x_detector.dimensions();
  int nsides                = x_dimensions.numsides();
  double inner_r            = x_dimensions.rmin();
  double dphi               = (2 * M_PI / nsides);
  double half_dphi          = dphi / 2;

  DetElement sdet(detector_name, detector_id);
  Volume mother_volume = desc.pickMotherVolume(sdet);

  Assembly detector_volume(detector_name);
  detector_volume.setAttributes(desc, x_detector.regionStr(), x_detector.limitsStr(),
                                x_detector.visStr());

  Transform3D detector_tr       = Translation3D(0, 0, offset) * RotationZ(0);
  PlacedVolume detector_physvol = mother_volume.placeVolume(detector_volume, detector_tr);
  sens.setType("calorimeter");

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_detector, sdet);

  detector_physvol.addPhysVolID("system", detector_id);
  sdet.setPlacement(detector_physvol);

  // Build the envelope assembly for layer 1
  DetElement envelope1_element("envelope1", detector_id);
  Assembly envelope1_volume("envelope1");
  auto& envelope1_params =
      DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(envelope1_element);
  if (x_detector.hasChild(_Unicode(envelope))) {
    xml_comp_t x_envelope = x_detector.child(_Unicode(envelope));
    envelope1_params.set<double>("envelope_r_min",
                                 getAttrOrDefault(x_envelope, _Unicode(envelope_r_min), 0.));
    envelope1_params.set<double>("envelope_r_max",
                                 getAttrOrDefault(x_envelope, _Unicode(envelope_r_max), 0.));
    envelope1_params.set<double>("envelope_z_min",
                                 getAttrOrDefault(x_envelope, _Unicode(envelope_z_min), 0.));
    envelope1_params.set<double>("envelope_z_max",
                                 getAttrOrDefault(x_envelope, _Unicode(envelope_z_max), 0.));
  }
  // Add the volume boundary material if configured
  for (xml_coll_t boundary_material(x_detector, _Unicode(boundary_material)); boundary_material;
       ++boundary_material) {
    xml_comp_t x_boundary_material = boundary_material;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, envelope1_params,
                                                    "boundary_material");
  }
  // Add the volume layer material if configured
  for (xml_coll_t layer_material(x_detector, _Unicode(layer_material)); layer_material;
       ++layer_material) {
    xml_comp_t x_layer_material = layer_material;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, envelope1_params,
                                                    "layer_material");
  }

  // Build the envelope assembly for layers 2 to 6 (no Acts support)
  DetElement envelope26_element("envelope26", detector_id);
  Assembly envelope26_volume("envelope26");

  // Build a sector for layer 1
  DetElement sector1_element("sector1", detector_id);
  Assembly sector1_volume("sector1");
  if (x_detector.hasChild(_Unicode(sectors))) {
    xml_comp_t x_sectors = x_detector.child(_Unicode(sectors));
    sector1_volume.setVisAttributes(desc.visAttributes(x_sectors.visStr()));
  }

  // Build a sector for layers 2 to 6
  DetElement sector26_element("sector26", detector_id);
  Assembly sector26_volume("sector26");
  if (x_detector.hasChild(_Unicode(sectors))) {
    xml_comp_t x_sectors = x_detector.child(_Unicode(sectors));
    sector26_volume.setVisAttributes(desc.visAttributes(x_sectors.visStr()));
  }

  // keep tracking of the total thickness
  double layer_pos_z = inner_r;

  // Parameters for computing the layer X dimension:
  double tan_half_dphi = std::tan(half_dphi);
  double layer_dim_y   = x_dimensions.z() / 2.;

  // Build the chips (one AstroPix sensor with its slices)
  using dd4hep::rec::VolPlane;
  std::map<std::string, Volume> volumes;
  std::map<std::string, std::vector<PlacedVolume>> sensitives;
  std::map<std::string, std::vector<VolPlane>> volplane_surfaces;
  std::map<std::string, std::array<double, 3>> chip_dims; // width, length, thickness
  for (xml_coll_t i_chip(x_detector, _Unicode(chip)); i_chip; ++i_chip) {
    xml_comp_t x_chip     = i_chip;
    std::string chip_name = x_chip.nameStr();

    if (volumes.find(chip_name) != volumes.end()) {
      printout(ERROR, "BarrelCalorimeterImaging", "Chip with name %s already exists",
               chip_name.c_str());
      throw std::runtime_error("Logics error in building chips.");
    }

    // Compute chip total thickness (and its footprint) from components
    double total_thickness = 0;
    double chip_width      = 0;
    double chip_length     = 0;
    xml_coll_t ci(x_chip, _U(component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
      chip_width  = std::max(chip_width, xml_comp_t(ci).width());
      chip_length = std::max(chip_length, xml_comp_t(ci).length());
    }
    chip_dims[chip_name] = {chip_width, chip_length, total_thickness};

    // Create the chip assembly volume
    Assembly chip_volume(chip_name);
    volumes[chip_name] = chip_volume;
    chip_volume.setVisAttributes(desc.visAttributes(x_chip.visStr()));

    // Add components
    int sensor_number       = 1;
    int ncomponents         = 0;
    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;
    for (xml_coll_t i_component(x_chip, _U(component)); i_component; ++i_component, ++ncomponents) {
      xml_comp_t x_component           = i_component;
      xml_comp_t x_component_pos       = x_component.position(false);
      xml_comp_t x_component_rot       = x_component.rotation(false);
      const std::string component_name = _toString(ncomponents, "component%d");
      Box component_box(x_component.width() / 2, x_component.length() / 2,
                        x_component.thickness() / 2);
      Volume component_volume(component_name, component_box,
                              desc.material(x_component.materialStr()));
      component_volume.setAttributes(desc, x_component.regionStr(), x_component.limitsStr(),
                                     x_component.visStr());

      // Loop over the slices for this component
      int slice_num      = 1;
      double slice_pos_z = -x_component.thickness() / 2;
      for (xml_coll_t i_slice(x_component, _U(slice)); i_slice; ++i_slice) {
        xml_comp_t x_slice     = i_slice;
        std::string slice_name = Form("slice%d", slice_num);
        double slice_dim_x     = x_component.width() / 2;
        double slice_dim_y     = x_component.length() / 2;
        double slice_dim_z     = x_slice.thickness() / 2;
        Box slice_shape(slice_dim_x, slice_dim_y, slice_dim_z);
        Volume slice_volume(slice_name, slice_shape, desc.material(x_slice.materialStr()));
        slice_volume.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(),
                                   x_slice.visStr());

        // Place slice
        PlacedVolume slice_physvol = component_volume.placeVolume(
            slice_volume, Position(0, 0, slice_pos_z + x_slice.thickness() / 2));

        // Set sensitive
        if (x_slice.isSensitive()) {
          slice_physvol.addPhysVolID("slice", slice_num);
          slice_volume.setSensitiveDetector(sens);
          sensitives[chip_name].push_back(slice_physvol);
        }

        // Increment Z position of slice
        slice_pos_z += x_slice.thickness();
        ++slice_num;
      }

      // Utility variable for the relative z-offset based off the previous components
      const double zoff = thickness_sum + x_component.thickness() / 2.0;
      PlacedVolume component_physvol;
      if (x_component_pos && x_component_rot) {
        Position component_pos(x_component_pos.x(0), x_component_pos.y(0),
                               x_component_pos.z(0) + zoff);
        RotationZYX component_rot(x_component_rot.z(0), x_component_rot.y(0), x_component_rot.x(0));
        component_physvol =
            chip_volume.placeVolume(component_volume, Transform3D(component_rot, component_pos));
      } else if (x_component_rot) {
        Position component_pos(0, 0, zoff);
        RotationZYX component_rot(x_component_rot.z(0), x_component_rot.y(0), x_component_rot.x(0));
        component_physvol =
            chip_volume.placeVolume(component_volume, Transform3D(component_rot, component_pos));
      } else if (x_component_pos) {
        Position component_pos(x_component_pos.x(0), x_component_pos.y(0),
                               x_component_pos.z(0) + zoff);
        component_physvol = chip_volume.placeVolume(component_volume, component_pos);
      } else {
        Position component_pos(0, 0, zoff);
        component_physvol = chip_volume.placeVolume(component_volume, component_pos);
      }

      if (x_component.isSensitive()) {
        component_physvol.addPhysVolID("sensor", sensor_number++);
        component_volume.setSensitiveDetector(sens);
        sensitives[chip_name].push_back(component_physvol);
      }

      thickness_sum += x_component.thickness();
      thickness_so_far += x_component.thickness();

      // apply relative offsets in z-position used to stack components side-by-side
      if (x_component_pos) {
        thickness_sum += x_component_pos.z(0);
        thickness_so_far += x_component_pos.z(0);
      }
    }
  }
  // Build the modules: a row of identical chips along the stave length (local y), with a gap
  // between neighbouring chips. Chip k is centred at y = -length/2 + (k + 1/2) chip length + k gap.
  struct ModuleInfo {
    Volume volume;
    std::string chip;
    int n_chips;
    double width, length, thickness;
    std::vector<PlacedVolume> chip_physvols;
  };
  std::map<std::string, ModuleInfo> modules;
  for (xml_coll_t i_mod(x_detector, _U(module)); i_mod; ++i_mod) {
    xml_comp_t x_mod      = i_mod;
    std::string mod_name  = x_mod.nameStr();
    std::string chip_name = x_mod.attr<std::string>(_Unicode(chip));
    const int n_chips     = x_mod.attr<int>(_Unicode(chips));
    const double chip_gap = x_mod.gap();
    if (modules.count(mod_name) || !volumes.count(chip_name) || n_chips < 1 || chip_gap < 0.) {
      printout(ERROR, "BarrelCalorimeterImaging",
               "module %s: duplicate name, unknown chip %s, %d chips or chip gap %g",
               mod_name.c_str(), chip_name.c_str(), n_chips, chip_gap);
      throw std::runtime_error("Logics error in building modules.");
    }
    const auto [chip_w, chip_l, chip_t] = chip_dims[chip_name];
    ModuleInfo& mod                     = modules[mod_name];
    mod.chip                            = chip_name;
    mod.n_chips                         = n_chips;
    mod.width                           = chip_w;
    mod.length                          = n_chips * chip_l + (n_chips - 1) * chip_gap;
    mod.thickness                       = chip_t;
    mod.volume =
        Volume(mod_name, Box(mod.width / 2., mod.length / 2., mod.thickness / 2.), desc.air());
    mod.volume.setAttributes(desc, x_mod.regionStr(), x_mod.limitsStr(), x_mod.visStr());
    for (int k = 0; k < n_chips; ++k) {
      const double y            = -mod.length / 2. + (k + 0.5) * chip_l + k * chip_gap;
      PlacedVolume chip_physvol = mod.volume.placeVolume(volumes[chip_name], Position(0, y, 0));
      chip_physvol.addPhysVolID("chip", k);
      mod.chip_physvols.push_back(chip_physvol);
    }
  }

  // Readout field ranges, to check that the stave, module and chip indices can be encoded
  const IDDescriptor id_spec = sens.readout().idSpec();
  auto check_id_range        = [&](int layer, const char* field, int max_index) {
    const auto* f = id_spec.field(field);
    // largest unsigned value (BitFieldElement::maxValue() returns 2^width for unsigned fields)
    const long max_value = (1L << f->width()) - 1;
    if (max_index > max_value) {
      printout(ERROR, "BarrelCalorimeterImaging",
               "layer %d: %s index %d exceeds the readout field %s (max %ld)", layer, field,
               max_index, f->name().c_str(), max_value);
      throw std::runtime_error("Imaging layer indices do not fit the readout.");
    }
  };

  // Half-stave volumes and their DetElement templates, shared between layers with identical staves
  struct HalfStaveInfo {
    Volume volume;
    DetElement element; // template, cloned for each placement
  };
  std::map<std::string, HalfStaveInfo> half_staves;

  if (x_detector.hasChild(_Unicode(support))) {
    buildSupport(desc, sector1_volume, x_detector.child(_Unicode(support)),
                 {inner_r, layer_pos_z, x_dimensions.z(), half_dphi});
  }
  //Loop over the sets of layer elements in the detector.
  int layer_num = 1;
  for (xml_coll_t i_layer(x_detector, _U(layer)); i_layer; ++i_layer) {
    xml_comp_t x_layer         = i_layer;
    int layer_repeat           = x_layer.repeat();
    double layer_thickness     = x_layer.thickness();
    bool layer_has_frame       = x_layer.hasChild(_Unicode(frame));
    double layer_space_between = getAttrOrDefault(x_layer, _Unicode(space_between), 0.);
    double layer_space_before  = getAttrOrDefault(x_layer, _Unicode(space_before), 0.);
    // Frame dimensions (also used to locate the staves); without a frame the "frame" spans the layer
    double frame_height    = layer_thickness;
    double frame_thickness = 0.;
    double frame_side      = 0.; // thickness of the side walls, 0 for none
    if (layer_has_frame) {
      xml_comp_t x_frame = x_layer.child(_Unicode(frame));
      frame_height       = x_frame.height();
      frame_thickness    = x_frame.thickness();
      frame_side         = getAttrOrDefault(x_frame, _Unicode(side_thickness), 0.);
    }
    // Horizontal width of a side wall, whose faces are parallel to the sector side
    const double frame_side_dx = frame_side / std::cos(half_dphi);
    layer_pos_z += layer_space_before;

    // Loop over number of repeats for this layer.
    for (int layer_j = 0; layer_j < layer_repeat; layer_j++) {

      // Make an envelope for this layer
      std::string layer_name = Form("layer%d", layer_num);
      double layer_dim_x     = tan_half_dphi * layer_pos_z;
      auto layer_mat         = desc.air();

      Position layer_pos(0, 0, layer_pos_z + layer_thickness / 2.);
      double layer_trd_x1 = layer_dim_x;
      double layer_trd_x2 = layer_dim_x + layer_thickness * tan_half_dphi;
      double layer_trd_y1 = layer_dim_y;
      double layer_trd_y2 = layer_trd_y1;
      double layer_trd_z  = layer_thickness / 2; // account for frame
      Trapezoid layer_shape(layer_trd_x1, layer_trd_x2, layer_trd_y1, layer_trd_y2, layer_trd_z);
      Volume layer_volume(layer_name, layer_shape, layer_mat);
      DetElement layer_element(layer_name, detector_id);
      if (layer_num == 1) {
        sector1_element.add(layer_element);
      } else {
        sector26_element.add(layer_element);
      }

      // Set region, limitset, and vis of layer.
      layer_volume.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Loop over the staves for this layer, if the tray is enabled
      int stave_index = 0; // readout index of the next half-stave
      for (xml_coll_t i_stave(x_layer, _U(stave)); i_stave; ++i_stave) {
        xml_comp_t x_stave = i_stave;
        // Check if we enabled this silicon tray (actual silicon detector in the slots) at the top level
        // allowing us to easily disable layers in the XML file without needing multiple
        // copies of the XML file
        if (getAttrOrDefault(x_stave, _Unicode(enable), 1) == 0) {
          // disabled
          continue;
        }

        int stave_repeat   = x_stave.repeat();
        double stave_thick = x_stave.thickness();
        double stave_dim_x = x_stave.width() / 2.0;
        double stave_dim_y = x_stave.length() / 2.0;
        double stave_dim_z = stave_thick / 2.0;
        // Staves are flat and sit on one of two rails, alternating low (even index) and high (odd)
        double rail_low  = getAttrOrDefault(x_stave, _Unicode(rail_low), 0.);
        double rail_high = getAttrOrDefault(x_stave, _Unicode(rail_high), 0.);
        double wall_gap  = getAttrOrDefault(x_stave, _Unicode(wall_gap), 0.);

        // Stave bottom sits at the rail height above the outer face of the inner carbon plate
        const double z_plate = -frame_height / 2. + frame_thickness;
        auto stave_pos_z     = [&](int j) {
          return z_plate + (j % 2 == 0 ? rail_low : rail_high) + stave_dim_z;
        };

        // Arrange staves symmetrically around center of layer. The wall inner face is slanted, so
        // the closest approach is the outer bottom corner of a low stave: place the outermost staves
        // such that the perpendicular distance from that corner to the inner face is wall_gap
        const double wall_x_low   = layer_dim_x +
                                    (z_plate + rail_low + layer_thickness / 2.) * tan_half_dphi -
                                    frame_side_dx; // |x| of the inner face at the low stave bottom
        const double layer_wall_x = wall_x_low - wall_gap / std::cos(half_dphi);
        double stave_pos_x        = -layer_wall_x + stave_dim_x;
        double stave_pitch        = 0.;
        if (stave_repeat == 1) {
          stave_pos_x = 0.; // a single stave is centered
        } else {
          stave_pitch = -2.0 * stave_pos_x / (stave_repeat - 1);
        }

        // Check that the staves fit: no overlap between staves, with the outer plate or the walls
        {
          const double pos_x_outer = -stave_pos_x; // |x| of the outermost staves
          // perpendicular distance from the corner of the outermost low stave to the inner face
          const double wall_clearance =
              (wall_x_low - (pos_x_outer + stave_dim_x)) * std::cos(half_dphi);
          const bool bad_pitch = stave_repeat > 1 && stave_pitch < stave_dim_x;
          const bool bad_tiers = stave_repeat > 1 && stave_pitch < 2. * stave_dim_x &&
                                 std::abs(rail_high - rail_low) < stave_thick;
          const bool bad_plate =
              std::max(rail_low, rail_high) + stave_thick > frame_height - 2. * frame_thickness;
          if (bad_pitch || bad_tiers || bad_plate || wall_clearance < 0.) {
            printout(ERROR, "BarrelCalorimeterImaging",
                     "layer %d: %d staves do not fit: pitch %g, stave width %g, rails %g/%g, "
                     "stave thickness %g, space between plates %g, wall clearance %g",
                     layer_num, stave_repeat, stave_pitch, 2. * stave_dim_x, rail_low, rail_high,
                     stave_thick, frame_height - 2. * frame_thickness, wall_clearance);
            throw std::runtime_error("Staves do not fit in the layer.");
          }
        }

        // Each stave position holds n_halves identical half-length staves, back to back along the
        // stave length (local y). Half h of stave position j (both from 0, h = 0 at -y) is readout
        // stave index n_halves * j + h. A half-stave carries n_modules modules at a pitch of module
        // length + module gap, centred in the half-stave (see EcalBarrel_Module_gap). Modules and
        // chips count from the half-stave's local -y end; halves in the +y half of the stave are
        // rotated by 180 degrees about the layer normal, so in every half module 0 and chip 0 are
        // on the outer side of the stave (chip-local x and y point the other way there).
        xml_comp_t x_halves          = x_stave.child(_Unicode(half_staves));
        const int n_halves           = x_halves.number();
        const int n_modules          = x_halves.attr<int>(_Unicode(modules));
        const double module_gap      = x_halves.gap();
        const std::string module_str = x_stave.moduleStr();
        if (!modules.count(module_str) || n_halves < 1 || n_modules < 1) {
          printout(ERROR, "BarrelCalorimeterImaging",
                   "layer %d: unknown module %s, or %d half-staves with %d modules", layer_num,
                   module_str.c_str(), n_halves, n_modules);
          throw std::runtime_error("Logics error in building staves.");
        }
        const ModuleInfo& mod     = modules[module_str];
        const double half_dim_y   = stave_dim_y / n_halves;
        const double module_span  = n_modules * mod.length + (n_modules - 1) * module_gap;
        const double half_end_gap = (2. * half_dim_y - module_span) / 2.;

        // Check that the modules fit the half-stave and the indices fit the readout
        if (module_gap < 0. || half_end_gap < 0. || mod.width > 2. * stave_dim_x ||
            mod.thickness > stave_thick) {
          printout(ERROR, "BarrelCalorimeterImaging",
                   "layer %d: %d modules of %g x %g x %g with gap %g do not fit a half-stave of "
                   "%g x %g x %g (end gap %g)",
                   layer_num, n_modules, mod.width, mod.length, mod.thickness, module_gap,
                   2. * stave_dim_x, 2. * half_dim_y, stave_thick, half_end_gap);
          throw std::runtime_error("Modules do not fit in the stave.");
        }
        check_id_range(layer_num, "stave", stave_index + n_halves * stave_repeat - 1);
        check_id_range(layer_num, "module", n_modules - 1);
        check_id_range(layer_num, "chip", mod.n_chips - 1);
        printout(DEBUG, "BarrelCalorimeterImaging",
                 "layer %d: %d x %d half-staves of %d modules, module gap %g, end gap %g",
                 layer_num, stave_repeat, n_halves, n_modules, module_gap, half_end_gap);

        // Make one half-stave (or reuse an identical one): modules with their chips, and the
        // DetElement tree half-stave > module > chip > sensitive slice
        const std::string half_key =
            Form("%s|%.9g|%.9g|%.9g|%d|%.9g|%s|%s|%s", module_str.c_str(), 2. * stave_dim_x,
                 2. * half_dim_y, stave_thick, n_modules, module_gap, x_stave.regionStr().c_str(),
                 x_stave.limitsStr().c_str(), x_stave.visStr().c_str());
        if (!half_staves.count(half_key)) {
          HalfStaveInfo& half = half_staves[half_key];
          half.volume         = Volume(Form("stave_half%zu", half_staves.size()),
                                       Box(stave_dim_x, half_dim_y, stave_dim_z), desc.air());
          half.volume.setAttributes(desc, x_stave.regionStr(), x_stave.limitsStr(),
                                    x_stave.visStr());
          // No parent, since it will be added on cloning and placement
          half.element = DetElement("stave", detector_id);
          for (int i = 0; i < n_modules; ++i) {
            const double y =
                -half_dim_y + half_end_gap + mod.length / 2. + i * (mod.length + module_gap);
            PlacedVolume module_physvol = half.volume.placeVolume(mod.volume, Position(0, y, 0));
            module_physvol.addPhysVolID("module", i);
            DetElement module_element(half.element, Form("module%d", i), i);
            module_element.setPlacement(module_physvol);
            for (int k = 0; k < mod.n_chips; ++k) {
              DetElement chip_element(module_element, Form("chip%d", k), k);
              chip_element.setPlacement(mod.chip_physvols[k]);
              // Add sensitive volumes
              for (auto& sensitive_physvol : sensitives[mod.chip]) {
                DetElement sensitive_element(chip_element, sensitive_physvol.volume().name(), k);
                sensitive_element.setPlacement(sensitive_physvol);

                auto& sensitive_element_params =
                    DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(
                        sensitive_element);
                sensitive_element_params.set<std::string>("axis_definitions", "XYZ");
              }
            }
          }
        }
        const HalfStaveInfo& half = half_staves[half_key];

        // Loop over number of repeats for this stave
        for (int stave_j = 0; stave_j < stave_repeat; stave_j++) {
          for (int half_h = 0; half_h < n_halves; ++half_h) {
            // Half-stave placement
            const int stave_id  = stave_index + n_halves * stave_j + half_h;
            const double half_y = -stave_dim_y + (2 * half_h + 1) * half_dim_y;
            Position half_pos(stave_pos_x, half_y, stave_pos_z(stave_j));
            RotationZ half_rot(half_y > 0. ? M_PI : 0.);
            PlacedVolume half_physvol =
                layer_volume.placeVolume(half.volume, Transform3D(half_rot, half_pos));
            half_physvol.addPhysVolID("stave", stave_id);
            DetElement half_element = half.element.clone(Form("stave%d", stave_id), stave_id);
            half_element.setPlacement(half_physvol);
            layer_element.add(half_element);
          }

          // Increment X position of stave
          stave_pos_x += stave_pitch;
        }
        stave_index += n_halves * stave_repeat;
      }

      // Place frame if defined.
      if (layer_has_frame) {
        xml_comp_t x_frame  = x_layer.child(_Unicode(frame));
        auto frame_material = desc.material(x_frame.materialStr());

        std::string frame1_name = Form("frame_inner%d", layer_num);
        double frame1_thick     = frame_thickness;
        double frame1_trd_x1 =
            layer_dim_x + (layer_thickness / 2 - frame_height / 2) * tan_half_dphi;
        double frame1_trd_x2 =
            layer_dim_x + (layer_thickness / 2 - frame_height / 2 + frame1_thick) * tan_half_dphi;
        double frame1_trd_y1 = layer_trd_y1;
        double frame1_trd_y2 = frame1_trd_y1;
        double frame1_trd_z  = frame1_thick / 2.;
        Trapezoid frame1_shape(frame1_trd_x1, frame1_trd_x2, frame1_trd_y1, frame1_trd_y2,
                               frame1_trd_z);
        Volume frame1_volume(frame1_name, frame1_shape, frame_material);
        layer_volume.placeVolume(frame1_volume,
                                 Position(0, 0, -frame_height / 2.0 + frame_thickness / 2.0));

        std::string frame2_name = Form("frame_outer%d", layer_num);
        double frame2_thick     = frame_thickness;
        double frame2_trd_x1 =
            layer_dim_x + (layer_thickness / 2 + frame_height / 2 - frame2_thick) * tan_half_dphi;
        double frame2_trd_x2 =
            layer_dim_x + (layer_thickness / 2 + frame_height / 2) * tan_half_dphi;
        double frame2_trd_y1 = layer_trd_y1;
        double frame2_trd_y2 = frame2_trd_y1;
        double frame2_trd_z  = frame2_thick / 2.;
        Trapezoid frame2_shape(frame2_trd_x1, frame2_trd_x2, frame2_trd_y1, frame2_trd_y2,
                               frame2_trd_z);
        Volume frame2_volume(frame2_name, frame2_shape, frame_material);
        layer_volume.placeVolume(frame2_volume,
                                 Position(0, 0, +frame_height / 2.0 - frame_thickness / 2.0));

        // Side walls between the plates: parallelogram slabs along the trapezoid sides, with the
        // outer face on the side and the full height between the plates
        if (frame_side > 0.) {
          const double wall_dz = (frame_height - 2. * frame_thickness) / 2.;
          const double wall_x =
              layer_dim_x + (layer_thickness / 2.) * tan_half_dphi - frame_side_dx / 2.;
          for (int sign : {-1, 1}) {
            // phi = 0 slants the wall towards +x with z, phi = pi towards -x
            Trap wall_shape(wall_dz, half_dphi, sign > 0 ? 0. : M_PI, layer_trd_y1,
                            frame_side_dx / 2., frame_side_dx / 2., 0., layer_trd_y1,
                            frame_side_dx / 2., frame_side_dx / 2., 0.);
            Volume wall_volume(Form("frame_side_%s%d", sign < 0 ? "left" : "right", layer_num),
                               wall_shape, frame_material);
            layer_volume.placeVolume(wall_volume, Position(sign * wall_x, 0, 0));
          }
        }
      }

      // Place layer into sector
      PlacedVolume layer_physvol = (layer_num == 1)
                                       ? sector1_volume.placeVolume(layer_volume, layer_pos)
                                       : sector26_volume.placeVolume(layer_volume, layer_pos);
      layer_physvol.addPhysVolID("layer", layer_num);
      layer_element.setPlacement(layer_physvol);

      // Increment to next layer Z position. Do not add space_between for the last layer
      layer_pos_z += layer_thickness;
      if (layer_j < layer_repeat - 1) {
        layer_pos_z += layer_space_between;
      }
      ++layer_num;
    }
  }

  // Phi start for a sector.
  double phi = M_PI / nsides;

  // Create nsides sectors.
  for (int sector_num = 0; sector_num < nsides; sector_num++, phi -= dphi) {
    // Compute the sector position
    Transform3D tr(RotationZYX(0, phi, M_PI * 0.5), Translation3D(0, 0, 0));

    // Place 1st layer
    PlacedVolume sector1_physvol = envelope1_volume.placeVolume(sector1_volume, tr);
    sector1_physvol.addPhysVolID("sector", sector_num + 1);
    DetElement sd1 =
        (sector_num == 0) ? sector1_element : sector1_element.clone(Form("sector1_%d", sector_num));
    sd1.setPlacement(sector1_physvol);
    envelope1_element.add(sd1);

    // Place layers 2 through 6
    PlacedVolume sector26_physvol = envelope26_volume.placeVolume(sector26_volume, tr);
    sector26_physvol.addPhysVolID("sector", sector_num + 1);
    DetElement sd26 = (sector_num == 0) ? sector26_element
                                        : sector26_element.clone(Form("sector26_%d", sector_num));
    sd26.setPlacement(sector26_physvol);
    envelope26_element.add(sd26);
  }

  // Place envelope
  PlacedVolume envelope1_physvol = detector_volume.placeVolume(envelope1_volume);
  envelope1_element.setPlacement(envelope1_physvol);
  sdet.add(envelope1_element);

  PlacedVolume envelope26_physvol = detector_volume.placeVolume(envelope26_volume);
  envelope26_element.setPlacement(envelope26_physvol);
  sdet.add(envelope26_element);

  return sdet;
}

// simple aluminum sheet cover
// dimensions: (inner r, position in z, length, phi)
static void buildSupport(Detector& desc, Volume& sector_volume, xml_comp_t x_support,
                         const std::tuple<double, double, double, double>& dimensions) {
  auto [inner_r, pos_z, sector_length, hphi] = dimensions;
  double support_thickness = getAttrOrDefault(x_support, _Unicode(thickness), 0.5 * cm);
  auto material            = desc.material(x_support.materialStr());
  double trd_y             = sector_length / 2.;
  double trd_x1_support    = std::tan(hphi) * pos_z;
  double trd_x2_support    = std::tan(hphi) * (pos_z + support_thickness);

  Trapezoid s_shape(trd_x1_support, trd_x2_support, trd_y, trd_y, support_thickness / 2.);
  Volume s_vol("support_layer", s_shape, material);
  s_vol.setVisAttributes(desc.visAttributes(x_support.visStr()));
  sector_volume.placeVolume(s_vol, Position(0.0, 0.0, pos_z + support_thickness / 2.));
}
DECLARE_DETELEMENT(epic_EcalBarrelImaging, create_detector)
