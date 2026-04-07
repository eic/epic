// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hepDetectorHelper.h"
#include "XML/Utilities.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <fmt/core.h>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

namespace {

using Placements = vector<PlacedVolume>;

// One entry in the hard-coded module stack definition.
struct ComponentTemplate {
  double thickness{0.0};
  string material;
  string vis;
  bool sensitive{false};
  double x_override{-1.0};
  double y_override{-1.0};
};

// One named RSU module type used by the tiled disk CSV.
struct ModuleTemplate {
  string name;
  string vis;
  vector<ComponentTemplate> components;
  double total_thickness{0.0};
  double x_size{0.0};
  double y_size{0.0};
};

// One placement row loaded from the CSV file.
struct TileRow {
  string disk_key;
  string module_name;
  double x_min{0.0};
  double y_min{0.0};
  double x_size{0.0};
  double y_size{0.0};
  double dz{0.0};
  bool facing_positive_z{true};
  bool enabled{true};
  int csv_line{0};
};

// The per-layer geometric boundary used for placement filtering and solid building.
struct DiskBoundary {
  string disk_key;
  int layer_id{0};
  double rmin{0.0};
  double rmax{0.0};
  double zstart{0.0};
  double length{0.0};
  double center_z{0.0};
  string vis;
  struct CircularOpening {
    double center_x{0.0};
    double center_y{0.0};
    double radius{0.0};
  };
  bool has_beampipe_opening{false};
  CircularOpening lepton_opening;
  CircularOpening hadron_opening;
};

// Cached built module volume and its sensitive surfaces, reused across many placements.
struct TilePrototype {
  Volume volume;
  Placements sensitives;
  vector<VolPlane> surfaces;
};

// Trim whitespace around CSV fields and free-form text values.
string trim(string value) {
  auto not_space = [](unsigned char ch) { return !std::isspace(ch); };
  auto begin_it  = std::find_if(value.begin(), value.end(), not_space);
  if (begin_it == value.end()) {
    return "";
  }
  auto end_it = std::find_if(value.rbegin(), value.rend(), not_space).base();
  return string(begin_it, end_it);
}

// Minimal CSV splitter for the placement table. The current workflow only needs
// plain comma-separated fields without embedded quoted commas.
vector<string> split_csv_line(const string& line) {
  vector<string> fields;
  string field;
  std::stringstream stream(line);
  while (std::getline(stream, field, ',')) {
    fields.push_back(trim(field));
  }
  if (!line.empty() && line.back() == ',') {
    fields.emplace_back("");
  }
  return fields;
}

// Shared boolean parser used for CSV enabled flags.
bool parse_bool(const string& value) {
  string lower;
  lower.reserve(value.size());
  for (char ch : value) {
    lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  }
  return !(lower == "0" || lower == "false" || lower == "no" || lower == "off");
}

// Interpret the optional CSV facing column into a local module orientation.
bool parse_facing(const string& value, bool& facing_positive_z) {
  string lower;
  lower.reserve(value.size());
  for (char ch : value) {
    lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  }

  if (lower.empty() || lower == "+z" || lower == "z+" || lower == "positive" ||
      lower == "positive_z" || lower == "plusz") {
    facing_positive_z = true;
    return true;
  }
  if (lower == "-z" || lower == "z-" || lower == "negative" || lower == "negative_z" ||
      lower == "minusz") {
    facing_positive_z = false;
    return true;
  }
  return false;
}

// Build the small set of supported tiled-disk module types.
map<string, ModuleTemplate> builtin_module_templates(Detector& description) {
  // Workflow note:
  // Tile footprints and stack-up for the SVT disk RSU modules are defined here in C++.
  // The tiled XML only names modules in the placement CSV; it intentionally does not
  // carry local <module> blocks for EIC_LAS_6RSU / EIC_LAS_5RSU.
  auto build = [&](const string& name, double x_size, double y_size) {
    ModuleTemplate module_template;
    module_template.name   = name;
    module_template.vis    = "TrackerModuleVis";
    module_template.x_size = x_size;
    module_template.y_size = y_size;

    const array<ComponentTemplate, 3> components{{
        {description.constant<double>("SiEndcapModuleFlatSheetCF_thickness"), "CarbonFiber",
         "TrackerSupportVis", false, -1.0, -1.0},
        {description.constant<double>("SiEndcapModuleGlue_thickness"), "SVT_Endcap_Glue",
         "TrackerServiceVis", false, -1.0, -1.0},
        {description.constant<double>("SiTrackerSensor_thickness"), "Silicon",
         "TrackerLayerVis", true, -1.0, -1.0},
    }};

    for (const auto& component : components) {
      module_template.total_thickness += component.thickness;
      module_template.components.push_back(component);
    }
    return module_template;
  };

  map<string, ModuleTemplate> module_templates;
  module_templates.emplace("EIC_LAS_6RSU", build("EIC_LAS_6RSU", 130.0 * mm, 30.0 * mm));
  module_templates.emplace("EIC_LAS_5RSU", build("EIC_LAS_5RSU", 105.0 * mm, 30.0 * mm));
  return module_templates;
}

// Resolve a CSV path either directly or relative to DETECTOR_PATH, matching the
// usual EPIC installed-layout workflow.
string resolve_input_file(const string& configured_path) {
  namespace fs = std::filesystem;
  fs::path direct(configured_path);
  if (fs::exists(direct)) {
    return direct.string();
  }

  auto detector_path_env = std::getenv("DETECTOR_PATH");
  if (detector_path_env != nullptr) {
    fs::path detector_root(detector_path_env);
    fs::path combined = detector_root / configured_path;
    if (fs::exists(combined)) {
      return combined.string();
    }
  }

  return configured_path;
}

// Check whether a tile rectangle fits inside the layer cross-section.
// For layers with beampipe openings, this rejects any rectangle intersecting
// either circular exclusion, not just tiles with corners inside the hole.
bool tile_inside_disk(const TileRow& row, const DiskBoundary& disk) {
  const array<pair<double, double>, 4> corners{{
      {row.x_min, row.y_min},
      {row.x_min + row.x_size, row.y_min},
      {row.x_min, row.y_min + row.y_size},
      {row.x_min + row.x_size, row.y_min + row.y_size},
  }};
  const double x_max = row.x_min + row.x_size;
  const double y_max = row.y_min + row.y_size;

  auto rectangle_intersects_opening = [&](const DiskBoundary::CircularOpening& opening) {
    const double closest_x = std::clamp(opening.center_x, row.x_min, x_max);
    const double closest_y = std::clamp(opening.center_y, row.y_min, y_max);
    const double dx        = closest_x - opening.center_x;
    const double dy        = closest_y - opening.center_y;
    return dx * dx + dy * dy < opening.radius * opening.radius;
  };

  for (const auto& [x, y] : corners) {
    double r = std::hypot(x, y);
    if (r > disk.rmax) {
      return false;
    }
    if (!disk.has_beampipe_opening && r < disk.rmin) {
      return false;
    }
  }

  if (disk.has_beampipe_opening &&
      (rectangle_intersects_opening(disk.lepton_opening) ||
       rectangle_intersects_opening(disk.hadron_opening))) {
    return false;
  }
  return true;
}

// Construct the mother volume solid for one layer.
// Without a beampipe opening this is the legacy tube. With an opening, build a
// full disk and subtract the two circular exclusions that approximate the pipes.
Solid build_disk_solid(const string& layer_name, const DiskBoundary& disk) {
  if (!disk.has_beampipe_opening) {
    return Tube(disk.rmin, disk.rmax, disk.length / 2.0);
  }

  // With explicit beampipe openings, build the layer as a solid disk and subtract the
  // two cylindrical exclusions. The XML opening constants are expected to be evaluated
  // at the disk face nearest the IP for thick layers.
  const double opening_half_length = disk.length;
  Solid layer_solid                = Tube(0.0, disk.rmax, disk.length / 2.0);

  const Tube lepton_hole(0.0, disk.lepton_opening.radius, opening_half_length);
  const Transform3D lepton_tf(RotationZYX(0.0, 0.0, 0.0),
                              Position(disk.lepton_opening.center_x, disk.lepton_opening.center_y, 0.0));
  layer_solid =
      SubtractionSolid(layer_name + "_minus_lepton", layer_solid, lepton_hole, lepton_tf);

  const Tube hadron_hole(0.0, disk.hadron_opening.radius, opening_half_length);
  const Transform3D hadron_tf(RotationZYX(0.0, 0.0, 0.0),
                              Position(disk.hadron_opening.center_x, disk.hadron_opening.center_y, 0.0));
  layer_solid =
      SubtractionSolid(layer_name + "_minus_hadron", layer_solid, hadron_hole, hadron_tf);

  return layer_solid;
}

// Ensure the placed module thickness fits within the XML layer thickness.
bool tile_inside_layer_z(const TileRow& row, const ModuleTemplate& module_template,
                         const DiskBoundary& disk) {
  const double module_half_thickness = module_template.total_thickness / 2.0;
  const double layer_half_thickness  = disk.length / 2.0;
  return (row.dz - module_half_thickness >= -layer_half_thickness) &&
         (row.dz + module_half_thickness <= layer_half_thickness);
}

// Load the placement CSV once for the detector and attach module dimensions from the
// built-in template map. This routine is intentionally tolerant: malformed or unknown
// rows are skipped with warnings instead of aborting the geometry build.
vector<TileRow> load_tile_rows(const string& file_name,
                               const map<string, ModuleTemplate>& module_templates) {
  vector<TileRow> rows;
  std::ifstream input(file_name);
  if (!input.is_open()) {
    printout(WARNING, "TileEndcapTracker",
             fmt::format("unable to open tile placement file '{}'", file_name));
    return rows;
  }

  string line;
  int line_number = 0;
  map<string, size_t> header_index;
  bool header_loaded = false;
  while (std::getline(input, line)) {
    ++line_number;
    string stripped = trim(line);
    if (stripped.empty()) {
      continue;
    }

    // Allow a commented header line such as:
    // # disk,module,x_min_mm,y_min_mm,dz_mm,facing
    if (!header_loaded && !stripped.empty() && stripped[0] == '#') {
      stripped = trim(stripped.substr(1));
    } else if (!stripped.empty() && stripped[0] == '#') {
      continue;
    }

    vector<string> fields = split_csv_line(stripped);
    if (!header_loaded) {
      for (size_t idx = 0; idx < fields.size(); ++idx) {
        header_index[fields[idx]] = idx;
      }
      const array<string, 4> required_headers{
          "disk", "module", "x_min_mm", "y_min_mm"};
      bool missing_header = false;
      for (const auto& header : required_headers) {
        if (!header_index.count(header)) {
          printout(WARNING, "TileEndcapTracker",
                   fmt::format("missing required CSV header '{}' in '{}'", header, file_name));
          missing_header = true;
        }
      }
      if (missing_header) {
        return {};
      }
      header_loaded = true;
      continue;
    }

    auto get_field = [&](const string& key) -> string {
      auto iter = header_index.find(key);
      if (iter == header_index.end() || iter->second >= fields.size()) {
        return "";
      }
      return fields[iter->second];
    };

    TileRow row;
    row.disk_key    = get_field("disk");
    row.module_name = get_field("module");
    row.csv_line    = line_number;
    if (row.disk_key.empty()) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("skipping CSV line {} with empty disk key in '{}'", line_number,
                           file_name));
      continue;
    }
    if (row.module_name.empty()) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("skipping CSV line {} with empty module name in '{}'", line_number,
                           file_name));
      continue;
    }
    auto module_iter = module_templates.find(row.module_name);
    if (module_iter == module_templates.end()) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("skipping CSV line {} with unknown module '{}' in '{}'", line_number,
                           row.module_name, file_name));
      continue;
    }

    try {
      // The CSV carries placement only. Geometry dimensions come from the built-in
      // module template selected by the row's module name.
      row.x_min  = std::stod(get_field("x_min_mm")) * mm;
      row.y_min  = std::stod(get_field("y_min_mm")) * mm;
      row.x_size = module_iter->second.x_size;
      row.y_size = module_iter->second.y_size;
      string dz_value = get_field("dz_mm");
      if (!dz_value.empty()) {
        row.dz = std::stod(dz_value) * mm;
      }
    } catch (const std::exception&) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("skipping malformed CSV line {} in '{}'", line_number, file_name));
      continue;
    }

    string facing_value = get_field("facing");
    if (!facing_value.empty() && !parse_facing(facing_value, row.facing_positive_z)) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("skipping CSV line {} with invalid facing '{}' in '{}'", line_number,
                           facing_value, file_name));
      continue;
    }

    string enabled_value = get_field("enabled");
    if (!enabled_value.empty()) {
      row.enabled = parse_bool(enabled_value);
    }

    if (!row.enabled) {
      continue;
    }
    if (row.x_size <= 0 || row.y_size <= 0) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("skipping CSV line {} with non-positive dimensions in '{}'",
                           line_number, file_name));
      continue;
    }
    rows.push_back(row);
  }

  return rows;
}

// Build one reusable module volume per module type and cache it. Individual placements
// later reuse this prototype volume rather than rebuilding the stack for every CSV row.
TilePrototype build_tile_prototype(Detector& description, SensitiveDetector& sens,
                                   const ModuleTemplate& module_template) {
  TilePrototype prototype;
  Material vacuum = description.vacuum();
  const double x_size = module_template.x_size;
  const double y_size = module_template.y_size;
  Box module_solid(x_size / 2.0, y_size / 2.0, module_template.total_thickness / 2.0);
  prototype.volume = Volume(module_template.name, module_solid, vacuum);
  prototype.volume.setVisAttributes(description.visAttributes(module_template.vis));

  double z_position          = -module_template.total_thickness / 2.0;
  double thickness_so_far    = 0.0;
  int component_id           = 0;
  int sensor_id              = 1;
  for (const auto& component : module_template.components) {
    double comp_x = component.x_override > 0.0 ? component.x_override : x_size;
    double comp_y = component.y_override > 0.0 ? component.y_override : y_size;

    Material material = description.material(component.material);
    string component_name = _toString(component_id, "component%d");
    Box comp_solid(comp_x / 2.0, comp_y / 2.0, component.thickness / 2.0);
    Volume component_volume(component_name, comp_solid, material);
    component_volume.setVisAttributes(description.visAttributes(component.vis));

    z_position += component.thickness / 2.0;
    PlacedVolume pv = prototype.volume.placeVolume(component_volume, Position(0, 0, z_position));
    if (component.sensitive) {
      pv.addPhysVolID("sensor", sensor_id);
      component_volume.setSensitiveDetector(sens);
      prototype.sensitives.push_back(pv);

      const double inner_thickness = thickness_so_far + component.thickness / 2.0;
      const double outer_thickness =
          module_template.total_thickness - thickness_so_far - component.thickness / 2.0;
      Vector3D u(1.0, 0.0, 0.0);
      Vector3D v(0.0, 1.0, 0.0);
      Vector3D n(0.0, 0.0, 1.0);
      Vector3D o(0.0, 0.0, 0.0);
      SurfaceType type(SurfaceType::Sensitive);
      prototype.surfaces.emplace_back(component_volume, type, inner_thickness, outer_thickness, u,
                                      v, n, o);
      ++sensor_id;
    }
    z_position += component.thickness / 2.0;
    thickness_so_far += component.thickness;
    ++component_id;
  }

  return prototype;
}

} // namespace

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  // Detector-level workflow:
  // 1. read the detector XML block and shared tile-placement CSV
  // 2. build each layer mother volume, including optional beampipe cutouts
  // 3. filter CSV rows against that layer geometry
  // 4. place cached module prototypes and attach sensitive surfaces
  xml_det_t x_det = e;
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  bool reflect    = x_det.reflect(false);

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);
  Volume motherVol = description.pickMotherVolume(sdet);
  Material air     = description.material("Air");
  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  // Detector-local XML provides layers/openings/CSV path; module types themselves are
  // shared and injected here from the built-in RSU template map.
  map<string, ModuleTemplate> module_templates = builtin_module_templates(description);

  // The tiled SVT workflow uses one shared placement CSV by default. A detector-local
  // <tile_placements .../> block can still override that file if needed later.
  string tile_file   = "compact/tracking/SVT_endcap_tiles.csv";
  string tile_format = "csv";
  xml_comp_t x_tile_placements(x_det.child("tile_placements", false));
  if (x_tile_placements) {
    tile_file   = getAttrOrDefault<string>(x_tile_placements, _Unicode(file), tile_file);
    tile_format = getAttrOrDefault<string>(x_tile_placements, _Unicode(format), tile_format);
  }
  if (tile_format != "csv") {
    printout(WARNING, "TileEndcapTracker",
             fmt::format("unsupported tile placement format '{}' for '{}'; attempting CSV",
                         tile_format, det_name));
  }
  vector<TileRow> tile_rows;
  if (!tile_file.empty()) {
    tile_rows = load_tile_rows(resolve_input_file(tile_file), module_templates);
  } else {
    printout(WARNING, "TileEndcapTracker",
             fmt::format("detector '{}' is missing a tile_placements file declaration", det_name));
  }

  std::map<string, TilePrototype> tile_cache;
  int mod_num = 1;
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    // Read one disk layer description from XML and convert it into the placement
    // boundary object used by the rest of the plugin.
    xml_comp_t x_layer(li);
    const int layer_id = x_layer.id();

    xml_comp_t x_env  = x_layer.child(_U(envelope));
    DiskBoundary disk;
    disk.layer_id   = layer_id;
    disk.rmin       = x_env.attr<double>(_Unicode(rmin));
    disk.rmax       = x_env.attr<double>(_Unicode(rmax));
    disk.zstart     = x_env.attr<double>(_Unicode(zstart));
    disk.length     = x_env.attr<double>(_Unicode(length));
    disk.center_z   = disk.zstart + disk.length / 2.0;
    disk.vis        = x_env.attr<string>(_Unicode(vis));
    disk.disk_key = x_layer.hasAttr(_Unicode(name))
                        ? x_layer.nameStr()
                        : det_name + "_disk" + std::to_string(layer_id);
    xml_comp_t x_beampipe_opening(x_layer.child(_Unicode(beampipe_opening), false));
    if (x_beampipe_opening) {
      // The opening describes the beampipe exclusion in the layer-local x-y plane as
      // two circles: the centered lepton pipe and the shifted hadron pipe.
      disk.has_beampipe_opening = true;
      disk.lepton_opening.center_x =
          getAttrOrDefault(x_beampipe_opening, _Unicode(lepton_center_x), 0.0);
      disk.lepton_opening.center_y =
          getAttrOrDefault(x_beampipe_opening, _Unicode(lepton_center_y), 0.0);
      disk.lepton_opening.radius =
          getAttrOrDefault(x_beampipe_opening, _Unicode(lepton_radius), 0.0);
      disk.hadron_opening.center_x =
          getAttrOrDefault(x_beampipe_opening, _Unicode(hadron_center_x),
                           disk.lepton_opening.center_x);
      disk.hadron_opening.center_y =
          getAttrOrDefault(x_beampipe_opening, _Unicode(hadron_center_y),
                           disk.lepton_opening.center_y);
      disk.hadron_opening.radius =
          getAttrOrDefault(x_beampipe_opening, _Unicode(hadron_radius), disk.lepton_opening.radius);
      if (disk.lepton_opening.radius <= 0.0 || disk.hadron_opening.radius <= 0.0) {
        printout(ERROR, "TileEndcapTracker",
                 fmt::format("detector '{}' disk '{}' has non-positive beampipe opening radii",
                             det_name, disk.disk_key));
        std::_Exit(EXIT_FAILURE);
      }
    }
    string layer_name = det_name + string("_layer") + to_string(layer_id);
    Volume layer_vol(layer_name, build_disk_solid(layer_name, disk), air);
    layer_vol.setVisAttributes(description.visAttributes(disk.vis));

    // The detector XML uses reflect="true" for the negative-z side. The layer solid
    // itself is built in local coordinates and then mirrored here during placement.
    PlacedVolume layer_pv;
    if (reflect) {
      layer_pv = assembly.placeVolume(
          layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -disk.center_z)));
      layer_name += "_N";
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, disk.center_z));
      layer_name += "_P";
    }
    layer_pv.addPhysVolID("layer", layer_id);

    DetElement layer_element(sdet, layer_name, layer_id);
    layer_element.setPlacement(layer_pv);
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layer_element);
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    int placed_in_layer = 0;
    for (const auto& row : tile_rows) {
      if (row.disk_key != disk.disk_key) {
        continue;
      }

      // Placement workflow:
      // 1. select rows targeting this disk
      // 2. reject rows that intersect the disk boundary or beampipe opening
      // 3. reject rows that do not fit in z
      // 4. place the requested module template
      if (!tile_inside_disk(row, disk)) {
        printout(WARNING, "TileEndcapTracker",
                 fmt::format("skipping tile line {} for '{}' (disk '{}'): x_min={} mm y_min={} mm "
                             "x_size={} mm y_size={} mm violates disk outer boundary/opening",
                             row.csv_line, det_name, row.disk_key, row.x_min / mm, row.y_min / mm,
                             row.x_size / mm, row.y_size / mm));
        continue;
      }
      auto module_template_it = module_templates.find(row.module_name);
      if (module_template_it == module_templates.end()) {
        printout(WARNING, "TileEndcapTracker",
                 fmt::format("skipping tile line {} for '{}' (disk '{}'): module '{}' not found",
                             row.csv_line, det_name, row.disk_key, row.module_name));
        continue;
      }
      const ModuleTemplate& module_template = module_template_it->second;

      if (!tile_inside_layer_z(row, module_template, disk)) {
        printout(WARNING, "TileEndcapTracker",
                 fmt::format("skipping tile line {} for '{}' (disk '{}'): dz={} mm with module "
                             "thickness {} mm exceeds layer half-thickness {} mm",
                             row.csv_line, det_name, row.disk_key, row.dz / mm,
                             module_template.total_thickness / mm, disk.length / (2.0 * mm)));
        continue;
      }

      string cache_key = row.module_name;
      auto cache_it    = tile_cache.find(cache_key);
      if (cache_it == tile_cache.end()) {
        cache_it =
            tile_cache.emplace(cache_key, build_tile_prototype(description, sens, module_template))
                .first;
      }
      TilePrototype& prototype = cache_it->second;

      // Convert the CSV rectangle origin into the module center used by DD4hep.
      const double x_center = row.x_min + row.x_size / 2.0;
      const double y_center = row.y_min + row.y_size / 2.0;
      string module_name = _toString(layer_id, "layer%d") + _toString(mod_num, "_module%d");
      module_name += reflect ? "_neg" : "_pos";
      DetElement module_element(layer_element, module_name, mod_num);
      // The optional facing column flips the module around local y so the sensor can
      // face either +z or -z within the layer volume.
      Transform3D module_transform = row.facing_positive_z
                                         ? Transform3D(RotationZYX(0.0, 0.0, 0.0),
                                                       Position(x_center, y_center, row.dz))
                                         : Transform3D(RotationZYX(0.0, M_PI, 0.0),
                                                       Position(x_center, y_center, row.dz));
      PlacedVolume module_pv = layer_vol.placeVolume(prototype.volume, module_transform);
      module_pv.addPhysVolID("module", mod_num);
      module_element.setPlacement(module_pv);

      // Reattach the cached sensitive surfaces to the concrete placed module instance.
      for (size_t idx = 0; idx < prototype.sensitives.size(); ++idx) {
        PlacedVolume sens_pv = prototype.sensitives[idx];
        DetElement comp_elt(module_element, sens_pv.volume().name(), mod_num);
        auto& comp_params =
            DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
        comp_params.set<string>("axis_definitions", "XYZ");
        comp_elt.setPlacement(sens_pv);
        volSurfaceList(comp_elt)->push_back(prototype.surfaces[idx]);
      }

      ++mod_num;
      ++placed_in_layer;
    }

    if (placed_in_layer == 0) {
      printout(WARNING, "TileEndcapTracker",
               fmt::format("no valid tiles placed for '{}' disk '{}'", det_name, disk.disk_key));
    }
  }

  // Place the full detector assembly in its mother volume. The tiny signed z offset is
  // kept from the original pattern to avoid coincident-placement edge cases.
  PlacedVolume pv = motherVol.placeVolume(assembly, Position(0, 0, reflect ? -1.0e-9 : 1.0e-9));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_TileEndcapTracker, create_detector)
