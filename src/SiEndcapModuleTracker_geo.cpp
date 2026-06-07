// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Shujie Li

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
#include <tuple>
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
  double x_offset{0.0};
  double y_offset{0.0};
  int x_repeat{1};
  bool rsu_twelve_tile_pattern{false};
  bool rsu_end_electronics{false};
  bool lec_after_rsu{false};
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
struct ModuleRow {
  string disk_key;
  string module_name;
  double x_min{0.0};
  double y_min{0.0};
  double x_size{0.0};
  double y_size{0.0};
  double dz{0.0};
  bool facing_positive_z{true};
  string handedness;
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
struct ModulePrototype {
  Volume volume;
  Placements sensitives;
  vector<VolPlane> surfaces;
};

// Layer-local corrugated carbon frame configuration from:
// <frame type="corrugated" material="..." thickness="..." h="..."
//        theta="..." d="..." y0="..." z_center_offset="..." vis="..."/>
// The frame is an insensitive layer support component. It reuses the same
// DiskBoundary as the tiled modules, so XML users do not duplicate disk geometry.
struct CorrugatedFrameConfig {
  string name{"frame"};
  string material{"CarbonFiber"};
  string vis{"SVTSupportVis"};
  double thickness{0.2 * mm};
  double height{4.2 * mm};
  double theta{36.0 * degree};
  double half_pitch{15.8 * mm};
  double y0{0.0};
  double z_center_offset{0.0};
  string rowwise_placement;
};

// One positive-y corrugation cell. The builder mirrors nonzero rows to negative y.
struct CorrugationRow {
  double row_y{0.0};
  double height{0.0};
  double theta{0.0};
  double half_pitch{0.0};
  int csv_line{0};
};

// Axis-aligned y footprint used to require full containment in the disk x-y cross-section.
struct FrameFootprint {
  double y_min{0.0};
  double y_max{0.0};
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

bool parse_handedness(const string& value, string& handedness) {
  string lower;
  lower.reserve(value.size());
  for (char ch : value) {
    lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  }

  if (lower.empty()) {
    handedness = "";
    return true;
  }
  if (lower == "left" || lower == "right") {
    handedness = lower;
    return true;
  }
  return false;
}

// to do: add RSU structure, as well as detailed of the module.
// Build the small set of supported tiled-disk module types.
map<string, ModuleTemplate> builtin_module_templates(Detector& description) {
  // Workflow note:
  // Module footprints and stack-up for the SVT disk RSU modules are defined here in C++.
  // The compact XML only names modules in the placement CSV; it intentionally does not
  // carry local <module> blocks for EIC_LAS_6RSU / EIC_LAS_5RSU.
  auto build = [&](const string& name, double x_size, double y_size) {
    ModuleTemplate module_template;
    module_template.name   = name;
    module_template.vis    = "SVTModuleVis";
    module_template.x_size = x_size;
    module_template.y_size = y_size;

    const array<ComponentTemplate, 3> components{{
        {description.constant<double>("SiEndcapModuleFlatSheetCF_thickness"), "CarbonFiber",
         "SVTSupportVis", false, -1.0, -1.0},
        {description.constant<double>("SiEndcapModuleGlue_thickness"), "SVT_Endcap_Glue",
         "SVTGlueVis", false, -1.0, -1.0},
        {description.constant<double>("SiTrackerSensor_thickness"), "Silicon", "SVTSensorVis", true,
         -1.0, -1.0},
    }};

    for (const auto& component : components) {
      module_template.total_thickness += component.thickness;
      module_template.components.push_back(component);
    }
    return module_template;
  };

  auto build_corrugated_6rsu = [&](const string& name) {
    ModuleTemplate module_template;
    module_template.name   = name;
    module_template.vis    = "TrackerModuleVis";
    module_template.x_size = description.constant<double>("SiEndcapModule6RSU_package_length");
    module_template.y_size = description.constant<double>("SiEndcapModule_width_corrugated");

    const double rsu_chain_length = 6.0 * description.constant<double>("SiEndcapRSU_length");
    const double sensor_x_offset =
        -module_template.x_size / 2.0 +
        description.constant<double>("SiEndcapModule6RSU_left_extension") +
        description.constant<double>("SiEndcapModule6RSU_sensor_left_margin") +
        rsu_chain_length / 2.0;
    const array<ComponentTemplate, 3> components{{
        {description.constant<double>("SiEndcapModuleCF_thickness"), "CarbonFiber", "SVTSupportVis",
         false, -1.0, -1.0},
        {description.constant<double>("SiEndcapAdhesive_thickness"), "SVT_Endcap_Glue",
         "SVTGlueVis", false, -1.0, -1.0},
        {description.constant<double>("SiEndcapSensor_thickness"), "Silicon", "SVTSensorVis", true,
         rsu_chain_length, description.constant<double>("SiEndcapRSU_width"), sensor_x_offset, 0.0,
         6, true, true},
    }};

    for (const auto& component : components) {
      module_template.total_thickness += component.thickness;
      module_template.components.push_back(component);
    }
    return module_template;
  };

  // hard-coded EIC-LAS module size.
  map<string, ModuleTemplate> module_templates;
  module_templates.emplace("EIC_LAS_6RSU", build("EIC_LAS_6RSU", 130.0 * mm, 30.0 * mm));
  module_templates.emplace("EIC_LAS_5RSU", build("EIC_LAS_5RSU", 105.0 * mm, 30.0 * mm));
  module_templates.emplace("EIC_LAS_6RSU_CORR", build_corrugated_6rsu("EIC_LAS_6RSU_CORR"));
  return module_templates;
}

double corrugated_6rsu_sensor_x_offset(Detector& description, const ModuleTemplate& module_template,
                                       const string& handedness) {
  const double rsu_chain_length = 6.0 * description.constant<double>("SiEndcapRSU_length");
  double end_extension          = description.constant<double>("SiEndcapModule6RSU_left_extension");
  double sensor_margin = description.constant<double>("SiEndcapModule6RSU_sensor_left_margin");
  if (handedness == "right") {
    end_extension = description.constant<double>("SiEndcapModule6RSU_right_extension");
    sensor_margin = description.constant<double>("SiEndcapModule6RSU_sensor_right_margin");
  }
  return -module_template.x_size / 2.0 + end_extension + sensor_margin + rsu_chain_length / 2.0;
}

ModuleTemplate with_corrugated_handedness(Detector& description, ModuleTemplate module_template,
                                          const string& handedness) {
  if (module_template.name != "EIC_LAS_6RSU_CORR" || handedness.empty()) {
    return module_template;
  }

  const double sensor_x_offset =
      corrugated_6rsu_sensor_x_offset(description, module_template, handedness);
  for (auto& component : module_template.components) {
    if (component.sensitive) {
      component.x_offset = sensor_x_offset;
    }
    if (component.rsu_end_electronics) {
      component.lec_after_rsu = handedness == "right";
    }
  }
  return module_template;
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

// Check whether a module rectangle fits inside the layer cross-section. For layers
// with explicit openings, reject any rectangle intersecting either beampipe
// exclusion circle, not just modules with corners inside the hole.
bool module_inside_disk(const ModuleRow& row, const DiskBoundary& disk) {
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

  if (disk.has_beampipe_opening && (rectangle_intersects_opening(disk.lepton_opening) ||
                                    rectangle_intersects_opening(disk.hadron_opening))) {
    return false;
  }
  return true;
}

// Remove an x interval occupied by a circular opening from the candidate row span.
void subtract_interval(vector<pair<double, double>>& intervals, double remove_min,
                       double remove_max) {
  vector<pair<double, double>> kept;
  for (const auto& [x_min, x_max] : intervals) {
    if (remove_max <= x_min || remove_min >= x_max) {
      kept.emplace_back(x_min, x_max);
      continue;
    }
    if (remove_min > x_min) {
      kept.emplace_back(x_min, std::min(remove_min, x_max));
    }
    if (remove_max < x_max) {
      kept.emplace_back(std::max(remove_max, x_min), x_max);
    }
  }
  intervals = kept;
}

// Compute x spans where the full row footprint remains inside the annulus and
// outside the beampipe openings. Returned spans may be shortened/split, but each
// placed box is fully contained by construction.
vector<pair<double, double>> contained_x_intervals(const DiskBoundary& disk,
                                                   const FrameFootprint& footprint) {
  const double tolerance = 1.0e-6 * mm;
  const double y_abs_max = std::max(std::abs(footprint.y_min), std::abs(footprint.y_max));
  if (y_abs_max >= disk.rmax - tolerance) {
    return {};
  }

  const double outer_x = std::sqrt(disk.rmax * disk.rmax - y_abs_max * y_abs_max) - tolerance;
  if (outer_x <= 0.0) {
    return {};
  }

  vector<pair<double, double>> intervals{{-outer_x, outer_x}};

  auto subtract_circle = [&](double center_x, double center_y, double radius) {
    if (radius <= 0.0) {
      return;
    }

    double dy = 0.0;
    if (center_y < footprint.y_min) {
      dy = footprint.y_min - center_y;
    } else if (center_y > footprint.y_max) {
      dy = center_y - footprint.y_max;
    }
    if (dy >= radius) {
      return;
    }

    const double dx = std::sqrt(radius * radius - dy * dy) + tolerance;
    subtract_interval(intervals, center_x - dx, center_x + dx);
  };

  subtract_circle(0.0, 0.0, disk.rmin);
  if (disk.has_beampipe_opening) {
    subtract_circle(disk.lepton_opening.center_x, disk.lepton_opening.center_y,
                    disk.lepton_opening.radius);
    subtract_circle(disk.hadron_opening.center_x, disk.hadron_opening.center_y,
                    disk.hadron_opening.radius);
  }

  intervals.erase(std::remove_if(intervals.begin(), intervals.end(),
                                 [](const auto& interval) {
                                   return interval.second - interval.first <= 2.0e-6 * mm;
                                 }),
                  intervals.end());
  return intervals;
}

// Parse the tiled-layer frame block. Only type="corrugated" is handled here; other
// frame types can be added later without changing the layer or module placement logic.
bool parse_corrugated_frame(xml_comp_t x_frame, CorrugatedFrameConfig& config) {
  const string frame_type = getAttrOrDefault<string>(x_frame, _Unicode(type), "");
  if (frame_type != "corrugated") {
    return false;
  }

  config.name       = getAttrOrDefault<string>(x_frame, _Unicode(name), config.name);
  config.material   = getAttrOrDefault<string>(x_frame, _Unicode(material), config.material);
  config.vis        = getAttrOrDefault<string>(x_frame, _Unicode(vis), config.vis);
  config.thickness  = getAttrOrDefault(x_frame, _Unicode(thickness), config.thickness);
  config.height     = getAttrOrDefault(x_frame, _Unicode(h), config.height);
  config.theta      = getAttrOrDefault(x_frame, _Unicode(theta), config.theta);
  config.half_pitch = getAttrOrDefault(x_frame, _Unicode(d), config.half_pitch);
  config.y0         = getAttrOrDefault(x_frame, _Unicode(y0), config.y0);
  config.z_center_offset =
      getAttrOrDefault(x_frame, _Unicode(z_center_offset), config.z_center_offset);
  config.rowwise_placement =
      getAttrOrDefault<string>(x_frame, _Unicode(rowwise_placement), config.rowwise_placement);
  return true;
}

vector<CorrugationRow> load_corrugation_rows(const string& file_name, const string& disk_key) {
  std::ifstream input(file_name);
  if (!input) {
    printout(WARNING, "SiEndcapModuleTracker",
             fmt::format("could not open corrugation row CSV '{}'", file_name));
    return {};
  }

  map<string, size_t> header_index;
  bool header_loaded = false;
  vector<CorrugationRow> shared_rows;
  vector<CorrugationRow> disk_rows;
  string line;
  int line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    string trimmed = trim(line);
    if (trimmed.empty() || trimmed[0] == '#') {
      continue;
    }

    vector<string> fields = split_csv_line(line);
    if (!header_loaded) {
      for (size_t idx = 0; idx < fields.size(); ++idx) {
        header_index[fields[idx]] = idx;
      }

      const vector<string> required_headers{"disk", "row_y_mm", "h_mm", "theta_deg"};
      bool missing_header = false;
      for (const auto& header : required_headers) {
        if (!header_index.count(header)) {
          printout(WARNING, "SiEndcapModuleTracker",
                   fmt::format("corrugation row CSV '{}' is missing required header '{}'",
                               file_name, header));
          missing_header = true;
        }
      }
      if (!header_index.count("d_mm") && !header_index.count("half_pitch_mm") &&
          !header_index.count("pitch_mm")) {
        printout(WARNING, "SiEndcapModuleTracker",
                 fmt::format("corrugation row CSV '{}' needs one of d_mm, half_pitch_mm, or "
                             "pitch_mm",
                             file_name));
        missing_header = true;
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

    string row_disk       = get_field("disk");
    const bool shared_row = row_disk == "*" || row_disk == "all";
    const bool disk_row   = row_disk == disk_key;
    if (!(shared_row || disk_row)) {
      continue;
    }
    string enabled_value = get_field("enabled");
    if (!enabled_value.empty() && !parse_bool(enabled_value)) {
      continue;
    }

    try {
      CorrugationRow row;
      row.csv_line            = line_number;
      row.row_y               = std::stod(get_field("row_y_mm")) * mm;
      row.height              = std::stod(get_field("h_mm")) * mm;
      row.theta               = std::stod(get_field("theta_deg")) * degree;
      string half_pitch_value = get_field("d_mm");
      if (half_pitch_value.empty()) {
        half_pitch_value = get_field("half_pitch_mm");
      }
      if (!half_pitch_value.empty()) {
        row.half_pitch = std::stod(half_pitch_value) * mm;
      } else {
        row.half_pitch = std::stod(get_field("pitch_mm")) * mm / 2.0;
      }

      if (row.row_y < -1.0e-9 * mm) {
        printout(WARNING, "SiEndcapModuleTracker",
                 fmt::format("skipping corrugation row CSV line {} in '{}': row_y_mm must be "
                             "non-negative for mirrored-row mode",
                             line_number, file_name));
        continue;
      }
      if (disk_row) {
        disk_rows.push_back(row);
      } else {
        shared_rows.push_back(row);
      }
    } catch (const std::exception&) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping malformed corrugation row CSV line {} in '{}'", line_number,
                           file_name));
    }
  }

  auto rows = disk_rows.empty() ? shared_rows : disk_rows;
  std::sort(rows.begin(), rows.end(), [](const CorrugationRow& lhs, const CorrugationRow& rhs) {
    return lhs.row_y < rhs.row_y;
  });
  return rows;
}

// Place the corrugated carbon frame as an insensitive layer component.
// Assumptions:
// - h is the surface-to-surface peak/valley envelope.
// - d is half the pitch between same-surface flat-strip centers.
// - theta is the web angle with respect to the x-y plane.
// - y0=0 keeps the central lower flat centered at y=0; it is split there so no
//   single frame volume crosses the xz plane.
int place_corrugated_frame(Detector& description, Volume& layer_vol, const DiskBoundary& disk,
                           const CorrugatedFrameConfig& config) {
  if (config.thickness <= 0.0 || config.height <= config.thickness || config.half_pitch <= 0.0 ||
      config.theta <= 0.0) {
    printout(
        WARNING, "SiEndcapModuleTracker",
        fmt::format("skipping corrugated frame for disk '{}': invalid dimensions", disk.disk_key));
    return 0;
  }

  const double y_step       = config.height / std::tan(config.theta);
  const double flat_length  = config.half_pitch - y_step;
  const double center_delta = config.height - config.thickness;
  if (config.rowwise_placement.empty() && (flat_length <= 0.0 || center_delta <= 0.0)) {
    printout(WARNING, "SiEndcapModuleTracker",
             fmt::format("skipping corrugated frame for disk '{}': flat length is non-positive",
                         disk.disk_key));
    return 0;
  }

  Material frame_material = description.material(config.material);

  int placed_count = 0;
  int frame_index  = 0;

  auto place_piece = [&](const string& kind, double y_center, double y_half_span, double z_center,
                         double local_y_half, double local_z_half, double angle) {
    const double y_min = y_center - y_half_span;
    const double y_max = y_center + y_half_span;
    vector<tuple<double, double, double>> segments{{y_min, y_max, local_y_half}};
    if (y_min < 0.0 && y_max > 0.0) {
      if (std::abs(angle) > 1.0e-12 || std::abs(y_half_span - local_y_half) > 1.0e-9 * mm) {
        // Sloped pieces are not split because their projected y span does not map
        // linearly to the local box axis; skip them instead of crossing y=0.
        return;
      }
      segments = {{y_min, 0.0, (0.0 - y_min) / 2.0}, {0.0, y_max, (y_max - 0.0) / 2.0}};
    }

    for (const auto& [segment_y_min, segment_y_max, segment_local_y_half] : segments) {
      if (segment_y_max - segment_y_min <= 1.0e-9 * mm) {
        continue;
      }
      const double segment_y_center = (segment_y_min + segment_y_max) / 2.0;
      const FrameFootprint footprint{segment_y_min, segment_y_max};
      const auto intervals = contained_x_intervals(disk, footprint);
      for (const auto& [x_min, x_max] : intervals) {
        const double x_half = (x_max - x_min) / 2.0;
        if (x_half <= 0.0) {
          continue;
        }

        const double x_center = (x_min + x_max) / 2.0;
        const string volume_name =
            fmt::format("{}_{}_{}_{}", disk.disk_key, config.name, kind, frame_index++);
        Box frame_solid(x_half, segment_local_y_half, local_z_half);
        Volume frame_vol(volume_name, frame_solid, frame_material);
        frame_vol.setVisAttributes(description.visAttributes(config.vis));
        layer_vol.placeVolume(frame_vol,
                              Transform3D(RotationZYX(0.0, 0.0, angle),
                                          Position(x_center, segment_y_center, z_center)));
        ++placed_count;
      }
    }
  };

  auto place_corrugation_cell = [&](const string& row_tag, double lower_y, double height,
                                    double half_pitch, double theta) {
    if (height <= config.thickness || half_pitch <= 0.0 || theta <= 0.0) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping corrugation row '{}' for disk '{}': invalid dimensions",
                           row_tag, disk.disk_key));
      return;
    }

    const double row_y_step       = height / std::tan(theta);
    const double row_flat_length  = half_pitch - row_y_step;
    const double row_center_delta = height - config.thickness;
    if (row_flat_length <= 0.0 || row_center_delta <= 0.0) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping corrugation row '{}' for disk '{}': flat length is "
                           "non-positive",
                           row_tag, disk.disk_key));
      return;
    }

    const double z_low       = config.z_center_offset - row_center_delta / 2.0;
    const double z_high      = config.z_center_offset + row_center_delta / 2.0;
    const double web_angle   = std::atan2(row_center_delta, row_y_step);
    const double web_length  = std::hypot(row_center_delta, row_y_step);
    const double flat_y_half = row_flat_length / 2.0;
    const double flat_z_half = config.thickness / 2.0;
    const double web_y_half  = web_length / 2.0;
    const double web_z_half  = config.thickness / 2.0;
    const double web_y_projection =
        std::abs(std::cos(web_angle)) * web_y_half + std::abs(std::sin(web_angle)) * web_z_half;
    const double upper_y = lower_y + half_pitch;

    place_piece(row_tag + "_lower_flat", lower_y, flat_y_half, z_low, flat_y_half, flat_z_half,
                0.0);
    place_piece(row_tag + "_upper_flat", upper_y, flat_y_half, z_high, flat_y_half, flat_z_half,
                0.0);

    const double positive_web_y = lower_y + row_flat_length / 2.0 + row_y_step / 2.0;
    const double negative_web_y = lower_y - row_flat_length / 2.0 - row_y_step / 2.0;
    place_piece(row_tag + "_web_pos", positive_web_y, web_y_projection, config.z_center_offset,
                web_y_half, web_z_half, web_angle);
    place_piece(row_tag + "_web_neg", negative_web_y, web_y_projection, config.z_center_offset,
                web_y_half, web_z_half, -web_angle);
  };

  if (!config.rowwise_placement.empty()) {
    const auto rows =
        load_corrugation_rows(resolve_input_file(config.rowwise_placement), disk.disk_key);
    if (!rows.empty()) {
      for (const auto& row : rows) {
        const string row_tag = fmt::format("row{}", row.csv_line);
        place_corrugation_cell(row_tag, config.y0 + row.row_y, row.height, row.half_pitch,
                               row.theta);
        if (std::abs(row.row_y) > 1.0e-9 * mm) {
          place_corrugation_cell(row_tag + "_mirror", config.y0 - row.row_y, row.height,
                                 row.half_pitch, row.theta);
        }
      }
      return placed_count;
    }

    printout(WARNING, "SiEndcapModuleTracker",
             fmt::format("no corrugation rows from '{}' matched disk '{}'; using XML fallback",
                         config.rowwise_placement, disk.disk_key));
  }

  const int n_min =
      static_cast<int>(std::floor((-disk.rmax - config.y0) / (2.0 * config.half_pitch))) - 2;
  const int n_max =
      static_cast<int>(std::ceil((disk.rmax - config.y0) / (2.0 * config.half_pitch))) + 2;

  for (int n = n_min; n <= n_max; ++n) {
    const double lower_y = config.y0 + 2.0 * n * config.half_pitch;
    place_corrugation_cell(_toString(n, "uniform%d"), lower_y, config.height, config.half_pitch,
                           config.theta);
  }

  return placed_count;
}

// Construct the Geant4/DD4hep truth envelope for one tiled layer. Layers with
// beampipe openings use a boolean solid so the mother volume follows both the
// lepton and hadron pipe cutouts; ACTS gets explicit envelope parameters later.
Solid build_disk_solid(const string&, const DiskBoundary& disk) {
  if (!disk.has_beampipe_opening) {
    return Tube(disk.rmin, disk.rmax, disk.length / 2.0);
  }

  const double cut_half_length = disk.length;
  Solid disk_solid             = Tube(0.0, disk.rmax, disk.length / 2.0);
  auto subtract_opening        = [&](const DiskBoundary::CircularOpening& opening) {
    Tube opening_cut(0.0, opening.radius, cut_half_length);
    disk_solid = SubtractionSolid(disk_solid, opening_cut,
                                  Position(opening.center_x, opening.center_y, 0.0));
  };

  subtract_opening(disk.lepton_opening);
  const bool same_opening =
      std::abs(disk.lepton_opening.center_x - disk.hadron_opening.center_x) < 1e-9 &&
      std::abs(disk.lepton_opening.center_y - disk.hadron_opening.center_y) < 1e-9 &&
      std::abs(disk.lepton_opening.radius - disk.hadron_opening.radius) < 1e-9;
  if (!same_opening) {
    subtract_opening(disk.hadron_opening);
  }
  return disk_solid;
}

// Ensure the placed module thickness fits within the XML layer thickness.
bool module_inside_layer_z(const ModuleRow& row, const ModuleTemplate& module_template,
                           const DiskBoundary& disk) {
  const double module_half_thickness = module_template.total_thickness / 2.0;
  const double layer_half_thickness  = disk.length / 2.0;
  return (row.dz - module_half_thickness >= -layer_half_thickness) &&
         (row.dz + module_half_thickness <= layer_half_thickness);
}

// Load the placement CSV once for the detector and attach module dimensions from the
// built-in template map. This routine is intentionally tolerant: malformed or unknown
// rows are skipped with warnings instead of aborting the geometry build.
vector<ModuleRow> load_module_rows(const string& file_name,
                                   const map<string, ModuleTemplate>& module_templates) {
  vector<ModuleRow> rows;
  std::ifstream input(file_name);
  if (!input.is_open()) {
    printout(WARNING, "SiEndcapModuleTracker",
             fmt::format("unable to open module placement file '{}'", file_name));
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
      const array<string, 4> required_headers{"disk", "module", "x_min_mm", "y_min_mm"};
      bool missing_header = false;
      for (const auto& header : required_headers) {
        if (!header_index.count(header)) {
          printout(WARNING, "SiEndcapModuleTracker",
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

    ModuleRow row;
    row.disk_key    = get_field("disk");
    row.module_name = get_field("module");
    row.csv_line    = line_number;
    if (row.disk_key.empty()) {
      printout(
          WARNING, "SiEndcapModuleTracker",
          fmt::format("skipping CSV line {} with empty disk key in '{}'", line_number, file_name));
      continue;
    }
    if (row.module_name.empty()) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping CSV line {} with empty module name in '{}'", line_number,
                           file_name));
      continue;
    }
    auto module_iter = module_templates.find(row.module_name);
    if (module_iter == module_templates.end()) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping CSV line {} with unknown module '{}' in '{}'", line_number,
                           row.module_name, file_name));
      continue;
    }

    try {
      // The CSV carries placement only. Geometry dimensions come from the built-in
      // module template selected by the row's module name.
      row.x_min       = std::stod(get_field("x_min_mm")) * mm;
      row.y_min       = std::stod(get_field("y_min_mm")) * mm;
      row.x_size      = module_iter->second.x_size;
      row.y_size      = module_iter->second.y_size;
      string dz_value = get_field("dz_mm");
      if (!dz_value.empty()) {
        row.dz = std::stod(dz_value) * mm;
      }
    } catch (const std::exception&) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping malformed CSV line {} in '{}'", line_number, file_name));
      continue;
    }

    string facing_value = get_field("facing");
    if (!facing_value.empty() && !parse_facing(facing_value, row.facing_positive_z)) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping CSV line {} with invalid facing '{}' in '{}'", line_number,
                           facing_value, file_name));
      continue;
    }

    string handedness_value = get_field("handedness");
    if (!parse_handedness(handedness_value, row.handedness)) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping CSV line {} with invalid handedness '{}' in '{}'", line_number,
                           handedness_value, file_name));
      continue;
    }
    const string corrugated_suffix = "_CORR";
    const bool corrugated_module =
        row.module_name.size() >= corrugated_suffix.size() &&
        row.module_name.compare(row.module_name.size() - corrugated_suffix.size(),
                                corrugated_suffix.size(), corrugated_suffix) == 0;
    if (corrugated_module && row.handedness.empty()) {
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping CSV line {} for corrugated module '{}' without explicit "
                           "left/right handedness in '{}'",
                           line_number, row.module_name, file_name));
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
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("skipping CSV line {} with non-positive dimensions in '{}'", line_number,
                           file_name));
      continue;
    }
    rows.push_back(row);
  }

  return rows;
}

// Build one reusable module volume per module type and cache it. Individual placements
// later reuse this prototype volume rather than rebuilding the stack for every CSV row.
ModulePrototype build_module_prototype(Detector& description, SensitiveDetector& sens,
                                       const ModuleTemplate& module_template) {
  ModulePrototype prototype;
  Material vacuum     = description.vacuum();
  const double x_size = module_template.x_size;
  const double y_size = module_template.y_size;
  Box module_solid(x_size / 2.0, y_size / 2.0, module_template.total_thickness / 2.0);
  prototype.volume = Volume(module_template.name, module_solid, vacuum);
  prototype.volume.setVisAttributes(description.visAttributes(module_template.vis));

  double z_position       = -module_template.total_thickness / 2.0;
  double thickness_so_far = 0.0;
  int component_id        = 0;
  int sensor_id           = 1;
  for (const auto& component : module_template.components) {
    double comp_x      = component.x_override > 0.0 ? component.x_override : x_size;
    double comp_y      = component.y_override > 0.0 ? component.y_override : y_size;
    const int x_repeat = std::max(1, component.x_repeat);

    Material material = description.material(component.material);
    z_position += component.thickness / 2.0;

    auto attach_sensitive = [&](Volume& component_volume, PlacedVolume& pv) {
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
    };

    if (component.rsu_twelve_tile_pattern) {
      const double rsu_pitch_x     = comp_x / x_repeat;
      const double half_rsu_x      = rsu_pitch_x / 2.0;
      const double half_rsu_y      = comp_y / 2.0;
      const double bias_width      = description.constant<double>("SiEndcapRSU_bias_width");
      const double periphery_width = description.constant<double>("SiEndcapRSU_periphery_width");
      const double backbone_width  = description.constant<double>("SiEndcapRSU_backbone_width");
      const double powerswitch_width =
          description.constant<double>("SiEndcapRSU_powerswitch_width");
      const double tile_x =
          std::max(0.0, (half_rsu_x - backbone_width - 3.0 * powerswitch_width) / 3.0);
      const double active_y = std::max(0.0, half_rsu_y - bias_width - periphery_width);

      auto place_box_with_thickness = [&](const string& name, double box_x, double box_y,
                                          double box_thickness, double pos_x, double pos_y,
                                          const string& vis, bool sensitive) {
        if (box_x <= 0.0 || box_y <= 0.0) {
          return;
        }
        Box box_solid(box_x / 2.0, box_y / 2.0, box_thickness / 2.0);
        Volume box_volume(name, box_solid, material);
        box_volume.setVisAttributes(description.visAttributes(vis));
        PlacedVolume pv = prototype.volume.placeVolume(
            box_volume,
            Position(component.x_offset + pos_x, component.y_offset + pos_y, z_position));
        if (sensitive) {
          attach_sensitive(box_volume, pv);
        }
      };
      auto place_box = [&](const string& name, double box_x, double box_y, double pos_x,
                           double pos_y, const string& vis, bool sensitive) {
        place_box_with_thickness(name, box_x, box_y, component.thickness, pos_x, pos_y, vis,
                                 sensitive);
      };

      if (component.rsu_end_electronics) {
        const double lec_length    = description.constant<double>("SiEndcapLEC_length");
        const double rec_length    = description.constant<double>("SiEndcapREC_length");
        const double lec_width     = description.constant<double>("SiEndcapLEC_width");
        const double rec_width     = description.constant<double>("SiEndcapREC_width");
        const double lec_thickness = description.constant<double>("SiEndcapLEC_thickness");
        const double rec_thickness = description.constant<double>("SiEndcapREC_thickness");
        const double lec_x         = component.lec_after_rsu ? comp_x / 2.0 + lec_length / 2.0
                                                             : -comp_x / 2.0 - lec_length / 2.0;
        const double rec_x         = component.lec_after_rsu ? -comp_x / 2.0 - rec_length / 2.0
                                                             : comp_x / 2.0 + rec_length / 2.0;

        place_box_with_thickness(_toString(component_id, "component%d_lec"), lec_length, lec_width,
                                 lec_thickness, lec_x, 0.0, "SVTElectronicsVis", false);
        place_box_with_thickness(_toString(component_id, "component%d_rec"), rec_length, rec_width,
                                 rec_thickness, rec_x, 0.0, "SVTElectronicsVis", false);
      }

      for (int repeat_idx = 0; repeat_idx < x_repeat; ++repeat_idx) {
        const double rsu_x_min = -comp_x / 2.0 + repeat_idx * rsu_pitch_x;
        const string base_name =
            _toString(component_id, "component%d") + _toString(repeat_idx + 1, "_rsu%d");

        for (int x_half = 0; x_half < 2; ++x_half) {
          const double half_x_min        = rsu_x_min + x_half * half_rsu_x;
          const double backbone_x_center = half_x_min + backbone_width / 2.0;

          place_box(base_name + _toString(x_half + 1, "_xhalf%d_backbone"), backbone_width, comp_y,
                    backbone_x_center, 0.0, "SVTReadoutVis", false);

          for (int tile_idx = 0; tile_idx < 3; ++tile_idx) {
            const double tile_x_min =
                half_x_min + backbone_width + tile_idx * (tile_x + powerswitch_width);
            const double tile_x_center        = tile_x_min + tile_x / 2.0;
            const double powerswitch_x_center = tile_x_min + tile_x + powerswitch_width / 2.0;

            place_box(base_name + _toString(x_half + 1, "_xhalf%d") +
                          _toString(tile_idx + 1, "_tile%d_powerswitch"),
                      powerswitch_width, comp_y, powerswitch_x_center, 0.0, "SVTReadoutVis", false);

            for (int y_half = 0; y_half < 2; ++y_half) {
              const double half_y_min             = -comp_y / 2.0 + y_half * half_rsu_y;
              const bool lower_half               = y_half == 0;
              const double first_passive_width    = lower_half ? periphery_width : bias_width;
              const double second_passive_width   = lower_half ? bias_width : periphery_width;
              const double first_passive_y_center = half_y_min + first_passive_width / 2.0;
              const double active_y_center = half_y_min + first_passive_width + active_y / 2.0;
              const double second_passive_y_center =
                  half_y_min + first_passive_width + active_y + second_passive_width / 2.0;
              const string tile_name = base_name + _toString(x_half + 1, "_xhalf%d") +
                                       _toString(tile_idx + 1, "_tile%d") +
                                       _toString(y_half + 1, "_yhalf%d");

              place_box(tile_name + (lower_half ? "_periphery" : "_bias"), tile_x,
                        first_passive_width, tile_x_center, first_passive_y_center, "SVTReadoutVis",
                        false);
              place_box(tile_name + "_sensor", tile_x, active_y, tile_x_center, active_y_center,
                        component.vis, true);
              place_box(tile_name + (lower_half ? "_bias" : "_periphery"), tile_x,
                        second_passive_width, tile_x_center, second_passive_y_center,
                        "SVTReadoutVis", false);
            }
          }
        }
      }
    } else {
      string component_name = _toString(component_id, "component%d");
      Box comp_solid(comp_x / 2.0, comp_y / 2.0, component.thickness / 2.0);
      Volume component_volume(component_name, comp_solid, material);
      component_volume.setVisAttributes(description.visAttributes(component.vis));

      PlacedVolume pv = prototype.volume.placeVolume(
          component_volume, Position(component.x_offset, component.y_offset, z_position));
      if (component.sensitive) {
        attach_sensitive(component_volume, pv);
      }
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
  // 1. read the detector XML block and shared module-placement CSV
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
  // <module_placements .../> block can still override that file if needed later.
  string module_file   = "compact/tracking/SVT_endcap_modules.csv";
  string module_format = "csv";
  xml_comp_t x_module_placements(x_det.child("module_placements", false));
  if (x_module_placements) {
    module_file   = getAttrOrDefault<string>(x_module_placements, _Unicode(file), module_file);
    module_format = getAttrOrDefault<string>(x_module_placements, _Unicode(format), module_format);
  }
  if (module_format != "csv") {
    printout(WARNING, "SiEndcapModuleTracker",
             fmt::format("unsupported module placement format '{}' for '{}'; attempting CSV",
                         module_format, det_name));
  }
  vector<ModuleRow> module_rows;
  if (!module_file.empty()) {
    module_rows = load_module_rows(resolve_input_file(module_file), module_templates);
  } else {
    printout(
        WARNING, "SiEndcapModuleTracker",
        fmt::format("detector '{}' is missing a module_placements file declaration", det_name));
  }

  std::map<string, ModulePrototype> module_cache;
  int mod_num = 1;
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    // Read one disk layer description from XML and convert it into the placement
    // boundary object used by the rest of the plugin.
    xml_comp_t x_layer(li);
    const int layer_id = x_layer.id();

    xml_comp_t x_env = x_layer.child(_U(envelope));
    DiskBoundary disk;
    disk.layer_id = layer_id;
    disk.rmin     = x_env.attr<double>(_Unicode(rmin));
    disk.rmax     = x_env.attr<double>(_Unicode(rmax));
    disk.zstart   = x_env.attr<double>(_Unicode(zstart));
    disk.length   = x_env.attr<double>(_Unicode(length));
    disk.center_z = disk.zstart + disk.length / 2.0;
    disk.vis      = x_env.attr<string>(_Unicode(vis));
    disk.disk_key = x_layer.hasAttr(_Unicode(name)) ? x_layer.nameStr()
                                                    : det_name + "_disk" + std::to_string(layer_id);
    xml_comp_t x_beampipe_opening(x_layer.child(_Unicode(beampipe_opening), false));
    if (x_beampipe_opening) {
      // The opening describes the beampipe exclusion in the layer-local x-y plane as
      // two circles: the centered lepton pipe and the shifted hadron pipe. Both
      // circles define the subtraction envelope and remain active for module rejection.
      disk.has_beampipe_opening = true;
      disk.lepton_opening.center_x =
          getAttrOrDefault(x_beampipe_opening, _Unicode(lepton_center_x), 0.0);
      disk.lepton_opening.center_y =
          getAttrOrDefault(x_beampipe_opening, _Unicode(lepton_center_y), 0.0);
      disk.lepton_opening.radius =
          getAttrOrDefault(x_beampipe_opening, _Unicode(lepton_radius), 0.0);
      disk.hadron_opening.center_x = getAttrOrDefault(x_beampipe_opening, _Unicode(hadron_center_x),
                                                      disk.lepton_opening.center_x);
      disk.hadron_opening.center_y = getAttrOrDefault(x_beampipe_opening, _Unicode(hadron_center_y),
                                                      disk.lepton_opening.center_y);
      disk.hadron_opening.radius =
          getAttrOrDefault(x_beampipe_opening, _Unicode(hadron_radius), disk.lepton_opening.radius);
      if (disk.lepton_opening.radius <= 0.0 || disk.hadron_opening.radius <= 0.0) {
        printout(ERROR, "SiEndcapModuleTracker",
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
    // Give ACTS the simple layer envelope tolerances used by B0-style assembly
    // layers while the DD4hep volume keeps the exact beampipe subtraction solid.
    layerParams.set<double>("envelope_r_min", 0 * dd4hep::mm / dd4hep::mm);
    layerParams.set<double>("envelope_r_max", 0 * dd4hep::mm / dd4hep::mm);
    layerParams.set<double>("envelope_z_min", dd4hep::mm / dd4hep::mm);
    layerParams.set<double>("envelope_z_max", dd4hep::mm / dd4hep::mm);

    for (xml_coll_t fi(x_layer, _U(frame)); fi; ++fi) {
      // Layer frame workflow:
      // A type="corrugated" frame is an insensitive carbon support component that
      // shares this tiled layer's envelope, z placement, reflection, and openings.
      xml_comp_t x_frame(fi);
      CorrugatedFrameConfig frame_config;
      if (!parse_corrugated_frame(x_frame, frame_config)) {
        continue;
      }
      const int placed_frame_count =
          place_corrugated_frame(description, layer_vol, disk, frame_config);
      if (placed_frame_count == 0) {
        printout(WARNING, "SiEndcapModuleTracker",
                 fmt::format("no corrugated frame pieces placed for '{}' disk '{}'", det_name,
                             disk.disk_key));
      }
    }

    int placed_in_layer = 0;
    for (const auto& row : module_rows) {
      if (row.disk_key != disk.disk_key) {
        continue;
      }

      // Placement workflow:
      // 1. select rows targeting this disk
      // 2. reject rows that intersect the disk boundary or beampipe opening
      // 3. reject rows that do not fit in z
      // 4. place the requested module template
      if (!module_inside_disk(row, disk)) {
        printout(
            WARNING, "SiEndcapModuleTracker",
            fmt::format("skipping module line {} for '{}' (disk '{}'): x_min={} mm y_min={} mm "
                        "x_size={} mm y_size={} mm violates disk outer boundary/opening",
                        row.csv_line, det_name, row.disk_key, row.x_min / mm, row.y_min / mm,
                        row.x_size / mm, row.y_size / mm));
        continue;
      }
      auto module_template_it = module_templates.find(row.module_name);
      if (module_template_it == module_templates.end()) {
        printout(WARNING, "SiEndcapModuleTracker",
                 fmt::format("skipping module line {} for '{}' (disk '{}'): module '{}' not found",
                             row.csv_line, det_name, row.disk_key, row.module_name));
        continue;
      }
      ModuleTemplate module_template =
          with_corrugated_handedness(description, module_template_it->second, row.handedness);

      if (!module_inside_layer_z(row, module_template, disk)) {
        printout(WARNING, "SiEndcapModuleTracker",
                 fmt::format("skipping module line {} for '{}' (disk '{}'): dz={} mm with module "
                             "thickness {} mm exceeds layer half-thickness {} mm",
                             row.csv_line, det_name, row.disk_key, row.dz / mm,
                             module_template.total_thickness / mm, disk.length / (2.0 * mm)));
        continue;
      }

      string cache_key =
          row.handedness.empty() ? row.module_name : row.module_name + "_" + row.handedness;
      auto cache_it = module_cache.find(cache_key);
      if (cache_it == module_cache.end()) {
        cache_it =
            module_cache
                .emplace(cache_key, build_module_prototype(description, sens, module_template))
                .first;
      }
      ModulePrototype& prototype = cache_it->second;

      // Convert the CSV rectangle origin into the module center used by DD4hep.
      const double x_center = row.x_min + row.x_size / 2.0;
      const double y_center = row.y_min + row.y_size / 2.0;
      string module_name    = _toString(layer_id, "layer%d") + _toString(mod_num, "_module%d");
      module_name += reflect ? "_neg" : "_pos";
      DetElement module_element(layer_element, module_name, mod_num);
      // The optional facing column flips the module around local y so the sensor can
      // face either +z or -z within the layer volume.
      Transform3D module_transform =
          row.facing_positive_z
              ? Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(x_center, y_center, row.dz))
              : Transform3D(RotationZYX(0.0, M_PI, 0.0), Position(x_center, y_center, row.dz));
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
      printout(WARNING, "SiEndcapModuleTracker",
               fmt::format("no valid modules placed for '{}' disk '{}'", det_name, disk.disk_key));
    }
  }

  // Place the full detector assembly in its mother volume. The tiny signed z offset is
  // kept from the original pattern to avoid coincident-placement edge cases.
  PlacedVolume pv = motherVol.placeVolume(assembly, Position(0, 0, reflect ? -1.0e-9 : 1.0e-9));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_SiEndcapModuleTracker, create_detector)
