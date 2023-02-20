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
#include "Math/Point2D.h"
#include "TGeoPolygon.h"
#include "XML/Layering.h"
#include <functional>

using namespace std;
using namespace dd4hep;

typedef ROOT::Math::XYPoint Point;
// fiber placement helpers, defined below
struct FiberGrid {
  int ix = 0, iy = 0;
  vector<Point> points;
  Point mean_centroid = Point(0., 0.);
  Assembly *assembly_ptr = nullptr;

  // initialize with grid id and points
  FiberGrid(int i, int j, const vector<Point> &pts) : ix(i), iy(j), points(pts) {
    if (pts.empty()) {
      return;
    }

    double mx = 0., my = 0.;
    for (auto &p : pts) {
      mx += p.x();
      my += p.y();
    }
    mx /= static_cast<double>(pts.size());
    my /= static_cast<double>(pts.size());
    mean_centroid = Point(mx, my);
  };
};

vector<Point> fiberPositions(double r, double sx, double sz, double trx, double trz, double phi, double stol = 1e-2);
std::pair<int, int>   getNdivisions(double x, double z, double dx, double dz);
vector<FiberGrid> gridPoints(int div_x, int div_z, double x, double z, double phi);

// geometry helpers
void buildFibers(Detector& desc, SensitiveDetector& sens, Volume& mother, xml_comp_t x_fiber,
                 const std::tuple<double, double, double, double>& dimensions);
void buildSupport(Detector& desc, Volume& mother, xml_comp_t x_support,
                  const std::tuple<double, double, double, double>& dimensions);

// barrel ecal layers contained in an assembly
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  Layering   layering(e);
  xml_det_t  x_det    = e;
  Material   air      = desc.air();
  int        det_id   = x_det.id();
  string     det_name = x_det.nameStr();
  double     offset   = x_det.attr<double>(_Unicode(offset));
  xml_comp_t x_dim    = x_det.dimensions();
  int        nsides   = x_dim.numsides();
  double     inner_r  = x_dim.rmin();
  double     dphi     = (2 * M_PI / nsides);
  double     hphi     = dphi / 2;

  DetElement sdet(det_name, det_id);
  Volume     motherVol = desc.pickMotherVolume(sdet);

  Assembly     envelope(det_name);
  Transform3D  tr_global = Translation3D(0, 0, offset) * RotationZ(hphi);
  PlacedVolume env_phv   = motherVol.placeVolume(envelope, tr_global);
  sens.setType("calorimeter");

  env_phv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_phv);

  // build a single stave
  DetElement stave_det("stave0", det_id);
  Assembly   mod_vol("stave");

  // keep tracking of the total thickness
  double l_pos_z = inner_r;
  { // =====  buildBarrelStave(desc, sens, module_volume) =====
    // Parameters for computing the layer X dimension:
    double tan_hphi = std::tan(hphi);
    double l_dim_y  = x_dim.z() / 2.;

    // Loop over the sets of layer elements in the detector.
    int l_num = 1;
    for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
      xml_comp_t x_layer         = li;
      int        repeat          = x_layer.repeat();
      double     l_space_between = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_between), 0.);
      double     l_space_before  = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_before), 0.);
      l_pos_z += l_space_before;
      // Loop over number of repeats for this layer.
      for (int j = 0; j < repeat; j++) {
        string l_name      = Form("layer%d", l_num);
        double l_thickness = layering.layer(l_num - 1)->thickness(); // Layer's thickness.
        double l_dim_x     = tan_hphi * l_pos_z;
        l_pos_z += l_thickness;

        Position   l_pos(0, 0, l_pos_z - l_thickness / 2.); // Position of the layer.
        double     l_trd_x1 = l_dim_x;
        double     l_trd_x2 = l_dim_x + l_thickness * tan_hphi;
        double     l_trd_y1 = l_dim_y;
        double     l_trd_y2 = l_trd_y1;
        double     l_trd_z  = l_thickness / 2;
        Trapezoid  l_shape(l_trd_x1, l_trd_x2, l_trd_y1, l_trd_y2, l_trd_z);
        Volume     l_vol(l_name, l_shape, air);
        DetElement layer(stave_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
        int    s_num   = 1;
        double s_pos_z = -(l_thickness / 2.);
        for (xml_coll_t si(x_layer, _U(slice)); si; ++si) {
          xml_comp_t x_slice  = si;
          string     s_name   = Form("slice%d", s_num);
          double     s_thick  = x_slice.thickness();
          double     s_trd_x1 = l_dim_x + (s_pos_z + l_thickness / 2) * tan_hphi;
          double     s_trd_x2 = l_dim_x + (s_pos_z + l_thickness / 2 + s_thick) * tan_hphi;
          double     s_trd_y1 = l_trd_y1;
          double     s_trd_y2 = s_trd_y1;
          double     s_trd_z  = s_thick / 2.;
          Trapezoid  s_shape(s_trd_x1, s_trd_x2, s_trd_y1, s_trd_y2, s_trd_z);
          Volume     s_vol(s_name, s_shape, desc.material(x_slice.materialStr()));
          DetElement slice(layer, s_name, det_id);

          // build fibers
          if (x_slice.hasChild(_Unicode(fiber))) {
            buildFibers(desc, sens, s_vol, x_slice.child(_Unicode(fiber)), {s_trd_x1, s_thick, l_dim_y, hphi});
          }

          if (x_slice.isSensitive()) {
            s_vol.setSensitiveDetector(sens);
          }
          s_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

          // Slice placement.
          PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick / 2));
          slice_phv.addPhysVolID("slice", s_num);
          slice.setPlacement(slice_phv);
          // Increment Z position of slice.
          s_pos_z += s_thick;
          ++s_num;
        }

        // Set region, limitset, and vis of layer.
        l_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

        PlacedVolume layer_phv = mod_vol.placeVolume(l_vol, l_pos);
        layer_phv.addPhysVolID("layer", l_num);
        layer.setPlacement(layer_phv);
        // Increment to next layer Z position. Do not add space_between for the last layer
        if (j < repeat - 1) {
          l_pos_z += l_space_between;
        }
        ++l_num;
      }
    }
  }
  // Phi start for a stave.
  double phi = M_PI / nsides;
  // Create nsides staves.
  for (int i = 0; i < nsides; i++, phi -= dphi) { // i is module number
    // Compute the stave position
    Transform3D  tr(RotationZYX(0, phi, M_PI * 0.5), Translation3D(0, 0, 0));
    PlacedVolume pv = envelope.placeVolume(mod_vol, tr);
    pv.addPhysVolID("module", i + 1);
    DetElement sd = (i == 0) ? stave_det : stave_det.clone(Form("stave%d", i));
    sd.setPlacement(pv);
    sdet.add(sd);
  }

  // optional stave support
  if (x_det.hasChild(_U(staves))) {
    xml_comp_t x_staves = x_det.staves();
    mod_vol.setVisAttributes(desc.visAttributes(x_staves.visStr()));
    if (x_staves.hasChild(_U(support))) {
      buildSupport(desc, mod_vol, x_staves.child(_U(support)), {inner_r, l_pos_z, x_dim.z(), hphi});
    }
  }

  // Set envelope volume attributes.
  envelope.setAttributes(desc, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  return sdet;
}

void buildFibers(Detector& desc, SensitiveDetector& sens, Volume& s_vol, xml_comp_t x_fiber,
                 const std::tuple<double, double, double, double>& dimensions)
{
  auto [s_trd_x1, s_thick, s_length, hphi] = dimensions;
  double      f_radius                     = getAttrOrDefault(x_fiber, _U(radius), 0.1 * cm);
  double      f_cladding_thickness         = getAttrOrDefault(x_fiber, _Unicode(cladding_thickness), 0.0 * cm);
  double      f_spacing_x                  = getAttrOrDefault(x_fiber, _Unicode(spacing_x), 0.122 * cm);
  double      f_spacing_z                  = getAttrOrDefault(x_fiber, _Unicode(spacing_z), 0.134 * cm);
  std::string f_id_grid                    = getAttrOrDefault<std::string>(x_fiber, _Unicode(identifier_grid), "grid");
  std::string f_id_fiber = getAttrOrDefault<std::string>(x_fiber, _Unicode(identifier_fiber), "fiber");

  // Set up the readout grid for the fiber layers
  // Trapezoid is divided into segments with equal dz and equal number of divisions in x
  // Every segment is a polygon that can be attached later to the lightguide
  // The grid size is assumed to be ~2x2 cm (starting values). This is to be larger than
  // SiPM chip (for GlueX 13mmx13mm: 4x4 grid 3mmx3mm with 3600 50×50 μm pixels each)
  // See, e.g., https://arxiv.org/abs/1801.03088 Fig. 2d

  // fiber and its cladding
  double f_radius_core = f_radius - f_cladding_thickness;
  Tube   f_tube_clad(f_radius_core, f_radius, s_length);
  Volume f_vol_clad("fiber_vol", f_tube_clad, desc.material(x_fiber.materialStr()));
  Tube   f_tube_core(0, f_radius_core, s_length);
  Volume f_vol_core("fiber_core_vol", f_tube_core, desc.material(x_fiber.materialStr()));
  if (x_fiber.isSensitive()) {
    f_vol_core.setSensitiveDetector(sens);
  }
  f_vol_core.setAttributes(desc, x_fiber.regionStr(), x_fiber.limitsStr(), x_fiber.visStr());


  // Calculate number of divisions
  auto        grid_div = getNdivisions(s_trd_x1, s_thick, 2.0 * cm, 2.0 * cm);
  // Calculate polygonal grid coordinates (vertices)
  auto        grids = gridPoints(grid_div.first, grid_div.second, s_trd_x1, s_thick, hphi);
  vector<int> f_id_count(grid_div.first * grid_div.second, 0);
  auto        f_pos = fiberPositions(f_radius, f_spacing_x, f_spacing_z, s_trd_x1, s_thick, hphi);
  // a helper struct to speed up searching
  struct Fiber {
    Point pos;
    bool assigned = false;
    Fiber (const Point &p) : pos(p) {};
  };
  std::vector<Fiber> fibers(f_pos.begin(), f_pos.end());

  // build assembly for each grid and put fibers in
  for (auto &gr : grids) {
    Assembly grid_vol(Form("fiber_grid_%i_%i", gr.ix, gr.iy));

    // loop over all fibers that are not assigned to a grid
    int f_id = 1;
    for (auto &fi : fibers) {
      if (fi.assigned) {
        continue;
      }

      // use TGeoPolygon to help check if fiber is inside a grid
      TGeoPolygon poly(gr.points.size());
      vector<double> vx, vy;
      transform(gr.points.begin(), gr.points.end(), back_inserter(vx), mem_fn(&Point::x));
      transform(gr.points.begin(), gr.points.end(), back_inserter(vy), mem_fn(&Point::y));
      poly.SetXY(vx.data(), vy.data());
      poly.FinishPolygon();

      double f_xy[2] = {fi.pos.x(), fi.pos.y()};
      if (not poly.Contains(f_xy)) {
        continue;
      }

      // place fiber in grid
      auto p = fi.pos - gr.mean_centroid;
      auto core_phv = grid_vol.placeVolume(f_vol_core, Position(p.x(), p.y(), 0.));
      core_phv.addPhysVolID(f_id_fiber, f_id);
      grid_vol.placeVolume(f_vol_clad, Position(p.x(), p.y(), 0.));
      fi.assigned = true;
      f_id ++;
    }

    // only add this if this grid has fibers
    if (f_id > 1) {
        // fiber is along y-axis of the layer volume, so grids are arranged on X-Z plane
        Transform3D gr_tr(RotationZYX(0., 0., M_PI*0.5), Position(gr.mean_centroid.x(), 0., gr.mean_centroid.y()));
        auto grid_phv = s_vol.placeVolume(grid_vol, gr_tr);
        grid_phv.addPhysVolID(f_id_grid, gr.ix + gr.iy * grid_div.first + 1);
        grid_vol.ptr()->Voxelize("");
    }
  }

  /*
  // sanity check
  size_t missing_fibers = 0;
  for (auto &fi : fibers) {
    if (not fi.assigned) {
      missing_fibers++;
    }
  }
  std::cout << "built " << fibers.size() << " fibers, "
            << missing_fibers << " of them failed to find a grid" << std::endl;
  */
}

// DAWN view seems to have some issue with overlapping solids even if they were unions
// The support is now built without overlapping
void buildSupport(Detector& desc, Volume& mod_vol, xml_comp_t x_support,
                  const std::tuple<double, double, double, double>& dimensions)
{
  auto [inner_r, l_pos_z, stave_length, hphi] = dimensions;

  double support_thickness = getAttrOrDefault(x_support, _Unicode(thickness), 5. * cm);
  double beam_thickness    = getAttrOrDefault(x_support, _Unicode(beam_thickness), support_thickness / 4.);
  // sanity check
  if (beam_thickness > support_thickness / 3.) {
    std::cerr << Form("beam_thickness (%.2f) cannot be greater than support_thickness/3 (%.2f), shrink it to fit",
                      beam_thickness, support_thickness / 3.)
              << std::endl;
    beam_thickness = support_thickness / 3.;
  }
  Assembly env_vol("support_envelope");
  double   trd_y          = stave_length / 2.;
  double   trd_x1_support = std::tan(hphi) * l_pos_z;
  // FIXME trd_x2_support is filled but unused
  // double   trd_x2_support = std::tan(hphi) * (l_pos_z + support_thickness);

  double grid_size        = getAttrOrDefault(x_support, _Unicode(grid_size), 25. * cm);
  int    n_cross_supports = std::floor(trd_y - beam_thickness) / grid_size;
  // number of "beams" running the length of the stave.
  // @TODO make it configurable
  int n_beams = getAttrOrDefault(x_support, _Unicode(n_beams), 3);
  ;
  double beam_width = 2. * trd_x1_support / (n_beams + 1); // quick hack to make some gap between T beams
  double beam_gap   = getAttrOrDefault(x_support, _Unicode(beam_gap), 3. * cm);

  // build T-shape beam
  double                  beam_space_x    = beam_width + beam_gap;
  [[maybe_unused]] double beam_space_z    = support_thickness - beam_thickness;
  double                  cross_thickness = support_thickness - beam_thickness;
  double                  beam_pos_z      = beam_thickness / 2.;
  [[maybe_unused]] double beam_center_z   = support_thickness / 2. - beam_pos_z;

  Box        beam_vert_s(beam_thickness / 2., trd_y, cross_thickness / 2.);
  Box        beam_hori_s(beam_width / 2., trd_y, beam_thickness / 2.);
  UnionSolid T_beam_s(beam_hori_s, beam_vert_s, Position(0., 0., support_thickness / 2.));
  Volume     H_beam_vol("H_beam", T_beam_s, desc.material(x_support.materialStr()));
  H_beam_vol.setVisAttributes(desc, x_support.visStr());
  // place H beams first
  double beam_start_x = -(n_beams - 1) * (beam_width + beam_gap) / 2.;
  for (int i = 0; i < n_beams; ++i) {
    Position beam_pos(beam_start_x + i * (beam_width + beam_gap), 0., -support_thickness / 2. + beam_pos_z);
    env_vol.placeVolume(H_beam_vol, beam_pos);
  }

  // place central crossing beams that connects the H beams
  double cross_x = beam_space_x - beam_thickness;
  Box    cross_s(cross_x / 2., beam_thickness / 2., cross_thickness / 2.);
  Volume cross_vol("cross_center_beam", cross_s, desc.material(x_support.materialStr()));
  cross_vol.setVisAttributes(desc, x_support.visStr());
  for (int i = 0; i < n_beams - 1; ++i) {
    env_vol.placeVolume(cross_vol, Position(beam_start_x + beam_space_x * (i + 0.5), 0., beam_pos_z));
    for (int j = 1; j < n_cross_supports; j++) {
      env_vol.placeVolume(cross_vol, Position(beam_start_x + beam_space_x * (i + 0.5), -j * grid_size, beam_pos_z));
      env_vol.placeVolume(cross_vol, Position(beam_start_x + beam_space_x * (i + 0.5), j * grid_size, beam_pos_z));
    }
  }

  // place edge crossing beams that connects the neighbour support
  // @TODO: connection part is still using boolean volumes, maybe problematic to DAWN
  double           cross_edge_x = trd_x1_support + beam_start_x - beam_thickness / 2.;
  double           cross_trd_x1 = cross_edge_x + std::tan(hphi) * beam_thickness;
  double           cross_trd_x2 = cross_trd_x1 + 2. * std::tan(hphi) * cross_thickness;
  double           edge_pos_x   = beam_start_x - cross_trd_x1 / 2. - beam_thickness / 2;
  Trapezoid        cross_s2_trd(cross_trd_x1 / 2., cross_trd_x2 / 2., beam_thickness / 2., beam_thickness / 2.,
                                cross_thickness / 2.);
  Box              cross_s2_box((cross_trd_x2 - cross_trd_x1) / 4., beam_thickness / 2., cross_thickness / 2.);
  SubtractionSolid cross_s2(cross_s2_trd, cross_s2_box, Position((cross_trd_x2 + cross_trd_x1) / 4., 0., 0.));
  Volume           cross_vol2("cross_edge_beam", cross_s2, desc.material(x_support.materialStr()));
  cross_vol2.setVisAttributes(desc, x_support.visStr());
  env_vol.placeVolume(cross_vol2, Position(edge_pos_x, 0., beam_pos_z));
  env_vol.placeVolume(cross_vol2, Transform3D(Translation3D(-edge_pos_x, 0., beam_pos_z) * RotationZ(M_PI)));
  for (int j = 1; j < n_cross_supports; j++) {
    env_vol.placeVolume(cross_vol2, Position(edge_pos_x, -j * grid_size, beam_pos_z));
    env_vol.placeVolume(cross_vol2, Position(edge_pos_x, j * grid_size, beam_pos_z));
    env_vol.placeVolume(cross_vol2,
                        Transform3D(Translation3D(-edge_pos_x, -j * grid_size, beam_pos_z) * RotationZ(M_PI)));
    env_vol.placeVolume(cross_vol2,
                        Transform3D(Translation3D(-edge_pos_x, j * grid_size, beam_pos_z) * RotationZ(M_PI)));
  }

  mod_vol.placeVolume(env_vol, Position(0.0, 0.0, l_pos_z + support_thickness / 2.));
}

// Fill fiber lattice into trapezoid starting from position (0,0) in x-z coordinate system
vector<Point> fiberPositions(double r, double sx, double sz, double trx, double trz, double phi, double stol)
{
  // r      - fiber radius
  // sx, sz - spacing between fibers in x, z
  // trx    - half-length of the shorter (bottom) base of the trapezoid
  // trz    - height of the trapezoid
  // phi    - angle between z and trapezoid arm
  // stol   - spacing tolerance

  vector<Point> positions;
  int z_layers = floor((trz / 2 - r - stol) / sz); // number of layers that fits in half trapezoid-z

  double px = 0., pz = 0.;

  for (int l = -z_layers; l < z_layers + 1; l++) {
    vector<Point> xline;
    pz           = l * sz;
    double x_max = trx + (trz / 2. + pz) * tan(phi) - stol; // calculate max x at particular z_pos
    (l % 2 == 0) ? px = 0. : px = sx / 2;                   // account for spacing/2 shift

    while (px < (x_max - r)) {
      xline.push_back(Point(px, pz));
      if (px != 0.)
        xline.push_back(Point(-px, pz)); // using symmetry around x=0
      px += sx;
    }

    // Sort fiber IDs for a better organization
    sort(xline.begin(), xline.end(), [](const Point& p1, const Point& p2) { return p1.x() < p2.x(); });
    positions.insert(positions.end(), xline.begin(), xline.end());
  }
  return positions;
}

// Calculate number of divisions for the readout grid for the fiber layers
std::pair<int, int> getNdivisions(double x, double z, double dx, double dz)
{
  // x and z defined as in vector<Point> fiberPositions
  // dx, dz - size of the grid in x and z we want to get close to with the polygons
  // See also descripltion when the function is called

  double SiPMsize = 13.0 * mm;
  double grid_min = SiPMsize + 3.0 * mm;

  if (dz < grid_min) {
    dz = grid_min;
  }

  if (dx < grid_min) {
    dx = grid_min;
  }

  int nfit_cells_z = floor(z / dz);
  int n_cells_z    = nfit_cells_z;

  if (nfit_cells_z == 0)
    n_cells_z++;

  int nfit_cells_x = floor((2 * x) / dx);
  int n_cells_x    = nfit_cells_x;

  if (nfit_cells_x == 0)
    n_cells_x++;

  return std::make_pair(n_cells_x, n_cells_z);
}

// Calculate dimensions of the polygonal grid in the cartesian coordinate system x-z
vector<FiberGrid> gridPoints(int div_x, int div_z, double x, double z, double phi)
{
  // x, z and phi defined as in vector<Point> fiberPositions
  // div_x, div_z - number of divisions in x and z
  vector<FiberGrid> results;
  double dz = z / div_z;

  for (int iz = 0; iz < div_z + 1; iz++) {
    for (int ix = 0; ix < div_x + 1; ix++) {
      double A_z = -z / 2 + iz * dz;
      double B_z = -z / 2 + (iz + 1) * dz;

      double len_x_for_z        = 2 * (x + iz * dz * tan(phi));
      double len_x_for_z_plus_1 = 2 * (x + (iz + 1) * dz * tan(phi));

      double dx_for_z        = len_x_for_z / div_x;
      double dx_for_z_plus_1 = len_x_for_z_plus_1 / div_x;

      double A_x = -len_x_for_z / 2. + ix * dx_for_z;
      double B_x = -len_x_for_z_plus_1 / 2. + ix * dx_for_z_plus_1;

      double C_z = B_z;
      double D_z = A_z;
      double C_x = B_x + dx_for_z_plus_1;
      double D_x = A_x + dx_for_z;

      auto A = Point(A_x, A_z);
      auto B = Point(B_x, B_z);
      auto C = Point(C_x, C_z);
      auto D = Point(D_x, D_z);

      // vertex points filled in the clock-wise direction
      results.emplace_back(FiberGrid(ix, iz, {A, B, C, D}));
    }
  }

  return results;
}

DECLARE_DETELEMENT(epic_EcalBarrelInterlayers, create_detector)
