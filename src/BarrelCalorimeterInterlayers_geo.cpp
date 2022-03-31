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
#include "XML/Layering.h"
#include "Math/Point2D.h"
#include "TGeoPolygon.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

typedef ROOT::Math::XYPoint Point;
// fiber placement helpers, defined below
vector<vector<Point>> fiberPositions(double radius, double x_spacing, double z_spacing,
                                     double x, double z, double phi, double spacing_tol = 1e-2);
std::pair<int, int> getNdivisions(double x, double z, double dx, double dz);
vector<tuple<int, Point, Point, Point, Point>> gridPoints(int div_x, int div_z, double x, double z, double phi);

// geometry helpers
void buildFibers(Detector& desc, SensitiveDetector &sens, Volume &mother, xml_comp_t x_fiber,
                 const std::tuple<double, double, double, double> &dimensions);
void buildSupport(Detector& desc, Volume &mother, xml_comp_t x_support,
                  const std::tuple<double, double, double, double> &dimensions);


// barrel ecal layers contained in an assembly
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)  {
  Layering      layering (e);
  xml_det_t     x_det     = e;
  Material      air       = desc.air();
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  double        offset    = x_det.attr<double>(_Unicode(offset));
  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();
  double        inner_r   = x_dim.rmin();
  double        dphi      = (2*M_PI/nsides);
  double        hphi      = dphi/2;

  DetElement    sdet      (det_name, det_id);
  Volume        motherVol = desc.pickMotherVolume(sdet);

  Assembly      envelope  (det_name);
  Transform3D   tr        = Translation3D(0, 0, offset) * RotationZ(hphi);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope, tr);
  sens.setType("calorimeter");

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);

  // build a single stave
  DetElement    stave_det("stave0", det_id);
  Assembly      mod_vol("stave");

  // keep tracking of the total thickness
  double l_pos_z = inner_r;
  { // =====  buildBarrelStave(desc, sens, module_volume) =====
    // Parameters for computing the layer X dimension:
    double tan_hphi = std::tan(hphi);
    double l_dim_y  = x_dim.z()/2.;

    // Loop over the sets of layer elements in the detector.
    int l_num = 1;
    for(xml_coll_t li(x_det, _U(layer)); li; ++li)  {
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();
      double l_space_between = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_between), 0.);
      double l_space_before = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_before), 0.);
      l_pos_z += l_space_before;
      // Loop over number of repeats for this layer.
      for (int j = 0; j < repeat; j++)    {
        string l_name = Form("layer%d", l_num);
        double l_thickness = layering.layer(l_num - 1)->thickness();  // Layer's thickness.
        double l_dim_x = tan_hphi* l_pos_z;
        l_pos_z += l_thickness;

        Position   l_pos(0, 0, l_pos_z - l_thickness/2.);      // Position of the layer.
	    double l_trd_x1 = l_dim_x;
	    double l_trd_x2 = l_dim_x + l_thickness*tan_hphi;
	    double l_trd_y1 = l_dim_y;
	    double l_trd_y2 = l_trd_y1;
	    double l_trd_z  = l_thickness/2;
        Trapezoid  l_shape(l_trd_x1, l_trd_x2, l_trd_y1, l_trd_y2, l_trd_z);
        Volume     l_vol(l_name, l_shape, air);
        DetElement layer(stave_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
        int s_num = 1;
        double s_pos_z = -(l_thickness / 2.);
        for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
          xml_comp_t x_slice = si;
          string     s_name  = Form("slice%d", s_num);
          double     s_thick = x_slice.thickness();
	      double s_trd_x1 = l_dim_x + (s_pos_z + l_thickness/2)*tan_hphi;
	      double s_trd_x2 = l_dim_x + (s_pos_z + l_thickness/2 + s_thick)*tan_hphi;
	      double s_trd_y1 = l_trd_y1;
	      double s_trd_y2 = s_trd_y1;
	      double s_trd_z  = s_thick/2.;
          Trapezoid  s_shape(s_trd_x1, s_trd_x2, s_trd_y1, s_trd_y2, s_trd_z);
          Volume     s_vol(s_name, s_shape, desc.material(x_slice.materialStr()));
          DetElement slice(layer, s_name, det_id);

          // build fibers
          if (x_slice.hasChild(_Unicode(fiber))) {
              buildFibers(desc, sens, s_vol, x_slice.child(_Unicode(fiber)), {s_trd_x1, s_thick, l_dim_y, hphi});
          }

          if ( x_slice.isSensitive() ) {
            s_vol.setSensitiveDetector(sens);
          }
          s_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

          // Slice placement.
          PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick/2));
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
  for (int i = 0; i < nsides; i++, phi -= dphi)      { // i is module number
    // Compute the stave position
    Transform3D tr(RotationZYX(0, phi, M_PI*0.5), Translation3D(0, 0, 0));
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

void buildFibers(Detector& desc, SensitiveDetector &sens, Volume &s_vol, xml_comp_t x_fiber,
                 const std::tuple<double, double, double, double> &dimensions)
{
  auto [s_trd_x1, s_thick, s_length, hphi] = dimensions;
  double f_radius = getAttrOrDefault(x_fiber, _U(radius), 0.1 * cm);
  double f_spacing_x = getAttrOrDefault(x_fiber, _Unicode(spacing_x), 0.122 * cm);
  double f_spacing_z = getAttrOrDefault(x_fiber, _Unicode(spacing_z), 0.134 * cm);
  std::string f_id_grid = getAttrOrDefault(x_fiber, _Unicode(identifier_grid), "grid");
  std::string f_id_fiber = getAttrOrDefault(x_fiber, _Unicode(identifier_fiber), "fiber");

  // Set up the readout grid for the fiber layers
  // Trapezoid is divided into segments with equal dz and equal number of divisions in x
  // Every segment is a polygon that can be attached later to the lightguide
  // The grid size is assumed to be ~2x2 cm (starting values). This is to be larger than
  // SiPM chip (for GlueX 13mmx13mm: 4x4 grid 3mmx3mm with 3600 50×50 μm pixels each)
  // See, e.g., https://arxiv.org/abs/1801.03088 Fig. 2d

  // Calculate number of divisions
  auto grid_div = getNdivisions(s_trd_x1, s_thick, 2.0*cm, 2.0*cm);
  // Calculate polygonal grid coordinates (vertices)
  auto grid_vtx = gridPoints(grid_div.first, grid_div.second, s_trd_x1, s_thick, hphi);
  Tube f_tube(0, f_radius, s_length);
  Volume f_vol("fiber_vol", f_tube, desc.material(x_fiber.materialStr()));

  vector<int> f_id_count(grid_div.first*grid_div.second, 0);
  auto f_pos = fiberPositions(f_radius, f_spacing_x, f_spacing_z, s_trd_x1, s_thick, hphi);
  // std::cout << f_pos.size() << " lines, ~" << f_pos.front().size() << " fibers each line" << std::endl;
  for (size_t il = 0; il < f_pos.size(); ++il) {
    auto &line = f_pos[il];
    if (line.empty()) {
      continue;
    }
    double l_pos_y = line.front().y();
    // use assembly as intermediate volume container to reduce number of daughter volumes
    Assembly lfibers(Form("fiber_array_line_%lu", il));
    for (auto &p : line) {
      int f_grid_id = -1;
      int f_id = -1;
      // Check to which grid fiber belongs to
      for (auto &poly_vtx : grid_vtx) {
        if (p.y() != l_pos_y) {
          std::cerr << Form("Expected the same y position from a same line: %.2f, but got %.2f", l_pos_y, p.y())
                    << std::endl;
          continue;
        }
        auto [grid_id, vtx_a, vtx_b, vtx_c, vtx_d] = poly_vtx;
        double poly_x[4] = {vtx_a.x(), vtx_b.x(), vtx_c.x(), vtx_d.x()};
        double poly_y[4] = {vtx_a.y(), vtx_b.y(), vtx_c.y(), vtx_d.y()};
        double f_xy[2] = {p.x(), p.y()};

        TGeoPolygon poly(4);
        poly.SetXY(poly_x, poly_y);
        poly.FinishPolygon();

        if(poly.Contains(f_xy)) {
          f_grid_id = grid_id;
          f_id = f_id_count[grid_id];
          f_id_count[grid_id]++;
        }
      }

      if ( x_fiber.isSensitive() ) {
        f_vol.setSensitiveDetector(sens);
      }
      f_vol.setAttributes(desc, x_fiber.regionStr(), x_fiber.limitsStr(), x_fiber.visStr());

      // Fiber placement
      // Transform3D f_tr(RotationZYX(0,0,M_PI*0.5),Position(p.x(), 0, p.y()));
      // PlacedVolume fiber_phv = s_vol.placeVolume(f_vol, Position(p.x(), 0., p.y()));
      PlacedVolume fiber_phv = lfibers.placeVolume(f_vol, Position(p.x(), 0., 0.));
      fiber_phv.addPhysVolID(f_id_grid, f_grid_id + 1).addPhysVolID(f_id_fiber, f_id + 1);
    }
    lfibers.ptr()->Voxelize("");
    Transform3D l_tr(RotationZYX(0,0,M_PI*0.5),Position(0., 0, l_pos_y));
    s_vol.placeVolume(lfibers, l_tr);
  }
}

// DAWN view seems to have some issue with overlapping solids even if they were unions
// The support is now built without overlapping
void buildSupport(Detector& desc, Volume &mod_vol, xml_comp_t x_support,
                  const std::tuple<double, double, double, double> &dimensions)
{
  auto [inner_r, l_pos_z, stave_length, hphi] = dimensions;

  double support_thickness = getAttrOrDefault(x_support, _Unicode(thickness), 5. * cm);
  double beam_thickness    = getAttrOrDefault(x_support, _Unicode(beam_thickness), support_thickness/4.);
  // sanity check
  if (beam_thickness > support_thickness/3.) {
    std::cerr << Form("beam_thickness (%.2f) cannot be greater than support_thickness/3 (%.2f), shrink it to fit",
                      beam_thickness, support_thickness/3.) << std::endl;
    beam_thickness = support_thickness/3.;
  }
  double   trd_x1_support = std::tan(hphi) * l_pos_z;
  double   trd_x2_support = std::tan(hphi) * (l_pos_z + support_thickness);
  double   trd_y          = stave_length / 2.;
  Assembly env_vol        ("support_envelope");

  double grid_size        = getAttrOrDefault(x_support, _Unicode(grid_size), 25. * cm);
  int    n_cross_supports = std::floor(trd_y - beam_thickness)/grid_size;
  // number of "beams" running the length of the stave.
  // @TODO make it configurable
  int    n_beams          = getAttrOrDefault(x_support, _Unicode(n_beams), 3);;
  double beam_width       = 2. * trd_x1_support / (n_beams + 1); // quick hack to make some gap between T beams
  double beam_gap         = getAttrOrDefault(x_support, _Unicode(beam_gap), 3.*cm);

  // build T-shape beam
  double beam_space_x    = beam_width + beam_gap;
  double beam_space_z    = support_thickness - beam_thickness;
  double cross_thickness = support_thickness - beam_thickness;
  double beam_pos_z      = beam_thickness / 2.;
  double beam_center_z   = support_thickness / 2. - beam_pos_z;

  Box        beam_vert_s(beam_thickness / 2., trd_y, cross_thickness / 2.);
  Box        beam_hori_s(beam_width / 2., trd_y, beam_thickness / 2.);
  UnionSolid T_beam_s(beam_hori_s, beam_vert_s, Position(0., 0., support_thickness / 2.));
  Volume H_beam_vol("H_beam", T_beam_s, desc.material(x_support.materialStr()));
  H_beam_vol.setVisAttributes(desc, x_support.visStr());
  // place H beams first
  double beam_start_x = - (n_beams - 1) * (beam_width + beam_gap) / 2.;
  for (int i = 0; i < n_beams; ++i) {
    Position beam_pos(beam_start_x + i * (beam_width + beam_gap), 0., - support_thickness / 2. + beam_pos_z);
    env_vol.placeVolume(H_beam_vol, beam_pos);
  }

  // place central crossing beams that connects the H beams
  double cross_x = beam_space_x - beam_thickness;
  Box cross_s(cross_x / 2., beam_thickness / 2., cross_thickness / 2.);
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
  double    cross_edge_x = trd_x1_support + beam_start_x - beam_thickness / 2.;
  double    cross_trd_x1 = cross_edge_x + std::tan(hphi) * beam_thickness;
  double    cross_trd_x2 = cross_trd_x1 + 2.* std::tan(hphi) * cross_thickness;
  double    edge_pos_x   = beam_start_x - cross_trd_x1 / 2. - beam_thickness / 2;
  Trapezoid cross_s2_trd (cross_trd_x1 / 2., cross_trd_x2 / 2.,
                          beam_thickness / 2., beam_thickness / 2., cross_thickness / 2.);
  Box       cross_s2_box ((cross_trd_x2 - cross_trd_x1)/4., beam_thickness / 2., cross_thickness / 2.);
  SubtractionSolid cross_s2(cross_s2_trd, cross_s2_box, Position((cross_trd_x2 + cross_trd_x1)/4., 0., 0.));
  Volume cross_vol2("cross_edge_beam", cross_s2, desc.material(x_support.materialStr()));
  cross_vol2.setVisAttributes(desc, x_support.visStr());
  env_vol.placeVolume(cross_vol2, Position(edge_pos_x, 0., beam_pos_z));
  env_vol.placeVolume(cross_vol2, Transform3D(Translation3D(-edge_pos_x, 0., beam_pos_z) * RotationZ(M_PI)));
  for (int j = 1; j < n_cross_supports; j++) {
    env_vol.placeVolume(cross_vol2, Position(edge_pos_x, -j * grid_size, beam_pos_z));
    env_vol.placeVolume(cross_vol2, Position(edge_pos_x, j * grid_size, beam_pos_z));
    env_vol.placeVolume(cross_vol2, Transform3D(Translation3D(-edge_pos_x, -j * grid_size, beam_pos_z) * RotationZ(M_PI)));
    env_vol.placeVolume(cross_vol2, Transform3D(Translation3D(-edge_pos_x, j * grid_size, beam_pos_z) * RotationZ(M_PI)));
  }

  mod_vol.placeVolume(env_vol, Position(0.0, 0.0, l_pos_z + support_thickness/2.));
}

// Fill fiber lattice into trapezoid starting from position (0,0) in x-z coordinate system
vector<vector<Point>> fiberPositions(double radius, double x_spacing, double z_spacing, double x, double z, double phi,
                                     double spacing_tol)
{
  // z_spacing - distance between fiber layers in z
  // x_spacing - distance between fiber centers in x
  // x - half-length of the shorter (bottom) base of the trapezoid
  // z - height of the trapezoid
  // phi - angle between z and trapezoid arm

  vector<vector<Point>> positions;
  int z_layers = floor((z / 2 - radius - spacing_tol) / z_spacing); // number of layers that fit in z/2

  double z_pos = 0.;
  double x_pos = 0.;

  for (int l = -z_layers; l < z_layers + 1; l++) {
    vector<Point> xline;
    z_pos                = l * z_spacing;
    double x_max         = x + (z / 2. + z_pos) * tan(phi) - spacing_tol; // calculate max x at particular z_pos
    (l % 2 == 0) ? x_pos = 0. : x_pos = x_spacing / 2;                    // account for spacing/2 shift

    while (x_pos < (x_max - radius)) {
      xline.push_back(Point(x_pos, z_pos));
      if (x_pos != 0.)
        xline.push_back(Point(-x_pos, z_pos)); // using symmetry around x=0
      x_pos += x_spacing;
    }
    // Sort fiber IDs for a better organization
    sort(xline.begin(), xline.end(), [](const Point& p1, const Point& p2) { return p1.x() < p2.x(); });
    positions.emplace_back(std::move(xline));
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
vector<tuple<int, Point, Point, Point, Point>> gridPoints(int div_x, int div_z, double x, double z, double phi)
{
  // x, z and phi defined as in vector<Point> fiberPositions
  // div_x, div_z - number of divisions in x and z
  double dz = z / div_z;

  std::vector<std::tuple<int, Point, Point, Point, Point>> points;

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

      int id = ix + div_x * iz;

      auto A = Point(A_x, A_z);
      auto B = Point(B_x, B_z);
      auto C = Point(C_x, C_z);
      auto D = Point(D_x, D_z);

      // vertex points filled in the clock-wise direction
      points.push_back(make_tuple(id, A, B, C, D));
    }
  }

  return points;
}

DECLARE_DETELEMENT(athena_EcalBarrelInterlayers, create_detector)

