// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Dmitry Kalinkin, Wouter Deconinck

#include <array>
#include <limits>
#include <string>

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Printout.h>

using std::string;
using namespace dd4hep;

static Ref_t create_detector(Detector& lcdd, xml_h handle, SensitiveDetector sens)
{
  xml_det_t  det_handle      = handle;
  xml_dim_t  envelope_handle = det_handle.dimensions();
  xml_dim_t  sectors_handle  = det_handle.child(_Unicode(sectors));
  xml_dim_t  rows_handle     = sectors_handle.child(_Unicode(rows));
  xml_dim_t  dim_handle      = rows_handle.dimensions();
  Material   tower_mat       = lcdd.material(det_handle.materialStr());
  double     row_rmin        = dim_handle.inner_r();
  DetElement det_element{det_handle.nameStr(), det_handle.id()};

  // envelope
  std::vector<double> v_rmin, v_rmax, v_z;

  {
    auto dim    = det_handle.dimensions();
    auto rmin   = dim.rmin();
    auto rmax   = dim.rmax();
    auto zmin   = dim.zmin();
    auto zmax   = dim.zmax();
    auto etamin = dd4hep::getAttrOrDefault<double>(det_handle.child(_U(dimensions)), _Unicode(etamin),
                                                   -std::numeric_limits<double>::max());
    auto etamax = dd4hep::getAttrOrDefault<double>(det_handle.child(_U(dimensions)), _Unicode(etamax),
                                                   +std::numeric_limits<double>::max());
    auto theta  = [](const auto eta) { return 2.0 * atan(exp(-eta)); };

    // backward nose cone
    printout(DEBUG, "SciGlassCalorimeter", "etamin cutout: etamin = %f, thetamin = %f", etamin, theta(etamin));
    if (zmin * tan(theta(etamin)) < rmin) {
      // no cutout: regular end face
      printout(DEBUG, "SciGlassCalorimeter", "no etamin cutout: zmin * tan(theta(etamin)) = %f < rmin = %f",
               zmin * tan(theta(etamin)), rmin);
      v_z.emplace_back(zmin);
      v_rmin.emplace_back(rmin);
      v_rmax.emplace_back(rmax);
    } else {
      // with cutout: first add zmin side
      printout(DEBUG, "SciGlassCalorimeter", "etamin cutout: zmin * tan(theta(etamin)) = %f > rmin = %f",
               zmin * tan(theta(etamin)), rmin);
      auto z = std::max(zmin, rmax / tan(theta(etamin))); // zmin or furthest backwards
      v_z.emplace_back(z);
      v_rmin.emplace_back(std::max(rmin, z * tan(theta(etamin))));
      v_rmax.emplace_back(rmax);
      // then where cutout starts
      v_z.emplace_back(rmin / tan(theta(etamin)));
      v_rmin.emplace_back(rmin);
      v_rmax.emplace_back(rmax);
    }

    // forward nose cone
    printout(DEBUG, "SciGlassCalorimeter", "etamax cutout: etamax = %f, thetamax = %f", etamax, theta(etamax));
    if (zmax * tan(theta(etamax)) < rmin) {
      // no cutout: regular end face
      printout(DEBUG, "SciGlassCalorimeter", "no etamax cutout: zmax * tan(theta(etamin)) = %f < rmin = %f",
               zmax * tan(theta(etamax)), rmin);
      v_z.emplace_back(zmax);
      v_rmin.emplace_back(rmin);
      v_rmax.emplace_back(rmax);
    } else {
      // with cutout: first add where cutout starts
      printout(DEBUG, "SciGlassCalorimeter", "etamax cutout: zmax * tan(theta(etamax)) = %f > rmin = %f",
               zmax * tan(theta(etamax)), rmin);
      v_z.emplace_back(rmin / tan(theta(etamax)));
      v_rmin.emplace_back(rmin);
      v_rmax.emplace_back(rmax);
      // then add zmax side
      auto z = std::min(zmax, rmax / tan(theta(etamax))); // zmax or furthest forward
      v_z.emplace_back(z);
      v_rmin.emplace_back(std::max(rmin, z * tan(theta(etamax))));
      v_rmax.emplace_back(rmax);
    }
  }

  // create polycone
  Polycone envelope_shape(0.0, 2 * M_PI, v_rmin, v_rmax, v_z);
  Volume   envelope_v{det_handle.nameStr(), envelope_shape, lcdd.material("Air")};

  PlacedVolume envelope_pv = lcdd.pickMotherVolume(det_element).placeVolume(envelope_v);
  envelope_pv.addPhysVolID("system", det_handle.id());
  det_element.setPlacement(envelope_pv);

  sens.setType("calorimeter");

  double support_inner_r = std::numeric_limits<double>::max();

  if (det_handle.hasChild(_Unicode(outer_supports))) {
    xml_comp_t outer_supports = det_handle.child(_Unicode(outer_supports));

    for (xml_coll_t layer_handle{outer_supports, _Unicode(layer)}; layer_handle; ++layer_handle) {
      xml_comp_t outer_support_handle = layer_handle;
      Material   outer_support_mat    = lcdd.material(outer_support_handle.materialStr());

      Tube              outer_support_tube_full_shape{outer_support_handle.inner_r(),
                                         outer_support_handle.inner_r() + outer_support_handle.thickness(),
                                         (envelope_handle.zmax() - envelope_handle.zmin()) / 2};
      IntersectionSolid outer_support_tube_shape{
          outer_support_tube_full_shape,
          envelope_shape,
      };

      support_inner_r = std::min(support_inner_r, outer_support_handle.inner_r());

      Volume outer_support_tube_v{"outer_support", outer_support_tube_shape, outer_support_mat};
      outer_support_tube_v.setVisAttributes(lcdd.visAttributes(outer_support_handle.visStr()));
      envelope_v.placeVolume(outer_support_tube_v);
    }
  }

  double sector_rmin = envelope_handle.rmin();

  if (det_handle.hasChild(_Unicode(wedge_box))) {
    xml_comp_t wedge_box_handle = det_handle.child(_Unicode(wedge_box));
    Material   wedge_box_mat    = lcdd.material(wedge_box_handle.materialStr());

    // Approximate bottoms of wedge boxes with a single tube common for all sectors
    // The actual wedge bottom may have to be a plane, not a tube section

    Tube              wedge_box_tube_full_shape{wedge_box_handle.inner_r(),
                                   wedge_box_handle.inner_r() + wedge_box_handle.thickness(),
                                   (envelope_handle.zmax() - envelope_handle.zmin()) / 2};
    IntersectionSolid wedge_box_tube_shape{
        wedge_box_tube_full_shape,
        envelope_shape,
    };

    Volume wedge_box_tube_v{"wedge_box_placeholder", wedge_box_tube_shape, wedge_box_mat};
    wedge_box_tube_v.setVisAttributes(lcdd.visAttributes(wedge_box_handle.visStr()));
    envelope_v.placeVolume(wedge_box_tube_v);

    const double side_rmin         = wedge_box_handle.inner_r() + wedge_box_handle.thickness();
    const double side_rmax_default = std::min(envelope_handle.rmax(),
                                              // subtract a unit of thickness to avoid an overlap
                                              support_inner_r - wedge_box_handle.thickness());
    const double side_rmax         = dd4hep::getAttrOrDefault<double>(wedge_box_handle, _U(outer_r), side_rmax_default);

    sector_rmin = side_rmin;

    Box wedge_box_side_box_shape{
        (side_rmax - side_rmin) / 2,
        wedge_box_handle.thickness() / 2, (envelope_handle.zmax() - envelope_handle.zmin()) / 2};
    Volume wedge_box_side_v[2];
    for (int side = -1; side <= 1; side += 2) {
      IntersectionSolid wedge_box_side_shape{envelope_shape, wedge_box_side_box_shape,
                                             Position{(side_rmax + side_rmin) / 2, side * wedge_box_handle.gap() / 2,
                                                      (envelope_handle.zmax() + envelope_handle.zmin()) / 2}};

      wedge_box_side_v[(side + 1) / 2] = Volume({"wedge_box_side", wedge_box_side_shape, wedge_box_mat});
      wedge_box_side_v[(side + 1) / 2].setVisAttributes(lcdd.visAttributes(wedge_box_handle.visStr()));
    }

    int    sector     = 0;
    double sector_phi = sectors_handle.phi0();
    for (; sector < sectors_handle.number(); sector++, sector_phi += sectors_handle.deltaphi()) {
      for (int side = -1; side <= 1; side += 2) {
        envelope_v.placeVolume(wedge_box_side_v[(side + 1) / 2],
                               Transform3D{RotationZ{sector_phi + side * sectors_handle.deltaphi() / 2}});
      }
    }

    // TODO: The endcap sides of the box are not implemented
  }

  int    row     = 0;
  double row_phi = -rows_handle.deltaphi() / 2 * (rows_handle.number() - 1);

  Tube              sector_tube_shape{sector_rmin, std::min(support_inner_r, envelope_handle.rmax()),
                         (envelope_handle.zmax() - envelope_handle.zmin()) / 2,
                         -rows_handle.deltaphi() / 2 * rows_handle.number(),
                         rows_handle.deltaphi() / 2 * rows_handle.number()};
  IntersectionSolid sector_shape{envelope_shape, sector_tube_shape,
                                 Position{0., 0., (envelope_handle.zmax() + envelope_handle.zmin()) / 2}};

  Volume sector_v{"sector", sector_shape, lcdd.material("Air")};
  sector_v.setVisAttributes(lcdd.visAttributes(det_handle.visStr()));

  for (; row < rows_handle.number(); row++, row_phi += rows_handle.deltaphi()) {

    const double tower_gap_longitudinal = dim_handle.gap();

    // negative rapidity towers will be counted backwards from -1
    std::array<int, 2>    tower_ids               = {-1, 0};
    std::array<double, 2> betas                   = {0., 0.};
    std::array<double, 2> beta_prevs              = {0., 0.};
    std::array<double, 2> dzs                     = {0., 0.};
    std::array<double, 2> flare_angle_polar_prevs = {0., 0.};

    for (xml_coll_t family_handle{rows_handle, _Unicode(family)}; family_handle; ++family_handle) {
      const int dir_sign               = family_handle.attr<double>(_Unicode(dir_sign));
      int&      tower_id               = tower_ids[(dir_sign > 0) ? 1 : 0];
      double&   beta                   = betas[(dir_sign > 0) ? 1 : 0];
      double&   beta_prev              = beta_prevs[(dir_sign > 0) ? 1 : 0];
      double&   dz                     = dzs[(dir_sign > 0) ? 1 : 0];
      double&   flare_angle_polar_prev = flare_angle_polar_prevs[(dir_sign > 0) ? 1 : 0];

      xml_dim_t          family_dim_handle   = family_handle;
      const double       length              = family_dim_handle.z_length();
      const auto         flare_angle_polar   = family_dim_handle.attr<double>(_Unicode(flare_angle_polar));
      const unsigned int number              = family_dim_handle.number();
      const auto         flare_angle_at_face = family_dim_handle.attr<double>(_Unicode(flare_angle_at_face));

      const double z  = length / 2;
      // Face parameters (see doc/sciglass_tower_front_view.svg for definitions)
      const double y1 = family_dim_handle.y1();
      const double y2 = y1 + length * tan(flare_angle_polar);
      double       x1 = family_dim_handle.x1();
      double       x2 = family_dim_handle.x1() + (2 * y1) * tan(flare_angle_at_face);
      double       x3, x4;
      if (family_dim_handle.hasAttr(_Unicode(flare_angle_azimuthal))) {
        // Azimuthal flaring independently defined
        const auto flare_angle_azimuthal = family_dim_handle.attr<double>(_Unicode(flare_angle_azimuthal));
        x3 = x1 + length * (tan(flare_angle_azimuthal) - tan(flare_angle_polar) * tan(flare_angle_at_face));
        x4 = x2 + length * (tan(flare_angle_azimuthal) + tan(flare_angle_polar) * tan(flare_angle_at_face));
      } else {
        // Pyramidal shape
        // corresponds to: tan(flare_angle_azimuthal) = (x1 / y1 + tan(flare_angle_at_face)) * tan(flare_angle_polar)
        x3 = x1 * (y2 / y1);
        x4 = x2 * (y2 / y1);
      }

      if (dir_sign < 0) {
        // Mirror towers going towards negative rapidities
        std::swap(x1, x2);
        std::swap(x3, x4);
      }
      const double theta  = 0.;
      const double phi    = 0.;
      const double alpha1 = 0.;
      const double alpha2 = 0.;

      for (unsigned int tower = 0; tower < number-1; tower++, tower_id += dir_sign) {
        // see https://github.com/eic/epic/blob/main/doc/sciglass_tower_stacking.svg
        beta += flare_angle_polar_prev + flare_angle_polar;
        const double gamma = M_PI_2 - flare_angle_polar_prev - beta_prev;

        const double dz_prev = dz;
        dz += (tower_gap_longitudinal / cos(flare_angle_polar) + 2 * y1) * sin(M_PI - gamma - beta) / sin(gamma);
        const string t_name = _toString(row, "_row%d") + _toString(tower_id, "_tower%d");
        sector_v
	  .placeVolume(
		       Volume{t_name, Trap{z, theta, phi, y1, x1, x2, alpha1, y2, x3, x4, alpha2}, tower_mat},
		       Transform3D{RotationZ{-M_PI_2 + row_phi}} * Transform3D{Position{0. * cm, row_rmin, dir_sign * dz}} *
		       Transform3D{RotationX{-M_PI / 2 + dir_sign * beta}} * Transform3D{Position{0, dir_sign * y1, z}})
	  .addPhysVolID("row", row)
	  .addPhysVolID("tower", tower_id)
	  .volume()
	  .setSensitiveDetector(sens)
	  .setVisAttributes(lcdd.visAttributes(family_dim_handle.visStr()));
	
        if (sectors_handle.hasChild(_Unicode(carbon_fiber_support))) {
          xml_comp_t carbon_fiber_support_handle = sectors_handle.child(_Unicode(carbon_fiber_support));
          xml_comp_t cut_out_handle              = carbon_fiber_support_handle.child(_Unicode(cut_out));
          Material   carbon_fiber_support_mat    = lcdd.material(carbon_fiber_support_handle.materialStr());

          const double margin_horizontal = cut_out_handle.attr<double>(_Unicode(margin_horizontal));
          const double margin_top        = cut_out_handle.attr<double>(_Unicode(margin_top));
          const double margin_bottom     = cut_out_handle.attr<double>(_Unicode(margin_bottom));
          const double overhang_top      = carbon_fiber_support_handle.attr<double>(_Unicode(overhang_top));
          const double overhang_bottom   = carbon_fiber_support_handle.attr<double>(_Unicode(overhang_bottom));

          const double y1_long =
              y1 + tower_gap_longitudinal / 2 / cos(flare_angle_polar) - overhang_bottom * tan(flare_angle_polar);
          const double y2_long =
              y2 + tower_gap_longitudinal / 2 / cos(flare_angle_polar) + overhang_top * tan(flare_angle_polar);

          Trap trap_long_1{z + (overhang_top + overhang_bottom) / 2,
                           theta,
                           phi,
                           y1_long,
                           carbon_fiber_support_handle.thickness() / 4,
                           carbon_fiber_support_handle.thickness() / 4,
                           alpha1,
                           y2_long,
                           carbon_fiber_support_handle.thickness() / 4,
                           carbon_fiber_support_handle.thickness() / 4,
                           alpha2};
          Trap trap_long_2{
              z - (margin_top + margin_bottom) / 2,
              theta,
              phi,
              y1_long - margin_horizontal /
                            2,  // FIXME: y1_long/y2_long do not account for reduction of z by the vertical margins
              carbon_fiber_support_handle
                  .thickness(), // no division by 2 to ensure subtrahend volume is thicker than the minuend
              carbon_fiber_support_handle.thickness(),
              alpha1,
              y2_long - margin_horizontal / 2,
              carbon_fiber_support_handle.thickness(),
              carbon_fiber_support_handle.thickness(),
              alpha2};
          SubtractionSolid trap_long{trap_long_1, trap_long_2, Position{0., 0., (margin_bottom - margin_top) / 2}};

          for (int side = -1; side <= 1; side += 2) {
            sector_v
                .placeVolume(Volume{"fiber_structure_longitudinal", trap_long, carbon_fiber_support_mat},
                             Transform3D{RotationZ{-M_PI_2 + row_phi + side * rows_handle.deltaphi() / 2}} *
                                 Transform3D{Position{0. * cm, row_rmin, dir_sign * dz}} *
                                 Transform3D{RotationX{-M_PI / 2 + dir_sign * beta}} *
                                 Transform3D{Position{side * carbon_fiber_support_handle.thickness() / 4, dir_sign * y1,
                                                      z + (overhang_top - overhang_bottom) / 2}})
                .volume()
                .setVisAttributes(lcdd.visAttributes(carbon_fiber_support_handle.visStr()));
          }

          const double dz_trans = dz_prev + (tower_gap_longitudinal / 2 + carbon_fiber_support_handle.thickness() / 2) /
                                                cos(flare_angle_polar) * sin(M_PI - gamma - beta) / sin(gamma);

          const double non_overlap_long =
              sin(beta) * (tower_gap_longitudinal / cos(flare_angle_polar) + 2 * y1) / sin(gamma);

          const double beta_trans = beta_prev + flare_angle_polar_prev;
          const double x1_trans   = tan(rows_handle.deltaphi() / 2) * (row_rmin - overhang_bottom * cos(beta_trans)) -
                                  carbon_fiber_support_handle.thickness() / 2 / cos(rows_handle.deltaphi() / 2);
          const double x2_trans = tan(rows_handle.deltaphi() / 2) *
                                      (row_rmin + (2 * z + overhang_top + non_overlap_long) * cos(beta_trans)) -
                                  carbon_fiber_support_handle.thickness() / 2 / cos(rows_handle.deltaphi() / 2);

          Trap trap_trans_1{z + (overhang_top + non_overlap_long + overhang_bottom) / 2,
                            theta,
                            phi,
                            carbon_fiber_support_handle.thickness() / 2,
                            x1_trans,
                            x1_trans,
                            alpha1,
                            carbon_fiber_support_handle.thickness() / 2,
                            x2_trans,
                            x2_trans,
                            alpha2};
          Trap trap_trans_2{
              z - (margin_top + non_overlap_long + margin_bottom) / 2,
              theta,
              phi,
              carbon_fiber_support_handle
                  .thickness(), // no division by 2 to ensure subtrahend volume is thicker than the minuend
              x1_trans - margin_horizontal /
                             2, // FIXME: x1_trans/x2_trans do not account for reduction of z by the vertical margins
              x1_trans - margin_horizontal / 2,
              alpha1,
              carbon_fiber_support_handle.thickness(),
              x2_trans - margin_horizontal / 2,
              x2_trans - margin_horizontal / 2,
              alpha2};
          SubtractionSolid trap_trans{
              trap_trans_1, trap_trans_2,
              Position{0., 0., (margin_bottom - margin_top) / 2 + (overhang_bottom - overhang_top) / 2}};

          sector_v
              .placeVolume(Volume{"fiber_structure_transverse" + t_name, trap_trans, carbon_fiber_support_mat},
                           Transform3D{RotationZ{-M_PI_2 + row_phi}} *
                               Transform3D{Position{0. * cm, row_rmin, dir_sign * dz_trans}} *
                               Transform3D{RotationX{-M_PI / 2 + dir_sign * beta_trans}} *
                               Transform3D{Position{0, dir_sign * carbon_fiber_support_handle.thickness() / 2,
                                                    z + (overhang_top + non_overlap_long - overhang_bottom) / 2}})
              .volume()
              .setVisAttributes(lcdd.visAttributes(carbon_fiber_support_handle.visStr()));
        }

        beta_prev              = beta;
        flare_angle_polar_prev = flare_angle_polar;
      }
    }
  }

  int    sector     = 0;
  double sector_phi = sectors_handle.phi0();
  for (; sector < sectors_handle.number(); sector++, sector_phi += sectors_handle.deltaphi()) {
    envelope_v.placeVolume(sector_v, Transform3D{RotationZ{sector_phi}}).addPhysVolID("sector", sector);
  }

  envelope_v.setVisAttributes(lcdd.visAttributes(det_handle.visStr()));

  return det_element;
}

DECLARE_DETELEMENT(epic_EcalBarrelSciGlass, create_detector)
