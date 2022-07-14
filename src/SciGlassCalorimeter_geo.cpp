#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/OpticalSurfaces.h>
#include <DD4hep/Printout.h>
#include <Math/AxisAngle.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <math.h>
#include <tuple>

#include "GeometryHelpers.h"

using namespace dd4hep;

// helper function to get x, y, z if defined in a xml component
template <class XmlComp>
Position get_xml_xyz(const XmlComp& comp, dd4hep::xml::Strng_t name)
{
  Position pos(0., 0., 0.);
  if (comp.hasChild(name)) {
    auto child = comp.child(name);
    pos.SetX(dd4hep::getAttrOrDefault<double>(child, _Unicode(x), 0.));
    pos.SetY(dd4hep::getAttrOrDefault<double>(child, _Unicode(y), 0.));
    pos.SetZ(dd4hep::getAttrOrDefault<double>(child, _Unicode(z), 0.));
  }
  return pos;
}

static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml_det_t  x_det    = handle;
  auto       det_name = x_det.nameStr();
  DetElement det(det_name, x_det.id());
  sens.setType("calorimeter");

  // detector position and rotation
  auto det_pos = get_xml_xyz(x_det, _U(position));
  auto det_rot = get_xml_xyz(x_det, _U(rotation));

  // envelop
  auto                dim    = x_det.dimensions();
  auto                rmin   = dim.rmin();
  auto                rmax   = dim.rmax();
  auto                zmin   = dim.zmin();
  auto                zmax   = dim.zmax();
  auto                etamin = dd4hep::getAttrOrDefault<double>(x_det.child(_U(dimensions)), _Unicode(etamin),
                                                 -std::numeric_limits<double>::max());
  auto                etamax = dd4hep::getAttrOrDefault<double>(x_det.child(_U(dimensions)), _Unicode(etamax),
                                                 +std::numeric_limits<double>::max());
  std::vector<double> v_rmin, v_rmax, v_z;
  auto                theta = [](const auto eta) { return 2.0 * atan(exp(-eta)); };

  // backward nose cone
  printout(DEBUG, "SciGlassCalorimeter", "etamin cutout: etamin = %f, thetamin = %f", etamin, theta(etamin));
  if (-zmin * sin(theta(etamin)) < rmin) {
    // no cutout: regular end face
    printout(DEBUG, "SciGlassCalorimeter", "no etamin cutout: zmin * sin(theta(etamin)) = %f < rmin = %f",
             zmin * sin(theta(etamin)), rmin);
    v_z.emplace_back(zmin);
    v_rmin.emplace_back(rmin);
    v_rmax.emplace_back(rmax);
  } else {
    // with cutout: first add zmin side
    printout(DEBUG, "SciGlassCalorimeter", "etamin cutout: zmin * sin(theta(etamin)) = %f > rmin = %f",
             zmin * sin(theta(etamin)), rmin);
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
  if (zmax * sin(theta(etamax)) < rmin) {
    // no cutout: regular end face
    printout(DEBUG, "SciGlassCalorimeter", "no etamax cutout: zmax * sin(theta(etamin)) = %f < rmin = %f",
             zmax * sin(theta(etamax)), rmin);
    v_z.emplace_back(zmax);
    v_rmin.emplace_back(rmin);
    v_rmax.emplace_back(rmax);
  } else {
    // with cutout: first add where cutout starts
    printout(DEBUG, "SciGlassCalorimeter", "etamax cutout: zmax * sin(theta(etamax)) = %f > rmin = %f",
             zmax * sin(theta(etamax)), rmin);
    v_z.emplace_back(rmin / tan(theta(etamax)));
    v_rmin.emplace_back(rmin);
    v_rmax.emplace_back(rmax);
    // then add zmax side
    auto z = std::min(zmax, rmax / tan(theta(etamax))); // zmax or furthest forward
    v_z.emplace_back(z);
    v_rmin.emplace_back(std::max(rmin, z * tan(theta(etamax))));
    v_rmax.emplace_back(rmax);
  }

  // create polycone
  Polycone envShape(0.0, 2 * M_PI, v_rmin, v_rmax, v_z);
  Volume   env(det_name + "_envelope", envShape, desc.material("Air"));
  env.setVisAttributes(desc.visAttributes(x_det.visStr()));

  // number of modules (default to 128 in phi, and 1 in eta so one ring is placed)
  xml_comp_t x_mod     = x_det.child(_U(module));
  auto       n_phi     = getAttrOrDefault<size_t>(x_mod, _Unicode(repeat_phi), 128);
  auto       n_eta_pos = getAttrOrDefault<size_t>(x_mod, _Unicode(repeat_eta_pos), 1);
  auto       n_eta_neg = getAttrOrDefault<size_t>(x_mod, _Unicode(repeat_eta_neg), 1);

  // angle to exscribed chord
  auto xcrd = [](const float& phi) { return 2.0 * tan(phi / 2.0); }; // unsafe for phi = pi, no exscribed chord

  // angular size in phi
  auto dphi      = 2.0 * M_PI / n_phi;
  auto xcrd_dphi = rmin * xcrd(dphi); // size of module around phi around inscribed circle
  printout(DEBUG, "SciGlassCalorimeter", "rmin = %f, dphi = %f, xcrd_dphi = %f", rmin, dphi, xcrd_dphi);

  // module parameters
  auto mod_name                    = x_mod.nameStr();
  auto mod_length                  = getAttrOrDefault<float>(x_mod, _Unicode(length), 0.0);
  auto mod_round_backface_to       = getAttrOrDefault<float>(x_mod, _Unicode(round_backface_to), 0.0);
  auto mod_phi_projectivity_tilt   = getAttrOrDefault<float>(x_mod, _Unicode(phi_projectivity_tilt), 0.0);
  auto mod_eta_projectivity_offset = getAttrOrDefault<float>(x_mod, _Unicode(eta_projectivity_offset), 0.0);
  printout(DEBUG, "SciGlassCalorimeter", "Proj: phi = %f, eta = %f", mod_phi_projectivity_tilt,
           mod_eta_projectivity_offset);

  // solver (secant's method)
  auto solve_secant = [](const auto& f, auto x0, auto x1, const float eps = 1e-4f, const size_t n = 10u) {
    size_t i  = 0;
    auto   x2 = x1;
    while (fabs(x1 - x0) > eps && ++i < n) {
      x2 -= f(x1) * (x1 - x0) / (f(x1) - f(x0));
      x0 = x1;
      x1 = x2;
    }
    return std::pair{x2, i != n};
  };

  // start at center and move outwards in eta
  auto theta_min = 0.0;
  auto theta_max = 0.0;
  for (unsigned int k_eta = 0; k_eta < std::max(n_eta_neg, n_eta_pos); k_eta++) {
    // previous theta_max is current theta_min
    theta_min = theta_max;

    // inner face (default to touching blocks with square face)
    const auto mod_x1 = getAttrOrDefault<float>(x_mod, _Unicode(x1), xcrd_dphi * cos(mod_phi_projectivity_tilt));
    const auto mod_y1 = getAttrOrDefault<float>(x_mod, _Unicode(y1), mod_x1);
    // FIXME base should not be square or leakage will be present even for phi_projectivity_tilt == 0
    // instead of Trd2 with x1 and x2, we may need a Trap with x1, x2, x3, x4

    // theta_max update
    // dtheta is the solution of sin(dtheta / 2.0) / cos(theta_min + dtheta) = mod_y1 / 2.0 / rmin
    const auto dtheta0 = mod_y1 / 2.0 / rmin;
    const auto dtheta1 = 2.0 * asin(dtheta0);
    const auto f       = [&theta_min, &dtheta0](const auto& dtheta) {
      return sin(dtheta / 2.0) / cos(theta_min + dtheta) - dtheta0;
    };
    const auto [dtheta, valid] = solve_secant(f, dtheta0, dtheta1);
    if (!valid) {
      printout(WARNING, "SciGlassCalorimeter", "cannot solve for dtheta");
    }
    theta_max = theta_min + dtheta;

    // outer face (default to radial extension to square face)
    auto mod_expansion = cos(theta_max) / rmin;
    auto mod_x2        = getAttrOrDefault<float>(x_mod, _Unicode(x2), mod_x1 * (1 + mod_expansion * mod_length));
    auto mod_y2        = getAttrOrDefault<float>(x_mod, _Unicode(y2), mod_y1 * (1 + mod_expansion * mod_length));
    // round down to nearest multiple for limited block families
    if (mod_round_backface_to > 0.0) {
      mod_x2        = std::floor(mod_x2 / mod_round_backface_to) * mod_round_backface_to;
      mod_y2        = std::floor(mod_y2 / mod_round_backface_to) * mod_round_backface_to;
      mod_expansion = (mod_x2 / mod_x1 - 1.0) / mod_length;
    }
    printout(DEBUG, "SciGlassCalorimeter", "Trd2: x1 = %f, x2 = %f, y1 = %f, y2 = %f", mod_x1, mod_x2, mod_y1, mod_y2);

    // create module envelope
    Trd2     mod_env_trd(mod_x1 / 2.0, mod_x2 / 2.0, mod_y1 / 2.0, mod_y2 / 2.0, mod_length / 2.0);
    Material mod_env_mat(desc.material(x_mod.materialStr()));
    Volume   mod_env(mod_name + _toString((signed)k_eta, "_sector%d"), mod_env_trd, mod_env_mat);
    mod_env.setVisAttributes(desc.visAttributes(x_mod.visStr()));

    // place slices in module
    auto s_num   = 1u;
    auto s_rmin  = rmin;
    auto s_x1    = mod_x1;
    auto s_y1    = mod_y1;
    auto s_x2    = mod_x2;
    auto s_y2    = mod_y2;
    auto s_pos_z = -mod_length / 2.0;
    for (xml_coll_t si(x_mod, _U(slice)); si; ++si) {
      xml_comp_t x_slice     = si;
      auto       s_name      = Form("slice%d", s_num);
      const auto s_material  = desc.material(x_slice.materialStr());
      const auto s_thickness = dd4hep::getAttrOrDefault<float>(x_slice, _Unicode(thickness), 0.);
      s_x2                   = s_x1 * (1 + mod_expansion * s_thickness);
      s_y2                   = s_y1 * (1 + mod_expansion * s_thickness);
      Trd2   s_trd(s_x1 / 2.0, s_x2 / 2.0, s_y1 / 2.0, s_y2 / 2.0, s_thickness / 2.0);
      Volume s_vol(s_name, s_trd, s_material);
      s_vol.setVisAttributes(desc.visAttributes(x_slice.visStr()));
      if (x_slice.isSensitive()) {
        s_vol.setSensitiveDetector(sens);
      }
      // DetElement slice(stave_det, s_name, det_id);
      s_pos_z += s_thickness / 2.0;
      mod_env.placeVolume(s_vol, Position(0, 0, s_pos_z));
      s_pos_z += s_thickness / 2.0;
      // set starting face for next slice
      s_rmin += s_thickness / cos(dtheta / 2.0) * cos(M_PI_2 - theta_max);
      s_x1 = s_x2;
      s_y1 = s_y2;
    }

    // place around phi
    for (size_t j_phi = 0; j_phi < n_phi; j_phi++) {
      // azimuthal and polar angles
      const auto phi       = dphi * j_phi;
      const auto avg_theta = 0.5 * (theta_min + theta_max);

      // module center position
      const auto r = rmin + 0.5 * mod_x1 * sin(mod_phi_projectivity_tilt) + 0.5 * mod_y1 * sin(avg_theta) +
                     0.5 * mod_length * cos(avg_theta);
      const auto x = r * cos(phi);
      const auto y = r * sin(phi);
      const auto z = r / tan(M_PI_2 - avg_theta);

      // place negative module
      if (k_eta < n_eta_neg) {
        Transform3D tr_neg = Translation3D(x, y, -det_pos.z() + mod_eta_projectivity_offset - z) *
                             RotationZ(phi + mod_phi_projectivity_tilt) * RotationY(M_PI_2 + avg_theta);
        auto pv_neg = env.placeVolume(mod_env, tr_neg);
        pv_neg.addPhysVolID("sector", n_eta_neg + n_eta_pos - k_eta - 1);
        pv_neg.addPhysVolID("module", j_phi);
      }

      // place positive module
      if (k_eta < n_eta_pos) {
        Transform3D tr_pos = Translation3D(x, y, -det_pos.z() + mod_eta_projectivity_offset + z) *
                             RotationZ(phi + mod_phi_projectivity_tilt) * RotationY(M_PI_2 - avg_theta);
        auto pv_pos = env.placeVolume(mod_env, tr_pos);
        pv_pos.addPhysVolID("sector", k_eta);
        pv_pos.addPhysVolID("module", j_phi);
      }
    }
  }

  // placement
  Volume      motherVol = desc.pickMotherVolume(det);
  Transform3D tr =
      Translation3D(det_pos.x(), det_pos.y(), det_pos.z()) * RotationZYX(det_rot.z(), det_rot.y(), det_rot.x());
  auto pv_env = motherVol.placeVolume(env, tr);
  pv_env.addPhysVolID("system", det.id());
  det.setPlacement(pv_env);
  return det;
}

DECLARE_DETELEMENT(ecce_SciGlassCalorimeter, create_detector)
