// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong, Sylvester Joosten

//==========================================================================
//
//      <detector name ="DetName" type="Beampipe" >
//      <layer id="#(int)" inner_r="#(double)" outer_z="#(double)" >
//      <slice material="string" thickness="#(double)" >
//      </layer>
//      </detector>
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>
#include "XML/Utilities.h"
#include "DD4hepDetectorHelper.h"
#include <tuple>

using namespace std;
using namespace dd4hep;

/** \addtogroup beamline Beamline Instrumentation
 */

/** \addtogroup IRChamber Interaction Region Vacuum Chamber.
 * \brief Type: **IRChamber**.
 * \ingroup beamline
 *
 *
 * \code
 *   <detector>
 *   </detector>
 * \endcode
 *
 */
static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) {
  using namespace ROOT::Math;
  xml_det_t x_det       = e;
  string det_name       = x_det.nameStr();
  xml_comp_t x_dettype  = x_det.child(dd4hep::xml::Strng_t("type_flags"));
  unsigned int typeFlag = x_dettype.type();
  DetElement sdet(det_name, x_det.id());
  Assembly assembly(det_name + "_assembly");

  xml::Component IP_pipe_c = x_det.child(_Unicode(IP_pipe));

  Material IP_beampipe_wall_material =
      det.material(IP_pipe_c.attr<string>(_Unicode(wall_material)));
  Material IP_beampipe_coating_material =
      det.material(IP_pipe_c.attr<string>(_Unicode(coating_material)));
  Material m_Vacuum  = det.material("Vacuum");
  Material m_Wall    = det.material("StainlessSteel");
  Material m_Coating = det.material("Copper");

  // IP Beampipe
  double IP_beampipe_ID                = IP_pipe_c.attr<double>(_Unicode(ID));
  double IP_beampipe_wall_thickness    = IP_pipe_c.attr<double>(_Unicode(wall_thickness));
  double IP_beampipe_coating_thickness = IP_pipe_c.attr<double>(_Unicode(coating_thickness));
  double IP_acts_beampipe_OD           = IP_beampipe_ID - 5.0 * mm;
  double IP_acts_beampipe_ID           = IP_acts_beampipe_OD - 1.0 * mm;

  double upstream_straight_length   = IP_pipe_c.attr<double>(_Unicode(upstream_straight_length));
  double downstream_straight_length = IP_pipe_c.attr<double>(_Unicode(downstream_straight_length));

  // visualization
  auto wallVis      = det.visAttributes(x_det.attr<std::string>(_Unicode(vis_wall)));
  auto coatingVis   = det.visAttributes(x_det.attr<std::string>(_Unicode(vis_coating)));
  auto IPwallVis    = det.visAttributes(x_det.attr<std::string>(_Unicode(vis_IPwall)));
  auto IPcoatingVis = det.visAttributes(x_det.attr<std::string>(_Unicode(vis_IPcoating)));

  // colors: (r, g, b, alpha)
  wallVis.setColor(0.0, 0.0, 1.0, 1.0);      // blue
  coatingVis.setColor(1.0, 0.0, 0.0, 1.0);   // red
  IPwallVis.setColor(0.0, 1.0, 0.0, 1.0);    // green
  IPcoatingVis.setColor(1.0, 1.0, 0.0, 1.0); // yellow

  // central acts beampipe volume
  Tube central_tube(0.5 * IP_acts_beampipe_ID, 0.5 * IP_acts_beampipe_OD,
                    0.5 * (upstream_straight_length + downstream_straight_length));
  Volume central_volume("acts_central_beampipe_vol", central_tube, m_Vacuum);
  const double central_offset = -.5 * (upstream_straight_length - downstream_straight_length);
  DetElement central_det(sdet, "acts_beampipe_central", 1);

  // Set dd4hep variant parameters for conversion to ACTS tracking geometry
  central_det.setTypeFlag(typeFlag);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(central_det);
  int nBinPhi  = 144; // fix later. Should take this from a xml tag
  int nBinZ    = 10;  // fix later. Should take this from a xml tag
  params.set<bool>("layer_material", true);
  params.set<bool>("layer_material_representing", true);
  params.set<int>("layer_material_representing_binPhi", nBinPhi);
  params.set<int>("layer_material_representing_binZ", nBinZ);

  // -----------------------------
  // IP beampipe
  //
  // setup for the central IP beampipe:
  //
  // /-------\ Be wall
  //  /-----\  Au coating
  //   /---\   Vacuum padding (5mm)
  //    /-\    Fake vacuum beampipe (1mm)
  //     -     Vacuum filled inner beampipe
  //
  Tube downstream_IP_tube(IP_beampipe_ID / 2.0, IP_beampipe_ID / 2.0 + IP_beampipe_wall_thickness,
                          downstream_straight_length / 2.0);
  Tube downstream_IP_coating(IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness,
                             IP_beampipe_ID / 2.0, downstream_straight_length / 2.0);
  Tube downstream_IP_vacuum_padding(IP_acts_beampipe_OD / 2.0,
                                    IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness,
                                    downstream_straight_length / 2.0);
  Tube downstream_IP_acts_beampipe(IP_acts_beampipe_ID / 2.0, IP_acts_beampipe_OD / 2.0,
                                   downstream_straight_length / 2.0);
  Tube downstream_IP_vacuum_fill(0.0, IP_acts_beampipe_ID / 2.0, downstream_straight_length / 2.0);

  Tube upstream_IP_tube(IP_beampipe_ID / 2.0, IP_beampipe_ID / 2.0 + IP_beampipe_wall_thickness,
                        upstream_straight_length / 2.0);
  Tube upstream_IP_coating(IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness,
                           IP_beampipe_ID / 2.0, upstream_straight_length / 2.0);
  Tube upstream_IP_vacuum_padding(IP_acts_beampipe_OD / 2.0,
                                  IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness,
                                  upstream_straight_length / 2.0);
  Tube upstream_IP_acts_beampipe(IP_acts_beampipe_ID / 2.0, IP_acts_beampipe_OD / 2.0,
                                 upstream_straight_length / 2.0);
  Tube upstream_IP_vacuum_fill(0.0, IP_acts_beampipe_ID / 2.0, upstream_straight_length / 2.0);

  // create volumes
  Volume v_downstream_IP_vacuum_fill("v_downstream_IP_vacuum_fill", downstream_IP_vacuum_fill,
                                     m_Vacuum);
  Volume v_downstream_IP_acts_beampipe("v_downstream_IP_acts_beampipe", downstream_IP_acts_beampipe,
                                       m_Vacuum);
  Volume v_downstream_IP_vacuum_padding("v_downstream_IP_vacuum_padding",
                                        downstream_IP_vacuum_padding, m_Vacuum);
  Volume v_downstream_IP_coating("v_downstream_IP_coating", downstream_IP_coating,
                                 IP_beampipe_coating_material);
  Volume v_downstream_IP_tube("v_downstream_IP_tube", downstream_IP_tube,
                              IP_beampipe_wall_material);
  Volume v_upstream_IP_vacuum_fill("v_upstream_IP_vacuum_fill", upstream_IP_vacuum_fill, m_Vacuum);
  Volume v_upstream_IP_acts_beampipe("v_upstream_IP_acts_beampipe", upstream_IP_acts_beampipe,
                                     m_Vacuum);
  Volume v_upstream_IP_vacuum_padding("v_upstream_IP_vacuum_padding", upstream_IP_vacuum_padding,
                                      m_Vacuum);
  Volume v_upstream_IP_coating("v_upstream_IP_coating", upstream_IP_coating,
                               IP_beampipe_coating_material);
  Volume v_upstream_IP_tube("v_upstream_IP_tube", upstream_IP_tube, IP_beampipe_wall_material);

  // set names
  sdet.setAttributes(det, v_upstream_IP_coating, x_det.regionStr(), x_det.limitsStr(),
                     x_det.attr<std::string>(_Unicode(vis_IPcoating)));
  sdet.setAttributes(det, v_upstream_IP_tube, x_det.regionStr(), x_det.limitsStr(),
                     x_det.attr<std::string>(_Unicode(vis_IPwall)));
  sdet.setAttributes(det, v_downstream_IP_coating, x_det.regionStr(), x_det.limitsStr(),
                     x_det.attr<std::string>(_Unicode(vis_IPcoating)));
  sdet.setAttributes(det, v_downstream_IP_tube, x_det.regionStr(), x_det.limitsStr(),
                     x_det.attr<std::string>(_Unicode(vis_IPwall)));

  // place volumes
  assembly.placeVolume(v_upstream_IP_vacuum_fill, Position(0, 0, -upstream_straight_length / 2.0));
  central_volume.placeVolume(v_upstream_IP_acts_beampipe,
                             Position(0, 0, -upstream_straight_length / 2.0 - central_offset));
  assembly.placeVolume(v_upstream_IP_vacuum_padding,
                       Position(0, 0, -upstream_straight_length / 2.0));
  assembly.placeVolume(v_upstream_IP_coating, Position(0, 0, -upstream_straight_length / 2.0));
  assembly.placeVolume(v_upstream_IP_tube, Position(0, 0, -upstream_straight_length / 2.0));

  assembly.placeVolume(v_downstream_IP_vacuum_fill,
                       Position(0, 0, downstream_straight_length / 2.0));
  central_volume.placeVolume(v_downstream_IP_acts_beampipe,
                             Position(0, 0, downstream_straight_length / 2.0 - central_offset));
  assembly.placeVolume(v_downstream_IP_vacuum_padding,
                       Position(0, 0, downstream_straight_length / 2.0));
  assembly.placeVolume(v_downstream_IP_coating, Position(0, 0, downstream_straight_length / 2.0));
  assembly.placeVolume(v_downstream_IP_tube, Position(0, 0, downstream_straight_length / 2.0));

  auto central_pv = assembly.placeVolume(central_volume, Position(0, 0, +central_offset));
  central_det.setPlacement(central_pv);

  //---------------------------------------------------------------------------------------------------------
  // Helper function to create polycone pairs (wall, coating, and vacuum)
  auto zplane_to_polycones = [](xml::Component& x_pipe) {
    std::vector<double> zero, z;
    std::vector<double> rmax_wall, rmax_coating, rmax_vacuum;
    std::vector<double> rmin_wall, rmin_coating, rmin_vacuum;
    // thickness
    auto wall_thickness    = getAttrOrDefault(x_pipe, _Unicode(wall_thickness), 1 * mm);
    auto coating_thickness = getAttrOrDefault(x_pipe, _Unicode(coating_thickness), 30 * um);

    for (xml_coll_t x_zplane_i(x_pipe, _Unicode(zplane)); x_zplane_i; ++x_zplane_i) {
      xml_comp_t x_zplane = x_zplane_i;
      // z position
      z.push_back(x_zplane.attr<double>(_Unicode(z)));
      // outer radius
      rmax_wall.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 + wall_thickness);
      rmax_coating.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0);
      rmax_vacuum.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 - coating_thickness);
      // inner radius
      rmin_wall.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0);
      rmin_coating.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 - coating_thickness);
      rmin_vacuum.push_back(0);
    }
    return std::tuple<Polycone, Polycone, Polycone>({0, 2.0 * M_PI, rmin_wall, rmax_wall, z},
                                                    {0, 2.0 * M_PI, rmin_coating, rmax_coating, z},
                                                    {0, 2.0 * M_PI, rmin_vacuum, rmax_vacuum, z});
  };
  //---------------------------------------------------------------------------------------------------------
  // Helper function to create a racetrack volumes
  auto create_racetrack_solids = [](xml::Component& x_racetrack) {
    // get geometry parameters
    double semiCircle_rmin   = getAttrOrDefault(x_racetrack, _Unicode(semiCircle_rmin), 0.0);
    double wall_thickness    = getAttrOrDefault(x_racetrack, _Unicode(wall_thickness), 0.0);
    double coating_thickness = getAttrOrDefault(x_racetrack, _Unicode(coating_thickness), 0.0);
    double length            = getAttrOrDefault(x_racetrack, _Unicode(length), 0.0);
    double rectangle_h       = getAttrOrDefault(x_racetrack, _Unicode(rectangle_h), 0.0);

    // create semicylinders and boxes
    Tube wall_semiCircle(0, semiCircle_rmin + wall_thickness, length / 2.);
    Box wall_rectangle(semiCircle_rmin + wall_thickness, rectangle_h / 2., length / 2.);

    Tube coating_semiCircle(0, semiCircle_rmin, length / 2.);
    Box coating_rectangle(semiCircle_rmin, rectangle_h / 2., length / 2.);

    Tube vacuum_semiCircle(0, semiCircle_rmin - coating_thickness, length / 2.);
    Box vacuum_rectangle(semiCircle_rmin - coating_thickness, rectangle_h / 2., length / 2.);

    // unite semicylinders and boxes
    Transform3D upShift_tf(RotationZYX(0, 0, 0), Position(0, rectangle_h / 2., 0));
    Transform3D downShift_tf(RotationZYX(0, 0, 0), Position(0, -rectangle_h / 2., 0));

    UnionSolid wall_solid, coating_solid, vacuum_solid;

    wall_solid = UnionSolid("wall_solid_tmp", wall_rectangle, wall_semiCircle, upShift_tf);
    wall_solid = UnionSolid("wall_solid", wall_solid, wall_semiCircle, downShift_tf);

    coating_solid =
        UnionSolid("coating_solid_tmp", coating_rectangle, coating_semiCircle, upShift_tf);
    coating_solid = UnionSolid("coating_solid", coating_solid, coating_semiCircle, downShift_tf);

    vacuum_solid = UnionSolid("vacuum_solid_tmp", vacuum_rectangle, vacuum_semiCircle, upShift_tf);
    vacuum_solid = UnionSolid("vacuum_solid", vacuum_solid, vacuum_semiCircle, downShift_tf);

    // subtract vacuum
    BooleanSolid wall_subtract, coating_subtract, vacuum_subtract;

    wall_subtract = SubtractionSolid("wall_subtract", wall_solid, coating_solid, Transform3D());
    coating_subtract =
        SubtractionSolid("coating_subtract", coating_solid, vacuum_solid, Transform3D());
    vacuum_subtract = vacuum_solid;

    return std::tuple<Solid, Solid, Solid>(wall_subtract, coating_subtract, vacuum_subtract);
  };
  //---------------------------------------------------------------------------------------------------------
  // Helper function to build cones for interface
  auto build_interface_cones = [](std::vector<double>& wall_interface_cone1_rmax,
                                  std::vector<double>& coating_interface_cone1_rmax,
                                  std::vector<double>& vacuum_interface_cone1_rmax,
                                  std::vector<double>& wall_interface_cone2_rmax,
                                  std::vector<double>& coating_interface_cone2_rmax,
                                  std::vector<double>& vacuum_interface_cone2_rmax,
                                  std::vector<double>& interface_rmin,
                                  std::vector<double>& interface_z, double offset1_y,
                                  double offset2_y) {
    // wall cones
    Polycone wall_interface_cone1(0, 2.0 * M_PI, interface_rmin, wall_interface_cone1_rmax,
                                  interface_z);
    Polycone wall_interface_cone2(0, 2.0 * M_PI, interface_rmin, wall_interface_cone2_rmax,
                                  interface_z);
    // coating cones
    Polycone coating_interface_cone1(0, 2.0 * M_PI, interface_rmin, coating_interface_cone1_rmax,
                                     interface_z);
    Polycone coating_interface_cone2(0, 2.0 * M_PI, interface_rmin, coating_interface_cone2_rmax,
                                     interface_z);
    // vacuum cones
    Polycone vacuum_interface_cone1(0, 2.0 * M_PI, interface_rmin, vacuum_interface_cone1_rmax,
                                    interface_z);
    Polycone vacuum_interface_cone2(0, 2.0 * M_PI, interface_rmin, vacuum_interface_cone2_rmax,
                                    interface_z);

    Transform3D tf1(RotationZYX(0, 0, 0), Position(0, offset1_y, 0));
    Transform3D tf2(RotationZYX(0, 0, 0), Position(0, offset2_y, 0));

    // ---- Create wall
    // unite wall cones
    auto wall_interface_final = UnionSolid(wall_interface_cone1, wall_interface_cone2, tf1);
    wall_interface_final      = UnionSolid(wall_interface_final, wall_interface_cone2, tf2);
    Solid wall_interface_vac_cut(
        wall_interface_final); // cut cut out the interface volume from lepton beam pipe vac
    // subtract coating cones
    wall_interface_final =
        SubtractionSolid(wall_interface_final, coating_interface_cone1, Transform3D());
    wall_interface_final = SubtractionSolid(wall_interface_final, coating_interface_cone2, tf1);
    wall_interface_final = SubtractionSolid(wall_interface_final, coating_interface_cone2, tf2);

    // --- Create coating
    // unite coating cones
    auto coating_interface_final =
        UnionSolid(coating_interface_cone1, coating_interface_cone2, tf1);
    coating_interface_final = UnionSolid(coating_interface_final, coating_interface_cone2, tf2);
    // subtract vacuum cones
    coating_interface_final =
        SubtractionSolid(coating_interface_final, vacuum_interface_cone1, Transform3D());
    coating_interface_final =
        SubtractionSolid(coating_interface_final, vacuum_interface_cone2, tf1);
    coating_interface_final =
        SubtractionSolid(coating_interface_final, vacuum_interface_cone2, tf2);

    // --- Create vacuum
    // unite vacuum cones
    auto vacuum_interface_final = UnionSolid(vacuum_interface_cone1, vacuum_interface_cone2, tf1);
    vacuum_interface_final      = UnionSolid(vacuum_interface_final, vacuum_interface_cone2, tf2);

    return std::tuple<Solid, Solid, Solid, Solid>(wall_interface_final, coating_interface_final,
                                                  vacuum_interface_final, wall_interface_vac_cut);
  };
  //---------------------------------------------------------------------------------------------------------
  // Helper function to create an interface between racetrack and circular pipes
  auto create_interface_solids = [&](xml::Component& x_racetrack) {
    // get geometry parameters
    double semiCircle_rmin   = getAttrOrDefault(x_racetrack, _Unicode(semiCircle_rmin), 2.3 * cm);
    double wall_thickness    = getAttrOrDefault(x_racetrack, _Unicode(wall_thickness), 1.0 * mm);
    double coating_thickness = getAttrOrDefault(x_racetrack, _Unicode(coating_thickness), 30 * um);
    double rectangle_h       = getAttrOrDefault(x_racetrack, _Unicode(rectangle_h), 1.6 * cm);
    double cylRadius1        = getAttrOrDefault(x_racetrack, _Unicode(cylRadius1), 6.2 / 2. * cm);
    double cylRadius2        = getAttrOrDefault(x_racetrack, _Unicode(cylRadius2), 2.6 / 2. * cm);
    double straight_pipe_endz =
        getAttrOrDefault(x_racetrack, _Unicode(straight_pipe_endz), 66.385 * cm);
    double offset_z = getAttrOrDefault(x_racetrack, _Unicode(offset_z), 72.385 * cm);
    double length   = getAttrOrDefault(x_racetrack, _Unicode(length), 125.420 * cm);
    double interface_length_1 =
        getAttrOrDefault(x_racetrack, _Unicode(interface_length_1), 6.0 * cm);
    double interface_length_2 =
        getAttrOrDefault(x_racetrack, _Unicode(interface_length_2), 13.495 * cm);

    if (cylRadius2 - rectangle_h / 2. < 0) {
      printout(ERROR, "IP6BeamPipe_geo", "Negative cone base size: (%f - %f/2) < 0", cylRadius2,
               rectangle_h);
      throw std::runtime_error("IP6BeamPipe failed to build a cone");
    }

    std::vector<double> interface_rmin = {0, 0};

    // interface-1
    std::vector<double> interface_z_1 = {straight_pipe_endz,
                                         straight_pipe_endz + interface_length_1};
    // interface-2
    std::vector<double> interface_z_2 = {offset_z + length, offset_z + length + interface_length_2};

    // --- Central cones
    // interface-1
    std::vector<double> wall_interface_cone1_rmax_1    = {cylRadius1 + wall_thickness,
                                                          semiCircle_rmin + wall_thickness};
    std::vector<double> coating_interface_cone1_rmax_1 = {cylRadius1, semiCircle_rmin};
    std::vector<double> vacuum_interface_cone1_rmax_1  = {cylRadius1 - coating_thickness,
                                                          semiCircle_rmin - coating_thickness};
    // interface-2
    std::vector<double> wall_interface_cone1_rmax_2    = {semiCircle_rmin + wall_thickness,
                                                          cylRadius2 + wall_thickness};
    std::vector<double> coating_interface_cone1_rmax_2 = {semiCircle_rmin, cylRadius2};
    std::vector<double> vacuum_interface_cone1_rmax_2  = {semiCircle_rmin - coating_thickness,
                                                          cylRadius2 - coating_thickness};

    // --- Shifted cones
    // interface-1
    std::vector<double> wall_interface_cone2_rmax_1    = {semiCircle_rmin + wall_thickness,
                                                          semiCircle_rmin + wall_thickness};
    std::vector<double> coating_interface_cone2_rmax_1 = {semiCircle_rmin, semiCircle_rmin};
    std::vector<double> vacuum_interface_cone2_rmax_1  = {semiCircle_rmin - coating_thickness,
                                                          semiCircle_rmin - coating_thickness};
    // interface-2
    std::vector<double> wall_interface_cone2_rmax_2 = {
        semiCircle_rmin + wall_thickness, (cylRadius2 - rectangle_h / 2.) + wall_thickness};
    std::vector<double> coating_interface_cone2_rmax_2 = {semiCircle_rmin,
                                                          (cylRadius2 - rectangle_h / 2.)};
    std::vector<double> vacuum_interface_cone2_rmax_2  = {
        semiCircle_rmin - coating_thickness, (cylRadius2 - rectangle_h / 2.) - coating_thickness};

    // vertical shifts
    double offset1_y = rectangle_h / 2.;
    double offset2_y = -rectangle_h / 2.;

    // -- Build cones
    // interface-1
    auto interface_cones_1 = build_interface_cones(
        wall_interface_cone1_rmax_1, coating_interface_cone1_rmax_1, vacuum_interface_cone1_rmax_1,
        wall_interface_cone2_rmax_1, coating_interface_cone2_rmax_1, vacuum_interface_cone2_rmax_1,
        interface_rmin, interface_z_1, offset1_y, offset2_y);
    // interface-2
    auto interface_cones_2 = build_interface_cones(
        wall_interface_cone1_rmax_2, coating_interface_cone1_rmax_2, vacuum_interface_cone1_rmax_2,
        wall_interface_cone2_rmax_2, coating_interface_cone2_rmax_2, vacuum_interface_cone2_rmax_2,
        interface_rmin, interface_z_2, offset1_y, offset2_y);

    return std::tuple<Solid, Solid, Solid, Solid, Solid, Solid, Solid, Solid>(
        std::get<0>(interface_cones_1), std::get<1>(interface_cones_1),
        std::get<2>(interface_cones_1), std::get<3>(interface_cones_1),
        std::get<0>(interface_cones_2), std::get<1>(interface_cones_2),
        std::get<2>(interface_cones_2), std::get<3>(interface_cones_2));
  };
  //---------------------------------------------------------------------------------------------------------
  // Helper function to create union of lepton and hadron volumes
  auto create_volumes = [&](const std::string& name, xml::Component& x_pipe1,
                            xml::Component& x_pipe2, xml_coll_t& x_additional_subtraction_i) {
    auto pipe1_polycones = zplane_to_polycones(x_pipe1);
    auto pipe2_polycones = zplane_to_polycones(x_pipe2);

    auto crossing_angle    = getAttrOrDefault(x_pipe2, _Unicode(crossing_angle), 0.0);
    auto axis_intersection = getAttrOrDefault(x_pipe2, _Unicode(axis_intersection), 0.0);

    auto tf = Transform3D(Position(0, 0, axis_intersection)) *
              Transform3D(RotationY(crossing_angle)) *
              Transform3D(Position(0, 0, -axis_intersection));

    // union of wall, coating, and vacuum
    BooleanSolid wall_union =
        UnionSolid(std::get<0>(pipe1_polycones), std::get<0>(pipe2_polycones), tf);
    BooleanSolid coating_union =
        UnionSolid(std::get<1>(pipe1_polycones), std::get<1>(pipe2_polycones), tf);
    BooleanSolid vacuum_union =
        UnionSolid(std::get<2>(pipe1_polycones), std::get<2>(pipe2_polycones), tf);

    BooleanSolid wall_racetrack_final, coating_racetrack_final;
    BooleanSolid wall_interface_final_1, coating_interface_final_1;
    BooleanSolid wall_interface_final_2, coating_interface_final_2;
    BooleanSolid vacuum_interface, lepton_pipe_vac;

    // --- Create a vacuum cylinder to place a lepton pipe inside and cut it out of the hadron cone vacuum
    // --- This helps to avoid problems with nested solid unions and subtractions for the hadron vacuum
    std::vector<double> posz_vac = {
        getAttrOrDefault(x_pipe1, _Unicode(lepton_pipe_vac_tube_startz), 66.10 * cm),
        getAttrOrDefault(x_pipe1, _Unicode(lepton_pipe_vac_tube_endz), 494.556 * cm)};
    std::vector<double> rmin_vac = {0, 0};
    std::vector<double> rmax_vac = {
        getAttrOrDefault(x_pipe1, _Unicode(ipBeampipe_ID), 6.2 / 2. * cm) / 2. +
            getAttrOrDefault(x_pipe1, _Unicode(wall_thickness), 1.0 * mm),
        getAttrOrDefault(x_pipe1, _Unicode(ipBeampipe_ID), 6.2 / 2. * cm) / 2. +
            getAttrOrDefault(x_pipe1, _Unicode(wall_thickness), 1.0 * mm)};
    Polycone lepton_pipe_vac_tube(0, 2.0 * M_PI, rmin_vac, rmax_vac, posz_vac);

    // downstream side is more complex and requires additional effort to create all the volumes
    if (name == "downstream") {
      // subtract the lepton beam pipe vacuum from the hadron beam pipe vacuum
      vacuum_union = SubtractionSolid("vacuum_union", vacuum_union, lepton_pipe_vac_tube);

      xml::Component racetrack_lepton_c = x_pipe1.child(_Unicode(racetrack_lepton));

      // ---- Read geometry parameters ----
      double cylRadius1 = // cylinder radius on the IP side
          getAttrOrDefault(racetrack_lepton_c, _Unicode(cylRadius1), 6.2 / 2. * cm);
      double rtRadius = // racetrack radius
          getAttrOrDefault(racetrack_lepton_c, _Unicode(semiCircle_rmin), 2.3 * cm);
      double wall_thickness = // wall thickness
          getAttrOrDefault(racetrack_lepton_c, _Unicode(wall_thickness), 1.0 * mm);
      double coating_thickness = // coating thickness
          getAttrOrDefault(racetrack_lepton_c, _Unicode(coating_thickness), 30.0 * um);
      double straight_pipe_startz = // straight pipe on the IP side, start position
          getAttrOrDefault(racetrack_lepton_c, _Unicode(straight_pipe_startz), 66.10 * cm);
      double straight_pipe_endz = // straight pipe on the IP side, end position
          getAttrOrDefault(racetrack_lepton_c, _Unicode(straight_pipe_endz), 66.385 * cm);
      double elliptical_cut_rx_1 = // elliptical cut (IP side) rX for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_rx_1), 0.305 * m);
      double elliptical_cut_ry_1 = // elliptical cut (IP side) rY for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_ry_1), 0.021 * m);
      double elliptical_cut_rx_2 = // elliptical cut (non-IP side) rX for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_rx_2), 0.152 * m);
      double elliptical_cut_ry_2 = // elliptical cut (non-IP side) rY for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_ry_2), 0.021 * m);
      double rectangular_cut_a = // rectangular cut A for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(rectangular_cut_a), 0.81 / 2. * m);
      double rectangular_cut_b = // rectangular cut B for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(rectangular_cut_b), 0.021 * m);
      double elliptical_cut_dz = // thickness of the cut for the hadron beam opening
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_dz), 0.7 * cm);
      double elliptical_cut_offset_z_1 = // offset of the elliptical cut (IP side)
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_offset_z_1), (0.976) * m);
      double elliptical_cut_offset_z_2 = // offset of the elliptical cut (non-IP side)
          getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_offset_z_2),
                           (0.976 + 0.810) * m);
      double rectangular_cut_offset_z = // offset of the rectangular cut
          getAttrOrDefault(racetrack_lepton_c, _Unicode(rectangular_cut_offset_z),
                           (0.976 + 0.810 / 2.) * m);

      // ---- Create racetrack solids ----
      auto racetrack_solids = create_racetrack_solids(racetrack_lepton_c);

      // ---- Create an interface between racetrack and cylindrical beam pipe ----
      auto interface_solids = create_interface_solids(racetrack_lepton_c);
      vacuum_interface      = // unite interface vacuum solids
          UnionSolid("vacuum_interface", std::get<2>(interface_solids),
                     std::get<6>(interface_solids), Transform3D());

      // --- Create a straight pipe on the IP side ---
      std::vector<double> z            = {straight_pipe_startz, straight_pipe_endz};
      std::vector<double> wall_rmin    = {cylRadius1, cylRadius1};
      std::vector<double> wall_rmax    = {cylRadius1 + wall_thickness, cylRadius1 + wall_thickness};
      std::vector<double> coating_rmin = {cylRadius1 - coating_thickness,
                                          cylRadius1 - coating_thickness};
      std::vector<double> coating_rmax = wall_rmin;

      Polycone wall_pipe(0, 2.0 * M_PI, wall_rmin, wall_rmax, z);
      Polycone coating_pipe(0, 2.0 * M_PI, coating_rmin, coating_rmax, z);

      Polycone wall_pipe_2(0, 2.0 * M_PI, wall_rmin, wall_rmax, z);

      // --- Unite racetrack and straight pipe with the rest of the pipe
      double offset_z = getAttrOrDefault(racetrack_lepton_c, _Unicode(offset_z), 0.0);
      double length   = getAttrOrDefault(racetrack_lepton_c, _Unicode(length), 0.0);
      Transform3D racetrack_tf(RotationZYX(0, 0, 0), Position(0, 0, offset_z + length / 2.));

      UnionSolid wall_racetrack_pipe =
          UnionSolid("wall_racetrack_pipe", wall_pipe, std::get<0>(racetrack_solids), racetrack_tf);
      UnionSolid coating_racetrack_pipe = UnionSolid("coating_racetrack_pipe", coating_pipe,
                                                     std::get<1>(racetrack_solids), racetrack_tf);

      // --- Subtract racetrack and interface sections from vacuum
      lepton_pipe_vac = SubtractionSolid("vacuum_racetrack_cut_1", lepton_pipe_vac_tube,
                                         wall_racetrack_pipe, Transform3D());
      lepton_pipe_vac = SubtractionSolid("vacuum_racetrack_cut_2", lepton_pipe_vac,
                                         coating_racetrack_pipe, Transform3D());
      lepton_pipe_vac = SubtractionSolid("vacuum_interface_cut_3", lepton_pipe_vac,
                                         std::get<3>(interface_solids), Transform3D());
      lepton_pipe_vac = SubtractionSolid("vacuum_interface_cut_4", lepton_pipe_vac,
                                         std::get<7>(interface_solids), Transform3D());
      // --- Subtract wall and coating from vacuum
      lepton_pipe_vac =
          SubtractionSolid("vacuum_interface_cut_5", lepton_pipe_vac, wall_union, Transform3D());
      lepton_pipe_vac =
          SubtractionSolid("vacuum_interface_cut_6", lepton_pipe_vac, coating_union, Transform3D());

      // create a cut volume - vacuum = two ellipses + rectangle
      EllipticalTube elliptical_cut_1("elliptical_cut_1", elliptical_cut_rx_1, elliptical_cut_ry_1,
                                      elliptical_cut_dz);
      EllipticalTube elliptical_cut_2("elliptical_cut_2", elliptical_cut_rx_2, elliptical_cut_ry_2,
                                      elliptical_cut_dz);
      Box rectangular_cut_3("rectangular_cut_3", rectangular_cut_a, rectangular_cut_b,
                            elliptical_cut_dz);
      Transform3D tf_cut_1(RotationZYX(0, M_PI_2, 0),
                           Position(-rtRadius, 0, elliptical_cut_offset_z_1));
      Transform3D tf_cut_2(RotationZYX(0, M_PI_2, 0),
                           Position(-rtRadius, 0, elliptical_cut_offset_z_2));
      Transform3D tf_cut_3(RotationZYX(0, M_PI_2, 0),
                           Position(-rtRadius, 0, rectangular_cut_offset_z));

      // subtract from racetrack wall
      SubtractionSolid wall_racetrack_cut_1 =
          SubtractionSolid("wall_racetrack_cut_1", wall_racetrack_pipe, elliptical_cut_1, tf_cut_1);
      SubtractionSolid wall_racetrack_cut_2 = SubtractionSolid(
          "wall_racetrack_cut_1", wall_racetrack_cut_1, elliptical_cut_2, tf_cut_2);
      wall_racetrack_final = SubtractionSolid("wall_racetrack_final", wall_racetrack_cut_2,
                                              rectangular_cut_3, tf_cut_3);

      // subtract from racetrack coating
      SubtractionSolid coating_racetrack_cut_1 = SubtractionSolid(
          "coating_racetrack_cut_1", coating_racetrack_pipe, elliptical_cut_1, tf_cut_1);
      SubtractionSolid coating_racetrack_cut_2 = SubtractionSolid(
          "coating_racetrack_cut_2", coating_racetrack_cut_1, elliptical_cut_2, tf_cut_2);
      coating_racetrack_final = SubtractionSolid("coating_racetrack_final", coating_racetrack_cut_2,
                                                 rectangular_cut_3, tf_cut_3);
      // unite with lepton pipe vacuum
      lepton_pipe_vac =
          UnionSolid("lepton_pipe_vac_uni_7", lepton_pipe_vac, elliptical_cut_1, tf_cut_1);
      lepton_pipe_vac =
          UnionSolid("lepton_pipe_vac_uni_8", lepton_pipe_vac, elliptical_cut_2, tf_cut_2);
      lepton_pipe_vac =
          UnionSolid("lepton_pipe_vac_uni_9", lepton_pipe_vac, rectangular_cut_3, tf_cut_3);

      // subtract from hadron cone pipe vacuum
      vacuum_union =
          SubtractionSolid("vacuum_union_sub_1", vacuum_union, elliptical_cut_1, tf_cut_1);
      vacuum_union =
          SubtractionSolid("vacuum_union_sub_2", vacuum_union, elliptical_cut_2, tf_cut_2);
      vacuum_union =
          SubtractionSolid("vacuum_union_sub_3", vacuum_union, rectangular_cut_3, tf_cut_3);

      // subtract from interface-1 wall
      wall_interface_final_1 = SubtractionSolid(
          "wall_interface_final_1", std::get<0>(interface_solids), elliptical_cut_1, tf_cut_1);

      // subtract from interface-1 coating
      coating_interface_final_1 = SubtractionSolid(
          "coating_interface_final_1", std::get<1>(interface_solids), elliptical_cut_1, tf_cut_1);

      // subtract from interface vacuum
      vacuum_interface =
          SubtractionSolid("vacuum_interface", vacuum_interface, elliptical_cut_1, tf_cut_1);

      // get  interface-2 w/o subtraction
      wall_interface_final_2    = std::get<4>(interface_solids);
      coating_interface_final_2 = std::get<5>(interface_solids);
    }

    Solid wall_racetrack(wall_racetrack_final);
    Solid coating_racetrack(coating_racetrack_final);
    Solid wall_interface_1(wall_interface_final_1);
    Solid coating_interface_1(coating_interface_final_1);
    Solid wall_interface_2(wall_interface_final_2);
    Solid coating_interface_2(coating_interface_final_2);

    // subtract additional vacuum from wall and coating
    for (; x_additional_subtraction_i; ++x_additional_subtraction_i) {
      xml_comp_t x_additional_subtraction = x_additional_subtraction_i;
      auto additional_polycones           = zplane_to_polycones(x_additional_subtraction);
      auto additional_crossing_angle =
          getAttrOrDefault(x_additional_subtraction, _Unicode(crossing_angle), 0.0);
      auto additional_axis_intersection =
          getAttrOrDefault(x_additional_subtraction, _Unicode(axis_intersection), 0.0);
      auto additional_tf = Transform3D(Position(0, 0, additional_axis_intersection)) *
                           Transform3D(RotationY(additional_crossing_angle)) *
                           Transform3D(Position(0, 0, -additional_axis_intersection));

      wall_union = SubtractionSolid(wall_union, std::get<2>(additional_polycones), additional_tf);
      coating_union =
          SubtractionSolid(coating_union, std::get<2>(additional_polycones), additional_tf);
      vacuum_union =
          SubtractionSolid(vacuum_union, std::get<2>(additional_polycones), additional_tf);
    }

    Solid wall, coating, vacuum, wall_ipflange, coating_ipflange, vacuum_ipflange;

    if (name == "upstream") // upstream
    {
      // subtract vacuum from coating
      coating = SubtractionSolid(coating_union, vacuum_union);

      // subtract vacuum from wall
      wall = SubtractionSolid(wall_union, vacuum_union);

      // get vacuum
      vacuum = vacuum_union;
    } else // downstream
    {
      // subtract walls and coatings from vacuum
      vacuum_union =
          SubtractionSolid("vacuum_union_sub_2", vacuum_union, wall_union, Transform3D());
      vacuum_union =
          SubtractionSolid("vacuum_union_sub_3", vacuum_union, coating_union, Transform3D());

      // get vacuum
      vacuum = vacuum_union;

      // get wall
      wall = wall_union;

      // get coating
      coating = coating_union;

      // get a flange between the FWD and IP beam pipes
      xml::Component fwdipflange_c = x_pipe1.child(_Unicode(fwdipflange));
      auto fwdipflange_polycones   = zplane_to_polycones(fwdipflange_c);

      wall_ipflange    = std::get<0>(fwdipflange_polycones);
      coating_ipflange = std::get<1>(fwdipflange_polycones);
      vacuum_ipflange  = std::get<2>(fwdipflange_polycones);
    }

    return std::tuple<Volume, Volume, Volume, Volume, Volume, Volume, Volume, Volume, Volume,
                      Volume, Volume, Volume, Volume, Volume>(
        {"v_" + name + "_wall", wall, m_Wall}, {"v_" + name + "_coating", coating, m_Coating},
        {"v_" + name + "_vacuum", vacuum, m_Vacuum},
        {"v_" + name + "_vacuum_ipflange", vacuum_ipflange, m_Vacuum},
        {"v_" + name + "_wall_racetrack", wall_racetrack, m_Wall},
        {"v_" + name + "_coating_racetrack", coating_racetrack, m_Coating},
        {"v_" + name + "_wall_ipflange", wall_ipflange, m_Wall},
        {"v_" + name + "_coating_ipflange", coating_ipflange, m_Coating},
        {"v_" + name + "_wall_interface_1", wall_interface_1, m_Wall},
        {"v_" + name + "_coating_interface_1", coating_interface_1, m_Coating},
        {"v_" + name + "_wall_interface_2", wall_interface_2, m_Wall},
        {"v_" + name + "_coating_interface_2", coating_interface_2, m_Coating},
        {"v_" + name + "_vacuum_interface", vacuum_interface, m_Vacuum},
        {"v_" + name + "_lepton_pipe_vac", lepton_pipe_vac, m_Vacuum});
  };

  // -----------------------------
  // Upstream/BWD/Rear Side:
  // - incoming hadron tube: straight
  // - outgoing electron tube: tapering

  xml::Component upstream_c        = x_det.child(_Unicode(upstream));
  xml::Component incoming_hadron_c = upstream_c.child(_Unicode(incoming_hadron));
  xml::Component outgoing_lepton_c = upstream_c.child(_Unicode(outgoing_lepton));
  xml_coll_t additional_subtractions_upstream(upstream_c, _Unicode(additional_subtraction));
  auto volumes_upstream = create_volumes("upstream", outgoing_lepton_c, incoming_hadron_c,
                                         additional_subtractions_upstream);

  std::get<0>(volumes_upstream).setVisAttributes(wallVis);
  std::get<1>(volumes_upstream).setVisAttributes(coatingVis);

  auto tf_upstream = Transform3D(RotationZYX(0, 0, 0));
  if (getAttrOrDefault<bool>(upstream_c, _Unicode(reflect), true)) {
    tf_upstream = Transform3D(RotationZYX(0, M_PI, 0));
  }
  // place wall
  assembly.placeVolume(std::get<0>(volumes_upstream), tf_upstream);
  // place coating
  assembly.placeVolume(std::get<1>(volumes_upstream), tf_upstream);
  // place vacuum
  if (getAttrOrDefault<bool>(upstream_c, _Unicode(place_vacuum), true)) {
    assembly.placeVolume(std::get<2>(volumes_upstream), tf_upstream);
  }

  // -----------------------------
  // Downstream/FWD Side:
  // - incoming electron tube: tube with tube cut out
  // - outgoing hadron tube: cone centered at scattering angle

  xml::Component downstream_c      = x_det.child(_Unicode(downstream));
  xml::Component incoming_lepton_c = downstream_c.child(_Unicode(incoming_lepton));
  xml::Component outgoing_hadron_c = downstream_c.child(_Unicode(outgoing_hadron));
  xml_coll_t additional_subtractions_downstream(downstream_c, _Unicode(additional_subtraction));
  auto volumes_downstream = create_volumes("downstream", incoming_lepton_c, outgoing_hadron_c,
                                           additional_subtractions_downstream);

  // add colors
  std::get<0>(volumes_downstream).setVisAttributes(wallVis);
  std::get<1>(volumes_downstream).setVisAttributes(coatingVis);
  std::get<4>(volumes_downstream).setVisAttributes(wallVis);
  std::get<5>(volumes_downstream).setVisAttributes(coatingVis);
  std::get<6>(volumes_downstream).setVisAttributes(wallVis);
  std::get<7>(volumes_downstream).setVisAttributes(coatingVis);
  std::get<8>(volumes_downstream).setVisAttributes(wallVis);
  std::get<9>(volumes_downstream).setVisAttributes(coatingVis);
  std::get<10>(volumes_downstream).setVisAttributes(wallVis);
  std::get<11>(volumes_downstream).setVisAttributes(coatingVis);

  auto tf_downstream = Transform3D(RotationZYX(0, 0, 0));
  if (getAttrOrDefault<bool>(downstream_c, _Unicode(reflect), true)) {
    tf_downstream = Transform3D(RotationZYX(0, M_PI, 0));
  }

  // place wall
  assembly.placeVolume(std::get<0>(volumes_downstream), tf_downstream);
  // place coating
  assembly.placeVolume(std::get<1>(volumes_downstream), tf_downstream);
  // place vacuum
  if (getAttrOrDefault<bool>(downstream_c, _Unicode(place_vacuum), true)) {
    assembly.placeVolume(std::get<2>(volumes_downstream), tf_downstream);
    assembly.placeVolume(std::get<3>(volumes_downstream), tf_downstream);
    assembly.placeVolume(std::get<12>(volumes_downstream), tf_downstream);
    assembly.placeVolume(std::get<13>(volumes_downstream), tf_downstream);
  }
  // place racetrack wall
  assembly.placeVolume(std::get<4>(volumes_downstream), tf_downstream);
  // place racetrack coating
  assembly.placeVolume(std::get<5>(volumes_downstream), tf_downstream);
  // place FWD IP flange wall
  assembly.placeVolume(std::get<6>(volumes_downstream), tf_downstream);
  // place FWD IP flange coating
  assembly.placeVolume(std::get<7>(volumes_downstream), tf_downstream);
  // place interface wall
  assembly.placeVolume(std::get<8>(volumes_downstream), tf_downstream);
  assembly.placeVolume(std::get<10>(volumes_downstream), tf_downstream);
  // place interface coating
  assembly.placeVolume(std::get<9>(volumes_downstream), tf_downstream);
  assembly.placeVolume(std::get<11>(volumes_downstream), tf_downstream);

  // delete temporaries
  wallVis.destroy();
  coatingVis.destroy();
  IPwallVis.destroy();
  IPcoatingVis.destroy();

  // -----------------------------
  // final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  pv_assembly.addPhysVolID("system", sdet.id());
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(IP6BeamPipe, create_detector)
