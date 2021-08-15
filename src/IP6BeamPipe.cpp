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
static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector sens)  {

  using namespace ROOT::Math;
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = det.air();
  DetElement sdet        (det_name,x_det.id());
  Assembly   assembly    (det_name+"_assembly");
  Material   m_Cu      = det.material("Copper");
  Material   m_Al      = det.material("Aluminum");
  Material   m_Be      = det.material("Beryllium");
  Material   m_Au      = det.material("Gold");
  Material   m_Vacuum  = det.material("Vacuum");
  string     vis_name  = x_det.visStr();

  xml::Component IP_pipe_c = x_det.child(_Unicode(IP_pipe));

  // IP
  double IP_beampipe_OD             = IP_pipe_c.attr<double>(_Unicode(OD));
  double IP_beampipe_wall_thickness = IP_pipe_c.attr<double>(_Unicode(wall_thickness));
  double IP_beampipe_gold_thickness = IP_pipe_c.attr<double>(_Unicode(gold_thickness));
  double IP_beampipe_ID             = IP_beampipe_OD - IP_beampipe_gold_thickness - IP_beampipe_wall_thickness;

  double upstream_straight_length   = IP_pipe_c.attr<double>(_Unicode(upstream_straight_length));
  double downstream_straight_length = IP_pipe_c.attr<double>(_Unicode(downstream_straight_length));

  // -----------------------------
  // IP beampipe
  Tube downstream_IP_vacuum(0.0, IP_beampipe_ID/2.0, downstream_straight_length/2.0);
  Tube downstream_IP_gold(IP_beampipe_ID/2.0, IP_beampipe_ID/2.0 + IP_beampipe_gold_thickness, downstream_straight_length/2.0);
  Tube downstream_IP_tube(IP_beampipe_ID/2.0 + IP_beampipe_gold_thickness, IP_beampipe_OD/2.0, downstream_straight_length/2.0);
  Tube upstream_IP_vacuum(0.0, IP_beampipe_ID/2.0, upstream_straight_length/2.0);
  Tube upstream_IP_gold(IP_beampipe_ID/2.0, IP_beampipe_ID/2.0 + IP_beampipe_gold_thickness, upstream_straight_length/2.0);
  Tube upstream_IP_tube(IP_beampipe_ID/2.0 + IP_beampipe_gold_thickness, IP_beampipe_OD/2.0, upstream_straight_length/2.0);

  Volume v_downstream_IP_vacuum("v_downstream_IP_vacuum", downstream_IP_vacuum, m_Vacuum);
  Volume v_downstream_IP_gold("v_downstream_IP_gold", downstream_IP_gold, m_Au);
  Volume v_downstream_IP_tube("v_downstream_IP_tube", downstream_IP_tube, m_Be);
  Volume v_upstream_IP_vacuum("v_upstream_IP_vacuum", upstream_IP_vacuum, m_Vacuum);
  Volume v_upstream_IP_gold("v_upstream_IP_gold", upstream_IP_gold, m_Au);
  Volume v_upstream_IP_tube("v_upstream_IP_tube", upstream_IP_tube, m_Be);

  sdet.setAttributes(det, v_upstream_IP_gold  , x_det.regionStr(), x_det.limitsStr(), vis_name);
  sdet.setAttributes(det, v_upstream_IP_tube  , x_det.regionStr(), x_det.limitsStr(), vis_name);
  sdet.setAttributes(det, v_downstream_IP_gold, x_det.regionStr(), x_det.limitsStr(), vis_name);
  sdet.setAttributes(det, v_downstream_IP_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  assembly.placeVolume(v_upstream_IP_vacuum, Position(0, 0, -upstream_straight_length / 2.0));
  assembly.placeVolume(v_upstream_IP_gold, Position(0, 0, -upstream_straight_length / 2.0));
  assembly.placeVolume(v_upstream_IP_tube, Position(0, 0, -upstream_straight_length / 2.0));

  assembly.placeVolume(v_downstream_IP_vacuum, Position(0, 0, downstream_straight_length / 2.0));
  assembly.placeVolume(v_downstream_IP_gold, Position(0, 0, downstream_straight_length / 2.0));
  assembly.placeVolume(v_downstream_IP_tube, Position(0, 0, downstream_straight_length / 2.0));


  // Helper function to create polycone pairs (shell and vacuum)
  auto zplane_to_polycones = [](xml::Component& x_pipe) {
    std::vector<double> zero, rmax, rmin, z;
    for (xml_coll_t x_zplane_i(x_pipe, _Unicode(zplane)); x_zplane_i; ++x_zplane_i) {
      xml_comp_t x_zplane = x_zplane_i;
      auto thickness = getAttrOrDefault(x_zplane, _U(thickness), x_pipe.thickness());
      zero.push_back(0);
      rmax.push_back(x_zplane.attr<double>(_Unicode(OD)) / 2.0);
      rmin.push_back(x_zplane.attr<double>(_Unicode(OD)) / 2.0 - thickness);
      z.push_back(x_zplane.attr<double>(_Unicode(z)));
    }
    return std::make_pair<Polycone,Polycone>(
      {0, 2.0 * M_PI, rmin, rmax, z},
      {0, 2.0 * M_PI, zero, rmin, z}
    );
  };

  // -----------------------------
  // Upstream:
  // - incoming hadron tube: straight section, tapered section, straight section
  // - outgoing electron tube: tapered section, straight section

  xml::Component upstream_c = x_det.child(_Unicode(upstream));
  xml::Component incoming_hadron_c = upstream_c.child(_Unicode(incoming_hadron));
  xml::Component outgoing_lepton_c = upstream_c.child(_Unicode(outgoing_lepton));
  auto outgoing_lepton_polycones = zplane_to_polycones(outgoing_lepton_c);
  auto incoming_hadron_polycones = zplane_to_polycones(incoming_hadron_c);

  auto incoming_crossing_angle = getAttrOrDefault(incoming_hadron_c, _Unicode(crossing_angle), 0.0);
  auto incoming_axis_intersection = getAttrOrDefault(incoming_hadron_c, _Unicode(axis_intersection), 0.0);


  auto upstream_tf = Transform3D(Position(0,0,incoming_axis_intersection)) *
                     Transform3D(RotationY(incoming_crossing_angle)) *
                     Transform3D(Position(0,0,-incoming_axis_intersection));

  UnionSolid upstream_matter(outgoing_lepton_polycones.first,
                             incoming_hadron_polycones.first,
                             upstream_tf);
  UnionSolid upstream_vacuum(outgoing_lepton_polycones.second,
                             incoming_hadron_polycones.second,
                             upstream_tf);
  SubtractionSolid upstream(upstream_matter, upstream_vacuum);

  Volume v_upstream("v_upstream", upstream, m_Al);
  Volume v_upstream_vacuum("v_upstream_vacuum", upstream_vacuum, m_Vacuum);

  if (upstream_c.reflect(true)) {
    auto tf = Transform3D(RotationZYX(0, M_PI, 0));
    assembly.placeVolume(v_upstream, tf);
    if (getAttrOrDefault<bool>(upstream_c, _Unicode(place_vacuum), true)) {
      assembly.placeVolume(v_upstream_vacuum, tf);
    }
  } else {
    auto tf = Transform3D(RotationZYX(0, 0, 0));
    assembly.placeVolume(v_upstream, tf);
    if (getAttrOrDefault<bool>(upstream_c, _Unicode(place_vacuum), true)) {
      assembly.placeVolume(v_upstream_vacuum, tf);
    }
  }

  // -----------------------------
  // downstream:
  // - incoming electron tube: tube with tube cut out
  // - outgoing hadron tube: cone centered at scattering angle
  // (incoming electron tube internal touching to outgoing hadron tube)

  xml::Component downstream_c = x_det.child(_Unicode(downstream));
  xml::Component incoming_lepton_c = downstream_c.child(_Unicode(incoming_lepton));
  xml::Component outgoing_hadron_c = downstream_c.child(_Unicode(outgoing_hadron));
  auto incoming_lepton_polycones = zplane_to_polycones(incoming_lepton_c);
  auto outgoing_hadron_polycones = zplane_to_polycones(outgoing_hadron_c);

  auto outgoing_crossing_angle = getAttrOrDefault(outgoing_hadron_c, _Unicode(crossing_angle), 0.0);
  auto outgoing_axis_intersection = getAttrOrDefault(outgoing_hadron_c, _Unicode(axis_intersection), 0.0);

  auto downstream_tf = Transform3D(Position(0,0,outgoing_axis_intersection)) *
                       Transform3D(RotationY(outgoing_crossing_angle)) *
                       Transform3D(Position(0,0,-outgoing_axis_intersection));

  UnionSolid downstream_matter(incoming_lepton_polycones.first,
                               outgoing_hadron_polycones.first,
                               downstream_tf);
  UnionSolid downstream_vacuum(incoming_lepton_polycones.second,
                               outgoing_hadron_polycones.second,
                               downstream_tf);

  // subtract vacuum
  BooleanSolid downstream;
  if (getAttrOrDefault<bool>(downstream_c, _Unicode(subtract_vacuum), true)) {
    downstream = SubtractionSolid (downstream_matter, downstream_vacuum);
  } else {
    downstream = downstream_matter;
  }

  // subtract additional tubes
  for (xml_coll_t x_additional_subtraction_i(downstream_c, _Unicode(additional_subtraction)); x_additional_subtraction_i; ++x_additional_subtraction_i) {
    xml_comp_t x_additional_subtraction = x_additional_subtraction_i;
    auto additional_subtraction_polycones = zplane_to_polycones(x_additional_subtraction);
    auto crossing_angle = getAttrOrDefault(x_additional_subtraction, _Unicode(crossing_angle), 0.0);
    auto axis_intersection = getAttrOrDefault(x_additional_subtraction, _Unicode(axis_intersection), 0.0);
    auto tf = Transform3D(Position(0,0,axis_intersection)) *
                         Transform3D(RotationY(crossing_angle)) *
                         Transform3D(Position(0,0,-axis_intersection));
    downstream = SubtractionSolid(downstream, additional_subtraction_polycones.second, tf);
  }

  Volume v_downstream("v_downstream", downstream, m_Al);
  Volume v_downstream_vacuum("v_downstream_vacuum", downstream_vacuum, m_Vacuum);

  if (downstream_c.reflect(true)) {
    auto tf = Transform3D(RotationZYX(0, M_PI, 0));
    assembly.placeVolume(v_downstream, tf);
    if (getAttrOrDefault<bool>(downstream_c, _Unicode(place_vacuum), true)) {
      assembly.placeVolume(v_downstream_vacuum, tf);
    }
  } else {
    auto tf = Transform3D(RotationZYX(0, 0, 0));
    assembly.placeVolume(v_downstream, tf);
    if (getAttrOrDefault<bool>(downstream_c, _Unicode(place_vacuum), true)) {
      assembly.placeVolume(v_downstream_vacuum, tf);
    }
  }

  // -----------------------------
  // final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  pv_assembly.addPhysVolID("system",sdet.id()).addPhysVolID("barrel",1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox() ;
  return sdet;
}

DECLARE_DETELEMENT(IP6BeamPipe,create_detector)
