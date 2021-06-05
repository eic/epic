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
  Material   m_Cu    = det.material("Copper");
  Material   m_Al    = det.material("Aluminum");
  Material   m_Be    = det.material("Beryllium");
  string     vis_name  = x_det.visStr();

  int n = 0;

  xml::Component IP_pipe_c = x_det.child(_Unicode(IP_pipe));
  xml::Component upstream_c = x_det.child(_Unicode(upstream_pipe));
  xml::Component downstream_c = x_det.child(_Unicode(downstream_pipe));

  // IP
  double IP_beampipe_OD        = IP_pipe_c.attr<double>(_Unicode(OD));
  double IP_beampipe_thickness = IP_pipe_c.attr<double>(_Unicode(wall_thickness));
  double IP_beampipe_ID        = IP_beampipe_OD - IP_beampipe_thickness;
  double crossing_angle        = IP_pipe_c.attr<double>(_Unicode(crossing_angle));

  // upstream parameters
  double upstream_straight_length  = upstream_c.attr<double>(_Unicode(straight_length));;
  double upstream_total_length     = upstream_c.attr<double>(_Unicode(length));
  double upstream_conic_length     = upstream_total_length - upstream_straight_length;
  double upstream_beampipe_exit_OD = 229.96 * mm;
  double upstream_beampipe_exit_ID = upstream_beampipe_exit_OD - IP_beampipe_thickness;
  double upstream_epipe_thickness  = 4.0*mm;

  // downstream parameters
  double downstream_straight_length = downstream_c.attr<double>(_Unicode(straight_length));
  double downstream_taper_length    = 310.0 * mm;
  double downstream_total_length    = downstream_c.attr<double>(_Unicode(length));
  double downstream_conic_length    = downstream_total_length - downstream_straight_length ;
  double downstream_cone_rmax       = 110.00*mm/2.0;
  double downstream_epipe_thickness  = 4.0*mm; // just a guess
  double downstream_beampipe_exit_OD = 62.0 * mm;
  double downstream_beampipe_exit_ID = downstream_beampipe_exit_OD - downstream_epipe_thickness;
  double downstream_hpipe_OD         = IP_beampipe_ID/2; // just a guess
  double downstream_hpipe_thickness  = 4.0*mm; // just a guess

  double upstream_delta_r   = upstream_conic_length*std::tan(crossing_angle/2.0);
  double downstream_delta_r = (downstream_conic_length+downstream_straight_length)*std::tan(crossing_angle);


  // -----------------------------
  // IP beampipe
  Tube downstream_IP_tube(IP_beampipe_ID/2.0, IP_beampipe_OD/2.0, downstream_straight_length/2.0);
  Tube downstream_IP_vacuum(0.0, IP_beampipe_ID/2.0, downstream_straight_length/2.0);
  Tube upstream_IP_tube(IP_beampipe_ID/2.0, IP_beampipe_OD/2.0, upstream_straight_length/2.0);
  Tube upstream_IP_vacuum(0.0, IP_beampipe_OD/2.0, upstream_straight_length/2.0);

  Volume v_upstream_IP_tube("v_upstream_IP_tube", upstream_IP_tube, m_Be);
  Volume v_downstream_IP_tube("v_downstream_IP_tube", downstream_IP_tube, m_Be);

  //v_upstream_IP_tube.setVisAttributes(det,"GrayVis");
  //v_downstream_IP_tube.setVisAttributes(det,"RedVis");
  sdet.setAttributes(det, v_upstream_IP_tube  , x_det.regionStr(), x_det.limitsStr(), vis_name);
  sdet.setAttributes(det, v_downstream_IP_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  auto pv_upstream_IP_tube = assembly.placeVolume( v_upstream_IP_tube, Position(0, 0, -upstream_straight_length / 2.0));

  auto pv_downstream_IP_tube = assembly.placeVolume(
      v_downstream_IP_tube, Position(0, 0, downstream_straight_length / 2.0));

  // -----------------------------
  // upstream
  Tube upstream_electron_tube(IP_beampipe_ID/2.0 - upstream_epipe_thickness,
                              IP_beampipe_ID/2.0-0.01*upstream_epipe_thickness,  // leave a gap 1% of pipe thickness
                              upstream_conic_length/2.0);
  Cone upstream_conic_section(upstream_conic_length / 2.0,
                              IP_beampipe_ID / 2.0, IP_beampipe_OD / 2.0,
                              IP_beampipe_ID / 2.0 + upstream_delta_r,
                              IP_beampipe_OD / 2.0+ upstream_delta_r  + IP_beampipe_thickness);
  Cone upstream_conic_section_vacuum(upstream_conic_length / 2.0,
                                     0.0, IP_beampipe_ID / 2.0,
                                     0.0, IP_beampipe_ID / 2.0 + upstream_delta_r);
  Volume v_upstream_conic_section("v_upstream_conic_section", upstream_conic_section, m_Be);
  Volume v_upstream_electron_tube("v_upstream_electron_tube", upstream_electron_tube, m_Be);
  sdet.setAttributes(det, v_upstream_conic_section, x_det.regionStr(), x_det.limitsStr(), vis_name);
  sdet.setAttributes(det, v_upstream_electron_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);
  //Volume v_upstream_conic_section_vacuum("v_upstream_conic_section_vacuum", upstream_electron_tube, m_Al);

  //auto pv_upstream_conic_section = assembly.placeVolume(
  //    v_upstream_conic_section,
  //    Transform3D(
  //        Position(-upstream_delta_r/2.0,0, -upstream_straight_length - upstream_conic_length / 2.0)) *
  //        RotationY(crossing_angle / 2.0) * RotationX(M_PI));
  auto pv_upstream_electron_tube =
      assembly.placeVolume(v_upstream_electron_tube,
                           Position(0, 0, -upstream_straight_length - upstream_conic_length / 2.0));
  //DetElement de_upstream_conic_section(sdet,"de_upstream_conic_section",1);
  //de_upstream_conic_section.setPlacement(pv_upstream_conic_section);

  // -----------------------------
  // downstream
  Tube downstream_hadron_tube(downstream_hpipe_OD/2.0 - downstream_epipe_thickness,
                                downstream_hpipe_OD/2.0, downstream_conic_length / 2.0);
  Cone downstream_conic_section(downstream_conic_length / 2.0, IP_beampipe_ID / 2.0,
                                IP_beampipe_OD / 2.0, downstream_cone_rmax,
                                downstream_cone_rmax + downstream_epipe_thickness);
  Cone downstream_taper_section(downstream_taper_length / 2.0, 
                                  downstream_cone_rmax, downstream_cone_rmax + downstream_epipe_thickness,
                                  downstream_beampipe_exit_ID/2.0, downstream_beampipe_exit_OD/2.0);

  UnionSolid downstream_pipe_split0(downstream_conic_section, downstream_hadron_tube,
                                    Transform3D(Position(downstream_delta_r / 2.0, 0.0, 0.0)) *
                                    RotationY(crossing_angle));

  // Vacuum
  Cone downstream_conic_section_vacuum(downstream_conic_length / 2.0, 
                                       0.0, IP_beampipe_ID / 2.0,
                                       0.0, downstream_cone_rmax);
  Tube downstream_hadron_vacuum(0.0, downstream_hpipe_OD / 2.0 - downstream_epipe_thickness,
                                downstream_conic_length / 1.9);
  UnionSolid downstream_pipe_vacuum_split0(downstream_conic_section_vacuum, downstream_hadron_vacuum,
                                           Transform3D(Position(downstream_delta_r / 2.0, 0.0, 0.0)) *
                                           RotationY(crossing_angle));
  SubtractionSolid downstream_pipe_split1(downstream_pipe_split0, 
                                          downstream_pipe_vacuum_split0);
  Volume v_downstream_pipe_split1("v_downstream_pipe_split1", downstream_pipe_split1, m_Be);
  sdet.setAttributes(det, v_downstream_pipe_split1, x_det.regionStr(), x_det.limitsStr(), vis_name);
  auto pv_downstream_pipe_split1 = assembly.placeVolume(v_downstream_pipe_split1,
                                                        Position(0, 0, 
                                                                 downstream_straight_length + 
                                                                 downstream_conic_length / 2.0));
  Volume v_downstream_taper_section("v_downstream_taper_section", downstream_taper_section, m_Be);
  sdet.setAttributes(det, v_downstream_taper_section, x_det.regionStr(), x_det.limitsStr(), vis_name);
  auto   pv_downstream_taper_section = assembly.placeVolume(
      v_downstream_taper_section, Position(0, 0,
                                           downstream_straight_length + downstream_conic_length +
                                               downstream_taper_length / 2.0));

  //Volume v_downstream_conic_section("v_downstream_conic_section", downstream_conic_section, m_Al);
  //Volume v_downstream_hadron_tube("v_downstream_hadron_tube", downstream_hadron_tube, m_Al);
  //auto pv_downstream_hadron_tube = assembly.placeVolume(
  //    v_downstream_hadron_tube,
  //    Transform3D(Position(downstream_delta_r / 2.0, 0,
  //                         downstream_straight_length + downstream_conic_length / 2.0)) *
  //        RotationY(crossing_angle) );

  //auto pv_downstream_conic_section = assembly.placeVolume(
  //    v_downstream_conic_section, Position(0, 0, downstream_straight_length + downstream_conic_length / 2.0));
  //DetElement de_downstream_conic_section(sdet,"de_downstream_conic_section",1);
  //de_downstream_conic_section.setPlacement(pv_downstream_conic_section);

  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  pv_assembly.addPhysVolID("system",sdet.id()).addPhysVolID("barrel",1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox() ;
  return sdet;
}

DECLARE_DETELEMENT(IP6BeamPipe,create_detector)
