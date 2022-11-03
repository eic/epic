//==========================================================================
//   Beampipe described by end points with different radii at either end
//==========================================================================
//
//      <detector id="Pipe_in" name ="DetName" type="BackwardsBeamPipe" >
//        <Pipe wall_thickness="pipe_thickness" outerD1="start_radius" outerD2="end_radius"
//        end1z="start_z" end2z="end_z" end1x="start_x" end2x="end_x"/>
//      </detector>
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Al     = det.material("Aluminum");
  Material   m_Vacuum = det.material("Vacuum");
  string     vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");

  xml::Component Pipe_c = x_det.child(_Unicode(Pipe));

  // Get pipe dimensions from xml
  double thickness = Pipe_c.attr<double>(_Unicode(wall_thickness));
  double innerD1   = Pipe_c.hasAttr(_Unicode(innerD1)) ? Pipe_c.attr<double>(_Unicode(innerD1))
                                                          : Pipe_c.attr<double>(_Unicode(outerD1)) - 2 * thickness;
  double innerD2   = Pipe_c.hasAttr(_Unicode(innerD2)) ? Pipe_c.attr<double>(_Unicode(innerD2))
                                                          : Pipe_c.attr<double>(_Unicode(outerD2)) - 2 * thickness;

  double end1z     = Pipe_c.attr<double>(_Unicode(end1z));
  double end2z     = Pipe_c.attr<double>(_Unicode(end2z));
  double end1x     = dd4hep::getAttrOrDefault<double>(Pipe_c, _Unicode(end1x), 0.0 );
  double end2x     = dd4hep::getAttrOrDefault<double>(Pipe_c, _Unicode(end2x), 0.0 );
  
  double length = sqrt((end1z-end2z)*(end1z-end2z)+(end1x-end2x)*(end1x-end2x));
  double yrot   = atan((end1x-end2x)/(end1z-end2z));

  // -----------------------------
  // Make beampipe shape
  ConeSegment tube_vacuum(length / 2.0, 0.0, innerD1 / 2.0, 0.0, innerD2 / 2.0);
  ConeSegment tube_tube(length / 2.0, innerD1 / 2.0, innerD1 / 2.0 + thickness, innerD2 / 2.0,
                        innerD2 / 2.0 + thickness);

  Volume v_vacuum("v_vacuum", tube_vacuum, m_Vacuum);
  Volume v_tube("v_tube", tube_tube, m_Al);

  sdet.setAttributes(det, v_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  assembly.placeVolume(v_vacuum, Position(0, 0, -length / 2.0));
  assembly.placeVolume(v_tube,   Position(0, 0, -length / 2.0));

  // -----------------------------
  // final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly, Transform3D(RotationY(yrot),Position(end1x, 0.0, end1z)));
  pv_assembly.addPhysVolID("system", sdet.id()).addPhysVolID("barrel", 1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(BackwardsBeamPipe, create_detector)
