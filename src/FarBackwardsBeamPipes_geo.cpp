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

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Al     = description.material("Aluminum");
  Material   m_Vacuum = description.material("Vacuum");
  string     vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");

  //cout<<m_Al.Z()<<endl;

  //double preX = 0;
  //double preZ = 0;
  //double preRot = 0;
  //double firstX = 0, firstZ = 0, firstRot = 0;

  //BooleanSolid full_tube;
  //BooleanSolid full_vacuum;

  for( xml_coll_t pipe_coll(x_det, _Unicode(pipe)); pipe_coll; pipe_coll++ ) { // pipes
   
    xml_comp_t pipe( pipe_coll );
    
    double x = getAttrOrDefault<double>(pipe, _Unicode(xcenter), 0);
    double z = getAttrOrDefault<double>(pipe, _Unicode(zcenter), 0);
    double length = getAttrOrDefault<double>(pipe, _Unicode(length), 0);
    double yrot = getAttrOrDefault<double>(pipe, _Unicode(theta), 0);
    double r1 = getAttrOrDefault<double>(pipe, _Unicode(rout1), 0);
    double r2 = getAttrOrDefault<double>(pipe, _Unicode(rout2), 0);
    double thickness = getAttrOrDefault<double>(x_det, _Unicode(wall_thickness), 0);

    ConeSegment s_tube( length / 2.0, r1 - thickness, r1, r2 - thickness, r2 );
    ConeSegment s_vacuum( length / 2.0, 0, r1 - thickness, 0, r2 - thickness );
    Volume v_tube("v_tube", s_tube, m_Al);
    Volume v_vacuum("v_vacuum", s_vacuum, m_Vacuum);
    
    //sdet.setAttributes(det, v_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);
    v_tube.setVisAttributes(description.visAttributes( vis_name ) );

    assembly.placeVolume(v_tube, Transform3D( RotationY(yrot), Position(x, 0, z)));
    assembly.placeVolume(v_vacuum, Transform3D( RotationY(yrot), Position(x, 0, z)));

    //if( ! full_tube ) {
    //   full_tube = IntersectionSolid( s_tube, s_tube );
    //   //full_vacuum = IntersectionSolid( s_vacuum, s_vacuum );
    //   //firstX = (x1 + x2)/2.;
    //   //firstZ = (z1 + z2)/2.;
    //   //firstX = x;
    //   //firstZ = z;
    //   //firstRot = yrot;
    //} else {
    //  //full_tube = UnionSolid( full_tube, s_tube, Transform3D( Translation3D( (x1+x2)/2. - preX, 0, (z1+z2)/2. - preZ ) * RotationY( preRot - yrot ) ) );
    //  full_tube = UnionSolid( full_tube, s_tube, Transform3D( RotationY( preRot - yrot ), Translation3D( x - preX, 0, z - preZ ) ) );
    //  //full_vacuum = UnionSolid( full_vacuum, s_vacuum, Transform3D( Translation3D( (x1+x2)/2., 0, (z1+z2)/2. ) * RotationY( yrot) ) );
    //
    //  //Volume v_vacuum("v_vacuum", full_vacuum, m_Vacuum);
    //  //Volume v_tube("v_tube", full_tube, m_Al);

    //  //assembly.placeVolume(v_vacuum, Transform3D( RotationZYX(0.0,0.0,0.0), Position(0, 0, -30*m)));
    //  //assembly.placeVolume(v_tube, Transform3D( RotationZYX(0.0,0.0,0.0), Position(0, 0, -30*m)));
    //}
    //
    ////preX = (x1+x2)/2.;
    ////preZ = (z1+z2)/2.;
    //preX = x;
    //preZ = z;
    //preRot = yrot;
    //cout<<preX<<"  "<<preZ<<"  "<<preRot<<endl;
    

  }
 

  //Volume v_vacuum("v_vacuum", full_vacuum, m_Vacuum);
  //Volume v_tube("v_tube", full_tube, m_Al);

  //sdet.setAttributes(det, v_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  //assembly.placeVolume(v_tube, Transform3D( RotationY(firstRot), Position(firstX, 0, firstZ)));

  // -----------------------------
  // final placement
  auto pv_assembly =
      description.pickMotherVolume(sdet).placeVolume( assembly, Position(0.0, 0.0, 0.0));
  
  sdet.setPlacement(pv_assembly);
  
  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(FarBackwardsBeamPipes, create_detector)
