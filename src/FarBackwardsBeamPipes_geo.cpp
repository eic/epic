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
  //Material   m_Vacuum = det.material("Vacuum");
  string     vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BeamPipeVis");

  //cout<<m_Al.Z()<<endl;

  double preX = 0;
  double preZ = 0;
  double preRot = 0;
  double firstX = 0, firstZ = 0, firstRot = 0;

  BooleanSolid full_tube;
  //BooleanSolid full_vacuum;

  for( xml_coll_t pipe_coll(x_det, _Unicode(pipe)); pipe_coll; pipe_coll++ ) { // pipes
   
    xml_comp_t pipe( pipe_coll );
    double x1 = getAttrOrDefault<double>(pipe, _Unicode(end1x), 0);
    double x2 = getAttrOrDefault<double>(pipe, _Unicode(end2x), 0);
    double z1 = getAttrOrDefault<double>(pipe, _Unicode(end1z), 0);
    double z2 = getAttrOrDefault<double>(pipe, _Unicode(end2z), 0);
    double r1 = getAttrOrDefault<double>(pipe, _Unicode(rout1), 0);
    double r2 = getAttrOrDefault<double>(pipe, _Unicode(rout2), 0);
    //cout<<"x1: "<<x1<<endl;
    //cout<<"x2: "<<x2<<endl;
    //cout<<"z1: "<<z1<<endl;
    //cout<<"z2: "<<z2<<endl;
    //cout<<setprecision(6)<<"r1: "<<r1<<endl;
    //cout<<setprecision(6)<<"r2: "<<r2<<endl;
  
    double thickness = getAttrOrDefault<double>(x_det, _Unicode(wall_thickness), 0);
    double length = sqrt( pow(z1 - z2, 2) + pow(x1 - x2, 2) );
    double yrot   = atan( (x1 - x2) / (z1 - z2) );
    cout<<"yrot: "<<yrot<<endl;

    //ConeSegment s_vacuum( length / 2.0, 0.0, r1, 0.0, r2 );
    ConeSegment s_tube( length / 2.0, r1, r1 + thickness, r2, r2 + thickness);

    if( ! full_tube ) {
       full_tube = IntersectionSolid( s_tube, s_tube );
       //full_vacuum = IntersectionSolid( s_vacuum, s_vacuum );
       firstX = (x1 + x2)/2.;
       firstZ = (z1 + z2)/2.;
       firstRot = yrot;
    } else {
      full_tube = UnionSolid( full_tube, s_tube, Transform3D( Translation3D( (x1+x2)/2. - preX, 0, (z1+z2)/2. - preZ ) * RotationY( preRot - yrot ) ) );
      //full_vacuum = UnionSolid( full_vacuum, s_vacuum, Transform3D( Translation3D( (x1+x2)/2., 0, (z1+z2)/2. ) * RotationY( yrot) ) );
    
      //Volume v_vacuum("v_vacuum", full_vacuum, m_Vacuum);
      //Volume v_tube("v_tube", full_tube, m_Al);

      //assembly.placeVolume(v_vacuum, Transform3D( RotationZYX(0.0,0.0,0.0), Position(0, 0, -30*m)));
      //assembly.placeVolume(v_tube, Transform3D( RotationZYX(0.0,0.0,0.0), Position(0, 0, -30*m)));
    }
    
    preX = (x1+x2)/2.;
    preZ = (z1+z2)/2.;
    preRot = yrot;
    cout<<preX<<"  "<<preZ<<"  "<<preRot<<endl;
    

  }
 

  //Volume v_vacuum("v_vacuum", full_vacuum, m_Vacuum);
  Volume v_tube("v_tube", full_tube, m_Al);

  sdet.setAttributes(det, v_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  //assembly.placeVolume(v_vacuum, Transform3D( RotationZYX(0.0,0.0,0.0), Position(0, 0, -30*m)));
  assembly.placeVolume(v_tube, Transform3D( RotationY(firstRot), Position(firstX, 0, firstZ)));

  // -----------------------------
  // final placement
  auto pv_assembly =
      det.pickMotherVolume(sdet).placeVolume( assembly, Position(0.0, 0.0, 0.0));
  
  sdet.setPlacement(pv_assembly);
  
  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(FarBackwardsBeamPipes, create_detector)
