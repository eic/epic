/** \addtogroup PID
 * \brief Type: **Barrel bar detector (rectangular geom.) with frames surrounding the perimeter.
 * \author S. Joosten
 * Modified by M. Posik
 *
 * \ingroup PID 
 * \ingroup Tracking
 *
 *
 * \code
 * \endcode
 *
 * @{
*/
  
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include <array>

using namespace std;
using namespace dd4hep;

/** Barrel Bar detector with optional frame
 *
 * - Optional "frame" tag within the module element (frame eats part of the module width and length
 *   and surrounds the rectangular perimeter) 
 * - Detector is setup as a "tracker" so we can use the hits
 *
 */
static Ref_t create_MPGDDIRC_geo(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume>            Placements;
  xml_det_t                               x_det    = e;
//  Material                                air      = description.air();
  int                                     det_id   = x_det.id();
  string                                  det_name = x_det.nameStr();
  DetElement                              sdet(det_name, det_id);
  map<string, Volume>                     volumes;
  map<string, Placements>                 sensitives;
  map<string, std::vector<rec::VolPlane>> volplane_surfaces;
  PlacedVolume                            pv;
  dd4hep::xml::Dimension                  dimensions(x_det.dimensions());
  xml_dim_t                               mpgd_dirc_pos = x_det.position();
  Assembly 				  assembly(det_name);

  map<string, std::array<double, 2>> module_thicknesses;
  sens.setType("tracker");

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "MPGDDIRC_geo",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int    ncomponents     = 0;
    double total_thickness = 0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }
    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description, x_mod.visStr());

    // Optional module frame.
    // frame is 4  bars around the perimeter of the rectangular module. The frame will eat the
    // overlapping module area
    //
    //  ___
    // |___| <-- example module cross section (x-y plane), frame is flush with the
    //           bottom of the module and protrudes on the top if needed
    
    // Get frame width, as it impacts the main module for being built. We
    // construct the actual frame structure later (once we know the module width)
    double frame_width = 0;
    if (x_mod.hasChild(_U(frame))) {
      xml_comp_t m_frame = x_mod.child(_U(frame));
      frame_width        = m_frame.width();
    }

    double thickness_so_far    = 0.0;
    double thickness_sum       = -total_thickness / 2.0;
    double max_component_width = 0;
    double max_component_length = 0;
    double gas_thickness = 0.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t x_comp = mci;
      string     c_nam  = _toString(ncomponents, "component%d");
      string  comp_name = x_comp.nameStr();

      double box_width = x_comp.width();
      double box_length = x_comp.length();
      Box c_box;
      //Since MPGD frames are layed over the MPGD foils, the foil material is pressent under the frame as well. 
      //The gas volumes are not present under the frames, so our frames must eat only the gas module areas
      //
      //  ------------------- MPGD foil 
      //  --               -- Frame
      //  --  gas volume   -- Frame
      //  --               -- Frame
      //  ------------------- MPGD foil

      //Look for gas modules to subtract frame thickness from
      //FIXME: these module names are hard coded for now. Should find 
      //a way to set a arribut via the moduel tag to flag what components 
      //need to have frame thickness subtracted.
      if( (comp_name == "DriftGap" || comp_name == "WindowGasGap") ){
	box_width = x_comp.width() - 2.0 * frame_width;
	box_length = x_comp.length()- 2.0 * frame_width;
	max_component_width = box_width;
	max_component_length = box_length;
	gas_thickness += x_comp.thickness();
        c_box = {box_width / 2, box_length / 2, x_comp.thickness() / 2};
	//cout << "gas: " << comp_name << " \n";
	//cout << "box_width: " << box_width << endl;
	//cout << "box_length: " << box_length << endl;
	//cout << "box_thickness: " << x_comp.thickness() << endl;
      }else{
        c_box = {x_comp.width() / 2, x_comp.length() / 2, x_comp.thickness() / 2};
	//cout << "Not gas: " << comp_name << " \n";
	//cout << "box_comp_width: " << x_comp.width() << endl;
	//cout << "box_comp_length: " << x_comp.length() << endl;
	//cout << "box_comp_thickness: " << x_comp.thickness() << endl;
      }
      Volume c_vol{c_nam, c_box, description.material(x_comp.materialStr())};
      pv = m_vol.placeVolume(c_vol, Position(0, 0, thickness_sum + x_comp.thickness() / 2.0));

      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                     total_thickness - thickness_so_far - x_comp.thickness() / 2.0};
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        rec::Vector3D u(0., 1., 0.);
        rec::Vector3D v(0., 0., 1.);
        rec::Vector3D n(1., 0., 0.);

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        rec::SurfaceType type(rec::SurfaceType::Sensitive);

        rec::VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[m_nam].push_back(surf);
      }
      thickness_sum += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
    }
    // Now add-on the frame
    if (x_mod.hasChild(_U(frame))) {
      xml_comp_t m_frame         = x_mod.child(_U(frame));
      double     frame_thickness = getAttrOrDefault<double>(m_frame, _U(thickness), total_thickness);

        Box lframe_box{m_frame.width() / 2.0, (max_component_length + 2.0*m_frame.width()) / 2.0, frame_thickness / 2.0};
        Box rframe_box{m_frame.width() / 2.0, (max_component_length + 2.0*m_frame.width()) / 2.0, frame_thickness / 2.0};
        Box tframe_box{max_component_width / 2.0, m_frame.width() / 2.0, frame_thickness / 2.0};
        Box bframe_box{max_component_width / 2.0, m_frame.width() / 2.0, frame_thickness / 2.0};
      
      // Keep track of frame with so we can adjust the module bars appropriately

        Volume lframe_vol{"left_frame", lframe_box, description.material(m_frame.materialStr())};
        Volume rframe_vol{"right_frame", rframe_box, description.material(m_frame.materialStr())};
        Volume tframe_vol{"top_frame",    tframe_box, description.material(m_frame.materialStr())};
        Volume bframe_vol{"bottom_frame", bframe_box, description.material(m_frame.materialStr())};

        lframe_vol.setVisAttributes(description, m_frame.visStr());
        rframe_vol.setVisAttributes(description, m_frame.visStr());
        tframe_vol.setVisAttributes(description, m_frame.visStr());
        bframe_vol.setVisAttributes(description, m_frame.visStr());

        cout << "frame_thickness: " << frame_thickness << endl;
        cout << "total_thickness: " << total_thickness << endl;
        cout << "frame_thickness + total_thickness: " << frame_thickness + total_thickness << endl;

        m_vol.placeVolume(lframe_vol, Position(frame_width / 2.0 + max_component_width / 2, 0.0,
                                               frame_thickness / 2.0 - total_thickness / 2.0 - gas_thickness/ 2.0));
        m_vol.placeVolume(rframe_vol, Position(-frame_width / 2.0 - max_component_width / 2.0, 0.0,
                                                frame_thickness / 2.0 - total_thickness / 2.0 - gas_thickness/2.0));
        m_vol.placeVolume(tframe_vol, Position(0.0, frame_width / 2.0 + max_component_length/2,
                                               frame_thickness / 2.0 - total_thickness / 2.0 - gas_thickness/2.0));
        m_vol.placeVolume(bframe_vol, Position(0.0,-frame_width / 2.0 - max_component_length/2.0,
                                               frame_thickness / 2.0 - total_thickness / 2.0 - gas_thickness/2.0));
					       
    }
    pv = assembly.placeVolume(m_vol);
  }

  //build layers
  for(xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer = li;
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));

    double phi0	    = x_layout.phi0();      //starting phi of first module
    double phi_tilt = x_layout.phi_tilt();  //Phi tilit of module
    double rc 	    = x_layout.rc();  	    //Radius of the module
    int nphi        = x_layout.nphi();      //Number of modules in phi
    double rphi_dr  = x_layout.dr();        //The delta radius of every other module
    double phi_incr = (2 * M_PI) / nphi;    //Phi increment for one module
    double phic = phi0;                     //Phi of the module

    //loop over phi modules
    for(int ii= 0; ii < nphi; ii++) {
      double xc =  rc * std::cos(phic);
      double yc =  rc * std::sin(phic);
      //place P-side moduels
      Transform3D tr(RotationZYX(0.0, ((M_PI/2) - phic -phi_tilt), -M_PI /2 ),  Position(xc, yc, mpgd_dirc_pos.z() )); //in x-y plane,
      sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
      //assembly.setVisAttributes(description.invisible());
      pv = description.pickMotherVolume(sdet).placeVolume(assembly, tr);
      pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
      //place N-side modules
      Transform3D tr2(RotationZYX(0.0, ((M_PI /2 ) - phic - phi_tilt), -M_PI /2),  Position(xc, yc,mpgd_dirc_pos.z() - dimensions.length())); 
      pv = description.pickMotherVolume(sdet).placeVolume(assembly, tr2);
      pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
      //increment counters
      phic += phi_incr;
      rc += rphi_dr;
    }
  }  
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_MPGDDIRC, create_MPGDDIRC_geo)
