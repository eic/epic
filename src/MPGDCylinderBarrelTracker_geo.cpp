// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author Nivedith Ramasubramanian
 *
 * \ingroup trackers
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>
#include "DD4hepDetectorHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Micromegas Barrel Tracker with space frame
 *
 * - Designed to process "mpgd_barrel.xml" ("mpgd_barrel_ver1" as of 2024/02).
 *
 * - Derived from "BarrelTrackerWithFrame_geo.cpp".
 *
 * - "support" tag not addressed.
 *
 * - "frame" tag within the module element.
 *
 * \code
 * \endcode
 *
 *
 * @author Whitney Armstrong
 */
static Ref_t create_MPGDCylinderBarrelTracker(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     air      = description.air();
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  DetElement                   sdet(det_name, det_id);

  map<string, Volume>                volumes;
  map<string, Placements>            sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;

  PlacedVolume pv;

  //#define DEBUG_MPGDCylinderBarrelTracker
#ifdef DEBUG_MPGDCylinderBarrelTracker
  // TEMPORARILY INCREASE VERBOSITY level for debugging pruposes
  PrintLevel priorPrintLevel = printLevel();
  setPrintLevel(DEBUG);
#endif

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                         "boundary_material");
  }

  dd4hep::xml::Dimension dimensions(x_det.dimensions());
  Assembly assembly(det_name);

  sens.setType("tracker");
 
  // loop over the modules
  double total_zlength=0;
  int Nmod = 0;
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    Nmod++;
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    printout(DEBUG, "MPGDCylinderBarrelTracker","New module #%d \"%s\"",
	     Nmod,m_nam.c_str());

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR,"MPGDCylinderBarrelTracker","Module \"%s\" already exists.",
	       m_nam.c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int    ncomponents =0;
    int    sensor_number=1;
    double thickness_so_far=0;
    double total_thickness=0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      xml_comp_t x_comp = ci;
      printout(DEBUG,"MPGDCylinderBarrelTracker",
	       "0x%x: \"%s\" \t total_thickness %.4f cm",
	       ci,x_comp.nameStr().c_str(),total_thickness/cm);
      total_thickness += x_comp.thickness();
    }
    printout(DEBUG,"MPGDCylinderBarrelTracker",
	     " => total_thickness %.4f cm",
	     total_thickness/cm);
    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    double det_rmin   = dimensions.rmin();
    double det_length = dimensions.length();
    int nmod = 8; // the value of 2piR/width is 6.9, want to round up to 8. So +2
    double phi_start = - M_PI/nmod - (1*3.14/180); //causing 1 degrees of overlap on each side
    double phi_end =     M_PI/nmod + (1*3.14/180) ;
    double zthickness_frame;

    // Module frame.
    if (x_mod.hasChild(_U(frame))) {
      xml_comp_t m_frame = x_mod.child(_U(frame));
      const string frame_nam  = _toString("frame");
      zthickness_frame=m_frame.width();
      total_zlength = det_length+2*zthickness_frame;
      
      Tube frame_tube_1(det_rmin, det_rmin + total_thickness,zthickness_frame/2 , phi_start, phi_end);
      Volume frame_vol_1(frame_nam, frame_tube_1, description.material(m_frame.materialStr()));
      pv=m_vol.placeVolume(frame_vol_1, Position(0, 0, -(det_length+zthickness_frame)/2));
      pv=m_vol.placeVolume(frame_vol_1, Position(0, 0,  (det_length+zthickness_frame)/2));

      m_vol.setVisAttributes(description, m_frame.visStr());
      double frame_phi_end=zthickness_frame/(det_rmin); //converting the thickness of the frame to angular radians.

      Tube frame_tube_2(det_rmin, det_rmin + total_thickness,(det_length+2*zthickness_frame)/2 , phi_start-frame_phi_end, phi_start);
      Volume frame_vol_2(frame_nam, frame_tube_2, description.material(m_frame.materialStr()));
      pv=m_vol.placeVolume(frame_vol_2);

      Tube frame_tube_3(det_rmin, det_rmin + total_thickness,(det_length+2*zthickness_frame)/2 , phi_end, phi_end+frame_phi_end);
      Volume frame_vol_3(frame_nam, frame_tube_3, description.material(m_frame.materialStr()));
      pv=m_vol.placeVolume(frame_vol_3);

      frame_vol_1.setVisAttributes(description, m_frame.visStr());
      frame_vol_2.setVisAttributes(description, m_frame.visStr());
      frame_vol_3.setVisAttributes(description, m_frame.visStr());
      
   }else{
      total_zlength = det_length;
   }//end of if on module frame
 
    double comp_rmin=det_rmin;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t   x_comp = mci;
      const string c_nam  = x_comp.nameStr();
      
      double comp_thickness = x_comp.thickness();
      Tube         c_tube(comp_rmin, comp_rmin + comp_thickness, det_length/2, phi_start, phi_end);
      Volume       c_vol(c_nam, c_tube, description.material(x_comp.materialStr()));
      PlacedVolume c_pv;
      //Use this if you need to bring the module to the center and the do the rotation and translation.
      //pv=m_vol.placeVolume(c_vol, Position(-(2*(comp_rmin + comp_thickness/2) - det_rmin - total_thickness/2), 0, 0));
     
      pv=m_vol.placeVolume(c_vol);
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        pv.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                     total_thickness - thickness_so_far - x_comp.thickness() / 2.0};

	double a = thickness_so_far + x_comp.thickness() / 2.0;
        double b1 = total_thickness - thickness_so_far - x_comp.thickness() / 2.0;
        double b2 = total_thickness - (thickness_so_far + x_comp.thickness()) - 0.05;
	printout(DEBUG,"MPGDCylinderBarrelTracker",
		 "Module \"%s\" thickness a,b1,b2 %.4f,%.4f,%.4f cm",
		 m_nam.c_str(),a/cm,b1/cm,b2/cm);

        // -------- create a measurement plane for the tracking surface attached to the sensitive volume -----
        Vector3D u(-1., 0., 0.);
        Vector3D v(0., -1., 0.);
        Vector3D n(0., 0., 1.);
      
        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

	SurfaceType type(SurfaceType::Sensitive);

        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[m_nam].push_back(surf);
      }
      comp_rmin += comp_thickness;
      thickness_so_far += x_comp.thickness();
      
    }//end of module component loop
  }//end of module loop

  // now build the layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int        lay_id   = x_layer.id();
    string     m_nam    = x_layer.moduleStr();
    string     lay_nam  = det_name + _toString(x_layer.id(), "_layer%d");
    //Tube       lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length()+2.5*cm); //2.5cm?
    //Tube       lay_tub(x_barrel.inner_r()-0.5*cm, x_barrel.outer_r()+0.5*cm, x_barrel.z_length()+0*cm);
    Tube       lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length()/2);
    Volume     lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position   lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    
    std::cout << "x_barrel.z_length() = " << x_barrel.z_length() << std::endl;
    double det_zmin = -x_barrel.zmin();
    double det_zmax = x_barrel.zmax();
       
    double phi0     = x_layout.phi0();     // Starting phi of first module.
    //double rc       = x_layout.rc();       // Radius of the module center. Use this if you moved the module to the center in the module component loop
    double rc       = 1;       
    int    nphi     = x_layout.nphi();     // Number of modules in phi.
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double phic     = phi0;                // Phi of the module center.
    double nz       = z_layout.nz();       // Number of modules to place in z.

    Volume      module_env = volumes[m_nam];
    DetElement  lay_elt(sdet, lay_nam, lay_id);
    Placements& sensVols = sensitives[m_nam];

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto &layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }

    double a = (det_zmax+det_zmin)/2;
    double z_off = 0.1*cm; //half the offset between the inner middle two modules.
    double modz_pos[4]={det_zmin+(total_zlength)/2 - a,-(total_zlength+z_off)/2,(total_zlength+z_off)/2,det_zmax-(total_zlength)/2-a}; //these are the 4 central values in z where the four modules will be placed
    //  double z_off = 4*cm; //half the offset between the inner middle two modules.
    //double modz_pos[4]={det_zmin+(total_zlength)/2 - a,det_zmin+total_zlength-z_off+(total_zlength)/2-a,det_zmax-total_zlength+z_off-(total_zlength)/2-a,det_zmax-(total_zlength)/2-a}; //these are the 4 central values in z where the four modules will be placed

    int    module   = 1;
    // Loop over the number of modules in phi.
    for (int ii = 0; ii < nphi; ii++) {
      double y1,x1;
      //loop over modules in z
      for (int j = 0; j < nz; j++) {
	if(j==1||j==2){
	  if(ii%2==0)
	  {
	    y1  = rc * std::cos(phic);              // Basic x module position.
	    x1  = rc * std::sin(phic);              // Basic y module position.
	  }
	else{
	  y1  = (rc+0.5) * std::cos(phic);              // Basic x module position.
	  x1  = (rc+0.5) * std::sin(phic);              // Basic y module position.
	}}
	else{
	   if(ii%2==0)
	  {
	    y1  = (rc+1) * std::cos(phic);              // Basic x module position.
	    x1  = (rc+1) * std::sin(phic);              // Basic y module position.
	  }
	else{
	  y1  = (rc+1.5) * std::cos(phic);              // Basic x module position.
	  x1  = (rc+1.5) * std::sin(phic);              // Basic y module position.
	}}

	string     module_name = _toString(module, "module%d");
        DetElement mod_elt(lay_elt, module_name, module);

	printout(DEBUG,"MPGDCylinderBarrelTracker",
		 "Layer \"%s\" sector %d,%d,\"%s\": \t x1,y1,r1: %7.4f,%7.4f,%7.4f cm",
		 lay_nam.c_str(),ii,j,module_name.c_str(),x1/cm,y1/cm,sqrt(x1*x1+y1*y1)/cm);

	Transform3D tr(RotationZYX((M_PI/2-phic),0,0),Position(x1,y1,modz_pos[j]));
	pv = lay_vol.placeVolume(module_env,tr);
	pv.addPhysVolID("module", module);
	mod_elt.setPlacement(pv);

        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement   comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
          comp_de.setPlacement(sens_pv);

          auto &comp_de_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
          comp_de_params.set<string>("axis_definitions", "XYZ");
	  volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
        }
	module++;
      }
      phic += phi_incr; // Increment the phi placement of module.
      rc += rphi_dr;    // Increment the center radius according to dr parameter.
      rphi_dr *= -1;    // Flip sign of dr parameter.
    }

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }
   
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);

#ifdef DEBUG_MPGDCylinderBarrelTracker
  // Reset initial print level before exiting
  setPrintLevel(priorPrintLevel);
#endif

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(MPGDCylinderBarrelTracker, create_MPGDCylinderBarrelTracker)
