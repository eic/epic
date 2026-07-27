// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
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
#include "DD4hepDetectorHelper.h"
#include <array>
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Endcap Trapezoidal Tracker.
 *
 * @author Whitney Armstrong
 *
 */
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material vacuum = description.vacuum();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  bool reflect    = x_det.reflect(false);
  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);

  Material air     = description.material("Air");
  Volume motherVol = description.pickMotherVolume(sdet);
  int m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume> modules;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume pv;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");


  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();
   // xml_comp_t tube   = x_mod.tubs();

    double posY;
    double rmin =0;
    double rmax =0;
    double z =0;
    double dx =0;
    double total_thickness = 0.;
    double thickness_so_far = 0.0;
    /*
    double rmin              = tube.rmin();
    double rmax              = tube.rmax();
    double z                 = tube.z();
    double total_thickness = 0.;
    double thickness_so_far = 0.0;

    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
    total_thickness += xml_comp_t(ci).thickness();

    */
    xml_coll_t ci(x_mod, _U(module_component));

    xml::Element mShape = x_mod.child("shape");

    if(mShape){
    	rmin = mShape.attr<double>(_U(rmin));
    	rmax = mShape.attr<double>(_U(rmax));
    	z    = mShape.attr<double>(_U(z));
    	dx   = mShape.attr<double>("d1");
    	total_thickness = z;
    }


    double y1               = total_thickness / 2;


    Tube m_solid_temp1(rmin, rmax, z/2.0, 0.0, M_PI/2.0);

    Box m_box1_temp(dx/2, (rmax-rmin)/2,z/2.0);

    Transform3D box1Pos(Rotation3D(), Position(-dx/2, rmin+(rmax-rmin)/2, 0.0));
    UnionSolid combinedSolid(m_solid_temp1, m_box1_temp, box1Pos);

    Box m_box2_temp((rmax-rmin)/2,dx/2,z/2.0);
    Transform3D box2Pos(Rotation3D(), Position(rmin+(rmax-rmin)/2, -dx/2, 0.0));
    UnionSolid m_solid(combinedSolid, m_box2_temp, box2Pos);

    Volume m_volume(m_nam, m_solid, vacuum);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
        xml_comp_t c = ci;
        string comp_name  = c.nameStr();
        xml::Element cShape = c.child("component_shape");

        double c_thick = c.thickness();
        double comp_rmin = 0;
        double comp_rmax = 0;
        double comp_d1 =0;

        if(cShape){
        	comp_rmin = cShape.attr<double>(_U(rmin));
        	comp_rmax = cShape.attr<double>(_U(rmax));
			comp_d1   = cShape.attr<double>("d1");

			cout<<comp_d1<<endl;
        }

        //auto comp_rmin = rmin;
        //auto comp_rmax = rmax;

        Material c_mat = description.material(c.materialStr());
        string c_name = _toString(c_id, "component%d");

        if (!c.isSensitive()) {

        	Solid comp_shape;
        	if (comp_name.find("Frame") != std::string::npos) {

        		c_name = _toString(c_id, "frame%d") ;
        		cout<<c_name<<endl;
        	//	for(int i=0; i<4; i++){
        		double comp_d2 = cShape.attr<double>("d2");
        		double comp_rmin2 = cShape.attr<double>("rmin2");
        		double comp_rmax2 = cShape.attr<double>("rmax2");

        		Tube comp_t1(comp_rmin, comp_rmin2, c_thick/2.0, 0.0, M_PI/2.0);

        		Tube comp_t2(comp_rmax2, comp_rmax, c_thick/2.0, 0.0, M_PI/2.0);

        		Transform3D CTubePos(Rotation3D(), Position(0.0,0.0, 0.0));
        		UnionSolid comp_tube(comp_t1, comp_t2, CTubePos);

    			Box m_Cbox1_temp(comp_d1/2, (comp_rmax-comp_rmin)/2,c_thick/2.0);
    			Box m_Cbox2_temp((comp_d1-comp_d2)/2, (comp_rmax2-comp_rmin2)/2,c_thick/2.0);
        		Transform3D Cbox2Pos(Rotation3D(), Position((comp_d2)/2,0,0));
        		SubtractionSolid comp_boxL(m_Cbox1_temp, m_Cbox2_temp, Cbox2Pos);

        		Transform3D CBoxLPos(Rotation3D(), Position(-comp_d2/2, comp_rmin+(comp_rmax-comp_rmin)/2, 0.0));
        		UnionSolid  comp_s1_temp(comp_tube,comp_boxL,CBoxLPos );

    			Box m_Cbox1_temp2((comp_rmax-comp_rmin)/2,comp_d2/2,c_thick/2.0);
    			Box m_Cbox2_temp2((comp_rmax2-comp_rmin2)/2,(comp_d1-comp_d2)/2,c_thick/2.0);
    			Transform3D Cbox2Pos2(Rotation3D(), Position(0,(comp_d2)/2,0));
    			SubtractionSolid comp_s2_temp(m_Cbox1_temp2, m_Cbox2_temp2, Cbox2Pos2);

        		Transform3D CBoxRPos(Rotation3D(), Position(comp_rmin+(comp_rmax-comp_rmin)/2,-comp_d2/2, 0.0));
        		UnionSolid  comp_s1(comp_s1_temp,comp_s2_temp,CBoxRPos );
        		comp_shape = comp_s1;

        		Volume c_vol(c_name, comp_shape, c_mat);



        		pv = m_volume.placeVolume(c_vol, Position(0, 0, posY + c_thick / 2));



        	}else{
        		Tube comp_s1_temp1(comp_rmin, comp_rmax, c_thick/2.0, 0.0, M_PI/2.0);

        		Box m_Cbox1_temp(comp_d1/2, (comp_rmax-comp_rmin)/2,c_thick/2.0);
        		Transform3D Cbox1Pos(Rotation3D(), Position(-comp_d1/2, comp_rmin+(comp_rmax-comp_rmin)/2, 0.0));
        		UnionSolid CcombinedSolid(comp_s1_temp1, m_Cbox1_temp, Cbox1Pos);

        		Box m_Cbox2_temp((comp_rmax-comp_rmin)/2,comp_d1/2,c_thick/2.0);
        		Transform3D Cbox2Pos(Rotation3D(), Position(comp_rmin+(comp_rmax-comp_rmin)/2, -comp_d1/2, 0.0));
        		UnionSolid comp_s1(CcombinedSolid, m_Cbox2_temp, Cbox2Pos);

        		 comp_shape = comp_s1;

        		Volume c_vol(c_name, comp_shape, c_mat);


        		c_vol.setVisAttributes(description.visAttributes(c.visStr()));
        		pv = m_volume.placeVolume(c_vol, Position(0, 0, posY + c_thick / 2));
        	}
        }else{
/*
        // Gestione del volume sensibile
        if (c.isSensitive()) {
        */

        	for(int i=0; i<3; i++){
        		Solid comp_shape;
        	  double x_pos=0;
        	  double y_pos=0;
        	  double z_pos=0;
        	  double thetax=0;
        	  double thetay=0;
        	  double thetaz=0;

        	  string cname_s = c_name + "_" + std::to_string(i);
        	  cout<<c_name<<endl;
        		if(i==0){
                cout<<comp_rmin<<" "<<comp_rmax<<endl;
        		Tube comp_s1_temp1(comp_rmin, comp_rmax, c_thick/2.0, 0.0, M_PI/2.0);


        	    x_pos = 0;
        	    y_pos = 0;
        	    z_pos = posY + c_thick / 2;

          	     thetax=0;
          	     thetay=0;
          	     thetaz=0;
          	  comp_shape = comp_s1_temp1;

        	}

        		if(i==1){
        		Box comp_s1_temp1(comp_d1/2, (comp_rmax-comp_rmin)/2,c_thick/2.0);

        	    x_pos = -comp_d1/2;
        	    y_pos = comp_rmin+(comp_rmax-comp_rmin)/2;
        	    z_pos = posY + c_thick / 2;


          	     thetax=0;
          	     thetay=0;
          	     thetaz=0;
          	  comp_shape = comp_s1_temp1;

        	}

        		if(i==2){
        		Box comp_s1_temp1((comp_rmax-comp_rmin)/2,comp_d1/2,c_thick/2.0);

        	    x_pos = comp_rmin+(comp_rmax-comp_rmin)/2;
        	    y_pos = -comp_d1/2;
        	    z_pos = posY + c_thick / 2;

          	     thetax=0;
          	     thetay=0;
          	     thetaz=0;
          	  comp_shape = comp_s1_temp1;

        	}


        	  Volume c_vol(cname_s, comp_shape, c_mat);
              c_vol.setVisAttributes(description.visAttributes(c.visStr()));
        	  pv = m_volume.placeVolume(c_vol, Transform3D(RotationZYX(thetax, thetay,thetaz),Position(x_pos, y_pos, z_pos)));

            module_thicknesses[m_nam] = {thickness_so_far + c_thick / 2.0,
                                         total_thickness - thickness_so_far - c_thick / 2.0};
          //  sdet.check(n_sensor > 2, "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");

            pv.addPhysVolID("sensor", n_sensor);
            c_vol.setSensitiveDetector(sens);
            sensitives[m_nam].push_back(pv);
            ++n_sensor;

            // Creazione di un piano di misura per la superficie di tracking associata al volume sensibile
            Vector3D u(0., 0., -1.);
            Vector3D v(-1., 0., 0.);
            Vector3D n(0., 1., 0.);

            // Calcolo degli spessori interni ed esterni da assegnare alla superficie di tracking
            double inner_thickness = module_thicknesses[m_nam][0];
            double outer_thickness = module_thicknesses[m_nam][1];

            SurfaceType type(SurfaceType::Sensitive);
            VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); // Associa la superficie di misura
            volplane_surfaces[m_nam].push_back(surf);
        }
        }

        if (comp_name.find("frame") == std::string::npos){
        posY += c_thick;
        thickness_so_far += c_thick;
        }
    }

    modules[m_nam] = m_volume;

    posY=0;

  }

  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int l_id    = x_layer.id();
    int mod_num = 1;

    xml_comp_t l_env  = x_layer.child(_U(envelope));
    string layer_name = det_name + std::string("_layer") + std::to_string(l_id);

    std::string layer_vis = l_env.attr<std::string>(_Unicode(vis));
    double layer_rmin     = l_env.attr<double>(_Unicode(rmin));
    double layer_rmax     = l_env.attr<double>(_Unicode(rmax));
    double layer_length   = l_env.attr<double>(_Unicode(length));
    double layer_zstart   = l_env.attr<double>(_Unicode(zstart));
    double layer_center_z = layer_zstart + layer_length / 2.0;
    // printout(INFO,"ROOTGDMLParse","+++ Read geometry from GDML file file:%s",input.c_str());


   //std::cout << "SiTracker Endcap layer " << l_id << " zstart = " << layer_zstart/dd4hep::mm << "mm ( " <<
    		//layer_rmax/dd4hep::mm << " mm thick )\n";

    // Assembly    layer_assembly(layer_name);
    // assembly.placeVolume(layer_assembly);



    Tube layer_tub(layer_rmin, layer_rmax, layer_length / 2);
    Volume layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    PlacedVolume layer_pv;
    if (reflect) {
     layer_pv = assembly.placeVolume(layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));

     //   layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, -layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_N";
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_P";
    }
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);



    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layer_element);

    layerParams.set<double>("envelope_r_min", 1 / dd4hep::mm);
    layerParams.set<double>("envelope_r_max", 1 / dd4hep::mm);
    layerParams.set<double>("envelope_z_min", 1 / dd4hep::mm);
    layerParams.set<double>("envelope_z_max", 1 / dd4hep::mm);


    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }



    for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
      xml_comp_t x_ring    = ri;
      double r             = x_ring.r();
      double phi0          = x_ring.phi0(0);
      double zstart        = x_ring.zstart();
      double dz            = x_ring.dz(0);
      int nmodules         = x_ring.nmodules();
      string m_nam         = x_ring.moduleStr();
      Volume m_vol         = modules[m_nam];
      double iphi          = 2 * M_PI / nmodules;
      double phi           = phi0;
      Placements& sensVols = sensitives[m_nam];

      for (int k = 0; k < nmodules; ++k) {
        string m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
        double x      = -r * std::cos(M_PI/4);
        double y      = -r * std::sin(M_PI/4);


        if (!reflect) {

          DetElement module(layer_element, m_base + "_pos", det_id);
          pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(M_PI / 2 + phi,0,0),
                                                        Position(x, y, zstart + dz)));
          pv.addPhysVolID("module", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        } else {

            // cout <<dz<<endl;

        	 double dz_shift = (k % 2 == 0) ? -dz : dz;

        	// cout<<k<<" "<<dz_shift<<" "<<dz<<endl;
        	// pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0,0,0), Position(0, 0, -zstart - dz-dz_shift)));
        	pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(M_PI / 2 + phi,0,0), Position(0, 0, dz_shift)));


        	//std::cout << "Rotation Matrix:\n" << RotationZYX(M_PI / 2 + phi, 0, 0) << std::endl;
        	// pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, -zstart - dz)));

          pv.addPhysVolID("module", mod_num);
          DetElement r_module(layer_element, m_base + "_neg", det_id);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XYZ");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        }

        phi += iphi;
        ++mod_num;
      }
    }


  }
  pv = motherVol.placeVolume(assembly, Position(0, 0, (reflect ? -1.0e-9 : 1.0e-9)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);



  return sdet;
}

DECLARE_DETELEMENT(epic_MPDGEndcapTracker, create_detector)
