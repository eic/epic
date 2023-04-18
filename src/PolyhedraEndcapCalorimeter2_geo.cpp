// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Modified for ATHENA detector
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

void buildTiles(Detector& desc, SensitiveDetector &sens, Volume &s_vol_outer, Volume &s_vol_middle, Volume &s_vol_inner, DetElement &slice, int det_id, xml_comp_t x_tiles, double sliceZ, vector<PlacedVolume> sens_volumes)
{
  //auto [s_dim_x, s_dim_y, s_dim_z] = dimensions;
  //double f_radius = s_pos;
  //double f_spacing_eta = getAttrOrDefault(x_tiles, _Unicode(spacing_eta), 0.1);
  double f_spacing_phi = getAttrOrDefault(x_tiles, _Unicode(phi_width), 0.05);
  //int f_nbins_eta = getAttrOrDefault(x_tiles, _Unicode(nbins_eta), 20);
  int f_nbins_phi = getAttrOrDefault(x_tiles, _Unicode(nbins_phi), 2);
  std::string f_id_tiles = getAttrOrDefault(x_tiles, _Unicode(identifier_tiles), "tiles");

  /*
  double f_trd_x1 = s_dim_x;
  double f_trd_x2 = s_dim_x;
  double f_trd_y1 = s_dim_y;
  double f_trd_y2 = f_trd_y1;
  double f_trd_z  = s_dim_z/2;
  */

  int f_num = 1;

  //double eta_lo = 0.0;
  //double eta_hi = f_spacing_eta;
  //eta_hi = f_spacing_eta;
  /*
  double eta_values[28] = { 3.9736842, 3.82929, 3.65542, 3.48947, 3.33106, 3.17987, 
            3.03556, 2.89782, 2.76635, 2.64086, 2.52109, 2.40676, 2.29764,
            2.19349, 2.09299, 1.99721, 1.90582, 1.81882, 1.7357, 1.65646, 
            1.58079, 1.50851, 1.43941, 1.3734, 1.3104, 1.25011, 1.19095, 1.1158023};
            */
 double eta_values[21] = { 3.9736842, 3.38219, 3.01399, 2.74652, 2.5368, 2.36466,
            2.21898, 2.09299, 1.99721, 1.90582, 1.81882, 1.7357, 1.65646, 
            1.58079, 1.50851, 1.43941, 1.3734, 1.3104, 1.25011, 1.19095, 1.1158023};        
  double z_slice = sliceZ;
  double eta_in;
  double eta_out;


  for (int i = 0; i < 1; ++i) {
  
  eta_in = -1.0*eta_values[i];
  eta_out = -1.0*eta_values[i+1];
  double R = z_slice*tan(2*atan(exp(-1*(eta_in))));
  double delR = z_slice*(tan(2*atan(exp(-1*(eta_out))))-tan(2*atan(exp(-1*(eta_in)))));

  double f_trd_x1 = 2*R*tan(4*f_spacing_phi/2.)/2.0;
  double f_trd_x2 = 2*(R+delR)*tan(4*f_spacing_phi/2.)/2.0;
  double f_trd_y = delR/cos(4*f_spacing_phi/2.)/2.0;
  //double f_trd_y2 = f_trd_y1;
  double f_trd_z = 0.2;

  ///cout << R << "   " << delR <<  "   " << eta_in << "   " << eta_out << endl;
  cout << 1 << "-" << 15 << "   " << delR << " " << "cm." << "  " << eta_in << "  " << eta_out <<  "     " << (R + (delR/2.0)) << " " << "cm." << "        "  << "(" << 2*f_trd_x1 << "," << 2*f_trd_x2 << "," << 2*f_trd_y << "," << 2*f_trd_y << ")" << " " << "cm." << endl ; 
  ///cout << "Height" << "  " << delR <<  "   " << "Inner eta edge" << "   " << eta_in << "   " <<  "Outer eta edge" << "   " << eta_out << endl;
  /////cout << "Tile no" << "  " << 1 << " " << "to" << " " << 15 << endl; 
  /////cout << "Height" << "  " << delR <<  "   " << "Inner eta edge" << "   " << eta_in << "   " <<  "Outer eta edge" << "   " << eta_out << endl;
  /////cout << " Edges" << "  " << 2*f_trd_x1 <<  "   " << 2*f_trd_x2 << "   " << 2*f_trd_y << "   " << 2*f_trd_y << endl;

	double phi_lo = -4*f_spacing_phi;
	double phi_hi = phi_lo + 4*f_spacing_phi;


	for (int j = 0; j < 15; ++j) {

		


		double xcent = (R + (delR/2.0))*cos(M_PI/2.0 - ((2*j)+1)*(4*f_spacing_phi/2.0));
		double ycent = (R + (delR/2.0))*sin(M_PI/2.0 - ((2*j)+1)*(4*f_spacing_phi/2.0));
    


        string     f_name  = Form("tile%d%s%d", i, " ", j);

		//EightPointSolid f_shape(f_trd_z, vert);


       // Box        f_shape(dx/2.0,dy/2.0,s_dim_z/2.0);
		Trapezoid  f_shape(f_trd_x1, f_trd_x2, f_trd_z, f_trd_z, f_trd_y);
		Volume     f_vol(f_name, f_shape, desc.material(x_tiles.materialStr()));
		DetElement tower(slice, f_name, det_id);

        
        f_vol.setAttributes(desc, x_tiles.regionStr(), x_tiles.limitsStr(), x_tiles.visStr());

        // Slice placement.

        
        PlacedVolume tower_phv = s_vol_inner.placeVolume(f_vol, Transform3D(RotationZYX(0, ((2*j)+1)*(4*f_spacing_phi/2.0), -M_PI / 2.0), Position(xcent, ycent, 0.0)));
        
        if ( x_tiles.isSensitive() ) {
          sens.setType("calorimeter");
          f_vol.setSensitiveDetector(sens);
          sens_volumes.push_back(tower_phv);
        }
        tower_phv.addPhysVolID("tile", f_num);
        tower.setPlacement(tower_phv);
        

		f_num++;

		phi_lo += 4*f_spacing_phi;
		phi_hi += 4*f_spacing_phi;
	}

  }

  
  

  
  //cout<< "eta_out" << eta_max << endl;
 
 
  for (int i = 1; i < 5; ++i)
   {

  eta_in = -1.0*eta_values[i];
  eta_out = -1.0*eta_values[i+1];
  double R = z_slice*tan(2*atan(exp(-1*(eta_in))));
  double delR = z_slice*(tan(2*atan(exp(-1*(eta_out))))-tan(2*atan(exp(-1*(eta_in)))));

  double f_trd_x1 = 2*R*tan(2*f_spacing_phi/2.)/2.0;
  double f_trd_x2 = 2*(R+delR)*tan(2*f_spacing_phi/2.)/2.0;
  double f_trd_y = delR/cos(2*f_spacing_phi/2.)/2.0;
  //double f_trd_y2 = f_trd_y1;
  double f_trd_z = 0.2;

  /////cout << R << "   " << delR <<  "   " << eta_in << "   " << eta_out << endl;
  /*cout << "Tile no" << "  " << 15+((30*(i-1))+1) << " " << "to" << " " << 15+((30*(i-1))+30) << endl; 
  cout << "Height" << "  " << delR <<  "   " << "Inner eta edge" << "   " << eta_in << "  " <<  "Outer eta edge" << "   " << eta_out << endl;
  cout << " Edges" << "  " << 2*f_trd_x1 <<  "   " << 2*f_trd_x2 << "   " << 2*f_trd_y << "   " << 2*f_trd_y << endl;*/
  
  cout << 15+((30*(i-1))+1) << "-" << 15+((30*(i-1))+30) << "   " << delR << " " << "cm." << "  " << eta_in << "  " << eta_out << "     " << (R + (delR/2.0)) << " " << "cm." <<  "        " <<  "(" << 2*f_trd_x1 << "," << 2*f_trd_x2 << "," << 2*f_trd_y << "," << 2*f_trd_y << ")" << " " << "cm." << endl ; 

	double phi_lo = -2*f_spacing_phi;
	double phi_hi = phi_lo + 2*f_spacing_phi;


	for (int j = 0; j < 30; ++j) {

		


		double xcent = (R + (delR/2.0))*cos(M_PI/2.0 - ((2*j)+1)*(2*f_spacing_phi/2.0));
		double ycent = (R + (delR/2.0))*sin(M_PI/2.0 - ((2*j)+1)*(2*f_spacing_phi/2.0));
    


        string  f_name  = Form("tile%d%s%d", i, " ", j);

		//EightPointSolid f_shape(f_trd_z, vert);


       // Box        f_shape(dx/2.0,dy/2.0,s_dim_z/2.0);
		Trapezoid  f_shape(f_trd_x1, f_trd_x2, f_trd_z, f_trd_z, f_trd_y);
		Volume     f_vol(f_name, f_shape, desc.material(x_tiles.materialStr()));
		DetElement tower(slice, f_name, det_id);

        
        f_vol.setAttributes(desc, x_tiles.regionStr(), x_tiles.limitsStr(), x_tiles.visStr());

        // Slice placement.

        
        PlacedVolume tower_phv = s_vol_middle.placeVolume(f_vol, Transform3D(RotationZYX(0, ((2*j)+1)*(2*f_spacing_phi/2.0), -M_PI / 2.0), Position(xcent, ycent, 0.0)));
        
        if ( x_tiles.isSensitive() ) {
          sens.setType("calorimeter");
          f_vol.setSensitiveDetector(sens);
          sens_volumes.push_back(tower_phv);
        }
        tower_phv.addPhysVolID("tile", f_num);
        tower.setPlacement(tower_phv);
        

		f_num++;

		phi_lo += 2*f_spacing_phi;
		phi_hi += 2*f_spacing_phi;
	}
   }

  for (int i = 5; i < 20; ++i) {

  eta_in = -1.0*eta_values[i];
  eta_out = -1.0*eta_values[i+1];
  double R = z_slice*tan(2*atan(exp(-1*(eta_in))));
  double delR = z_slice*(tan(2*atan(exp(-1*(eta_out))))-tan(2*atan(exp(-1*(eta_in)))));

  double f_trd_x1 = 2*R*tan(f_spacing_phi/2.)/2.0;
  double f_trd_x2 = 2*(R+delR)*tan(f_spacing_phi/2.)/2.0;
  double f_trd_y = delR/cos(f_spacing_phi/2.)/2.0;
  //double f_trd_y2 = f_trd_y1;
  double f_trd_z = 0.2;

  /////cout << R << "   " << delR <<  "   " << eta_in << "   " << eta_out << endl;

  /*cout << "Tile no" << "  " << 135+((f_nbins_phi*(i-5))+1) << " " << "to" << " " << 135+((f_nbins_phi*(i-5))+f_nbins_phi) << endl; 
  cout << "Height" << "  " << delR <<  "   " << "Inner eta edge" << "   " << eta_in << "  " << "Outer eta edge" << "   " << eta_out << endl;
  cout << " Edges" << "  " << 2*f_trd_x1 <<  "   " << 2*f_trd_x2 << "   " << 2*f_trd_y << "   " << 2*f_trd_y << endl;
 */

  cout << 135+((f_nbins_phi*(i-5))+1) << "-" << 135+((f_nbins_phi*(i-5))+f_nbins_phi) << "   " << delR << " " << "cm." << "  " << eta_in << "  " << eta_out << "     " << (R + (delR/2.0)) << " " << "cm." <<  "        " << "(" << 2*f_trd_x1 << "," << 2*f_trd_x2 << "," << 2*f_trd_y << "," << 2*f_trd_y << ")" << " " << "cm." << endl; 

	double phi_lo = -f_spacing_phi;
	double phi_hi = phi_lo + f_spacing_phi;


	for (int j = 0; j < f_nbins_phi; ++j) {

		


		double xcent = (R + (delR/2.0))*cos(M_PI/2.0 - ((2*j)+1)*(f_spacing_phi/2.0));
		double ycent = (R + (delR/2.0))*sin(M_PI/2.0 - ((2*j)+1)*(f_spacing_phi/2.0));
    


       string     f_name  = Form("tile%d%s%d", i, " ", j);

		//EightPointSolid f_shape(f_trd_z, vert);


       // Box        f_shape(dx/2.0,dy/2.0,s_dim_z/2.0);
		Trapezoid  f_shape(f_trd_x1, f_trd_x2, f_trd_z, f_trd_z, f_trd_y);
		Volume     f_vol(f_name, f_shape, desc.material(x_tiles.materialStr()));
		DetElement tower(slice, f_name, det_id);

        
        f_vol.setAttributes(desc, x_tiles.regionStr(), x_tiles.limitsStr(), x_tiles.visStr());

        // Slice placement.

        
        PlacedVolume tower_phv = s_vol_outer.placeVolume(f_vol, Transform3D(RotationZYX(0, ((2*j)+1)*(f_spacing_phi/2.0), -M_PI / 2.0), Position(xcent, ycent, 0.0)));
        
        if ( x_tiles.isSensitive() ) {
          sens.setType("calorimeter");
          f_vol.setSensitiveDetector(sens);
          sens_volumes.push_back(tower_phv);
        }
        tower_phv.addPhysVolID("tile", f_num);
        tower.setPlacement(tower_phv);
        

		f_num++;

		phi_lo += f_spacing_phi;
		phi_hi += f_spacing_phi;
	}
  }

   }
  

  




static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  xml_det_t      x_det    = e;
  xml_dim_t      dim      = x_det.dimensions();
  int            det_id   = x_det.id();
  bool           reflect  = x_det.reflect(true);
  string         det_name = x_det.nameStr();
  Material       air      = description.air();
  int            numsides = dim.numsides();
  xml::Component pos      = x_det.position();
  double         rmin     = dim.rmin();
  double         rmax     = dim.rmax();
  double         zmin     = dim.zmin();
  Layering       layering(x_det);
  double         totalThickness = layering.totalThickness();

  double         zmax     = zmin + totalThickness;

  double       eta_max = 1*(log(tan(0.5*(atan(rmin/zmin)))));
  double       eta_min = 1*(log(tan(0.5*(atan(rmax/zmax)))));

  cout << "eta_max" << " " << eta_max << endl;
  cout << "eta_min" << " " << eta_min << endl;
  cout << "zmin" << " " << zmin << endl;
  cout << "rmax" << " " << rmax << endl;

  const std::vector< double > z = { totalThickness/2.0, -totalThickness/2.0 };

  const std::vector< double > r_min_outer = { -1*zmin*tan(2*atan(exp(-1*(-2.36466)))) , -1*zmax*tan(2*atan(exp(-1*(-2.36466))))};
  const std::vector< double > r_max_outer = { -1*zmin*tan(2*atan(exp(-1*(eta_min)))) , rmax  };

  const std::vector< double > r_min_middle = { -1*zmin*tan(2*atan(exp(-1*(-3.38219)))) , -1*zmax*tan(2*atan(exp(-1*(-3.38219))))};
  const std::vector< double > r_max_middle = { -1*zmin*tan(2*atan(exp(-1*(-2.36466)))) , -1*zmax*tan(2*atan(exp(-1*(-2.36466))))};

  const std::vector< double > r_min_inner = { rmin , -1*zmax*tan(2*atan(exp(-1*(eta_max))))};
  const std::vector< double > r_max_inner = { -1*zmin*tan(2*atan(exp(-1*(-3.38219)))) , -1*zmax*tan(2*atan(exp(-1*(-3.38219))))};

  Volume         endcapVol_outer("endcap_outer", Polyhedra(4*numsides, 0.0, 2*M_PI, z , r_min_outer, r_max_outer), air);
  Volume         endcapVol_middle("endcap_middle", Polyhedra(2*numsides, 0.0, 2*M_PI, z , r_min_middle, r_max_middle), air);
  Volume         endcapVol_inner("endcap_inner", Polyhedra(numsides, 0.0, 2*M_PI, z , r_min_inner, r_max_inner), air);

  Volume         endcapVol_forward("endcap", PolyhedraRegular(numsides, rmin, rmax, totalThickness), air);


  if (det_name == "PassiveSteelRingEndcapP")
  {
    endcapVol_outer = endcapVol_forward;
  }


  //Volume         endcapVol("endcap", PolyhedraRegular(numsides, rmin, rmax, totalThickness), air);
  DetElement     endcap("endcap", det_id);

  

  std::cout << z[0] << "  " << z[1] << "endcap" << endl;
  //std::cout << r_min[0] << "  " << r_min[1] << "endcap" << endl;
  //std::cout << r_max[0] << "  " << r_max[1] << "endcap" << endl;

  // std::cout << "rmin = " << rmin << "\n";
  // std::cout << "rmax = " << rmax << "\n";
  // std::cout << "nlayers = " << std::size(layering.layers()) << "\n";
  int    l_num     = 1;
  int    layerType = 0;
  double layerZ    = -totalThickness / 2;

 ////////////////start
  DetElement    stave_det("stave0",det_id);
  

  //double  dphi = (2*M_PI/numsides);
  //double  hphi = dphi/2;
  //double  dx = 0.0;
  //double trd_x1 = (2 * std::tan(hphi) * rmin + dx)/2 ;
  //double l_dim_x  = trd_x1;

  //double trd_y1 = dim.z()/2 ;
  //double stave_z  = trd_y1;

  //double l_pos_z = -(layering.totalThickness() / 2);

 ///////////////////end



  endcapVol_outer.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  endcapVol_middle.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  endcapVol_inner.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
 
  ///endcapVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

    xml_coll_t xc(x_det, _U(layer));
    xml_comp_t x_layer = xc;
     int     l_repeat = x_layer.repeat();

     int    s_num  = 1;
     if (l_repeat <= 0)  {
      throw std::runtime_error(x_det.nameStr() + "> Invalid repeat value");
     }
      for (int j = 0; j < l_repeat; ++j) {


        for (xml_coll_t xf(x_det, _U(layer)); xf; ++xf) {
        std::cout << "Layer No. = " << l_num << endl;
        cout << "Tile no" << "    " << "Height" << "    " << "eta_in" << "    " << "eta_out" << "   " << "Radial Coordinate" <<  "                  " << "Edges" << endl; 
        //////cout << "Height" << "  " << delR <<  "   " << "Inner eta edge" << "   " << eta_in << "   " <<  "Outer eta edge" << "   " << eta_out << endl;
        //////cout << " Edges" << "  " << 2*f_trd_x1 <<  "   " << 2*f_trd_x2 << "   " << 2*f_trd_y << "   " << 2*f_trd_y << endl;

       // std::cout << "xc = " << xc << "\n";
       x_layer = xf;
        double     l_thick = layering.layer(l_num - 1)->thickness();
        // std::cout << "xc = " << xc << "\n";
        string               l_name   = _toString(layerType, "layer%d");
        l_repeat = x_layer.repeat();
    
        const std::vector< double > z_l = { l_thick/2.0,  -l_thick/2.0 };

        const std::vector< double > r_min_l_outer = { -( zmin + totalThickness - l_num*l_thick )*tan(2*atan(exp(-1*(-2.36466)))), -( zmin + totalThickness - (l_num - 1.0)*l_thick )*tan(2*atan(exp(-1*(-2.36466)))) };
        const std::vector< double > r_max_l_outer = { -( zmin + totalThickness - l_num*l_thick )*tan(2*atan(exp(-1*(eta_min)))) , -( zmin + totalThickness - (l_num - 1.0)*l_thick )*tan(2*atan(exp(-1*(eta_min)))) };

        const std::vector< double > r_min_l_middle = { -( zmin + totalThickness - l_num*l_thick )*tan(2*atan(exp(-1*(-3.38219)))), -( zmin + totalThickness - (l_num - 1.0)*l_thick )*tan(2*atan(exp(-1*(-3.38219)))) };
        const std::vector< double > r_max_l_middle = { -( zmin + totalThickness - l_num*l_thick )*tan(2*atan(exp(-1*(-2.36466)))) , -( zmin + totalThickness - (l_num - 1.0)*l_thick )*tan(2*atan(exp(-1*(-2.36466)))) };

        const std::vector< double > r_min_l_inner = { -( zmin + totalThickness - l_num*l_thick )*tan(2*atan(exp(-1*(eta_max)))), -( zmin + totalThickness - (l_num - 1.0)*l_thick )*tan(2*atan(exp(-1*(eta_max)))) };
        const std::vector< double > r_max_l_inner = { -( zmin + totalThickness - l_num*l_thick )*tan(2*atan(exp(-1*(-3.38219)))) , -( zmin + totalThickness - (l_num - 1.0)*l_thick )*tan(2*atan(exp(-1*(-3.38219)))) };


        Volume               l_vol_outer(l_name, Polyhedra(4*numsides, 0.0, 2*M_PI, z_l , r_min_l_outer, r_max_l_outer), air);
        Volume               l_vol_middle(l_name, Polyhedra(2*numsides, 0.0, 2*M_PI, z_l , r_min_l_middle, r_max_l_middle), air);
        Volume               l_vol_inner(l_name, Polyhedra(numsides, 0.0, 2*M_PI, z_l , r_min_l_inner, r_max_l_inner), air);
       
        Volume               l_vol_forward(l_name, PolyhedraRegular(numsides, rmin, rmax, l_thick), air);

        
        if (det_name == "PassiveSteelRingEndcapP")
       {
          l_vol_outer = l_vol_forward;
       }
       

        vector<PlacedVolume> sensitives;

        //////std::cout << z_l[0] << "  " << z_l[1] << "layer" << "   " << l_num << endl;
        //std::cout << r_min_l[0] << "  " << r_min_l[1] << "layer" << l_num << endl;
        //std::cout << r_max_l[0] << "  " << r_max_l[1] << "layer" << l_num << endl;


        ////////////start
        DetElement layer(stave_det, l_name, det_id);
    
    
       ////////////end

       
       double sliceZ = -l_thick / 2;
       double z_slice_end = zmin + totalThickness - (l_num - 1.0)*l_thick;

        for (xml_coll_t xs(x_layer, _U(slice)); xs; ++xs) {
       xml_comp_t x_slice = xs;
       string     s_name  = _toString(s_num, "slice%d");
       double     s_thick = x_slice.thickness();
       Material   s_mat   = description.material(x_slice.materialStr());

       double z_slice_start = z_slice_end - s_thick;


       const std::vector< double > z_s = { s_thick/2.0,  -s_thick/2.0  };

       const std::vector< double > r_min_s_outer = { -1*z_slice_start*(tan(2*atan(exp(-1*(-2.36466))))), -1*z_slice_end*(tan(2*atan(exp(-1*(-2.36466))))) };
       const std::vector< double > r_max_s_outer = { -1*z_slice_start*tan(2*atan(exp(-1*(eta_min)))) , -1*z_slice_end*tan(2*atan(exp(-1*(eta_min))))  };

       const std::vector< double > r_min_s_middle = { -1*z_slice_start*(tan(2*atan(exp(-1*(-3.38219))))), -1*z_slice_end*(tan(2*atan(exp(-1*(-3.38219))))) };
       const std::vector< double > r_max_s_middle = { -1*z_slice_start*tan(2*atan(exp(-1*(-2.36466)))) , -1*z_slice_end*tan(2*atan(exp(-1*(-2.36466))))  };

       const std::vector< double > r_min_s_inner = { -1*z_slice_start*(tan(2*atan(exp(-1*(eta_max))))), -1*z_slice_end*(tan(2*atan(exp(-1*(eta_max))))) };
       const std::vector< double > r_max_s_inner = { -1*z_slice_start*tan(2*atan(exp(-1*(-3.38219)))) , -1*z_slice_end*tan(2*atan(exp(-1*(-3.38219))))  };
       
       //double z_slice = -0.5*(z_slice_start + z_slice_end);
       double z_slice = -1.0*z_slice_start;
       
       Volume     s_vol_outer(s_name, Polyhedra(4*numsides, 0.0, 2*M_PI , z_s , r_min_s_outer, r_max_s_outer), s_mat);
       Volume     s_vol_middle(s_name, Polyhedra(2*numsides, 0.0, 2*M_PI , z_s , r_min_s_middle, r_max_s_middle), s_mat);
       Volume     s_vol_inner(s_name, Polyhedra(numsides, 0.0, 2*M_PI , z_s , r_min_s_inner, r_max_s_inner), s_mat);

       Volume     s_vol_forward(s_name, PolyhedraRegular(numsides, rmin, rmax, s_thick), s_mat);

       
       if (det_name == "PassiveSteelRingEndcapP")
       {
          s_vol_outer = s_vol_forward;
       }
       


       /////////////start

       
       DetElement   slice(layer, s_name, det_id);
       //DetElement   slice_middle(layer_middle, s_name, det_id);
       //DetElement   slice_inner(layer_inner, s_name, det_id);


       /////std::cout << z_s[0] << "  " << z_s[1] << "slice" << s_num << endl;
       //std::cout << r_min_s[0] << "  " << r_min_s[1] << "slice" << s_num << endl ;
       //std::cout << r_max_s[0] << "  " << r_max_s[1] << "slice" << s_num << endl;

       ///////////end


       s_vol_outer.setVisAttributes(description.visAttributes(x_slice.visStr()));
       s_vol_middle.setVisAttributes(description.visAttributes(x_slice.visStr()));
       s_vol_inner.setVisAttributes(description.visAttributes(x_slice.visStr()));

       ////s_vol.setVisAttributes(description.visAttributes(x_slice.visStr()));

       sliceZ += s_thick / 2;

       /////////////////start

       //double layer_glob_rad = rmin+layering.totalThickness()/2.0 +l_pos_z+totalThickness/2.0+l_pos_z+s_thick/2.0;
      

       ////////////////end

       if (x_slice.hasChild(_Unicode(tiles))) {
              buildTiles(description, sens, s_vol_outer, s_vol_middle, s_vol_inner, slice, det_id, x_slice.child(_Unicode(tiles)), z_slice, sensitives);
          }

       PlacedVolume s_phv_outer = l_vol_outer.placeVolume(s_vol_outer, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, sliceZ)));
       PlacedVolume s_phv_middle;
       PlacedVolume s_phv_inner;

       /*
       
       PlacedVolume s_phv_outer = l_vol_outer.placeVolume(s_vol_outer, Position(0, 0, sliceZ));
       PlacedVolume s_phv_middle = l_vol_middle.placeVolume(s_vol_middle, Position(0, 0, sliceZ));
       PlacedVolume s_phv_inner = l_vol_inner.placeVolume(s_vol_inner, Position(0, 0, sliceZ));

       */

       if (s_num % 2 == 1)

       {
        s_phv_middle = l_vol_middle.placeVolume(s_vol_middle, Transform3D(RotationZYX(M_PI / 2*numsides, 0, 0), Position(0, 0, sliceZ)));
        s_phv_inner = l_vol_inner.placeVolume(s_vol_inner, Transform3D(RotationZYX(M_PI / 2*numsides, 0, 0), Position(0, 0, sliceZ)));
       }
       
       if (s_num % 2 == 0)

       {
        s_phv_middle = l_vol_middle.placeVolume(s_vol_middle, Transform3D(RotationZYX(M_PI / numsides, 0, 0), Position(0, 0, sliceZ)));
        s_phv_inner = l_vol_inner.placeVolume(s_vol_inner, Transform3D(RotationZYX(M_PI / numsides, 0, 0), Position(0, 0, sliceZ)));
       }
       ///PlacedVolume s_phv = l_vol.placeVolume(s_vol, Position(0, 0, sliceZ));


       s_phv_outer.addPhysVolID("slice", s_num);
       s_phv_middle.addPhysVolID("slice", s_num);
       s_phv_inner.addPhysVolID("slice", s_num);

       ///s_phv.addPhysVolID("slice", s_num);

       if (x_slice.isSensitive()) {
        sens.setType("calorimeter");

        s_vol_outer.setSensitiveDetector(sens);
        sensitives.push_back(s_phv_outer);

        s_vol_middle.setSensitiveDetector(sens);
        sensitives.push_back(s_phv_middle);

        s_vol_inner.setSensitiveDetector(sens);
        sensitives.push_back(s_phv_inner);

       }
       sliceZ += s_thick / 2;
       z_slice_end = z_slice_start;
       s_num++;
       }

       l_vol_outer.setVisAttributes(description.visAttributes(x_layer.visStr()));
       l_vol_middle.setVisAttributes(description.visAttributes(x_layer.visStr()));
       l_vol_inner.setVisAttributes(description.visAttributes(x_layer.visStr()));

       ////l_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

       string phys_lay = _toString(l_num, "layer%d");
       layerZ += l_thick / 2;

       DetElement   layer_elt(endcap, phys_lay, l_num);


       PlacedVolume pv_outer = endcapVol_outer.placeVolume(l_vol_outer, Position(0, 0, layerZ));
       PlacedVolume pv_middle = endcapVol_middle.placeVolume(l_vol_middle, Position(0, 0, layerZ));
       PlacedVolume pv_inner = endcapVol_inner.placeVolume(l_vol_inner, Position(0, 0, layerZ));

       /////PlacedVolume pv = endcapVol.placeVolume(l_vol, Position(0, 0, layerZ));

       pv_outer.addPhysVolID("layer", l_num);
       pv_middle.addPhysVolID("layer", l_num);
       pv_inner.addPhysVolID("layer", l_num);

       ///pv.addPhysVolID("layer", l_num);

       layer_elt.setPlacement(pv_outer);
       layer_elt.setPlacement(pv_middle);
       layer_elt.setPlacement(pv_inner);

       ////layer_elt.setPlacement(pv);


        for (size_t ic = 0; ic < sensitives.size(); ++ic) {
         PlacedVolume sens_pv = sensitives[ic];
         DetElement   comp_elt(layer_elt, sens_pv.volume().name(), l_num);
         comp_elt.setPlacement(sens_pv);
         }
       layerZ += l_thick / 2;
       ++l_num;
       }

    
    
       ////////start
     
    

       ////////end
    
       ++layerType;
      }

  double       z_pos = zmin + totalThickness / 2;

  PlacedVolume pv_outer;
  PlacedVolume pv_middle;
  PlacedVolume pv_inner;

  PlacedVolume pv;

  // Reflect it.
  Assembly   assembly(det_name);
  DetElement endcapAssyDE(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(endcapAssyDE);
  if (reflect) {
    pv_outer = assembly.placeVolume(endcapVol_outer, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, -z_pos)));
    pv_outer.addPhysVolID("barrel", 2);
    endcap.setPlacement(pv_outer);

    pv_middle = assembly.placeVolume(endcapVol_middle, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, -z_pos)));
    pv_middle.addPhysVolID("barrel", 2);
    endcap.setPlacement(pv_middle);

    pv_inner = assembly.placeVolume(endcapVol_inner, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, -z_pos)));
    pv_inner.addPhysVolID("barrel", 2);
    endcap.setPlacement(pv_inner);


    Ref_t(endcap)->SetName((det_name + "_backward").c_str());
    
  } else {
    pv_outer = assembly.placeVolume(endcapVol_outer, Transform3D(RotationZYX(M_PI / numsides, 0, 0), Position(0, 0, z_pos)));
    pv_outer.addPhysVolID("barrel", 1);
    endcap.setPlacement(pv_outer);

    Ref_t(endcap)->SetName((det_name + "_forward").c_str());
  }
  endcapAssyDE.add(endcap);

  pv = motherVol.placeVolume(assembly, Position(pos.x(), pos.y(), pos.z()));
  pv.addPhysVolID("system", det_id);
  endcapAssyDE.setPlacement(pv);


  return endcapAssyDE;
}

// clang-format off
DECLARE_DETELEMENT(epic_PolyhedraEndcapCalorimeter2, create_detector)
DECLARE_DETELEMENT(epic_PolyhedraEndcapCalorimeter, create_detector)
