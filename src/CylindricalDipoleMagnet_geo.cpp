// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Alex Jentsch, Wouter Deconinck, Whitney Armstrong

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "TMath.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace ROOT::Math;

void buildTubeElement(dd4hep::DetElement &sdet,
                  dd4hep::Assembly &assembly,
                  dd4hep::Detector &dtor,
                  dd4hep::xml::DetElement x_det,
                  string element,
                  dd4hep::Material material); 
void buildPolyElement(dd4hep::DetElement &sdet,
                  dd4hep::Assembly &assembly,
                  dd4hep::Detector &dtor,
                  dd4hep::xml::DetElement x_det,
                  string element,
                  dd4hep::Material material);

static Ref_t build_magnet(Detector& dtor, xml_h e, SensitiveDetector /* sens */) 
{
	xml_det_t x_det	= e;
	int det_id	= x_det.id();
	string det_name	= x_det.nameStr();

	DetElement sdet(det_name, det_id);
	Assembly assembly(det_name + "_assembly");

	// get materials
	Material iron		= dtor.material("Iron");
	Material nbti		= dtor.material("SolenoidCoil");
	Material steel_304l	= dtor.material("StainlessSteel");
	Material alum		= dtor.material("Al6061T6");
	Material steel_a53  	= dtor.material("StainlessSteelA53");

	// build magnet components
	buildTubeElement(sdet,assembly,dtor,x_det,"yoke",iron);
	buildTubeElement(sdet,assembly,dtor,x_det,"coil",nbti);
	buildTubeElement(sdet,assembly,dtor,x_det,"tube",steel_304l);
	buildPolyElement(sdet,assembly,dtor,x_det,"endplate",steel_304l);
	buildTubeElement(sdet,assembly,dtor,x_det,"yokeshield",steel_304l);
	buildTubeElement(sdet,assembly,dtor,x_det,"heatshieldbarrel",alum);
	buildPolyElement(sdet,assembly,dtor,x_det,"heatshieldend",alum);
	buildTubeElement(sdet,assembly,dtor,x_det,"cryobarrel",steel_a53);
	buildPolyElement(sdet,assembly,dtor,x_det,"cryoend",steel_a53);

	// final placement
	auto pv_assembly = dtor.pickMotherVolume(sdet).placeVolume(assembly);
	pv_assembly.addPhysVolID("system", det_id);
	sdet.setPlacement(pv_assembly);

	// update bounding box
	assembly->GetShape()->ComputeBBox();

	return sdet;
}

void buildPolyElement(dd4hep::DetElement &sdet,
                  dd4hep::Assembly &assembly,
                  dd4hep::Detector &dtor,
                  dd4hep::xml::DetElement x_det,
                  string element,
                  dd4hep::Material material) 
{
	// get all elems
	xml_coll_t elems_c(x_det, element.c_str());
	int elem_id = 1;

	// loop over elems
	for (; elems_c; ++elems_c) 
	{
		// get one element
      		xml_comp_t elem_c 	= elems_c;
		string elem_name	= elem_c.nameStr();
		// get placement coordinates
		xml_dim_t elem_pos	= elem_c.child(_U(placement));
		double elem_theta	= elem_pos.attr<double>("theta");
		std::vector<double> z;
		std::vector<double> rmax;
		std::vector<double> rmin;
		// loop over z-planes
		xml_coll_t zplanes_c(elem_c, _Unicode(zplane));
		for (; zplanes_c; ++zplanes_c)
		{
			// get z-plane
			xml_comp_t zplane_c = zplanes_c;
			z.push_back(zplane_c.attr<double>(_Unicode(z)));
			rmin.push_back(zplane_c.attr<double>(_Unicode(rmin)));
			rmax.push_back(zplane_c.attr<double>(_Unicode(rmax)));
		}

		// set attributes
		const string elem_vis = dd4hep::getAttrOrDefault<std::string>(elem_c, _Unicode(vis), "FFMagnetVis");
		sdet.setAttributes(dtor, assembly, x_det.regionStr(), x_det.limitsStr(), elem_vis);

		// build solid
		Polycone elem_tube(0, 2.0 * M_PI, rmin, rmax, z);
		Solid elem_final(elem_tube);
		
		// get all adds
		xml_coll_t adds_c(elem_c, _Unicode(add));
		// loop over adds
		for (; adds_c; ++adds_c)
		{
			Solid add_elem;
			// get one cut
			xml_comp_t add_c	= adds_c;
			// get shape
			string add_shape	= add_c.attr<string>("shape");
			// get placement coordinates
			xml_dim_t add_pos	= add_c.child(_U(placement));
			double add_theta	= add_pos.attr<double>("theta");
			// get dimentions
			xml_dim_t add_dim	= add_c.child(_U(dimensions));

			Transform3D tf_tmp(
				RotationZYX(0, add_theta, 0),
				Position(add_pos.x(),add_pos.y(),add_pos.z()));

			if(add_shape == "Prism")
			{
				double add_pdx1	= add_dim.attr<double>("pdx1");
				double add_pdx2	= add_dim.attr<double>("pdx2");
				double add_pdy1	= add_dim.attr<double>("pdy1");
				double add_pdy2	= add_dim.attr<double>("pdy2");
				double add_pdz	= add_dim.attr<double>("pdz");

				// build a solid
				Trapezoid add_prism(add_pdx1,add_pdx2,add_pdy1,add_pdy2,add_pdz);
				add_elem = add_prism;
			}
			else if(add_shape == "Tube")
			{
				double add_rmin		= add_dim.attr<double>("rmin");
				double add_rmax		= add_dim.attr<double>("rmax");
				double add_half_l	= add_dim.attr<double>("half_length");

				// build a solid
				Tube add_tube(add_rmin,add_rmax,add_half_l);
				add_elem = add_tube;
			}

			// unite the add with the element solid
			elem_final = UnionSolid("elem_final",elem_final,add_elem,tf_tmp);
		}

		// get all cuts
		xml_coll_t cuts_c(elem_c, _Unicode(cut));

		Solid cut_final;
		// loop over cuts
		for (; cuts_c; ++cuts_c)
		{
			// get one cut
			xml_comp_t cut_c = cuts_c;
			// get shape
			string cut_shape = dd4hep::getAttrOrDefault<std::string>(cut_c, _Unicode(shape), "Tube"); 

			// get placement coordinates
			xml_dim_t cut_pos	= cut_c.child(_U(placement));
			double cut_rotX		= cut_pos.attr<double>("rotX");
			double cut_rotY		= cut_pos.attr<double>("rotY");
			double cut_rotZ		= cut_pos.attr<double>("rotZ");
			// get rotation coordinates
			xml_dim_t cut_rot       = cut_c.child(_U(rotation));
			int cut_rot_num		= cut_rot.attr<int>("num");
			double cut_rot_step	= cut_rot.attr<double>("step");
			double cut_rot_start	= cut_rot.attr<double>("start");
			string cut_rot_axis 	= cut_rot.attr<string>("axis");

			if(cut_shape == "Cone")
			{
				// get dimentions
				xml_dim_t cut_dim	= cut_c.child(_U(dimensions));
				double cut_rmin1	= cut_dim.attr<double>("rmin1");
				double cut_rmax1	= cut_dim.attr<double>("rmax1");
				double cut_rmin2	= cut_dim.attr<double>("rmin2");
				double cut_rmax2	= cut_dim.attr<double>("rmax2");
				double cut_dz		= cut_dim.attr<double>("dz");

				// build a solid
				Cone cut_cone(cut_dz,cut_rmin1,cut_rmax1,cut_rmin2,cut_rmax2);
				cut_final = cut_cone;
			}
			else if(cut_shape == "Tube")
			{
				// get dimentions
				xml_dim_t cut_dim	= cut_c.child(_U(dimensions));
				double cut_rmin		= cut_dim.attr<double>("rmin");
				double cut_rmax		= cut_dim.attr<double>("rmax");
				double cut_half_l	= cut_dim.attr<double>("half_length");

				// build a solid
				Tube cut_tube(cut_rmin,cut_rmax,cut_half_l);
				cut_final = cut_tube;
			}

			// loop over rot steps
			for(int i = 0; i < cut_rot_num; i++)
			{
				Position pos_tmp(cut_pos.x(),cut_pos.y(),cut_pos.z());
				double ang_tmp = cut_rot_start + i * cut_rot_step;
				Rotation3D rot_tmp;
				if 	(cut_rot_axis == "X")  	{rot_tmp = RotationX(ang_tmp);}
				else if (cut_rot_axis == "Y")  	{rot_tmp = RotationY(ang_tmp);}
				else  				{rot_tmp = RotationZ(ang_tmp);}
				pos_tmp = rot_tmp * pos_tmp;

				Transform3D tf_tmp(RotationZYX(cut_rotZ,cut_rotY,cut_rotX),pos_tmp);
				// subtract the cut from the element solid
				elem_final = SubtractionSolid("elem_final",elem_final,cut_final,tf_tmp);
			}
		}

		// create volume
		Volume elem_vol(elem_name, elem_final, material);

		// placement 
		auto elem_pv = assembly.placeVolume(elem_vol,
			Transform3D(	Translation3D(elem_pos.x(),elem_pos.y(),elem_pos.z()) * 
					RotationY(elem_theta)));
		elem_pv.addPhysVolID(element, elem_id);
		DetElement elem_de(sdet, elem_name, elem_id);
		elem_de.setPlacement(elem_pv);
		elem_de.setAttributes(dtor, elem_vol, x_det.regionStr(), x_det.limitsStr(), elem_vis);
		elem_id++;
	}

	return;
}

void buildTubeElement(dd4hep::DetElement &sdet,
                  dd4hep::Assembly &assembly,
                  dd4hep::Detector &dtor,
                  dd4hep::xml::DetElement x_det,
                  string element,
                  dd4hep::Material material) 
{
	// get all elems
	xml_coll_t elems_c(x_det, element.c_str());
	int elem_id = 1;

	// loop over elems
	for (; elems_c; ++elems_c) 
	{
		// get one element
      		xml_comp_t elem_c 	= elems_c;
		string elem_name	= elem_c.nameStr();
		// get placement coordinates
		xml_dim_t elem_pos	= elem_c.child(_U(placement));
		double elem_theta	= elem_pos.attr<double>("theta");

		// set attributes
		const string elem_vis = 
			dd4hep::getAttrOrDefault<std::string>(elem_c, _Unicode(vis), "FFMagnetVis");
		sdet.setAttributes(dtor, assembly, x_det.regionStr(), x_det.limitsStr(), elem_vis);

		// get shape
		string elem_shape = dd4hep::getAttrOrDefault<std::string>(elem_c, _Unicode(shape), "Tube"); 

		Solid elem_sub;
		if(elem_shape == "Cone")
		{
			// get dimentions
			xml_dim_t elem_dim	= elem_c.child(_U(dimensions));
			double elem_rmin1	= elem_dim.attr<double>("rmin1");
			double elem_rmax1	= elem_dim.attr<double>("rmax1");
			double elem_rmin2	= elem_dim.attr<double>("rmin2");
			double elem_rmax2	= elem_dim.attr<double>("rmax2");
			double elem_dz		= elem_dim.attr<double>("dz");

			// build a solid
			Cone elem_cone(elem_dz,elem_rmin1,elem_rmax1,elem_rmin2,elem_rmax2);
			elem_sub = elem_cone;
		}
		else if(elem_shape == "ConeSegment")
		{
			// get dimentions
			xml_dim_t elem_dim	= elem_c.child(_U(dimensions));
			double elem_rmin1	= elem_dim.attr<double>("rmin1");
			double elem_rmax1	= elem_dim.attr<double>("rmax1");
			double elem_rmin2	= elem_dim.attr<double>("rmin2");
			double elem_rmax2	= elem_dim.attr<double>("rmax2");
			double elem_dz		= elem_dim.attr<double>("dz");
			double elem_sphi	= elem_dim.attr<double>("sphi");
			double elem_dphi	= elem_dim.attr<double>("dphi");

			// build a solid
			ConeSegment elem_conesegment(
				elem_dz,elem_rmin1,elem_rmax1,elem_rmin2,elem_rmax2,elem_sphi,
				elem_sphi+elem_dphi);
			elem_sub = elem_conesegment;
		}
		else if(elem_shape == "Tube")
		{
			// get dimentions
			xml_dim_t elem_dim	= elem_c.child(_U(dimensions));
			double elem_rmin	= elem_dim.attr<double>("rmin");
			double elem_rmax	= elem_dim.attr<double>("rmax");
			double elem_half_l	= elem_dim.attr<double>("half_length");
			double elem_sphi	= elem_dim.attr<double>("sphi"); 
			double elem_dphi	= elem_dim.attr<double>("dphi"); 

			// build solid
			Tube elem_tube(elem_rmin,elem_rmax,elem_half_l,elem_sphi,elem_sphi+elem_dphi);
			elem_sub = elem_tube;
		}

		Solid elem_final(elem_sub);
		
		// combine sub-elements
		if(elem_pos.hasAttr(_Unicode(phiNum)))
		{
			int phi_num 		= elem_pos.attr<int>("phiNum");
			double phi_step		= elem_pos.attr<double>("phiStep");
			double phi_start	= elem_pos.attr<double>("phiStart");

			// loop over steps
			for(int i = 0; i < phi_num; i++)
			{
				double phi_tmp = phi_start + i * phi_step;
      				Transform3D tf_tmp(RotationZ(phi_tmp),Position(0,0,0));
				// unite sub-elements
				elem_final = UnionSolid("elem_final",elem_final,elem_sub,tf_tmp);
			}
		}

		// get all cuts
		xml_coll_t cuts_c(elem_c, _Unicode(cut));
		// loop over cuts
		for (; cuts_c; ++cuts_c)
		{
			// get one cut
			xml_comp_t cut_c	= cuts_c;
			// get shape
			string add_shape	= dd4hep::getAttrOrDefault<std::string>(
							cut_c, _Unicode(shape), "Tube"); 
			// get placement coordinates
			xml_dim_t cut_pos	= cut_c.child(_U(placement));
			double cut_rotX		= cut_pos.attr<double>("rotX");
			double cut_rotY		= cut_pos.attr<double>("rotY");
			double cut_rotZ		= cut_pos.attr<double>("rotZ");
			// get rotation coordinates
			xml_dim_t cut_rot       = cut_c.child(_U(rotation));
			int cut_rot_num		= cut_rot.attr<int>("num");
			double cut_rot_step	= cut_rot.attr<double>("step");
			double cut_rot_start	= cut_rot.attr<double>("start");
			string cut_rot_axis 	= cut_rot.attr<string>("axis");

			Solid cut_elem;
			if(add_shape == "Box")
			{
				// get dimentions
				xml_dim_t cut_dim	= cut_c.child(_U(dimensions));
				double cut_dx		= cut_dim.attr<double>("dx");
				double cut_dy		= cut_dim.attr<double>("dy");
				double cut_dz		= cut_dim.attr<double>("dz");

				// build a solid
				Box cut_box(cut_dx,cut_dy,cut_dz);
				cut_elem = cut_box;
			}
			else if(add_shape == "Tube")
			{
				// get dimentions
				xml_dim_t cut_dim	= cut_c.child(_U(dimensions));
				double cut_rmin		= cut_dim.attr<double>("rmin");
				double cut_rmax		= cut_dim.attr<double>("rmax");
				double cut_half_l	= cut_dim.attr<double>("half_length");

				// build a solid
				Tube cut_tube(cut_rmin,cut_rmax,cut_half_l);
				cut_elem = cut_tube;
			}
			else if(add_shape == "Cone")
			{
				// get dimentions
				xml_dim_t cut_dim	= cut_c.child(_U(dimensions));
				double cut_rmin1	= cut_dim.attr<double>("rmin1");
				double cut_rmax1	= cut_dim.attr<double>("rmax1");
				double cut_rmin2	= cut_dim.attr<double>("rmin2");
				double cut_rmax2	= cut_dim.attr<double>("rmax2");
				double cut_dz		= cut_dim.attr<double>("dz");

				// build a solid
				Cone cut_cone(cut_dz,cut_rmin1,cut_rmax1,cut_rmin2,cut_rmax2);
				cut_elem = cut_cone;
			}
			// loop over rot steps
			for(int i = 0; i < cut_rot_num; i++)
			{
				Position pos_tmp(cut_pos.x(),cut_pos.y(),cut_pos.z());
				double ang_tmp = cut_rot_start + i * cut_rot_step;
				Rotation3D rot_tmp;
				if 	(cut_rot_axis == "X")  	{rot_tmp = RotationX(ang_tmp);}
				else if (cut_rot_axis == "Y")  	{rot_tmp = RotationY(ang_tmp);}
				else  				{rot_tmp = RotationZ(ang_tmp);}
				pos_tmp = rot_tmp * pos_tmp;

				Transform3D tf_tmp(RotationZYX(cut_rotZ,cut_rotY,cut_rotX),pos_tmp);
				// subtract the cut from the element solid
				elem_final = SubtractionSolid("elem_final",elem_final,cut_elem,tf_tmp);
			}
		}

		// create volume
		Volume elem_vol(elem_name, elem_final, material);

		// placement 
		auto elem_pv = assembly.placeVolume(elem_vol,
			Transform3D(	Translation3D(elem_pos.x(),elem_pos.y(),elem_pos.z()) * 
					RotationY(elem_theta)));
		elem_pv.addPhysVolID(element, elem_id);
		DetElement elem_de(sdet, elem_name, elem_id);
		elem_de.setPlacement(elem_pv);
		elem_de.setAttributes(dtor, elem_vol, x_det.regionStr(), x_det.limitsStr(), elem_vis);
		elem_id++;
	}

	return;
}

DECLARE_DETELEMENT(ip6_CylindricalDipoleMagnet, build_magnet)
