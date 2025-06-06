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
		// loop over cuts
		for (; cuts_c; ++cuts_c)
		{
			// get one cut
			xml_comp_t cut_c	= cuts_c;
			// get placement coordinates
			xml_dim_t cut_pos	= cut_c.child(_U(placement));
			double cut_rotX		= cut_pos.attr<double>("rotX");
			double cut_rotY		= cut_pos.attr<double>("rotY");
			double cut_rotZ		= cut_pos.attr<double>("rotZ");
			// get dimentions
			xml_dim_t cut_dim	= cut_c.child(_U(dimensions));
			double cut_rmin		= cut_dim.attr<double>("rmin");
			double cut_rmax		= cut_dim.attr<double>("rmax");
			double cut_half_l	= cut_dim.attr<double>("half_length");

			// build a solid
			Tube cut_tube(cut_rmin,cut_rmax,cut_half_l);

     			Transform3D tf_tmp(
				RotationZYX(cut_rotX,cut_rotY,cut_rotZ),	
				Position(cut_pos.x(),cut_pos.y(),cut_pos.z()));
			// subtract the cut from the element solid
			elem_final = SubtractionSolid("elem_final",elem_final,cut_tube,tf_tmp);
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
		// get dimentions
		xml_dim_t elem_dim	= elem_c.child(_U(dimensions));
		double elem_rmin	= elem_dim.attr<double>("rmin");
		double elem_rmax	= elem_dim.attr<double>("rmax");
		double elem_half_l	= elem_dim.attr<double>("half_length");
		double elem_sphi	= elem_dim.attr<double>("sphi"); 
		double elem_dphi	= elem_dim.attr<double>("dphi"); 

		// set attributes
		const string elem_vis = dd4hep::getAttrOrDefault<std::string>(elem_c, _Unicode(vis), "FFMagnetVis");
		sdet.setAttributes(dtor, assembly, x_det.regionStr(), x_det.limitsStr(), elem_vis);

		// build solid
		Tube elem_tube(elem_rmin,elem_rmax,elem_half_l,elem_sphi,elem_sphi+elem_dphi);
		Solid elem_final(elem_tube);
		
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
				elem_final = UnionSolid("elem_final",elem_final,elem_tube,tf_tmp);
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

/* from compact/far_forward/ion_beamline.xml  
    <detector id="B0PF_ID" name="B0PF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="B0PF_XPosition" y="0*m" z="B0PF_CenterPosition" theta="B0PF_RotationAngle" />
      <dimensions x="B0PF_InnerRadius*4" y="B0PF_InnerRadius*4" z="B0PF_Length" r="B0PF_InnerRadius*2.0" />
      <apperture x="B0PF_InnerRadius" y="B0PF_InnerRadius" r="B0PF_InnerRadius" />
      <coil vis="FFMagnetCoilVis" dx="2*cm" dy="1.5*cm" />
    </detector>
      <detector id="B0APF_ID" name="B0APF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="B0APF_XPosition" y="0*m" z="B0APF_CenterPosition" theta="B0APF_RotationAngle" />
      <dimensions x="B0APF_InnerRadius*4" y="B0APF_InnerRadius*4" z="B0APF_Length" r="B0APF_InnerRadius*2.0" />
      <apperture x="B0APF_InnerRadius" y="B0APF_InnerRadius" r="B0APF_InnerRadius" />
      <coil vis="FFMagnetCoilVis" dx="2*cm" dy="1.5*cm" />
    </detector>
      <detector id="Q1APF_ID" name="Q1APF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="Q1APF_XPosition" y="0*m" z="Q1APF_CenterPosition" theta="Q1APF_RotationAngle" />
      <dimensions x="Q1APF_InnerRadius*4" y="Q1APF_InnerRadius*4" z="Q1APF_Length" r="2.0*Q1APF_InnerRadius"/>
      <apperture x="Q1APF_InnerRadius*2" y="Q1APF_InnerRadius*2"  r="Q1APF_InnerRadius"/>
      <coil vis="FFMagnetCoilVis" dx="2*cm" dy="1.5*cm" />
    </detector>
    <detector id="Q1BPF_ID" name="Q1BPF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="Q1BPF_XPosition" y="0*m" z="Q1BPF_CenterPosition" theta="Q1BPF_RotationAngle" />
      <dimensions x="Q1BPF_InnerRadius*4" y="Q1BPF_InnerRadius*4" z="Q1BPF_Length"  r="2.0*Q1BPF_InnerRadius"/>
      <apperture x="Q1BPF_InnerRadius*2" y="Q1BPF_InnerRadius*2"  r="Q1BPF_InnerRadius"/>
      <coil vis="FFMagnetCoilVis" dx="1*cm" dy="0.5*cm" />
    </detector>
    <detector id="Q2PF_ID" name="Q2PF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="Q2PF_XPosition" y="0*m" z="Q2PF_CenterPosition" theta="Q2PF_RotationAngle" />
      <dimensions x="Q2PF_InnerRadius*4" y="Q2PF_InnerRadius*4" z="Q2PF_Length" r="2.0*Q2PF_InnerRadius" />
      <apperture x="Q2PF_InnerRadius*2" y="Q2PF_InnerRadius*2"  r="Q2PF_InnerRadius"/>
      <coil vis="FFMagnetCoilVis" dx="1*cm" dy="0.5*cm" />
    </detector>
    <detector id="B1PF_ID" name="B1PF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="B1PF_XPosition" y="0*m" z="B1PF_CenterPosition" theta="B1PF_RotationAngle" />
      <dimensions x="B1PF_InnerRadius*4" y="B1PF_InnerRadius*4" z="B1PF_Length" r="2.0*B1PF_InnerRadius"  />
      <apperture x="B1PF_InnerRadius*2" y="B1PF_InnerRadius*2"  r="B1PF_InnerRadius" />
      <coil vis="FFMagnetCoilVis" dx="1*cm" dy="0.5*cm" />
    </detector>
    <detector id="B1APF_ID" name="B1APF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="B1APF_XPosition" y="0*m" z="B1APF_CenterPosition" theta="B1APF_RotationAngle" />
      <dimensions x="B1APF_InnerRadius*4" y="B1APF_InnerRadius*4" z="B1APF_Length" r="2.0*B1APF_InnerRadius" />
      <apperture x="B1APF_InnerRadius*2" y="B1APF_InnerRadius*2"  r="B1APF_InnerRadius"/>
      <coil vis="FFMagnetCoilVis" dx="1*cm" dy="0.5*cm" />
    </detector>
    <detector id="B2PF_ID" name="B2PF_BeamlineMagnet" vis="FFMagnetVis" type="ip6_CylindricalDipoleMagnet">
      <placement x="B2PF_XPosition" y="0*m" z="B2PF_CenterPosition" theta="B2PF_RotationAngle" />
      <dimensions x="B2PF_InnerRadius*4" y="B2PF_InnerRadius*4" z="B2PF_Length" r="2.0*B2PF_InnerRadius" />
      <apperture x="B2PF_InnerRadius*2" y="B2PF_InnerRadius*2"  r="B2PF_InnerRadius"/>
      <coil vis="FFMagnetCoilVis" dx="1*cm" dy="0.5*cm" />
    </detector>
*/


/* from compact/far_forward/electron_beamline.xml  
    <!-- Q0eF magnet -->
    <detector name="Q0EF" type="ip6_CylindricalDipoleMagnet" vis="RedVis">
      <placement  x="0" y="0" z="(Q0EF_StartZ+Q0EF_EndZ)/2." theta="0"/>
      <dimensions x="Q0EF_InnerRadius*4" y="Q0EF_InnerRadius*4" z="Q0EF_StartZ-Q0EF_EndZ" r="1.9*Q0EF_InnerRadius" />
      <apperture  x="Q0EF_InnerRadius*2" y="Q0EF_InnerRadius*2" r="Q0EF_InnerRadius" />
      <coil dx="2*cm" dy="1.5*cm" />!--unchecked--
    </detector>

    <!-- inner vacuum for Q0eF -->
    <detector name="Q0EF_vac" type="DD4hep_TubeSegment" vis="VisFwElInvisible">
      <material name="Vacuum"/>
      <tubs rmin="0" rmax="Q0EF_InnerRadius" zhalf="(Q0EF_StartZ-Q0EF_EndZ)/2."/>
      <position x="0" y="0" z="(Q0EF_StartZ+Q0EF_EndZ)/2."/>
      <rotation x="0" y="0" z="0"/>
    </detector>

    <!-- Q1eF magnet -->
    <detector name="Q1EF" type="ip6_CylindricalDipoleMagnet" vis="RedVis">
      <placement  x="0" y="0" z="(Q1EF_StartZ+Q1EF_EndZ)/2." theta="0"/>
      <dimensions x="Q1EF_InnerRadius*4" y="Q1EF_InnerRadius*4" z="Q1EF_StartZ-Q1EF_EndZ" r="1.9*Q1EF_InnerRadius" />
      <apperture  x="Q1EF_InnerRadius*2" y="Q1EF_InnerRadius*2" r="Q1EF_InnerRadius" />
      <coil dx="2*cm" dy="1.5*cm" />!--unchecked--
    </detector>

    <!-- inner vacuum for Q1eF -->
    <detector name="Q1EF_vac" type="DD4hep_TubeSegment" vis="VisFwElInvisible">
      <material name="Vacuum"/>
      <tubs rmin="0" rmax="Q1EF_InnerRadius" zhalf="(Q1EF_StartZ-Q1EF_EndZ)/2."/>
      <position x="0" y="0" z="(Q1EF_StartZ+Q1EF_EndZ)/2."/>
      <rotation x="0" y="0" z="0"/>
    </detector>

*/

/* from compact/far_backward/magnets.xml 
    <detector name="Magnet_Q1eR" type="ip6_CylindricalDipoleMagnet" vis="FFMagnetVis">
      <placement  x="0" y="0" z="Q1eR_CenterPosition" theta="0*rad"/>
      <dimensions x="Q1eR_InnerRadius*4" y="Q1eR_InnerRadius*4" z="Q1eR_Length" r="1.5*Q1eR_InnerRadius" />
      <apperture  x="Q1eR_InnerRadius*2" y="Q1eR_InnerRadius*2" r="Q1eR_InnerRadius" />
      <coil dx="2*cm" dy="1.5*cm" /><!--unchecked-->
    </detector>

    <detector name="Magnet_Q2eR" type="ip6_CylindricalDipoleMagnet" vis="FFMagnetVis">
      <placement  x="0" y="0" z="Q2eR_CenterPosition" theta="0*rad"/>
      <dimensions x="Q2eR_InnerRadius*4" y="Q2eR_InnerRadius*4" z="Q2eR_Length" r="1.5*Q2eR_InnerRadius"/>
      <apperture  x="Q2eR_InnerRadius*2" y="Q2eR_InnerRadius*2" r="Q2eR_InnerRadius"/>
      <coil dx="2*cm" dy="1.5*cm" /><!--unchecked-->
    </detector>

    <detector name="Magnet_B2AeR" type="ip6_CylindricalDipoleMagnet" vis="FFMagnetVis">
      <placement  x="0" y="0" z="B2AeR_CenterPosition" theta="0*rad"/>
      <dimensions x="B2AeR_InnerRadius*4" y="B2AeR_InnerRadius*4" z="B2AeR_Length" r="1.5*B2AeR_InnerRadius"/>
      <apperture  x="B2AeR_InnerRadius*2" y="B2AeR_InnerRadius*2" r="B2AeR_InnerRadius"/>
      <coil dx="2*cm" dy="1.5*cm" /><!--unchecked-->
    </detector>

    <detector name="Magnet_B2BeR" type="ip6_CylindricalDipoleMagnet" vis="FFMagnetVis">
      <placement  x="0" y="0" z="B2BeR_CenterPosition" theta="0*rad"/>
      <dimensions x="B2BeR_InnerRadius*4" y="B2BeR_InnerRadius*4" z="B2BeR_Length" r="1.5*B2BeR_InnerRadius"/>
      <apperture  x="B2BeR_InnerRadius*2" y="B2BeR_InnerRadius*2" r="B2BeR_InnerRadius"/>
      <coil dx="2*cm" dy="1.5*cm" /><!--unchecked-->
    </detector>

*/

/* from compact/far_backward/beamline_extension_hadron.xml
    <detector name="Magnet_Q2PR" type="ip6_CylindricalDipoleMagnet" vis="RedVis">
      <placement  x="(Q2PR_StartX+Q2PR_EndX)/2" y="0" z="(Q2PR_StartZ+Q2PR_EndZ)/2" theta="Q1BPR_Theta"/>
      <dimensions x="Q2PR_InnerRadius*4" y="Q2PR_InnerRadius*4" z="Q2PR_Length" r="2.0*Q2PR_InnerRadius"/>
      <apperture  x="Q2PR_InnerRadius*2" y="Q2PR_InnerRadius*2" r="Q2PR_InnerRadius"/>
      <coil dx="2*cm" dy="1.5*cm" />
    </detector>
*/
