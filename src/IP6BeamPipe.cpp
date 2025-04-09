// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong, Sylvester Joosten

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
#include "XML/Utilities.h"
#include "DD4hepDetectorHelper.h"
#include <tuple>
#include "TGeoTessellated.h"
#include "DD4hep/Shapes.h"
#include "TGDMLParse.h"

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
static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) 
{
	using namespace ROOT::Math;
	xml_det_t x_det       = e;
	string det_name       = x_det.nameStr();
	xml_comp_t x_dettype  = x_det.child(dd4hep::xml::Strng_t("type_flags"));
	unsigned int typeFlag = x_dettype.type();
	DetElement sdet(det_name, x_det.id());
	Assembly assembly(det_name + "_assembly");

	xml::Component IP_pipe_c = x_det.child(_Unicode(IP_pipe));

	Material IP_beampipe_wall_material     	= det.material(IP_pipe_c.attr<string>(_Unicode(wall_material)));
	Material IP_beampipe_coating_material	= det.material(IP_pipe_c.attr<string>(_Unicode(coating_material)));
	Material m_Vacuum 			= det.material("Vacuum");
	Material m_Wall     			= det.material("StainlessSteel");
	Material m_Coating  	   		= det.material("Copper");
	string vis_name   			= x_det.visStr();

	// IP Beampipe
	double IP_beampipe_ID             	= IP_pipe_c.attr<double>(_Unicode(ID));
	double IP_beampipe_wall_thickness 	= IP_pipe_c.attr<double>(_Unicode(wall_thickness));
	double IP_beampipe_coating_thickness 	= IP_pipe_c.attr<double>(_Unicode(coating_thickness));
	double IP_acts_beampipe_OD 		= IP_beampipe_ID - 5.0 * mm;
	double IP_acts_beampipe_ID 		= IP_acts_beampipe_OD - 1.0 * mm;

	double upstream_straight_length   = IP_pipe_c.attr<double>(_Unicode(upstream_straight_length));
	double downstream_straight_length = IP_pipe_c.attr<double>(_Unicode(downstream_straight_length));

	// visualization
	VisAttr wallVis("wall");
	VisAttr coatingVis("coating");
	VisAttr IPwallVis("IPwall");
	VisAttr IPcoatingVis("IPcoating");

	// colors: (r, g, b, alpha)
    	wallVis.setColor(0.0, 0.0, 1.0, 1.0);  // blue
    	coatingVis.setColor(1.0, 0.0, 0.0, 1.0); // red
    	IPwallVis.setColor(0.0, 1.0, 0.0, 1.0); // green
    	IPcoatingVis.setColor(1.0, 1.0, 0.0, 1.0); // yellow

	// central acts beampipe volume
	Tube central_tube(
		0.5 * IP_acts_beampipe_ID, 
		0.5 * IP_acts_beampipe_OD,
		0.5 * (upstream_straight_length + downstream_straight_length));
	Volume central_volume("acts_central_beampipe_vol", central_tube, m_Vacuum);
	const double central_offset = -.5 * (upstream_straight_length - downstream_straight_length);
	DetElement central_det(sdet, "acts_beampipe_central", 1);

	// Set dd4hep variant parameters for conversion to ACTS tracking geometry
	central_det.setTypeFlag(typeFlag);
	auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(central_det);
	int nBinPhi  = 144; // fix later. Should take this from a xml tag
	int nBinZ    = 10;  // fix later. Should take this from a xml tag
	params.set<bool>("layer_material", true);
	params.set<bool>("layer_material_representing", true);
	params.set<int>("layer_material_representing_binPhi", nBinPhi);
	params.set<int>("layer_material_representing_binZ", nBinZ);

	// -----------------------------
	// IP beampipe
	//
	// setup for the central IP beampipe:
	//
	// /-------\ Be wall
	//  /-----\  Au coating
	//   /---\   Vacuum padding (5mm)
	//    /-\    Fake vacuum beampipe (1mm)
	//     -     Vacuum filled inner beampipe
	//
	Tube downstream_IP_tube(
		IP_beampipe_ID / 2.0, IP_beampipe_ID / 2.0 + IP_beampipe_wall_thickness, downstream_straight_length / 2.0);
	Tube downstream_IP_coating(
		IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness, IP_beampipe_ID / 2.0, downstream_straight_length / 2.0);
	Tube downstream_IP_vacuum_padding(
		IP_acts_beampipe_OD / 2.0, IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness, downstream_straight_length / 2.0);
  	Tube downstream_IP_acts_beampipe(
		IP_acts_beampipe_ID / 2.0, IP_acts_beampipe_OD / 2.0, downstream_straight_length / 2.0);
  	Tube downstream_IP_vacuum_fill(
		0.0, IP_acts_beampipe_ID / 2.0, downstream_straight_length / 2.0);

	Tube upstream_IP_tube(
		IP_beampipe_ID / 2.0, IP_beampipe_ID / 2.0 + IP_beampipe_wall_thickness, upstream_straight_length / 2.0);
	Tube upstream_IP_coating(
		IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness, IP_beampipe_ID / 2.0, upstream_straight_length / 2.0);
	Tube upstream_IP_vacuum_padding(
		IP_acts_beampipe_OD / 2.0, IP_beampipe_ID / 2.0 - IP_beampipe_coating_thickness, upstream_straight_length / 2.0);
	Tube upstream_IP_acts_beampipe(
		IP_acts_beampipe_ID / 2.0, IP_acts_beampipe_OD / 2.0, upstream_straight_length / 2.0);
	Tube upstream_IP_vacuum_fill(
		0.0, IP_acts_beampipe_ID / 2.0, upstream_straight_length / 2.0);

	// create volumes
	Volume v_downstream_IP_vacuum_fill("v_downstream_IP_vacuum_fill", downstream_IP_vacuum_fill, m_Vacuum);
	Volume v_downstream_IP_acts_beampipe("v_downstream_IP_acts_beampipe", downstream_IP_acts_beampipe, m_Vacuum);
	Volume v_downstream_IP_vacuum_padding("v_downstream_IP_vacuum_padding", downstream_IP_vacuum_padding, m_Vacuum);
	Volume v_downstream_IP_coating("v_downstream_IP_coating", downstream_IP_coating, IP_beampipe_coating_material);
	Volume v_downstream_IP_tube("v_downstream_IP_tube", downstream_IP_tube, IP_beampipe_wall_material);
	Volume v_upstream_IP_vacuum_fill("v_upstream_IP_vacuum_fill", upstream_IP_vacuum_fill, m_Vacuum);
	Volume v_upstream_IP_acts_beampipe("v_upstream_IP_acts_beampipe", upstream_IP_acts_beampipe, m_Vacuum);
	Volume v_upstream_IP_vacuum_padding("v_upstream_IP_vacuum_padding", upstream_IP_vacuum_padding, m_Vacuum);
	Volume v_upstream_IP_coating("v_upstream_IP_coating", upstream_IP_coating, IP_beampipe_coating_material);
	Volume v_upstream_IP_tube("v_upstream_IP_tube", upstream_IP_tube, IP_beampipe_wall_material);

	// set vis colors
	v_downstream_IP_tube.setVisAttributes(IPwallVis);
	v_upstream_IP_tube.setVisAttributes(IPwallVis);
	v_downstream_IP_coating.setVisAttributes(IPcoatingVis);
	v_upstream_IP_coating.setVisAttributes(IPcoatingVis);

	// set names
	sdet.setAttributes(det, v_upstream_IP_coating, x_det.regionStr(), x_det.limitsStr(), vis_name);
	sdet.setAttributes(det, v_upstream_IP_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);
	sdet.setAttributes(det, v_downstream_IP_coating, x_det.regionStr(), x_det.limitsStr(), vis_name);
	sdet.setAttributes(det, v_downstream_IP_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

	// place volumes
	assembly.placeVolume(v_upstream_IP_vacuum_fill, Position(0, 0, -upstream_straight_length / 2.0));
	central_volume.placeVolume(v_upstream_IP_acts_beampipe, Position(0, 0, -upstream_straight_length / 2.0 - central_offset));
	assembly.placeVolume(v_upstream_IP_vacuum_padding, Position(0, 0, -upstream_straight_length / 2.0));
	assembly.placeVolume(v_upstream_IP_coating, Position(0, 0, -upstream_straight_length / 2.0));
	assembly.placeVolume(v_upstream_IP_tube, Position(0, 0, -upstream_straight_length / 2.0));

	assembly.placeVolume(v_downstream_IP_vacuum_fill, Position(0, 0, downstream_straight_length / 2.0));
	central_volume.placeVolume(v_downstream_IP_acts_beampipe, Position(0,0,downstream_straight_length / 2.0 - central_offset));
	assembly.placeVolume(v_downstream_IP_vacuum_padding, Position(0, 0, downstream_straight_length / 2.0));
	assembly.placeVolume(v_downstream_IP_coating, Position(0, 0, downstream_straight_length / 2.0));
	assembly.placeVolume(v_downstream_IP_tube, Position(0, 0, downstream_straight_length / 2.0));

	auto central_pv = assembly.placeVolume(central_volume, Position(0, 0, +central_offset));
	central_det.setPlacement(central_pv);

	//---------------------------------------------------------------------------------------------------------
	// Helper function to create polycone pairs (wall, coating, and vacuum)
	auto zplane_to_polycones = [](xml::Component& x_pipe) 
	{
		std::vector<double> zero, z;
		std::vector<double> rmax_wall, rmax_coating, rmax_vacuum;
		std::vector<double> rmin_wall, rmin_coating, rmin_vacuum;
		// thickness
		auto wall_thickness = 
			getAttrOrDefault(x_pipe, _Unicode(wall_thickness), 1 * mm);
		auto coating_thickness = 
			getAttrOrDefault(x_pipe, _Unicode(coating_thickness), 30 * um);
	
		for (xml_coll_t x_zplane_i(x_pipe, _Unicode(zplane)); x_zplane_i; ++x_zplane_i) 
		{
			xml_comp_t x_zplane = x_zplane_i;
			// z position
			z.push_back(x_zplane.attr<double>(_Unicode(z)));
			// outer radius
			rmax_wall.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 + wall_thickness);
			rmax_coating.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0);
			rmax_vacuum.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 - coating_thickness);
			// inner radius
			rmin_wall.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0);
			rmin_coating.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 - coating_thickness);
			rmin_vacuum.push_back(0);
		}
		return std::tuple<Polycone, Polycone, Polycone>(
			{0, 2.0 * M_PI, rmin_wall, rmax_wall, z},
			{0, 2.0 * M_PI, rmin_coating, rmax_coating, z},
			{0, 2.0 * M_PI, rmin_vacuum, rmax_vacuum, z}
		);
	};
	//---------------------------------------------------------------------------------------------------------
	// Helper function to create a racetrack volumes 
  	auto create_racetrack_solids = [](xml::Component& x_racetrack)
	{
		// get geometry parameters
		double semiCircle_rmin = getAttrOrDefault(x_racetrack, _Unicode(semiCircle_rmin), 0.0);
		double wall_thickness = getAttrOrDefault(x_racetrack, _Unicode(wall_thickness), 0.0);
		double length = getAttrOrDefault(x_racetrack, _Unicode(length), 0.0);
		double rectangle_h = getAttrOrDefault(x_racetrack, _Unicode(rectangle_h), 0.0);
		
		// create semicylinders and boxes
		Tube wall_semiCircle(0, semiCircle_rmin +  wall_thickness, length/2.);
		Box wall_rectangle(semiCircle_rmin + wall_thickness, rectangle_h/2., length/2.);

		Tube coating_semiCircle(0, semiCircle_rmin, length/2.);
		Box coating_rectangle(semiCircle_rmin, rectangle_h/2., length/2.);

		Tube vacuum_semiCircle(0, semiCircle_rmin -  wall_thickness, length/2.);
		Box vacuum_rectangle(semiCircle_rmin - wall_thickness, rectangle_h/2., length/2.);

		// unite semicylinders and boxes
		Transform3D upShift_tf(RotationZYX(0, 0, 0), Position(0, rectangle_h/2., 0));
		Transform3D downShift_tf(RotationZYX(0, 0, 0), Position(0, -rectangle_h/2., 0));

		UnionSolid wall_solid, coating_solid, vacuum_solid;

		wall_solid = 
			UnionSolid("wall_solid_tmp",wall_rectangle, wall_semiCircle, upShift_tf);
		wall_solid = 
			UnionSolid("wall_solid",wall_solid, wall_semiCircle, downShift_tf);

		coating_solid = 
			UnionSolid("coating_solid_tmp",coating_rectangle, coating_semiCircle, upShift_tf);
		coating_solid = 
			UnionSolid("coating_solid",coating_solid, coating_semiCircle, downShift_tf);

		vacuum_solid = 
			UnionSolid("vacuum_solid_tmp",vacuum_rectangle, vacuum_semiCircle, upShift_tf);
		vacuum_solid = 
			UnionSolid("vacuum_solid",vacuum_solid, vacuum_semiCircle, downShift_tf);

		// subtract vacuum 
		BooleanSolid wall_subtract, coating_subtract, vacuum_subtract;

		wall_subtract = 
			SubtractionSolid("wall_subtract",wall_solid,coating_solid,Transform3D());
		coating_subtract =
			SubtractionSolid("coating_subtract",coating_solid,vacuum_solid,Transform3D());
		vacuum_subtract = vacuum_solid;
	
		return std::tuple<Solid,Solid,Solid>(wall_subtract,coating_subtract,vacuum_subtract);
	};
	//---------------------------------------------------------------------------------------------------------
	// Helper function to create union of lepton and hadron volumes
  	auto create_volumes = [&](
		const std::string& name, 
		xml::Component& x_pipe1,
		xml::Component& x_pipe2, 
		xml_coll_t& x_additional_subtraction_i)
	{
		auto pipe1_polycones = zplane_to_polycones(x_pipe1);
		auto pipe2_polycones = zplane_to_polycones(x_pipe2);

		auto crossing_angle    = getAttrOrDefault(x_pipe2, _Unicode(crossing_angle), 0.0);
		auto axis_intersection = getAttrOrDefault(x_pipe2, _Unicode(axis_intersection), 0.0);

		auto tf = 	Transform3D(Position(0, 0, axis_intersection)) *
				Transform3D(RotationY(crossing_angle)) *
				Transform3D(Position(0, 0, -axis_intersection));

		// union of wall, coating, and vacuum
		BooleanSolid wall_union = UnionSolid(
			std::get<0>(pipe1_polycones), std::get<0>(pipe2_polycones), tf);
		BooleanSolid coating_union = UnionSolid(
			std::get<1>(pipe1_polycones), std::get<1>(pipe2_polycones), tf);
		BooleanSolid vacuum_union = UnionSolid(
			std::get<2>(pipe1_polycones), std::get<2>(pipe2_polycones), tf);

		BooleanSolid wall_racetrack_final, coating_racetrack_final; 

		// downstream side is more complex and requires additional effort to create all the volumes 
		if(name == "downstream")
		{
			xml::Component racetrack_lepton_c = x_pipe1.child(_Unicode(racetrack_lepton));

			// ---- Read geometry parameters ----
			double cylRadius_1 = // cylinder radius on the IP side
				getAttrOrDefault(racetrack_lepton_c, _Unicode(cylRadius_1), 6.2/2. * cm);
			double rtRadius = // racetrack radius
				getAttrOrDefault(racetrack_lepton_c, _Unicode(semiCircle_rmin), 2.3 * cm);
			double wall_thickness = // wall thickness
				getAttrOrDefault(racetrack_lepton_c, _Unicode(wall_thickness), 1.0 * mm);
			double coating_thickness = // coating thickness
				getAttrOrDefault(racetrack_lepton_c, _Unicode(coating_thickness), 30.0 * um);
			double straight_pipe_startz = // straight pipe on the IP side, start position
				getAttrOrDefault(racetrack_lepton_c, _Unicode(straight_pipe_startz), 66.10 * cm);
			double straight_pipe_endz = // straight pipe on the IP side, end position
				getAttrOrDefault(racetrack_lepton_c, _Unicode(straight_pipe_endz), 66.385 * cm);
			double elliptical_cut_rx_1 = // elliptical cut (IP side) rX for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_rx_1), 0.305 * m);
			double elliptical_cut_ry_1 = // elliptical cut (IP side) rY for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_ry_1), 0.021 * m);
			double elliptical_cut_rx_2 = // elliptical cut (non-IP side) rX for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_rx_2), 0.152 * m);
			double elliptical_cut_ry_2 = // elliptical cut (non-IP side) rY for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_ry_2), 0.021 * m);
			double  rectangular_cut_a = // rectangular cut A for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(rectangular_cut_a), 0.81/2. * m);
			double  rectangular_cut_b = // rectangular cut B for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(rectangular_cut_b), 0.021 * m);
			double elliptical_cut_dz = // thickness of the cut for the hadron beam opening
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_dz), 1.0 * cm); 
			double elliptical_cut_offset_z_1 = // offset of the elliptical cut (IP side)
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_offset_z_1), (0.976) * m);
			double elliptical_cut_offset_z_2 = // offset of the elliptical cut (non-IP side)
				getAttrOrDefault(racetrack_lepton_c, _Unicode(elliptical_cut_offset_z_2), (0.976 + 0.810) * m);
			double rectangular_cut_offset_z = // offset of the rectangular cut
				getAttrOrDefault(racetrack_lepton_c, _Unicode(rectangular_cut_offset_z), (0.976 + 0.810/2.) * m);

			// ---- Create racetrack solids ----
			auto racetrack_solids = create_racetrack_solids(racetrack_lepton_c);

/*
// ---- Create an interface between racetrack and cylindrical beam pipe ----
TGDMLParse parser;
string fileName_1 = "/Users/andriinatochii/Desktop/EIC_IR6_BeamPipeModel_in_DD4Hep/interface_vol1.gdml";
string fileName_2 = "/Users/andriinatochii/Desktop/EIC_IR6_BeamPipeModel_in_DD4Hep/interface_vol2.gdml";
Volume interface_vol_1 = parser.GDMLReadFile(fileName_1.c_str());
Volume interface_vol_2 = parser.GDMLReadFile(fileName_2.c_str());
if (!interface_vol_1.isValid()) 
{
	printout(WARNING, "IP6BeamPipe", "%s", fileName_1.c_str());
	printout(WARNING, "IP6BeamPipe", "interface_vol_1 invalid, GDML parser failed!");
	std::_Exit(EXIT_FAILURE);
}
if (!interface_vol_2.isValid()) 
{
	printout(WARNING, "IP6BeamPipe", "%s", fileName_2.c_str());
	printout(WARNING, "IP6BeamPipe", "interface_vol_2 invalid, GDML parser failed!");
	std::_Exit(EXIT_FAILURE);
}

interface_vol_1.import();
TessellatedSolid interface_solid_1 = interface_vol_1.solid();
interface_solid_1->CloseShape(true, true, true); // tesselated solid not closed by import!
interface_vol_1.setMaterial(m_Wall);

interface_vol_2.import();
TessellatedSolid interface_solid_2 = interface_vol_2.solid();
interface_solid_2->CloseShape(true, true, true); // tesselated solid not closed by import!
interface_vol_2.setMaterial(m_Wall);
*/
			// --- Create a straight pipe on the IP side ---
			std::vector<double> z = {straight_pipe_startz, straight_pipe_endz};
			std::vector<double> wall_rmin = {cylRadius_1, cylRadius_1};
			std::vector<double> wall_rmax = {cylRadius_1 + wall_thickness, cylRadius_1 + wall_thickness};
			std::vector<double> coating_rmin = {cylRadius_1 - coating_thickness, cylRadius_1 - coating_thickness};
			std::vector<double> coating_rmax = wall_rmin;

			Polycone wall_pipe(0, 2.0 * M_PI, wall_rmin, wall_rmax, z);
			Polycone coating_pipe(0, 2.0 * M_PI, coating_rmin, coating_rmax, z);

			Polycone wall_pipe_2(0, 2.0 * M_PI, wall_rmin, wall_rmax, z);

			// --- Unite racetrack and straight pipe with the rest of the pipe
			double offset_z = getAttrOrDefault(racetrack_lepton_c, _Unicode(offset_z), 0.0); 
			double length = getAttrOrDefault(racetrack_lepton_c, _Unicode(length), 0.0);
			Transform3D racetrack_tf(RotationZYX(0, 0, 0), Position(0, 0, offset_z + length/2.));

			UnionSolid wall_racetrack_pipe = 
				UnionSolid("wall_racetrack_pipe",wall_pipe,std::get<0>(racetrack_solids), racetrack_tf);
			UnionSolid coating_racetrack_pipe = 
				UnionSolid("coating_racetrack_pipe",coating_pipe,std::get<1>(racetrack_solids), racetrack_tf);

			// --- Subtract racetrack sections from vacuum
			vacuum_union = 
				SubtractionSolid("vacuum_racetrack_cut_1",vacuum_union,wall_racetrack_pipe,Transform3D());
			vacuum_union = 
				SubtractionSolid("vacuum_racetrack_cut_2",vacuum_union,coating_racetrack_pipe,Transform3D());

			// create a cut volume - vacuum = two ellipses + rectangle
			EllipticalTube elliptical_cut_1(
				"elliptical_cut_1", elliptical_cut_rx_1, elliptical_cut_ry_1, elliptical_cut_dz);
			EllipticalTube elliptical_cut_2(
				"elliptical_cut_2", elliptical_cut_rx_2, elliptical_cut_ry_2, elliptical_cut_dz);
			Box rectangular_cut_3(
				"rectangular_cut_3", rectangular_cut_a, rectangular_cut_b, elliptical_cut_dz);
			Transform3D tf_cut_1(
				RotationZYX(0, M_PI_2, 0), Position(-rtRadius, 0, elliptical_cut_offset_z_1));
			Transform3D tf_cut_2(
				RotationZYX(0, M_PI_2, 0), Position(-rtRadius, 0, elliptical_cut_offset_z_2));
			Transform3D tf_cut_3(
				RotationZYX(0, M_PI_2, 0), Position(-rtRadius, 0, rectangular_cut_offset_z));

			// subtract from racetrack wall
			SubtractionSolid wall_racetrack_cut_1 = 
				SubtractionSolid("wall_racetrack_cut_1",wall_racetrack_pipe,elliptical_cut_1,tf_cut_1);
			SubtractionSolid wall_racetrack_cut_2 = 
				SubtractionSolid("wall_racetrack_cut_1",wall_racetrack_cut_1,elliptical_cut_2,tf_cut_2);
			wall_racetrack_final = 
				SubtractionSolid("wall_racetrack_final",wall_racetrack_cut_2,rectangular_cut_3,tf_cut_3);

			// subtract from racetrack coating
			SubtractionSolid coating_racetrack_cut_1 = 
				SubtractionSolid("coating_racetrack_cut_1",coating_racetrack_pipe,elliptical_cut_1,tf_cut_1);
			SubtractionSolid coating_racetrack_cut_2 = 
				SubtractionSolid("coating_racetrack_cut_2",coating_racetrack_cut_1,elliptical_cut_2,tf_cut_2);
			coating_racetrack_final = 
				SubtractionSolid("coating_racetrack_final",coating_racetrack_cut_2,rectangular_cut_3,tf_cut_3);
		}

		Solid wall_racetrack(wall_racetrack_final);
		Solid coating_racetrack(coating_racetrack_final);

		// subtract additional vacuum from wall and coating
		for (; x_additional_subtraction_i; ++x_additional_subtraction_i) 
		{
			xml_comp_t x_additional_subtraction = x_additional_subtraction_i;
			auto additional_polycones = zplane_to_polycones(x_additional_subtraction);
			auto additional_crossing_angle =
				getAttrOrDefault(x_additional_subtraction, _Unicode(crossing_angle), 0.0);
			auto additional_axis_intersection =
				getAttrOrDefault(x_additional_subtraction, _Unicode(axis_intersection), 0.0);
			auto additional_tf = 	Transform3D(Position(0, 0, additional_axis_intersection)) *
						Transform3D(RotationY(additional_crossing_angle)) *
						Transform3D(Position(0, 0, -additional_axis_intersection));
		
			wall_union = SubtractionSolid(
				wall_union, std::get<2>(additional_polycones), additional_tf);
			coating_union = SubtractionSolid(
				coating_union, std::get<2>(additional_polycones), additional_tf);
			vacuum_union = SubtractionSolid(
				vacuum_union, std::get<2>(additional_polycones), additional_tf);
		}

		Solid wall, coating, vacuum, wall_ipflange, coating_ipflange, vacuum_ipflange;

		if(name == "upstream") // upstream
		{
			// subtract vacuum from coating
			coating = SubtractionSolid(coating_union, vacuum_union);

			// subtract vacuum from wall
			wall = SubtractionSolid(wall_union, vacuum_union);

			// get vacuum
			vacuum = vacuum_union;
		}
		else // downstream
		{
			// subtract walls and coatings from vacuum
			auto main = vacuum_union;	main.setName("main"); 
			auto sub1 = wall_union; 	sub1.setName("sub1"); 
			auto sub2 = coating_union; 	sub2.setName("sub2"); 

			gGeoManager->GetListOfShapes()->Add((TGeoShape*)main.ptr());
			gGeoManager->GetListOfShapes()->Add((TGeoShape*)sub1.ptr());
			gGeoManager->GetListOfShapes()->Add((TGeoShape*)sub2.ptr());

			TGeoCompositeShape* composite = new TGeoCompositeShape("composite","main - sub1 - sub2"); 

			// get vacuum
			vacuum = (TGeoCompositeShape*)composite->Clone("vacuum_composite");
			delete composite;

			// get wall
			wall = wall_union;

			// get coating
			coating = coating_union;

			// get a flange between the FWD and IP beam pipes
			xml::Component fwdipflange_c = x_pipe1.child(_Unicode(fwdipflange));
			auto fwdipflange_polycones = zplane_to_polycones(fwdipflange_c);

			wall_ipflange = std::get<0>(fwdipflange_polycones);
			coating_ipflange = std::get<1>(fwdipflange_polycones);
			vacuum_ipflange = std::get<2>(fwdipflange_polycones);
		}

		return std::tuple<Volume, Volume, Volume, Volume, Volume, Volume, Volume, Volume>(
			{"v_" + name + "_wall", wall, m_Wall},
			{"v_" + name + "_coating", coating, m_Coating},
			{"v_" + name + "_vacuum", vacuum, m_Vacuum},
			{"v_" + name + "_vacuum_ipflange", vacuum_ipflange, m_Vacuum},
			{"v_" + name + "_wall_racetrack", wall_racetrack, m_Wall},
			{"v_" + name + "_coating_racetrack", coating_racetrack, m_Coating},
			{"v_" + name + "_wall_ipflange", wall_ipflange, m_Wall},
			{"v_" + name + "_coating_ipflange", coating_ipflange, m_Coating}
		);
	};

	// -----------------------------
	// Upstream/BWD/Rear Side:
	// - incoming hadron tube: straight
	// - outgoing electron tube: tapering

	xml::Component upstream_c		= x_det.child(_Unicode(upstream));
	xml::Component incoming_hadron_c 	= upstream_c.child(_Unicode(incoming_hadron));
	xml::Component outgoing_lepton_c 	= upstream_c.child(_Unicode(outgoing_lepton));
	xml_coll_t additional_subtractions_upstream(upstream_c, _Unicode(additional_subtraction));
	auto volumes_upstream = create_volumes(
		"upstream", 
		outgoing_lepton_c, 
		incoming_hadron_c,
		additional_subtractions_upstream);
    
	std::get<0>(volumes_upstream).setVisAttributes(wallVis);
	std::get<1>(volumes_upstream).setVisAttributes(coatingVis);

	auto tf_upstream = Transform3D(RotationZYX(0, 0, 0));
	if (getAttrOrDefault<bool>(upstream_c, _Unicode(reflect), true)) 
	{
		tf_upstream = Transform3D(RotationZYX(0, M_PI, 0));
	}
	// place wall
	assembly.placeVolume(std::get<0>(volumes_upstream), tf_upstream);
	// place coating
	assembly.placeVolume(std::get<1>(volumes_upstream), tf_upstream);
	// place vacuum
	if (getAttrOrDefault<bool>(upstream_c, _Unicode(place_vacuum), true)) 
	{
		assembly.placeVolume(std::get<2>(volumes_upstream), tf_upstream);
	}

	// -----------------------------
	// Downstream/FWD Side:
	// - incoming electron tube: tube with tube cut out
	// - outgoing hadron tube: cone centered at scattering angle

	xml::Component downstream_c      = x_det.child(_Unicode(downstream));
	xml::Component incoming_lepton_c = downstream_c.child(_Unicode(incoming_lepton));
	xml::Component outgoing_hadron_c = downstream_c.child(_Unicode(outgoing_hadron));
	xml_coll_t additional_subtractions_downstream(downstream_c, _Unicode(additional_subtraction));
	auto volumes_downstream = create_volumes(
		"downstream", 
		incoming_lepton_c, 
		outgoing_hadron_c,
		additional_subtractions_downstream);

	// add colors
	std::get<0>(volumes_downstream).setVisAttributes(wallVis);
	std::get<1>(volumes_downstream).setVisAttributes(coatingVis);
	std::get<3>(volumes_downstream).setVisAttributes(wallVis);
	std::get<4>(volumes_downstream).setVisAttributes(wallVis);
	std::get<5>(volumes_downstream).setVisAttributes(coatingVis);
	std::get<6>(volumes_downstream).setVisAttributes(wallVis);
	std::get<7>(volumes_downstream).setVisAttributes(coatingVis);

	auto tf_downstream = Transform3D(RotationZYX(0, 0, 0));
	if (getAttrOrDefault<bool>(downstream_c, _Unicode(reflect), true)) 
	{
		tf_downstream = Transform3D(RotationZYX(0, M_PI, 0));
	}

	// place wall
	assembly.placeVolume(std::get<0>(volumes_downstream), tf_downstream);
	// place coating
	assembly.placeVolume(std::get<1>(volumes_downstream), tf_downstream);
	// place vacuum
	if (getAttrOrDefault<bool>(downstream_c, _Unicode(place_vacuum), true)) 
	{
		assembly.placeVolume(std::get<2>(volumes_downstream), tf_downstream);
		assembly.placeVolume(std::get<3>(volumes_downstream), tf_downstream);
	}
	// place racetrack wall
	assembly.placeVolume(std::get<4>(volumes_downstream), tf_downstream);
	// place racetrack coating
	assembly.placeVolume(std::get<5>(volumes_downstream), tf_downstream);
	// place FWD IP flange wall
	assembly.placeVolume(std::get<6>(volumes_downstream), tf_downstream);
	// place FWD IP flange coating
	assembly.placeVolume(std::get<7>(volumes_downstream), tf_downstream);

	// -----------------------------
	// final placement
	auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
	pv_assembly.addPhysVolID("system", sdet.id());
	sdet.setPlacement(pv_assembly);
	assembly->GetShape()->ComputeBBox();

	return sdet;
}

DECLARE_DETELEMENT(IP6BeamPipe, create_detector)
