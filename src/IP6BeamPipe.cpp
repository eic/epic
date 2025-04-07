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
		for (xml_coll_t x_zplane_i(x_pipe, _Unicode(zplane)); x_zplane_i; ++x_zplane_i) 
		{
			xml_comp_t x_zplane = x_zplane_i;
			// thickness
			auto wall_thickness = 
				getAttrOrDefault(x_zplane, _Unicode(wall_thickness), 1 * mm);
			auto coating_thickness = 
				getAttrOrDefault(x_zplane, _Unicode(coating_thickness), 30 * um);
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
  	auto create_racetrack_solids = [&](
		xml::Component& x_racetrack)
	{
		double semiCircle_rmin = getAttrOrDefault(x_racetrack, _Unicode(semiCircle_rmin), 0.0);
		double wall_thickness = getAttrOrDefault(x_racetrack, _Unicode(wall_thickness), 0.0);
		double length = getAttrOrDefault(x_racetrack, _Unicode(length), 0.0);
		double rectangle_h = getAttrOrDefault(x_racetrack, _Unicode(rectangle_h), 0.0);
		
		Tube wall_semiCircle(0, semiCircle_rmin +  wall_thickness, length/2.);
		Box wall_rectangle(semiCircle_rmin + wall_thickness, rectangle_h/2., length/2.);

		Tube coating_semiCircle(0, semiCircle_rmin, length/2.);
		Box coating_rectangle(semiCircle_rmin, rectangle_h/2., length/2.);

		Tube vacuum_semiCircle(0, semiCircle_rmin -  wall_thickness, length/2.);
		Box vacuum_rectangle(semiCircle_rmin - wall_thickness, rectangle_h/2., length/2.);

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

		Transform3D zero_tf(RotationZYX(0, 0, 0), Position(0, 0, 0));
		BooleanSolid wall_subtract, coating_subtract, vacuum_subtract;

		wall_subtract = 
			SubtractionSolid("wall_subtract",wall_solid,vacuum_solid,zero_tf);
		coating_subtract = 
			SubtractionSolid("coating_subtract",coating_solid,vacuum_solid,zero_tf);
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

		if(name == "downstream")
		{
			xml::Component racetrack_lepton_c = x_pipe1.child(_Unicode(racetrack_lepton));
			auto racetrack_solids = create_racetrack_solids(racetrack_lepton_c);

			double offset_z = getAttrOrDefault(racetrack_lepton_c, _Unicode(offset_z), 0.0); 
			double length = getAttrOrDefault(racetrack_lepton_c, _Unicode(length), 0.0);
			Transform3D racetrack_tf(RotationZYX(0, 0, 0), Position(0, 0, offset_z + length/2.));

			wall_union = UnionSolid("wall_union",wall_union,std::get<0>(racetrack_solids),racetrack_tf);
			coating_union = UnionSolid("coating_union",coating_union,std::get<1>(racetrack_solids),racetrack_tf);
			vacuum_union = UnionSolid("vacuum_union",vacuum_union,std::get<2>(racetrack_solids),racetrack_tf);
		}

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
			// unite vacuums
			if(getAttrOrDefault<bool>(x_additional_subtraction, _Unicode(unite_vacuum), false))
			{
				vacuum_union = UnionSolid(
					vacuum_union, std::get<2>(additional_polycones), additional_tf);
			}
		}

		BooleanSolid wall, coating, vacuum;

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
			// subtract wall from vacuum
			vacuum = SubtractionSolid(vacuum_union, wall_union);
			// subtract coating from vacuum
			vacuum = SubtractionSolid(vacuum, coating_union);
			// get wall
			wall = wall_union;
			// get coating
			coating = coating_union;
		}

		return std::tuple<Volume, Volume, Volume>(
			{"v_" + name + "_wall", wall, m_Wall},
			{"v_" + name + "_coating", coating, m_Coating},
			{"v_" + name + "_vacuum", vacuum, m_Vacuum}
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

	std::get<0>(volumes_downstream).setVisAttributes(wallVis);
	std::get<1>(volumes_downstream).setVisAttributes(coatingVis);

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
	}

	// -----------------------------
	// final placement
	auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
	pv_assembly.addPhysVolID("system", sdet.id());
	sdet.setPlacement(pv_assembly);
	assembly->GetShape()->ComputeBBox();

	return sdet;
}

DECLARE_DETELEMENT(IP6BeamPipe, create_detector)
