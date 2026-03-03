// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Justin Chan

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;


static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) 
{
	using namespace ROOT::Math;
	xml_det_t x_det = e;
	string det_name = x_det.nameStr();
	DetElement sdet(det_name, x_det.id());
	Assembly assembly(det_name + "_assembly");
	Material m_Iron   = det.material("Iron");
	Material m_Copper = det.material("Copper");
	Material m_Tungsten = det.material("TungstenDens25");
	const string vis1 = getAttrOrDefault<string>(x_det, _Unicode(vis_name1), "AnlGreen");
	const string vis2 = getAttrOrDefault<string>(x_det, _Unicode(vis_name2), "AnlRed");
	const string vis3 = getAttrOrDefault<string>(x_det, _Unicode(vis_name3), "AnlGray");
	const string vis4 = getAttrOrDefault<string>(x_det, _Unicode(vis_name4), "AnlBlue");

	// Get position
	xml::Component pos = x_det.child(_Unicode(position));
	double x = pos.attr<double>(_Unicode(x));
	double y = pos.attr<double>(_Unicode(y));
	double z = pos.attr<double>(_Unicode(z));

	//--------------------------------------------------------------------------------//
	// Create the outer yoke box
	xml::Component dim_1 = x_det.child(_Unicode(dimensions_yoke_box));
	double height_1	= dim_1.attr<double>(_Unicode(y));
	double width_1	= dim_1.attr<double>(_Unicode(x));
	double depth_1	= dim_1.attr<double>(_Unicode(z));

	Box yoke_box (width_1 / 2., height_1 / 2., depth_1 / 2.);

  	// Create the yoke aperture
	xml::Component dim_2 = x_det.child(_Unicode(dimensions_yoke_aperture));
	std::vector<double> px, py;
	for (xml_coll_t p(dim_2, _U(point)); p; ++p) {
		xml_dim_t pt(p);
		px.push_back(pt.x());
		py.push_back(pt.y());
	}

	// Extrude polygon along Z with two sections (-dz, +dz), no offsets, scale=1
	std::vector<double> sec_z  = { -depth_1, +depth_1 };
	std::vector<double> sec_x  = {  0.0, 0.0 };
	std::vector<double> sec_y  = {  0.0, 0.0 };
	std::vector<double> zscale = {  1.0, 1.0 };

	Solid yoke_aperture = ExtrudedPolygon("yoke_aperture", px, py, sec_z, sec_x, sec_y, zscale);

	// Subtract the volume of the yoke aperture from the yoke box
	BooleanSolid yoke_body = SubtractionSolid(yoke_box, yoke_aperture);

	Volume v_yoke_body(det_name + "_vol_yoke_body", yoke_body, m_Iron);
	sdet.setAttributes(det, v_yoke_body, x_det.regionStr(), x_det.limitsStr(), vis1);
	assembly.placeVolume(v_yoke_body, Position(x, y, z));

	//--------------------------------------------------------------------------------//
	// Create a coil as a union of rectangle and cylinders
	xml::Component dim_3 = x_det.child(_Unicode(dimensions_coil_box));
	double height_3	= dim_3.attr<double>(_Unicode(y));
	double width_3	= dim_3.attr<double>(_Unicode(x));
	double depth_3	= dim_3.attr<double>(_Unicode(z));
	double round_3	= dim_3.attr<double>(_Unicode(r));
	double offset_3	= dim_3.attr<double>(_Unicode(dx));

	// Two overlapping boxes (half-lengths!)
	// boxA: W  x (H-2R) x T
	Solid boxA_3 = Box(width_3/2., height_3/2. - round_3, depth_3/2.);
	// boxB: (W-2R) x H x T
	Solid boxB_3 = Box(width_3/2. - round_3, height_3/2., depth_3/2.);
	
	Solid coil_box = UnionSolid("coil_box", boxA_3, boxB_3); 

	Solid coil_box_corner = Tube(0.0, round_3, depth_3/2.);

	const double cx = width_3/2. - round_3;
	const double cy = height_3/2. - round_3;

	// Unite with rounded corners
	coil_box = UnionSolid("coil_box_u1", coil_box, coil_box_corner, Position(+cx, +cy, 0.0));
	coil_box = UnionSolid("coil_box_u2", coil_box, coil_box_corner, Position(-cx, +cy, 0.0));
	coil_box = UnionSolid("coil_box_u3", coil_box, coil_box_corner, Position(-cx, -cy, 0.0));
	coil_box = UnionSolid("coil_box_u4", coil_box, coil_box_corner, Position(+cx, -cy, 0.0));

	// Create the coil aperture
	xml::Component dim_4 = x_det.child(_Unicode(dimensions_coil_aperture));
	double height_4	= dim_4.attr<double>(_Unicode(y));
	double width_4	= dim_4.attr<double>(_Unicode(x));
	double depth_4	= dim_4.attr<double>(_Unicode(z));

	Box coil_aperture (width_4 / 2., height_4 / 2., depth_4 / 2.);

	// Subtract the volume of the coil aperture from the coil box
	BooleanSolid coil_body = SubtractionSolid(coil_box, coil_aperture);

	Volume v_coil1_body(det_name + "_vol_coil1_body", coil_body, m_Copper);
	sdet.setAttributes(det, v_coil1_body, x_det.regionStr(), x_det.limitsStr(), vis2);
	assembly.placeVolume(v_coil1_body, Transform3D(RotationZYX(0, M_PI_2, 0), Position(x-offset_3, y, z)));

	Volume v_coil2_body(det_name + "_vol_coil2_body", coil_body, m_Copper);
	sdet.setAttributes(det, v_coil2_body, x_det.regionStr(), x_det.limitsStr(), vis2);
	assembly.placeVolume(v_coil2_body, Transform3D(RotationZYX(0, M_PI_2, 0), Position(x+offset_3, y, z)));

	//--------------------------------------------------------------------------------//
	// Create shielding
	if (x_det.hasChild(_Unicode(shielding)))
	{
		xml::Component shielding = x_det.child(_Unicode(shielding));
		double shielding_DR = shielding.attr<double>(_Unicode(DR));
		double shielding_IR = shielding.attr<double>(_Unicode(IR));
		double shielding_DZ = shielding.attr<double>(_Unicode(DZ));

		Solid shielding_tube = Tube(shielding_IR, shielding_IR + shielding_DR, shielding_DZ/2.);

		Volume v_shielding_body(det_name + "_vol_shielding_body", shielding_tube, m_Tungsten);
		sdet.setAttributes(det, v_shielding_body, x_det.regionStr(), x_det.limitsStr(), vis4);
		assembly.placeVolume(v_shielding_body, Transform3D(RotationZYX(0, 0, 0), Position(x, y, z)));
	}

	//--------------------------------------------------------------------------------//
	// Final placement
	auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(
	assembly, Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0, 0, 0)));

	sdet.setPlacement(pv_assembly);

	assembly->GetShape()->ComputeBBox();

	return sdet;
}

DECLARE_DETELEMENT(LumiMagnets, create_detector)
