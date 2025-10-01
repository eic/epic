// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Dhevan Gangadharan, Simon Gardner

//==========================================================================
//
// Places a chain of beam pipe segments within and between the beamline magnets.
//
// Approximation used for beam pipes in between magnets.
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>
#include "DD4hep/Shapes.h"
#include "DD4hep/Objects.h"
#include "TGeoTessellated.h"
#include "TGeoManager.h"
#include "TVector3.h"

using namespace std;
using namespace dd4hep;

//- post_d1ef
const double s_post_d1ef_start = 3000 * cm;
const double s_post_d1ef_stop = 2350 * cm;
const double rx_post_d1ef_min = 100/2. * mm;
const double ry_post_d1ef_min = 100/2. * mm;

//- pre_q1ef
const double s_pre_q1ef_start = 2300 * cm;
const double s_pre_q1ef_stop = 1237 * cm;
const double rx_pre_q1ef_min = 126/2. * mm;
const double ry_pre_q1ef_min = 100/2. * mm;

//- q1ef
const double s_q1ef_start = s_pre_q1ef_stop;
const double s_q1ef_stop = 1076 * cm;
const double rx_q1ef_min = rx_pre_q1ef_min;
const double ry_q1ef_min = ry_pre_q1ef_min;

//- post_q1ef
const double s_post_q1ef_start = s_q1ef_stop;
const double s_post_q1ef_stop = 1050 * cm;
const double rx_post_q1ef_min = rx_pre_q1ef_min;
const double ry_post_q1ef_min = ry_pre_q1ef_min;

//- pre_q0ef
const double s_pre_q0ef_start = 720 * cm;
const double s_pre_q0ef_stop = 700 * cm;
const double r_pre_q0ef_min = 54/2. * mm;

//- q0ef
const double s_q0ef_start = s_pre_q0ef_stop;
const double s_q0ef_stop = 580 * cm;
const double r_q0ef_min = r_pre_q0ef_min;

//- post_q0ef
const double s_post_q0ef_start = s_q0ef_stop;
const double s_post_q0ef_stop = 494.556 * cm;
const double r_post_q0ef_min = r_pre_q0ef_min;

//- tapper0
const double s_tapper0_start = s_post_d1ef_stop;
const double s_tapper0_stop = s_pre_q1ef_start;
const double rx_tapper0_min_start = rx_post_d1ef_min;
const double ry_tapper0_min_start = ry_post_d1ef_min;
const double rx_tapper0_min_stop = rx_pre_q1ef_min;
const double ry_tapper0_min_stop = ry_pre_q1ef_min;

//- tapper1
const double s_tapper1_start = s_post_q1ef_stop;
const double s_tapper1_middle = 1000 * cm;
const double s_tapper1_stop = s_pre_q0ef_start;
const double rx_tapper1_min_start = rx_pre_q1ef_min;
const double ry_tapper1_min_start = ry_pre_q1ef_min;
const double rx_tapper1_min_middle = 5.75 * cm;
const double ry_tapper1_min_middle = r_pre_q0ef_min;
const double rx_tapper1_min_stop = r_pre_q0ef_min;
const double ry_tapper1_min_stop = r_pre_q0ef_min;

//- bwd_polycone
const int n_bwd 		= 11; 
const double s_bwd[n_bwd] 	= {
	-20.0 * m, -16.0 * m,  -15.0 * m, -11.4140 * m, -11.3640 * m, -9.3640 * m, 
	-9.2880 * m, -7.4050 * m, -7.3420 * m,  -5.0 * m, -454.449 * cm};
const double d_bwd_min[n_bwd]  	= {  
	0.1 * m,   0.1 * m, 0.1960 * m,   0.1960 * m,   0.1540 * m,  0.1540 * m,  
	0.1286 * m,  0.1286 * m,  0.1114 * m, 0.062 * m, 0.062 * m};

const double screen_gap = 100 * um; // to avaid overlaps
const double coating_thickness = 30 * um;
const double screen_thickness = 1 * mm;

void BuildSectionSolid(	Detector& det,
			std::string name,
			double rx_min, double ry_min,
			double s_start, double s_stop,
			Assembly& assembly)
{
	double halfLength = std::abs(s_start - s_stop) / 2.0;
	double zpos       = s_start - (s_start - s_stop) / 2.0;

	// vacuum-1
	EllipticalTube cut_solid(rx_min, ry_min, halfLength + 1.0*dd4hep::cm);

	EllipticalTube vacuum_solid(rx_min, ry_min, halfLength);
	Volume vacuum_vol(
		name + "_vacuum",
		vacuum_solid,
		det.material("Vacuum"));
	assembly.placeVolume(vacuum_vol, Position(0,0,zpos));

	// coating (subtract vacuum-1)
	EllipticalTube coating_outer(
		rx_min + coating_thickness,
		ry_min + coating_thickness,
		halfLength);
	SubtractionSolid coating_solid(coating_outer, cut_solid);
	Volume coating_vol(
		name + "_coating",
		coating_solid,
		det.material("Copper"));
	assembly.placeVolume(coating_vol, Position(0,0,zpos));
	coating_vol.setVisAttributes(det.visAttributes("BrownVis"));

	// vacuum-2
	EllipticalTube vacuum2_solid(
		rx_min + coating_thickness + screen_gap,
		ry_min + coating_thickness + screen_gap,
		halfLength + 1.0*dd4hep::cm);

	// screen (subtract vacuum-2)
	EllipticalTube screen_outer(
		rx_min + coating_thickness + screen_gap + screen_thickness,
		ry_min + coating_thickness + screen_gap + screen_thickness,
		halfLength);
	SubtractionSolid screen_solid(screen_outer, vacuum2_solid);
	Volume screen_vol(
		name + "_screen",
		screen_solid,
		det.material("StainlessSteelP506"));
	screen_vol.setVisAttributes(det.visAttributes("YellowVis"));
	assembly.placeVolume(screen_vol, Position(0,0,zpos));

	return;
}

// Discretize ellipse into polygon points
std::vector<Position> MakeEllipse(double rx, double ry, double z, int nseg) 
{
	std::vector<Position> pts;
	pts.reserve(nseg);
	for (int i = 0; i < nseg; i++) 
	{
		double phi = 2.0 * M_PI * i / nseg;
		pts.emplace_back(rx * cos(phi), ry * sin(phi), z);
	}
	return pts;
}

TGeoTessellated* CreatePolygonHollowTube2(
	double rx1_in, double ry1_in, double z1,
	double rx2_in, double ry2_in, double z2,
	double thickness)
{
	auto* tess = new TGeoTessellated();
	const int nSegments = 360;

	std::vector<TVector3> outer1, outer2, inner1, inner2;

	// Generate vertices for inner and outer ellipses
	for (int i = 0; i < nSegments; ++i) 
	{
		double phi = 2.0 * M_PI * i / nSegments;

		// Inner ellipses
		inner1.emplace_back(rx1_in * cos(phi), ry1_in * sin(phi), z1);
		inner2.emplace_back(rx2_in * cos(phi), ry2_in * sin(phi), z2);

		// Outer ellipses = inner + thickness
		outer1.emplace_back(
			(rx1_in + thickness) * cos(phi),
			(ry1_in + thickness) * sin(phi), z1);
		outer2.emplace_back((
			rx2_in + thickness) * cos(phi),
			(ry2_in + thickness) * sin(phi), z2);
	}

	auto makeVertex = [](const TVector3& v) 
	{
		return TGeoTessellated::Vertex_t{v.X(), v.Y(), v.Z()};
	};

	// --- Outer wall ---
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(makeVertex(outer1[i]), makeVertex(outer2[i]), makeVertex(outer2[j]));
		tess->AddFacet(makeVertex(outer1[i]), makeVertex(outer2[j]), makeVertex(outer1[j]));
	}

	// --- Inner wall ---
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(inner2[j]), makeVertex(inner2[i]));
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(inner1[j]), makeVertex(inner2[j]));
	}

	// --- End cap at z1 ---
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(outer1[i]), makeVertex(outer1[j]));
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(outer1[j]), makeVertex(inner1[j]));
	}

	// --- End cap at z2 ---
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(makeVertex(outer2[i]), makeVertex(inner2[i]), makeVertex(inner2[j]));
		tess->AddFacet(makeVertex(outer2[i]), makeVertex(inner2[j]), makeVertex(outer2[j]));
	}

	tess->CloseShape();
	return tess;
}

TGeoTessellated* CreatePolygonFilledTube2(
	double rx1, double ry1, double z1,
	double rx2, double ry2, double z2)
{
	auto* tess = new TGeoTessellated();
	const int nSegments = 360;

	std::vector<TVector3> bottom, top;

	// Generate vertices for bottom and top ellipses
	for (int i = 0; i < nSegments; ++i) 
	{
		double phi = 2.0 * M_PI * i / nSegments;

		bottom.emplace_back(rx1 * cos(phi), ry1 * sin(phi), z1);
		top.emplace_back(rx2 * cos(phi), ry2 * sin(phi), z2);
	}

	auto makeVertex = [](const TVector3& v) 
	{
		return TGeoTessellated::Vertex_t{v.X(), v.Y(), v.Z()};
	};

	// --- Side walls ---
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;

		tess->AddFacet(makeVertex(bottom[i]), makeVertex(top[i]), makeVertex(top[j]));
		tess->AddFacet(makeVertex(bottom[i]), makeVertex(top[j]), makeVertex(bottom[j]));
	}

	// --- End cap at z1 (bottom) ---
	TVector3 centerBottom(0, 0, z1);
	auto centerB = makeVertex(centerBottom);
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(centerB, makeVertex(bottom[i]), makeVertex(bottom[j]));
	}

	// --- End cap at z2 (top) ---
	TVector3 centerTop(0, 0, z2);
	auto centerT = makeVertex(centerTop);
	for (int i = 0; i < nSegments; ++i) 
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(centerT, makeVertex(top[j]), makeVertex(top[i]));
	}

	tess->CloseShape();
    	return tess;
}

TGeoTessellated* CreatePolygonHollowTube3(
        double rx1_in, double ry1_in, double z1,
        double rx2_in, double ry2_in, double z2,
        double rx3_in, double ry3_in, double z3,
        double thickness)
{
	auto* tess = new TGeoTessellated();
	const int nSegments = 360;

	std::vector<TVector3> outer1, outer2, outer3;
	std::vector<TVector3> inner1, inner2, inner3;

	// Generate vertices for inner and outer ellipses
	for (int i = 0; i < nSegments; ++i)
	{
		double phi = 2.0 * M_PI * i / nSegments;

		// Inner ellipses
		inner1.emplace_back(rx1_in * cos(phi), ry1_in * sin(phi), z1);
		inner2.emplace_back(rx2_in * cos(phi), ry2_in * sin(phi), z2);
		inner3.emplace_back(rx3_in * cos(phi), ry3_in * sin(phi), z3);

		// Outer ellipses = inner + thickness
		outer1.emplace_back((rx1_in + thickness) * cos(phi),
			(ry1_in + thickness) * sin(phi), z1);
		outer2.emplace_back((rx2_in + thickness) * cos(phi),
			(ry2_in + thickness) * sin(phi), z2);
		outer3.emplace_back((rx3_in + thickness) * cos(phi),
			(ry3_in + thickness) * sin(phi), z3);
	}

	auto makeVertex = [](const TVector3& v)
	{
		return TGeoTessellated::Vertex_t{v.X(), v.Y(), v.Z()};
	};

	// --- Outer walls ---
	for (int i = 0; i < nSegments; ++i)
	{
		int j = (i + 1) % nSegments;
		// z1 → z2
		tess->AddFacet(makeVertex(outer1[i]), makeVertex(outer2[i]), makeVertex(outer2[j]));
		tess->AddFacet(makeVertex(outer1[i]), makeVertex(outer2[j]), makeVertex(outer1[j]));

		// z2 → z3
		tess->AddFacet(makeVertex(outer2[i]), makeVertex(outer3[i]), makeVertex(outer3[j]));
		tess->AddFacet(makeVertex(outer2[i]), makeVertex(outer3[j]), makeVertex(outer2[j]));
	}

	// --- Inner walls ---
	for (int i = 0; i < nSegments; ++i)
	{
		int j = (i + 1) % nSegments;

		// z1 → z2
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(inner2[j]), makeVertex(inner2[i]));
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(inner1[j]), makeVertex(inner2[j]));

		// z2 → z3
		tess->AddFacet(makeVertex(inner2[i]), makeVertex(inner3[j]), makeVertex(inner3[i]));
		tess->AddFacet(makeVertex(inner2[i]), makeVertex(inner2[j]), makeVertex(inner3[j]));
	}

	// --- End caps ---
	for (int i = 0; i < nSegments; ++i)
	{
		int j = (i + 1) % nSegments;
		// Cap at z1
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(outer1[i]), makeVertex(outer1[j]));
		tess->AddFacet(makeVertex(inner1[i]), makeVertex(outer1[j]), makeVertex(inner1[j]));
		// Cap at z3
		tess->AddFacet(makeVertex(outer3[i]), makeVertex(inner3[i]), makeVertex(inner3[j]));
		tess->AddFacet(makeVertex(outer3[i]), makeVertex(inner3[j]), makeVertex(outer3[j]));
	}

	tess->CloseShape();
	return tess;
}

TGeoTessellated* CreatePolygonFilledTube3(
	double rx1, double ry1, double z1,
	double rx2, double ry2, double z2,
	double rx3, double ry3, double z3)
{
	auto* tess = new TGeoTessellated();
	const int nSegments = 360;

	std::vector<TVector3> bottom, middle, top;

	// Generate vertices for three elliptical planes
	for (int i = 0; i < nSegments; ++i)
	{
		double phi = 2.0 * M_PI * i / nSegments;

		bottom.emplace_back(rx1 * cos(phi), ry1 * sin(phi), z1);
		middle.emplace_back(rx2 * cos(phi), ry2 * sin(phi), z2);
		top.emplace_back(rx3 * cos(phi), ry3 * sin(phi), z3);
	}

	auto makeVertex = [](const TVector3& v)
	{
		return TGeoTessellated::Vertex_t{v.X(), v.Y(), v.Z()};
	};

	// --- Side walls between bottom → middle and middle → top ---
	for (int i = 0; i < nSegments; ++i)
	{
		int j = (i + 1) % nSegments;

		// bottom → middle
		tess->AddFacet(makeVertex(bottom[i]), makeVertex(middle[i]), makeVertex(middle[j]));
		tess->AddFacet(makeVertex(bottom[i]), makeVertex(middle[j]), makeVertex(bottom[j]));

		// middle → top
		tess->AddFacet(makeVertex(middle[i]), makeVertex(top[i]), makeVertex(top[j]));
		tess->AddFacet(makeVertex(middle[i]), makeVertex(top[j]), makeVertex(middle[j]));
	}

	// --- End cap at z1 (bottom) ---
	TVector3 centerBottom(0, 0, z1);
	auto centerB = makeVertex(centerBottom);
	for (int i = 0; i < nSegments; ++i)
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(centerB, makeVertex(bottom[i]), makeVertex(bottom[j]));
	}

	// --- End cap at z3 (top) ---
	TVector3 centerTop(0, 0, z3);
	auto centerT = makeVertex(centerTop);
	for (int i = 0; i < nSegments; ++i)
	{
		int j = (i + 1) % nSegments;
		tess->AddFacet(centerT, makeVertex(top[j]), makeVertex(top[i]));
	}

	tess->CloseShape();
	return tess;
}

void BuildTapper1(Detector& det, Assembly& assembly) 
{
	// --- vacuum ---
	TGeoTessellated* rawShape_vacuum = CreatePolygonFilledTube3(
		rx_tapper1_min_start - screen_gap,
		ry_tapper1_min_start - screen_gap, s_tapper1_start,
		rx_tapper1_min_middle - screen_gap,
		ry_tapper1_min_middle - screen_gap, s_tapper1_middle,
		rx_tapper1_min_stop - screen_gap,
		ry_tapper1_min_stop - screen_gap, s_tapper1_stop);
	Solid tapper1_vacuum_solid(rawShape_vacuum);

	Volume tapper1_vacuum_vol("tapper1_vacuum", tapper1_vacuum_solid, det.material("Vacuum"));
	assembly.placeVolume(tapper1_vacuum_vol, Position(0,0,0));

	// --- coating ---
	TGeoTessellated* rawShape_coating = CreatePolygonHollowTube3(
		rx_tapper1_min_start,
		ry_tapper1_min_start, s_tapper1_start,
		rx_tapper1_min_middle,
		ry_tapper1_min_middle, s_tapper1_middle,
		rx_tapper1_min_stop,
		ry_tapper1_min_stop, s_tapper1_stop,
		coating_thickness);
	Solid tapper1_coating_solid(rawShape_coating);

	Volume tapper1_coating_vol("tapper1_coating", tapper1_coating_solid, det.material("Copper"));
	tapper1_coating_vol.setVisAttributes(det.visAttributes("BrownVis"));
	assembly.placeVolume(tapper1_coating_vol, Position(0,0,0));

	// --- screen ---
	TGeoTessellated* rawShape_screen = CreatePolygonHollowTube3(
		rx_tapper1_min_start + coating_thickness + screen_gap,
		ry_tapper1_min_start + coating_thickness + screen_gap, s_tapper1_start,
		rx_tapper1_min_middle + coating_thickness + screen_gap,
		ry_tapper1_min_middle + coating_thickness + screen_gap, s_tapper1_middle,
		rx_tapper1_min_stop + coating_thickness + screen_gap,
		ry_tapper1_min_stop + coating_thickness + screen_gap, s_tapper1_stop,
		screen_thickness);
	Solid tapper1_screen_solid(rawShape_screen);

	Volume tapper1_screen_vol("tapper1_screen", tapper1_screen_solid, det.material("StainlessSteelP506"));
	tapper1_screen_vol.setVisAttributes(det.visAttributes("YellowVis"));
	assembly.placeVolume(tapper1_screen_vol, Position(0,0,0));

	return;
}

void BuildTapper0(Detector& det, Assembly& assembly) 
{
	// --- vacuum ---
	TGeoTessellated* rawShape_vacuum = CreatePolygonFilledTube2(
		rx_tapper0_min_start - screen_gap,
		ry_tapper0_min_start - screen_gap, s_tapper0_start,
		rx_tapper0_min_stop - screen_gap,
		ry_tapper0_min_stop - screen_gap, s_tapper0_stop);

	Solid tapper0_vacuum_solid(rawShape_vacuum);

	Volume tapper0_vacuum_vol("tapper0_vacuum", tapper0_vacuum_solid, det.material("Vacuum"));
	assembly.placeVolume(tapper0_vacuum_vol, Position(0,0,0));

	// --- coating ---
	TGeoTessellated* rawShape_coating = CreatePolygonHollowTube2(
		rx_tapper0_min_start,
		ry_tapper0_min_start, s_tapper0_start,
		rx_tapper0_min_stop,
		ry_tapper0_min_stop, s_tapper0_stop,
		coating_thickness);

	Solid tapper0_coating_solid(rawShape_coating);

	Volume tapper0_coating_vol("tapper0_coating", tapper0_coating_solid, det.material("Copper"));
	tapper0_coating_vol.setVisAttributes(det.visAttributes("BrownVis"));
	assembly.placeVolume(tapper0_coating_vol, Position(0,0,0));

	// --- screen ---
	TGeoTessellated* rawShape_screen = CreatePolygonHollowTube2(
		rx_tapper0_min_start + coating_thickness + screen_gap,
		ry_tapper0_min_start + coating_thickness + screen_gap, s_tapper0_start,
		rx_tapper0_min_stop + coating_thickness + screen_gap,
		ry_tapper0_min_stop + coating_thickness + screen_gap, s_tapper0_stop,
		screen_thickness);
	Solid tapper0_screen_solid(rawShape_screen);

	Volume tapper0_screen_vol("tapper0_screen", tapper0_screen_solid, det.material("StainlessSteelP506"));
	tapper0_screen_vol.setVisAttributes(det.visAttributes("YellowVis"));
	assembly.placeVolume(tapper0_screen_vol, Position(0,0,0));

	return;
}

void BuildPolyconeSection(
	Detector& det, 
	std::string name,
	const double* pos, 
	const double* dia,
	std::size_t num,
	Assembly& assembly)
{
	std::vector<double> r_min(num), r_max(num), z_vec(num);

	// --- Vacuum ---
	for (size_t i=0; i<num; ++i) 
	{
		r_min[i] = 0;
		r_max[i] = dia[i]/2.;
		z_vec[i] = pos[i];
	}
	Solid vacuum_solid = Polycone(0, 2*M_PI, r_min, r_max, z_vec);
	Volume vacuum_vol(name + "_vacuum", vacuum_solid, det.material("Vacuum"));
	assembly.placeVolume(vacuum_vol, Position(0,0,0));

	// --- Coating ---
	for (size_t i=0; i<num; ++i) 
	{
		r_min[i] = dia[i]/2.;
		r_max[i] = dia[i]/2. + coating_thickness;
	}
	Solid coating_solid = Polycone(0, 2*M_PI, r_min, r_max, z_vec);
	Volume coating_vol(name + "_coating", coating_solid, det.material("Copper"));
	coating_vol.setVisAttributes(det.visAttributes("BrownVis"));
	assembly.placeVolume(coating_vol, Position(0,0,0));

	// --- Screen ---
	for (size_t i=0; i<num; ++i) 
	{
		r_min[i] = dia[i]/2. + coating_thickness;
		r_max[i] = r_min[i] + screen_thickness;
	}
	Solid screen_solid = Polycone(0, 2*M_PI, r_min, r_max, z_vec);
	Volume screen_vol(name + "_screen", screen_solid, det.material("StainlessSteelP506"));
	screen_vol.setVisAttributes(det.visAttributes("YellowVis"));
	assembly.placeVolume(screen_vol, Position(0,0,0));

	return;
}

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) 
{
	using namespace ROOT::Math;
	xml_det_t x_det = e;
	string det_name = x_det.nameStr();
	DetElement sdet(det_name, x_det.id());
	Assembly assembly(det_name + "_assembly");

	if(det_name == "Pipe_cen_to_pos") // forward side
	{
		//- post_d1ef
		BuildSectionSolid(det,"post_d1ef",rx_post_d1ef_min,ry_post_d1ef_min,s_post_d1ef_start,s_post_d1ef_stop,assembly);

		//- tapper0
		BuildTapper0(det,assembly);

		//- pre_q1ef
		BuildSectionSolid(det,"pre_q1ef",rx_pre_q1ef_min,ry_pre_q1ef_min,s_pre_q1ef_start,s_pre_q1ef_stop,assembly);

		//- q1ef
		BuildSectionSolid(det,"q1ef",rx_q1ef_min,ry_q1ef_min,s_q1ef_start,s_q1ef_stop,assembly);

		//- post_q1ef
		BuildSectionSolid(det,"post_q1ef",rx_post_q1ef_min,ry_post_q1ef_min,s_post_q1ef_start,s_post_q1ef_stop,assembly);

		//- pre_q0ef
		BuildSectionSolid(det,"pre_q0ef",r_pre_q0ef_min,r_pre_q0ef_min,s_pre_q0ef_start,s_pre_q0ef_stop,assembly);

		//- q0ef
		BuildSectionSolid(det,"q0ef",r_q0ef_min,r_q0ef_min,s_q0ef_start,s_q0ef_stop,assembly);

		//- post_q0ef
		BuildSectionSolid(det,"post_q0ef",r_post_q0ef_min,r_post_q0ef_min,s_post_q0ef_start,s_post_q0ef_stop,assembly);

		//- tapper1
		BuildTapper1(det,assembly);
	}
	else if(det_name == "Pipe_cen_to_neg") // rear side
	{
		//- pre_q1er_post_q2er
		BuildPolyconeSection(det,"pre_q1er_post_q2er",s_bwd,d_bwd_min,n_bwd,assembly);
	}

	// Final placement
	auto pv_assembly =
		det.pickMotherVolume(sdet).placeVolume(assembly, Position(0.0, 0.0, 0.0));
	sdet.setPlacement(pv_assembly);
	assembly->GetShape()->ComputeBBox();

	return sdet;
}

DECLARE_DETELEMENT(ElectronBeamPipeChain, create_detector)
