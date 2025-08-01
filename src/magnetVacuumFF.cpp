// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Alex Jentsch

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <XML/Helper.h>

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

static double getRotatedZ(double z, double x, double angle);
static double getRotatedX(double z, double x, double angle);

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) {

  xml_det_t x_det = e;
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly assembly(det_name + "_assembly");
  Material m_Vac  = det.material("Vacuum");
  string vis_name = x_det.visStr();

  PlacedVolume pv_assembly;

  //----------------------------------------------
  // Starting point is only the magnet centers,
  // lengths, rotations, and radii -->
  // everything else calculated internally to
  // make it easier to update later.
  //----------------------------------------------

  bool makeIP_B0pfVacuum   = true; //This is for the special gap location between IP and b0pf
  bool make_B2pf_EW_Vacuum = true; //This is for the gap after b2pf

  //information for actual FF magnets, with magnet centers as reference
  vector<double> radii_magnet;
  vector<double> lengths_magnet;
  vector<double> rotation_magnet;
  vector<double> x_elem_magnet;
  vector<double> y_elem_magnet;
  vector<double> z_elem_magnet;

  //calculated entrance/exit points of FF magnet
  vector<double> x_beg;
  vector<double> z_beg;
  vector<double> x_end;
  vector<double> z_end;

  //calculated center of gap regions between magnets, rotation, and length
  vector<double> angle_elem_gap;
  vector<double> z_gap;
  vector<double> x_gap;
  vector<double> length_gap;

  //storage elements for CutTube geometry element used for gaps
  vector<double> inRadius;
  vector<double> outRadius;
  vector<double> nxLow;
  vector<double> nyLow;
  vector<double> nzLow;
  vector<double> nxHigh;
  vector<double> nyHigh;
  vector<double> nzHigh;
  vector<double> phi_initial;
  vector<double> phi_final;

  for (xml_coll_t c(x_det, _U(element)); c; ++c) {

    xml_dim_t pos      = c.child(_U(placement));
    double pos_x       = pos.x();
    double pos_y       = pos.y();
    double pos_z       = pos.z();
    double pos_theta   = pos.attr<double>(_U(theta));
    xml_dim_t dims     = c.child(_U(dimensions)); //dimensions();
    double dim_z       = dims.z();
    xml_dim_t aperture = c.child(_Unicode(aperture));
    double app_r       = aperture.r();

    radii_magnet.push_back(app_r);        // cm
    lengths_magnet.push_back(dim_z);      //cm
    rotation_magnet.push_back(pos_theta); // radians
    x_elem_magnet.push_back(pos_x * dd4hep::cm);
    y_elem_magnet.push_back(pos_y * dd4hep::cm);
    z_elem_magnet.push_back(pos_z * dd4hep::cm);
  }

  int numMagnets = radii_magnet.size(); //number of actual FF magnets between IP and FF detectors
  int numGaps =
      numMagnets -
      2; //number of gaps between magnets (excluding the IP to B0pf transition -- special case, and the gao after B1apf)

  //-------------------------------------------
  // override numbers for the first element -->
  // doesn't use the actual B0pf geometry!!!
  // -->it's based on the B0 beam pipe
  // this needs to be fixed later to read-in
  // that beam pipe geometry
  //-------------------------------------------

  radii_magnet[0]    = 2.9;                 // cm
  lengths_magnet[0]  = 120.0;               // cm
  rotation_magnet[0] = -0.025;              // radians
  x_elem_magnet[0]   = 640.0 * sin(-0.025); // cm
  y_elem_magnet[0]   = 0.0;                 // cm
  z_elem_magnet[0]   = 640.0 * cos(-0.025); // cm

  //-------------------------------------------
  //calculate entrance/exit points of magnets
  //-------------------------------------------

  for (int i = 0; i < numMagnets; i++) {

    // need to use the common coordinate system -->
    // use x = z, and y = x to make things easier

    z_beg.push_back(getRotatedZ(-0.5 * lengths_magnet[i], 0.0, rotation_magnet[i]) +
                    z_elem_magnet[i]);
    z_end.push_back(getRotatedZ(0.5 * lengths_magnet[i], 0.0, rotation_magnet[i]) +
                    z_elem_magnet[i]);
    x_beg.push_back(getRotatedX(-0.5 * lengths_magnet[i], 0.0, rotation_magnet[i]) +
                    x_elem_magnet[i]);
    x_end.push_back(getRotatedX(0.5 * lengths_magnet[i], 0.0, rotation_magnet[i]) +
                    x_elem_magnet[i]);
  }

  //------------------------------------------
  // this part is a bit ugly for now -
  // it's to make the vacuum volume between the
  // end of the IP beam pipe and the beginning of
  // beginning of the B0pf magnet
  //
  // -->the volume will be calculated at the end
  //-------------------------------------------

  double endOfCentralBeamPipe_z =
      494.556 * dd4hep::cm +
      2. * dd4hep::mm; //extracted from central_beampipe.xml, line 112 + offset to avoid overlaps
  double diameterReduce = 11.0 * dd4hep::cm; //size reduction to avoid overlap with electron pipe
  double vacuumDiameterEntrance =
      25.792 * dd4hep::cm - diameterReduce; //extracted from central_beampipe.xml, line 64
  double vacuumDiameterExit =
      17.4 * dd4hep::cm; //15mrad @ entrance to magnet to not overlap electron magnet
  double crossingAngle          = -0.025; //radians
  double endOfCentralBeamPipe_x = endOfCentralBeamPipe_z * crossingAngle;

  //-----------------------------------------------
  //calculate gap region center, length, and angle
  //-----------------------------------------------

  for (int i = 1; i < numMagnets; i++) {

    angle_elem_gap.push_back((x_beg[i] - x_end[i - 1]) / (z_beg[i] - z_end[i - 1]));
    length_gap.push_back(sqrt(pow(z_beg[i] - z_end[i - 1], 2) + pow(x_beg[i] - x_end[i - 1], 2)));
    z_gap.push_back(z_end[i - 1] + 0.5 * length_gap[i - 1] * cos(angle_elem_gap[i - 1]));
    x_gap.push_back(x_end[i - 1] + 0.5 * length_gap[i - 1] * sin(angle_elem_gap[i - 1]));
  }

  //-----------------------------------------------
  // fill CutTube storage elements
  //-----------------------------------------------

  for (int gapIdx = 0; gapIdx < numGaps; gapIdx++) {

    inRadius.push_back(0.0);
    outRadius.push_back(radii_magnet[gapIdx]);
    phi_initial.push_back(0.0);
    phi_final.push_back(2 * M_PI);
    nxLow.push_back(-(length_gap[gapIdx] / 2.0) *
                    sin(rotation_magnet[gapIdx] - angle_elem_gap[gapIdx]));
    nyLow.push_back(0.0);
    nzLow.push_back(-(length_gap[gapIdx] / 2.0) *
                    cos(rotation_magnet[gapIdx] - angle_elem_gap[gapIdx]));
    nxHigh.push_back((length_gap[gapIdx] / 2.0) *
                     sin(rotation_magnet[gapIdx + 1] - angle_elem_gap[gapIdx]));
    nyHigh.push_back(0.0);
    nzHigh.push_back((length_gap[gapIdx] / 2.0) *
                     cos(rotation_magnet[gapIdx + 1] - angle_elem_gap[gapIdx]));
  }

  //-----------------------
  // inside magnets
  //-----------------------

  for (int pieceIdx = 0; pieceIdx < numMagnets; pieceIdx++) {

    std::string piece_name = Form("MagnetVacuum%d", pieceIdx);

    Tube magnetPiece(piece_name, 0.0, radii_magnet[pieceIdx], lengths_magnet[pieceIdx] / 2);
    Volume vpiece(piece_name, magnetPiece, m_Vac);
    sdet.setAttributes(det, vpiece, x_det.regionStr(), x_det.limitsStr(), vis_name);

    assembly.placeVolume(vpiece,
                         Transform3D(RotationY(rotation_magnet[pieceIdx]),
                                     Position(x_elem_magnet[pieceIdx], y_elem_magnet[pieceIdx],
                                              z_elem_magnet[pieceIdx])));
  }

  //--------------------------
  //between magnets
  //--------------------------

  for (int pieceIdx = numMagnets; pieceIdx < numGaps + numMagnets; pieceIdx++) {

    int correctIdx = pieceIdx - numMagnets;

    std::string piece_name = Form("GapVacuum%d", correctIdx);

    CutTube gapPiece(piece_name, inRadius[correctIdx], outRadius[correctIdx],
                     length_gap[correctIdx] / 2, phi_initial[correctIdx], phi_final[correctIdx],
                     nxLow[correctIdx], nyLow[correctIdx], nzLow[correctIdx], nxHigh[correctIdx],
                     nyHigh[correctIdx], nzHigh[correctIdx]);

    Volume vpiece(piece_name, gapPiece, m_Vac);
    sdet.setAttributes(det, vpiece, x_det.regionStr(), x_det.limitsStr(), vis_name);

    assembly.placeVolume(vpiece, Transform3D(RotationY(angle_elem_gap[correctIdx]),
                                             Position(x_gap[correctIdx], 0.0, z_gap[correctIdx])));
  }

  //--------------------------------------------------------------
  //make and place vacuum volume to connect IP beam pipe to B0pf
  //--------------------------------------------------------------

  if (makeIP_B0pfVacuum) {

    double specialGapLength = sqrt(pow(z_beg[0] - endOfCentralBeamPipe_z, 2) +
                                   pow(x_beg[0] - endOfCentralBeamPipe_x, 2)) -
                              0.1;
    double specialGap_z = 0.5 * specialGapLength * cos(crossingAngle) + endOfCentralBeamPipe_z;
    double specialGap_x = 0.5 * specialGapLength * sin(crossingAngle) + endOfCentralBeamPipe_x;

    std::string piece_name = Form("GapVacuum%d", numGaps + numMagnets);

    ConeSegment specialGap(piece_name, specialGapLength / 2, 0.0, vacuumDiameterEntrance / 2, 0.0,
                           vacuumDiameterExit / 2, 40 * deg, (360 - 40) * deg);

    Volume specialGap_v(piece_name, specialGap, m_Vac);
    sdet.setAttributes(det, specialGap_v, x_det.regionStr(), x_det.limitsStr(), vis_name);

    assembly.placeVolume(specialGap_v, Transform3D(RotationY(crossingAngle),
                                                   Position(specialGap_x, 0.0, specialGap_z)));
  }

  //----------------------------------------------------

  //--------------------------------------------------------------
  //make and place vacuum volume after the FF detector array up to end of the world
  //--------------------------------------------------------------
  if (make_B2pf_EW_Vacuum) {

    int pieceIdx           = numMagnets - 1; // last B2PF magnet
    std::string piece_name = Form("GapVacuum%d", numGaps + numMagnets + 1);
    double endGapLength    = (10000.0 - z_end[pieceIdx]) / cos(rotation_magnet[pieceIdx]);
    endGapLength =
        endGapLength -
        4 * radii_magnet[pieceIdx] *
            tan(-rotation_magnet[pieceIdx]); // shift to keep the tube inside the physical volume
    double endGap_z = 0.5 * endGapLength * cos(rotation_magnet[pieceIdx]) + z_end[pieceIdx];
    double endGap_x = 0.5 * endGapLength * sin(rotation_magnet[pieceIdx]) + x_end[pieceIdx];

    Tube vacuum_endWorld(piece_name, 0.0, 4 * radii_magnet[pieceIdx],
                         endGapLength / 2); // make larger tube than inner magnet radius
    Volume vpiece(piece_name, vacuum_endWorld, m_Vac);
    sdet.setAttributes(det, vpiece, x_det.regionStr(), x_det.limitsStr(),
                       "InvisibleNoDaughters"); // make invisible instead of AnlBlue

    assembly.placeVolume(vpiece, Transform3D(RotationY(rotation_magnet[pieceIdx]),
                                             Position(endGap_x, 0.0, endGap_z)));
  }
  //----------------------------------------------------

  pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  pv_assembly.addPhysVolID("system", x_det.id()).addPhysVolID("barrel", 1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();

  return sdet;
}

double getRotatedZ(double z, double x, double angle) { return z * cos(angle) - x * sin(angle); }

double getRotatedX(double z, double x, double angle) { return z * sin(angle) + x * cos(angle); }

DECLARE_DETELEMENT(magnetElementInnerVacuum, create_detector)
