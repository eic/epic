// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

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
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Steel  = det.material("StainlessSteel");
  Material   m_Vacuum = det.material("Vacuum");
  string     vis_name = x_det.visStr();

  xml::Component box_dim = x_det.child(_Unicode(dimensions));
  double         height  = box_dim.attr<double>(_Unicode(height));
  double         width   = box_dim.attr<double>(_Unicode(width));
  double         depth   = box_dim.attr<double>(_Unicode(depth));

  xml::Component col_XS    = x_det.child(_Unicode(collimator));
  double         colHeight = col_XS.attr<double>(_Unicode(height));
  double         colWidth  = col_XS.attr<double>(_Unicode(width));
  double         colXOff   = col_XS.attr<double>(_Unicode(xOff));

  xml::Component box_place = x_det.child(_Unicode(placement));
  double         zOff      = box_place.attr<double>(_Unicode(z));

  Box box_steel(width, height, depth);
  Box box_vacuum(colWidth, colHeight, depth);

  Volume v_steel("v_steel", box_steel, m_Steel);
  Volume v_vacuum("v_vacuum", box_vacuum, m_Vacuum);

  sdet.setAttributes(det, v_steel, x_det.regionStr(), x_det.limitsStr(), vis_name);

  assembly.placeVolume(v_steel, Position(0, 0, 0));
  v_steel.placeVolume(v_vacuum, Position(colXOff, 0, 0));

  // -----------------------------
  // final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(
      assembly, Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, zOff)));
  pv_assembly.addPhysVolID("system", sdet.id()).addPhysVolID("barrel", 1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(BackwardsCollimator, create_detector)
