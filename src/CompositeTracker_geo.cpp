//==========================================================================
//  Based off DD4hep_SubDetectorAssembly
//
//  This is a simple plugin to allow compositing different detectors
//  into a single barrel or endcap for ACTS
//  Note: positive/negative position strings differentiate between
//        positive/negative endcaps
//--------------------------------------------------------------------------

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"

#if defined(USE_ACTSDD4HEP)
#include "ActsDD4hep/ActsExtension.hpp"
#else
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#endif

using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_element(Detector& description, xml_h e, Ref_t)
{
  xml_det_t         x_det(e);
  const std::string det_name = x_det.nameStr();
  DetElement        sdet(det_name, x_det.id());
  Volume            vol;
  Position          pos;

  const bool usePos = x_det.hasChild(_U(position));

  sdet.setType("compound");
  xml::setDetectorTypeFlag(e, sdet);

  const std::string actsType = getAttrOrDefault(x_det, _Unicode(actsType), "endcap");
  printout(DEBUG, det_name, "+++ Creating composite tracking detector (type: " + actsType + ")");
  assert(actsType == "barrel" || actsType == "endcap");

  // ACTS extension
  {
    Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();
    detWorldExt->addType(actsType, "detector");
    sdet.addExtension<Acts::ActsExtension>(detWorldExt);
  }

  if (usePos) {
    pos = Position(x_det.position().x(), x_det.position().y(), x_det.position().z());
  }
  vol = Assembly(det_name);
  vol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  Volume       mother = description.pickMotherVolume(sdet);
  PlacedVolume pv;
  if (usePos) {
    pv = mother.placeVolume(vol, pos);
  } else {
    pv = mother.placeVolume(vol);
  }
  sdet.setPlacement(pv);
  for (xml_coll_t c(x_det, _U(composite)); c; ++c) {
    xml_dim_t         component = c;
    const std::string nam       = component.nameStr();
    description.declareParent(nam, sdet);
  }
  return sdet;
}

DECLARE_DETELEMENT(athena_CompositeTracker, create_element)
