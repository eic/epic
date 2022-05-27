#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "TMath.h"

// #if defined(USE_ACTSDD4HEP)
// #include "ActsDD4hep/ActsExtension.hpp"
// #else
// #include "Acts/Plugins/DD4hep/ActsExtension.hpp"
// #endif

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

/** \addtogroup trackers Tracking Detectors
 */

/** \addtogroup GEMdisc GEM Disc Tracker
 * \brief Type: **GEMTrackerDisc**.
 * \ingroup trackers
 *
 *  A simple GEM disc tracker.
 *
 * \code
 *   <detector>
 *   </detector>
 * \endcode
 *
 */
static Ref_t create_detector(Detector& lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t    x_det    = e;
  Material     carbon   = lcdd.material("CarbonFiber");
  int          det_id   = x_det.id();
  string       det_name = x_det.nameStr();
  PlacedVolume pv;

  DetElement sdet(det_name, det_id);
  Assembly   assembly(det_name + "_assembly");
  // Acts::ActsExtension* ecapDetExt = new Acts::ActsExtension();
  // ecapDetExt->addType("barrel", "detector");
  // sdet.addExtension<Acts::ActsExtension>(ecapDetExt);

  sens.setType("tracker");
  string module_name = "Timepix";

  int N_layers = 0;

  for (xml_coll_t lay(x_det, _U(layer)); lay; ++lay, ++N_layers) {

    xml_comp_t x_layer = lay;
    //    double radius = dd4hep::getAttrOrDefault(x_det, _Unicode(radius), 5.0*dd4hep::cm)
    double sensor_thickness = dd4hep::getAttrOrDefault(x_layer, _Unicode(sensor_thickness), 200 * dd4hep::um);
    double x_size   = dd4hep::getAttrOrDefault(x_layer, _Unicode(x_pix), 512) * lcdd.constantAsDouble("pixel_size");
    double y_size   = dd4hep::getAttrOrDefault(x_layer, _Unicode(y_pix), 448) * lcdd.constantAsDouble("pixel_size");
    double z        = dd4hep::getAttrOrDefault(x_layer, _Unicode(z), -10 * dd4hep::cm);
    int    layer_id = x_layer.id(); // attr<double>(  _Unicode(z) ) ;

    string layer_name = std::string("timepix_layer") + std::to_string(layer_id);

    Box    timepix_layer(x_size / 2, y_size / 2, sensor_thickness / 2);
    Volume timepix_layer_vol("timepix_layer_vol", timepix_layer, carbon);

    // DD4hep surfaces
    // Vector3D u( 1. , 0. , 0. ) ;
    // Vector3D v( 0. , 1. , 0. ) ;
    // Vector3D n( 0. , 0. , 1. ) ;
    // Vector3D o( 0. , 0. , 0. ) ;
    // double inner_thickness = thickness/2.0;
    // double outer_thickness = thickness/2.0;
    // SurfaceType type( SurfaceType::Sensitive ) ;
    // VolPlane    surf( gem_layer_vol, type, inner_thickness , outer_thickness , u,v,n,o ) ;

    timepix_layer_vol.setSensitiveDetector(sens);

    DetElement layer_DE(sdet, _toString(layer_id, "layer%d"), layer_id);

    // Acts::ActsExtension* detlayer = new Acts::ActsExtension();
    // detlayer->addType("sensitive disk", "layer");
    //  the local coordinate systems of modules in dd4hep and acts differ
    //  see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    // detlayer->addType("axes", "definitions", "XYZ");
    // layer_DE.addExtension<Acts::ActsExtension>(detlayer);

    // Assembly   layer_assembly( layer_name+"_assembly" );
    pv = assembly.placeVolume(timepix_layer_vol, Transform3D(RotationZ(0), Position(0.0, 0.0, z)));
    pv.addPhysVolID("layer", layer_id);
    layer_DE.setPlacement(pv);
    // layer_DE.setAttributes(lcdd, layer_assembly, "", "", "SiVertexLayerVis");
  }

  sdet.setAttributes(lcdd, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  // assembly.setVisAttributes(lcdd.invisible());

  pv = lcdd.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);

  assembly->GetShape()->ComputeBBox();
  return sdet;
}
DECLARE_DETELEMENT(my_Timepix, create_detector)
