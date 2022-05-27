#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Objects.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "TCanvas.h"
#include "TChain.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

/** Example loading ACTs.
 *
 *
 */
void test_ACTS(const char* compact = "ecce.xml")
{

  using namespace ROOT::Math;
  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compact);
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);

  // std::unique_ptr<const Acts::TrackingGeometry>
  auto acts_tracking_geometry = Acts::convertDD4hepDetector(detector.world(), Acts::Logging::Level::VERBOSE);
  // acts_tracking_geometry  = Acts::convertDD4hepDetector (detector.world(),Acts::Logging::Level::INFO);

  if (acts_tracking_geometry) {
    std::cout << "success?\n";
  }
  //  if(acts_tracking_geometry->highestTrackingVolume()) {
  //    std::cout << " volume name \n ";
  //    std::cout << acts_tracking_geometry->highestTrackingVolume()->volumeName() << std::endl;
  //  } else {
  //    std::cout << "derp\n";
  //  }
  //}
}
