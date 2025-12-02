// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Whitney Armstrong, Wouter Deconinck

#include <DD4hep/Detector.h>

#include <ActsExamples/DD4hepDetector/ePICDetector.hpp>

/**
 * Example: Loading and testing ACTS Gen3 geometry from a DD4hep compact file.
 */
void test_ACTS_gen3(const char* compact = "epic.xml") {
  Acts::GeometryContext trackingGeoCtx;
  ActsExamples::ePICDetector::Config detectorConfig;
  detectorConfig.logLevel       = Acts::Logging::Level::VERBOSE;
  detectorConfig.dd4hepLogLevel = Acts::Logging::Level::WARNING;
  detectorConfig.xmlFileNames.push_back(compact);
  detectorConfig.name = compact;
  ActsExamples::ePICDetector epicDetector(detectorConfig, trackingGeoCtx);
}
