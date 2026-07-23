// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Wouter Deconinck, Copilot

#pragma once

#include <string>
#include <vector>

struct G4OCCTSectorVertex {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct G4OCCTSectorTriangle {
  G4OCCTSectorVertex v[3];
};

struct G4OCCTSectorGeometry {
  std::vector<G4OCCTSectorTriangle> triangles;
};

G4OCCTSectorGeometry ImportSectorSTEPWithG4OCCT(const std::string& name, const std::string& path);
