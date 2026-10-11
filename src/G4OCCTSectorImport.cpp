// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Wouter Deconinck, Copilot

#include "G4OCCTSectorImport.h"

#ifdef WITH_G4OCCT_SUPPORT

#include "G4OCCT/G4OCCTSolid.hh"

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <stdexcept>
#include <utility>

G4OCCTSectorGeometry ImportSectorSTEPWithG4OCCT(const std::string& name, const std::string& path) {
  G4OCCTSolid* solid = nullptr;
  try {
    solid = G4OCCTSolid::FromSTEP(name, path);
  } catch (const std::exception& ex) {
    throw std::runtime_error("BarrelHCal G4OCCT import failed for '" + path + "' (" + ex.what() +
                             ")");
  }

  const TopoDS_Shape shape = solid->GetOCCTShape();
  delete solid;
  solid = nullptr;

  static constexpr double kRelativeDeflection = 0.01;
  BRepMesh_IncrementalMesh mesher(shape, kRelativeDeflection,
                                  /*isRelative=*/true);
  (void)mesher;

  G4OCCTSectorGeometry result;
  for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
    const TopoDS_Face& face = TopoDS::Face(explorer.Current());
    TopLoc_Location location;
    const Handle(Poly_Triangulation) & tri = BRep_Tool::Triangulation(face, location);
    if (tri.IsNull() || tri->NbTriangles() == 0) {
      continue;
    }

    const gp_Trsf& transform   = location.Transformation();
    const bool reverse_winding = (face.Orientation() == TopAbs_REVERSED);

    for (int i = 1; i <= tri->NbTriangles(); ++i) {
      int n1, n2, n3;
      tri->Triangle(i).Get(n1, n2, n3);
      if (reverse_winding) {
        std::swap(n2, n3);
      }

      const gp_Pnt p1 = tri->Node(n1).Transformed(transform);
      const gp_Pnt p2 = tri->Node(n2).Transformed(transform);
      const gp_Pnt p3 = tri->Node(n3).Transformed(transform);

      G4OCCTSectorTriangle facet;
      facet.v[0] = {p1.X(), p1.Y(), p1.Z()};
      facet.v[1] = {p2.X(), p2.Y(), p2.Z()};
      facet.v[2] = {p3.X(), p3.Y(), p3.Z()};
      result.triangles.push_back(facet);
    }
  }

  if (result.triangles.empty()) {
    throw std::runtime_error("BarrelHCal G4OCCT tessellation produced no "
                             "triangles for '" +
                             path + "'");
  }

  return result;
}

#else

#include <stdexcept>

G4OCCTSectorGeometry ImportSectorSTEPWithG4OCCT(const std::string&, const std::string&) {
  throw std::runtime_error("This epic build was configured without G4OCCT support");
}

#endif
