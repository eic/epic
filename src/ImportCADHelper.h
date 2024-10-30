#pragma once
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <array>
#include "DD4hepDetectorHelper.h"
#include <cmath>


#include "TVector3.h"
//#include "TGDMLParse.h"



//typedef Object::Vertex_t Vertex;

//declare the conversion function
dd4hep::rec::Vector3D vertex_to_vector3D(dd4hep::TessellatedSolid::Vertex vertex);