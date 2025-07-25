#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoBBox.h"
#include <vector>
#include <iostream>
#include <fstream>

// Plot the x,y coordinates of all vertices extracted from detector_geometry.root file (in TGeoTessellated objects)
// Used to test the SVT OB geometry constructed using imported gdml files
// draws circle of the OB radii to check correct placement of staves
// https://github.com/eic/epic/tree/SVTOB_UK
// Sam Henry 25 April 2025

ofstream outputFile;
int level = 0;
TGraph* g1 = new TGraph();
int np=0;
double L3min=500, L3max=0, L4min=500, L4max=0;
double L3Zmin=500, L3Zmax=-500, L4Zmin=500, L4Zmax=-500;


   // Function to recursively extract vertices from a node
   void extractVertices(TGeoNode* node, TGeoHMatrix parentMatrix) {
        auto matrix = node->GetMatrix();
        const Double_t* trans = matrix->GetTranslation();
        const Double_t* rot = matrix->GetRotationMatrix();
        if(trans[0]>0.0 || trans[1]>0.0){
        
        for(int j=0; j<level; ++j) cout << "\t";
    	cout << node->GetName() << " (" << trans[0] << ", " << trans[1] << ", " << trans[2] << ")  (" << rot[0] ;
    	for(int k=1; k<8; ++k)  cout << "," << rot[k] ;
    	cout << ", " << rot[8] << ")"  << endl;
                     
                     }
        TGeoVolume* nodeVolume = node->GetVolume();
        TGeoShape* shape = nodeVolume->GetShape();
        
        TGeoHMatrix combinedMatrix(parentMatrix);
        combinedMatrix.Multiply(node->GetMatrix());
        
        if (shape->IsA() == TGeoTessellated::Class()) {
            TGeoTessellated* tes = (TGeoTessellated*)shape;
            int Nvert = tes->GetNvertices();
            for(int i=0; i<Nvert; ++i){
               auto v = tes->GetVertex(i);
               
               double local[3] = {v.x(), v.y(), v.z() };
               double global[3];
               combinedMatrix.LocalToMaster(local, global);
               
//               outputFile << v.x() << "\t" << v.y() << "\t" << v.z()   <<  endl;
               double r = sqrt(global[0]*global[0] + global[1]*global[1]);
//               g1->SetPoint(np,global[2], r);
               g1->SetPoint(np,global[0], global[1]);
               np++;
               if (r>L4max) L4max=r;
               if(r<L4min && r>35) L4min=r;
               if(r<L3min) L3min=r;
               if(r>L3max && r<35) L3max=r;
               double z=global[2];
               if (z>L4Zmax && r>35) L4Zmax=z;
               if(z<L4Zmin && r>35) L4Zmin=z;
               if(z<L3Zmin && r<35) L3Zmin=z;
               if(z>L3Zmax && r<35) L3Zmax=z;

            }
        }
        // Recursively process child nodes
        for (int i = 0; i < nodeVolume->GetNdaughters(); ++i) {
            level += 1;
            extractVertices(nodeVolume->GetNode(i), combinedMatrix);
            level -= 1;
        }
    }

void ExtractVertexCoordinates(TString filename="detector_geometry.root") {
   
    auto geom = TGeoManager::Import(filename);
    TGeoVolume* volume = geom->FindVolumeFast("world_volume");
    outputFile.open("vertex_output.txt");
    for (int i =0; i<volume->GetNdaughters(); ++i){
            TGeoNode* topNode = volume->GetNode(i);
            if (topNode) {
                TGeoHMatrix identityMatrix;
                extractVertices(topNode, identityMatrix);
            }
    }
    outputFile.close();
TCanvas *c1 = new TCanvas("c1","hist",200,10,1400,1000);
  gStyle->SetOptStat(0);
g1->Draw("AP");

    TEllipse* circle1 = new TEllipse(0.0,0.0,26.9);
    circle1->SetLineColor(kRed);  
    circle1->SetFillStyle(0); // Transparent fill
    circle1->Draw();
    
    TEllipse* circle2 = new TEllipse(0.0,0.0,27.5);
    circle2->SetLineColor(kRed);  
    circle2->SetFillStyle(0); // Transparent fill
    circle2->Draw();
   TEllipse* circle3 = new TEllipse(0.0,0.0,42.1);
    circle3->SetLineColor(kRed);  
    circle3->SetFillStyle(0); // Transparent fill
    circle3->Draw();
   TEllipse* circle4 = new TEllipse(0.0,0.0,42.7);
    circle4->SetLineColor(kRed);  
circle4->SetFillStyle(0); // Transparent fill
    circle4->Draw();

   c1->Update();
   cout << "L3 minimum radius: " << L3min << " = 26.9 - " << 26.9-L3min << endl;
   cout << "L3 maximum radius: " << L3max << " = 26.9 + " << L3max-26.9 << endl;
   cout << "L4 minimum radius: " << L4min << " = 42.1 - " << 42.1-L4min << endl;
   cout << "L4 maximum radius: " << L4max << " = 42.1 + " << L4max-42.1 << endl;
   cout << "L3 minimum Z: " << L3Zmin  << " -> length " << 2.0*abs(L3Zmin) << endl;
   cout << "L3 maximum Z: " << L3Zmax << " -> length " << 2.0*L3Zmax << endl;
   cout << "L4 minimum Z: " << L4Zmin <<" -> length " << 2.0*abs(L4Zmin)<<  endl;
   cout << "L4 maximum Z: " << L4Zmax<< " -> length " << 2.0*L4Zmax  << endl;
   
}


