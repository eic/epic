#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
 
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"
#include "Math/GenVector/RotationZYX.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/BoostY.h"
#include "Math/GenVector/BoostZ.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"
 
using namespace ROOT::Math;

  // Read edm4hep.root file and plot x,y of all hits
// Used to test the SVT OB geometry constructed using imported gdml files
// draws circle of the OB radii to check correct placement of staves
// https://github.com/eic/epic/tree/SVTOB_UK
// updated Sam Henry 11 July 2025


void read_edm4hep_SiBarrel(TString infile="TestHitMap.edm4hep.root") {
  
  TGraph* g1 = new TGraph();

  TFile* eicFile = new TFile(infile);
  TTree* eicTree = (TTree*)eicFile->Get("events");
  int Nevents = eicTree->GetEntries();
  TTreeReader tr(eicTree);
  TTreeReaderArray<Double_t> SiB_x(tr,"SiBarrelHits.position.x");
  TTreeReaderArray<Double_t> SiB_y(tr,"SiBarrelHits.position.y");
  TTreeReaderArray<Double_t> SiB_z(tr,"SiBarrelHits.position.z");
  int EventCount=0, Nhits=0, L3hits=0, L4hits=0;
  while (tr.Next()) {
    for(int i = 0; i < SiB_x.GetSize(); i++) {
  	double R = sqrt(SiB_x[i]*SiB_x[i] + SiB_y[i]*SiB_y[i] );
  	double phi=atan(SiB_y[i]/SiB_x[i]);
	if (SiB_x[i] < 0) 
  	    phi += TMath::Pi(); 
  	if (phi < 0.0)
  	    phi += 2.0*TMath::Pi();
        if(R>300){
  //       g1->SetPoint(Nhits,SiB_x[i],SiB_y[i]);
//         g1->SetPoint(Nhits,SiB_z[i],R);
           g1->SetPoint(Nhits,SiB_z[i],phi);
           Nhits++;
        }
         if (R>300) L4hits++; else L3hits++;
         }
    EventCount++;
  }
  cout << "Events: " << Nevents << endl;
  cout << "Counted: " << EventCount << endl;
  cout << "Hits: " << Nhits << endl;
  cout << "L3 hits: " << L3hits << endl;
  cout << "L4 hits " << L4hits << endl;
  TCanvas *c1 = new TCanvas("c1","hist",200,10,1400,1000);
  gStyle->SetOptStat(0);

  g1->Draw("AP");
//  c1->Print("SiBarrel.pdf");
}
