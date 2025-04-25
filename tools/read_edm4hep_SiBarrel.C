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
// Sam Henry 25 April 2025


void read_edm4hep_SiBarrel(TString infile="TestHitMap.edm4hep.root") {
  
  int nbins=100;
  double maxR=500.; 
  double binwidth = 2.0 * maxR / (double)nbins;

  TH1* h0 = new TH1D("h0","Radius [mm]",nbins,0.0,maxR);
  TGraph* g1 = new TGraph();

  TFile* eicFile = new TFile(infile);
  TTree* eicTree = (TTree*)eicFile->Get("events");

  int Nevents = eicTree->GetEntries();

  TTreeReader tr(eicTree);

  TTreeReaderArray<Double_t> SiB_x(tr,"SiBarrelHits.position.x");
  TTreeReaderArray<Double_t> SiB_y(tr,"SiBarrelHits.position.y");
  TTreeReaderArray<Double_t> SiB_z(tr,"SiBarrelHits.position.z");
  int EventCount=0, Nhits=0;
  while (tr.Next()) {

    for(int i = 0; i < SiB_x.GetSize(); i++) {
    if(SiB_z[i]>0.0){
  	double R = sqrt(SiB_x[i]*SiB_x[i] + SiB_y[i]*SiB_y[i] );
  // double R = SiB_y[i];
         h0->Fill(R);
         g1->SetPoint(Nhits,SiB_x[i],SiB_y[i]);
         Nhits++;
         }
    }
    EventCount++;
  }
  cout << "Events: " << Nevents << endl;
  cout << "Counted: " << EventCount << endl;
  cout << "Hits: " << Nhits << endl;
  TCanvas *c1 = new TCanvas("c1","hist",200,10,1400,1000);
  gStyle->SetOptStat(0);

  h0->GetXaxis()->SetTitle("Radius");
  h0->GetYaxis()->SetTitle("Hits");
  //h0->Draw();
  g1->Draw("AP");
//  c1->Print("Hist_SiBarrel.pdf");
}
