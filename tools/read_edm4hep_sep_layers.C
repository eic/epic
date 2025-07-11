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
// plots upper, lower, layers of two radii in different colours

// draws circle of the OB radii to check correct placement of staves
// https://github.com/eic/epic/tree/SVTOB_UK
// Sam Henry 3 July 2025


void read_edm4hep_sep_layers(TString infile="TestHitMap.edm4hep.root") {

  bool XY = true; //true = XY plot, false = Z-phi plot

  double minZ=-500, maxZ=500;
  double minR=300, maxR= 500;
  double Rx1 = 0.0, Rx2=0, Ro=427, Ri=421;
  double dTh = 2.0 * TMath::Pi() / 70.0;

  TGraph* g1 = new TGraph(); // upper layer, outer radius
  TGraph* g2 = new TGraph(); // lower layer, inner radius
  TGraph* g3 = new TGraph(); // upper layer, inner radius
  TGraph* g4 = new TGraph(); // lower layer, outer radius
  TGraph* gc = new TGraph(); // radius of cut line
  
  int n1=0, n2=0, n3=0, n4=0, nc=0;
  
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
  	double phi = atan(SiB_y[i]/SiB_x[i]);
  	if (SiB_x[i] < 0) 
  	    phi += TMath::Pi(); 
  	if (phi < 0.0)
  	    phi += 2.0*TMath::Pi();
  	   
        if(SiB_z[i]>minZ && SiB_z[i]<maxZ && R>minR && R<maxR){

          if(R>Ro){
            if (XY)  g1->SetPoint(n1,SiB_x[i],SiB_y[i]);
            else g1->SetPoint(n1,SiB_z[i],phi);
	   n1++;
          }    else   if(R<Ri){
 	     if (XY) g2->SetPoint(n2,SiB_x[i],SiB_y[i]);
             else g2->SetPoint(n2,SiB_z[i],phi);
	   n2++;
           } else {
             double phi0 = (int)(35.0*(phi/(2.*TMath::Pi())))* 2.0*dTh + dTh; 
             double dphi = abs(phi - phi0);
             double rc = (Ri-Rx2) + ((Ro+Rx1)-(Ri-Rx2))*(dphi/(dTh));
             gc->SetPoint(nc, rc*cos(phi), rc*sin(phi)); nc++;
             if (R<rc){
              if (XY) g3->SetPoint(n3,SiB_x[i],SiB_y[i]);
              else g3->SetPoint(n3,SiB_z[i],phi);
	       n3++;
             } else{
               if (XY) g4->SetPoint(n4,SiB_x[i],SiB_y[i]);
               else g4->SetPoint(n4,SiB_z[i],phi);
	       n4++;
             }
          }
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
  
  
if (!XY){    
  g1->GetXaxis()->SetLimits(-500,500);
g1->GetHistogram()->SetMinimum(-0.5);
g1->GetHistogram()->SetMaximum(6.78);
g1->SetTitle("Hit map 1M muons; Z; phi");

g2->GetXaxis()->SetLimits(-500,500);
g2->GetHistogram()->SetMinimum(-0.5);
g2->GetHistogram()->SetMaximum(6.78);
g2->SetTitle("Hit map 1M muons; Z; phi");

g3->GetXaxis()->SetLimits(-500,500);
g3->GetHistogram()->SetMinimum(-0.5);
g3->GetHistogram()->SetMaximum(6.78);
g3->SetTitle("Hit map 1M muons; Z; phi");

g4->GetXaxis()->SetLimits(-500,500);
g4->GetHistogram()->SetMinimum(-0.5);
g4->GetHistogram()->SetMaximum(6.78);
g4->SetTitle("Hit map 1M muons; Z; phi"); 
 } else{
 g1->SetTitle("Hit map 1M muons; X; Y");
 }

g1->SetMarkerColor(4);
g1->Draw("AP");
g2->SetMarkerColor(3);
g2->Draw("P");
g3->SetMarkerColor(6);
g3->Draw("P");  
g4->SetMarkerColor(7);
g4->Draw("P");
gc->SetMarkerColor(2);  
if (XY) gc->Draw("P");
  
  if(XY){
   TEllipse* circle3 = new TEllipse(0.0,0.0,421); 
    circle3->SetLineColor(kRed);  
    circle3->SetFillStyle(0); // Transparent fill
    circle3->Draw();
    TEllipse* circle4 = new TEllipse(0.0,0.0,427); 
     circle4->SetLineColor(kRed);  
     circle4->SetFillStyle(0); // Transparent fill
    circle4->Draw();
    }
   c1->Update();
 
//  c1->Print("SiBarrel.pdf");
}
