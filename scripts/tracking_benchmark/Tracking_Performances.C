// Code to extract the Tracking Performances
// Shyam Kumar; INFN Bari, Italy
// shyam.kumar@ba.infn.it; shyam.kumar@cern.ch

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#define mpi 0.139  // 1.864 GeV/c^2

void Tracking_Performances(TString filename="tracking_output",TString particle="pi-", double mom=0.1, Double_t pTcut = 0.15, bool truth_seeding=false)
{

  // style of the plot
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleOffset(.85,"X");gStyle->SetTitleOffset(.85,"Y");
   gStyle->SetTitleSize(.05,"X");gStyle->SetTitleSize(.05,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(1);
   
   bool debug=true;	  
  // Tree with reconstructed tracks
   const int nbins_eta = 5;
   int theta_val[nbins_eta+1] ={3,50,45,135,130,177};
   int nfiles = 100; 
   double eta[nbins_eta+1]={-3.5,-2.5,-1.0,1.0,2.5,3.5};
   double pt[nbins_eta+1]={0.5,1.0,2.0,5.0,10.0,20.1};
   TH1D *histp[nbins_eta]; 
   
    TH3D *h_d0xy_3d= new TH3D("h_d0xy_3d","Transverse Pointing Resolution",500,-0.1,0.1,70,-3.5,3.5,201,0.,20.1);
    TH3D *h_d0z_3d= new TH3D("h_d0z_3d","Longitudinal Pointing Resolution",500,-0.1,0.1,70,-3.5,3.5,201,0.,20.1);
   
   for (int i=0; i<nbins_eta; i++){
   histp[i] = new TH1D(Form("hist_etabin%d",i),Form("hist_etabin%d",i),600,-0.3,0.3);
   histp[i]->SetTitle(Form("%1.1f < #eta < %1.1f && p = %1.1f ",eta[i],eta[i+1],mom));
   histp[i]->SetName(Form("hist_mom_%1.1f_%1.1f_pmax_%1.1f",mom,eta[i],eta[i+1]));
   }
   
   TFile* file = TFile::Open(filename.Data());
   if (!file) {printf("file not found !!!"); return;}
   TTreeReader myReader("events", file); // name of tree and file
   if (debug) cout<<"Filename: "<<file->GetName()<<"\t NEvents: "<<myReader.GetEntries()<<endl;
  
   TTree *tree = dynamic_cast<TTree*>(file->Get("events"));
   if (!(tree->GetBranch("CentralCKFSeededTrackParameters") || tree->GetBranch("CentralCKFTruthSeededTrackParameters"))) {
     cerr << "Found neither CentralCKFSeededTrackParameters nor CentralCKFTruthSeededTrackParameters!" << endl;
     return;
   }
   bool is_old_style = tree->GetBranch("CentralCKFSeededTrackParameters");

   TString dir = "";
   TString dist_dir_mom = ""; TString dist_dir_dca = "";
   TString tag = "";
   if (truth_seeding) {
      dist_dir_mom = "mom_resol_truth"; dist_dir_dca = "dca_resol_truth"; dir = "truthseed";
      if (is_old_style) {
        tag = "";
      } else {
        tag = "TruthSeeded";
      }
   } else {
      dist_dir_mom = "mom_resol_realseed"; dist_dir_dca = "dca_resol_realseed"; dir = "realseed";
      if (is_old_style) {
        tag = "Seeded";
      } else {
        tag = "";
      }
   }

   // MC and Reco information 
   TTreeReaderArray<Float_t> charge(myReader, "MCParticles.charge"); 
   TTreeReaderArray<Double_t> vx_mc(myReader, "MCParticles.vertex.x"); 
   TTreeReaderArray<Double_t> vy_mc(myReader, "MCParticles.vertex.y"); 
   TTreeReaderArray<Double_t> vz_mc(myReader, "MCParticles.vertex.z"); 
   TTreeReaderArray<Double_t> px_mc(myReader, "MCParticles.momentum.x"); 
   TTreeReaderArray<Double_t> py_mc(myReader, "MCParticles.momentum.y"); 
   TTreeReaderArray<Double_t> pz_mc(myReader, "MCParticles.momentum.z"); 
   TTreeReaderArray<Int_t> status(myReader, "MCParticles.generatorStatus"); 
   TTreeReaderArray<Int_t> pdg(myReader, "MCParticles.PDG"); 
   TTreeReaderArray<Int_t> match_flag(myReader, Form("CentralCKF%sTrackParameters.type",tag.Data()));
   TTreeReaderArray<Float_t> d0xy(myReader, Form("CentralCKF%sTrackParameters.loc.a",tag.Data()));
   TTreeReaderArray<Float_t> d0z(myReader, Form("CentralCKF%sTrackParameters.loc.b",tag.Data()));
   TTreeReaderArray<Float_t> theta(myReader, Form("CentralCKF%sTrackParameters.theta",tag.Data()));
   TTreeReaderArray<Float_t> phi(myReader, Form("CentralCKF%sTrackParameters.phi",tag.Data()));
   TTreeReaderArray<Float_t> qoverp(myReader, Form("CentralCKF%sTrackParameters.qOverP",tag.Data()));

  int count =0;
  int matchId = 1; // Always matched track assigned the index 0 
  while (myReader.Next()) 
  {
   if (match_flag.GetSize()==0) continue;  // Events with no reco tracks skip them
   for (int i = 0; i < matchId; ++i){
   
     for (int j = 0; j < pdg.GetSize(); ++j){
     	
      if (status[j] !=1 && pdg.GetSize()!=1) continue;
      Double_t ptmc = sqrt(px_mc[j]*px_mc[j]+py_mc[j]*py_mc[j]); 
      
      if (fabs(ptmc) < pTcut) continue;

      Double_t pmc = (1./charge[j])*sqrt(px_mc[j]*px_mc[j]+py_mc[j]*py_mc[j]+pz_mc[j]*pz_mc[j]); // 1./(q/p); similar to prec
      Double_t prec = 1./qoverp[j]; 

      Double_t pzrec = prec*TMath::Cos(theta[j]);  Double_t pt_rec = sqrt(prec*prec-pzrec*pzrec);  
      Double_t pzmc = pz_mc[j];  
      
      Double_t etamc = -1.0*TMath::Log(TMath::Tan((TMath::ACos(pzmc/fabs(pmc)))/2));
      Double_t p_resol = (prec-pmc)/pmc;
      
      for (int ibin=0; ibin<nbins_eta; ++ibin){ 
      if(etamc>eta[ibin] && etamc<eta[ibin+1]) histp[ibin]->Fill(p_resol); 
      }
      h_d0xy_3d->Fill(d0xy[j]*0.1, etamc, ptmc); // cm
      h_d0z_3d->Fill(d0z[j]*0.1, etamc, ptmc); // cm
      } // Generated Tracks  
     } // Reco Tracks
       
   }// event loop ends    
  
   TFile *fout_mom = new TFile(Form("%s/%s/mom/Performances_mom_%1.1f_%s_%s.root",dir.Data(),particle.Data(),mom,dist_dir_mom.Data(),particle.Data()),"recreate");
   fout_mom->cd();
   for (int ibin=0; ibin<nbins_eta; ++ibin) histp[ibin]->Write();
   fout_mom->Close();

   TFile *fout_dca = new TFile(Form("%s/%s/dca/Performances_dca_%1.1f_%s_%s.root",dir.Data(),particle.Data(),mom,dist_dir_dca.Data(),particle.Data()),"recreate");
   fout_dca->cd();
    h_d0xy_3d->Write();
    h_d0z_3d->Write();
    fout_dca->Close();
}




