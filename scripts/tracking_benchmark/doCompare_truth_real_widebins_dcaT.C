// Code to compare the tracking performances: Truth seeding vs real seeding
// Shyam Kumar; shyam.kumar@ba.infn.it; shyam055119@gmail.com

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

void draw_req_DCA(double etamin, double etamax, double xmin=0., double xmax=0.);
void doCompare_truth_real_widebins_dcaT(TString particle = "pi-",double etamin=-1.0, double etamax=1.0, Bool_t drawreq=1, TString epic ="24.06.0", TString eicrecon = "v1.14.0") // name = p, pt for getting p or pt dependence fitted results
{

//=== style of the plot=========
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleOffset(1.0,"XY");
   gStyle->SetTitleSize(.04,"XY");
   gStyle->SetLabelSize(.04,"XY");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(1);
  
   const Int_t nptbins = 10;
   double pt[nptbins] ={0.2, 0.3, 0.5,1.0, 1.5, 2.0, 5.0, 8.0, 10., 15.0};
   Double_t variation = 0.1; // 10 % variation 
   std::vector<double> momV_truth, momV_real;
   std::vector<double> dcaTresolV_truth, err_dcaTresolV_truth, dcaTresolV_real, err_dcaTresolV_real;
   momV_truth.clear(); momV_real.clear(); dcaTresolV_truth.clear(); err_dcaTresolV_truth.clear(); dcaTresolV_real.clear(); err_dcaTresolV_real.clear();
   
   
   TString symbolname = "";
   if (particle == "pi-") symbolname = "#pi^{-}"; 
   else symbolname = particle;
   // Write the parameters in text file (for comparison)
   ofstream outfile;
   outfile.open ("DCAT_resol.txt",ios_base::app);  
   
   
   TF1 *f1=new TF1("f1","FitPointingAngle",0.,30.0,2);
   f1->SetParLimits(0,0.,50000);	
   f1->SetParLimits(1,0.,50000);	
  
   
   TCanvas *c_dcaxy = new TCanvas("c_dcaxy","c_dcaxy",1400,1000);
   c_dcaxy->SetMargin(0.10, 0.05 ,0.1,0.05);
   c_dcaxy->SetGridy();
 
   //Reading the root file
    TFile *fDCA_truth, *fDCA_real;
    TGraphErrors *gr_dcaT_truth, *gr_dcaT_real, *gr_dcaZ_truth, *gr_dcaZ_real;
	 
    TMultiGraph *mgDCAT; 
    TLegend *lDCAT; 
     mgDCAT = new TMultiGraph("mgDCAT",";p_{T} (GeV/c); #sigma_{DCA_{T}} (#mum)");
     lDCAT = new TLegend(0.65,0.80,0.90,0.93);
     lDCAT->SetTextSize(0.03);
     lDCAT->SetBorderSize(0);
     lDCAT->SetHeader(Form("%s ePIC(%s/%s): %1.1f < #eta < %1.1f",symbolname.Data(),epic.Data(),eicrecon.Data(),etamin,etamax),"C");
      
      fDCA_truth = TFile::Open(Form("truthseed/%s/dca/final_hist_dca_truthseed.root",particle.Data()));
      fDCA_real = TFile::Open(Form("realseed/%s/dca/final_hist_dca_realseed.root",particle.Data()));
	 
	 // Truth seeding histograms
      TH3D *hist_d0xy_truth = (TH3D*) fDCA_truth->Get("h_d0xy_3d");
      TH3D *hist_d0xy_real = (TH3D*) fDCA_real->Get("h_d0xy_3d");
      
      // d0xy calculation for truth/real (binning are same in both cases)
     Int_t etamin_bin = hist_d0xy_truth->GetYaxis()->FindBin(etamin+0.0001);
     Int_t etamax_bin = hist_d0xy_truth->GetYaxis()->FindBin(etamax-0.0001);
	
     TF1 *func_truth = new TF1("func_truth","gaus",-0.5,0.5);
     TF1 *func_real = new TF1("func_real","gaus",-0.5,0.5);
	
    for(int iptbin=0; iptbin<nptbins; ++iptbin){
    	
   TCanvas *cp = new TCanvas("cp","cp",1400,1000);
   cp->SetMargin(0.10, 0.05 ,0.1,0.07);
     
    double ptmin = (1.0-variation)*pt[iptbin]; // 10% range  
    double ptmax = (1.0+variation)*pt[iptbin]; // 10% range 
  
    Int_t ptmin_bin = hist_d0xy_truth->GetZaxis()->FindBin(ptmin+0.0001);
    Int_t ptmax_bin = hist_d0xy_truth->GetZaxis()->FindBin(ptmax-0.0001);
		
    TH1D *histd0xy_truth_1d = (TH1D*)hist_d0xy_truth->ProjectionX(Form("histd0xy_truth_eta%1.1f_%1.1f_pt%1.1f_%1.1f",etamin,etamax,ptmin,ptmax),etamin_bin,etamax_bin,ptmin_bin,ptmax_bin,"o");
    histd0xy_truth_1d->SetTitle(Form("d0_{xy} (truth): %1.1f <#eta< %1.1f && %1.2f <p_{T}< %1.2f",etamin,etamax,ptmin,ptmax));   
    histd0xy_truth_1d->SetName(Form("eta_%1.1f_%1.1f_d0xy_truth_pt_%1.1f",etamin,etamax,pt[iptbin]));  
   
   //if (histd0xy_truth_1d->GetEntries()<100) continue;   
   double mu_truth = histd0xy_truth_1d->GetMean(); 
   double sigma_truth = histd0xy_truth_1d->GetStdDev();
   func_truth->SetRange(mu_truth-2.0*sigma_truth,mu_truth+2.0*sigma_truth); // fit with in 2 sigma range
   histd0xy_truth_1d->Fit(func_truth,"NR+");
   mu_truth = func_truth->GetParameter(2);
   sigma_truth = func_truth->GetParError(2);
   func_truth->SetRange(mu_truth-2.0*sigma_truth,mu_truth+2.0*sigma_truth); // fit with in 2 sigma range
   histd0xy_truth_1d->Fit(func_truth,"R+");
   float truth_par2 = func_truth->GetParameter(2)*10000; // cm to mum 10000 factor
   float truth_par2_err = func_truth->GetParError(2)*10000;
   momV_truth.push_back(pt[iptbin]);
   dcaTresolV_truth.push_back(truth_par2);
   err_dcaTresolV_truth.push_back(truth_par2_err);
 
    TH1D *histd0xy_real_1d = (TH1D*)hist_d0xy_real->ProjectionX(Form("histd0xy_real_eta%1.1f_%1.1f_pt%1.1f_%1.1f",etamin,etamax,ptmin,ptmax),etamin_bin,etamax_bin,ptmin_bin,ptmax_bin,"o");
    histd0xy_real_1d->SetTitle(Form("d0_{xy} (real): %1.1f <#eta< %1.1f && %1.2f <p_{T}< %1.2f",etamin,etamax,ptmin,ptmax));   
    histd0xy_real_1d->SetName(Form("eta_%1.1f_%1.1f_d0xy_real_pt_%1.1f",etamin,etamax,pt[iptbin])); 
   
   //if (histd0xy_real_1d->GetEntries()<100) continue;    
   double mu_real = histd0xy_real_1d->GetMean(); 
   double sigma_real = histd0xy_real_1d->GetStdDev();
   func_real->SetRange(mu_real-2.0*sigma_real,mu_real+2.0*sigma_real); // fit with in 2 sigma range
   histd0xy_real_1d->Fit(func_real,"NR+");
   mu_real = func_real->GetParameter(1); 
   sigma_real = func_real->GetParameter(2);
   func_real->SetRange(mu_real-2.0*sigma_real,mu_real+2.0*sigma_real); // fit with in 2 sigma range
   histd0xy_real_1d->Fit(func_real,"R+");
   float real_par2 = func_real->GetParameter(2)*10000; // cm to mum 10000 factor
   float real_par2_err = func_real->GetParError(2)*10000;
   momV_real.push_back(pt[iptbin]);
   dcaTresolV_real.push_back(real_par2);
   err_dcaTresolV_real.push_back(real_par2_err);
   
   cp->cd();
   histd0xy_truth_1d->Draw();
   cp->SaveAs(Form("Debug_Plots/truth/%s/dca/truth_dcaxy_resol_mom%1.1f_%1.1f_eta_%1.1f.png",particle.Data(),pt[iptbin],etamin,etamax));
   cp->Clear();
   cp->cd();
   histd0xy_real_1d->Draw();
   cp->SaveAs(Form("Debug_Plots/real/%s/dca/real_dcaxy_resol_mom%1.1f_%1.1f_eta_%1.1f.png",particle.Data(),pt[iptbin],etamin,etamax));
   delete cp;
   }   // ptbin
      
	const int size_truth = momV_truth.size();
	double pt_truth[size_truth], err_pt_truth[size_truth], sigma_dcaxy_truth[size_truth], err_sigma_dcaxy_truth[size_truth]; 
	
	for (int i=0; i<size_truth; i++){
	pt_truth[i] = momV_truth.at(i);
	sigma_dcaxy_truth[i] = dcaTresolV_truth.at(i);
	err_sigma_dcaxy_truth[i] = err_dcaTresolV_truth.at(i);
	err_pt_truth[i] = 0.;
	}
	
     const int size_real = momV_real.size();
	double pt_real[size_real], err_pt_real[size_real], sigma_dcaxy_real[size_real], err_sigma_dcaxy_real[size_real]; 
	
	for (int i=0; i<size_real; i++){
	pt_real[i] = momV_real.at(i);
     sigma_dcaxy_real[i] = dcaTresolV_real.at(i);
	err_sigma_dcaxy_real[i] = err_dcaTresolV_real.at(i);
	err_pt_real[i] = 0.;
	}

     TFile *fout = new TFile(Form("Final_Results/%s/dca/dcaxy_resol_%1.1f_eta_%1.1f.root",particle.Data(),etamin,etamax),"recreate");
	TGraphErrors *gr1 = new TGraphErrors(size_truth,pt_truth,sigma_dcaxy_truth,err_pt_truth,err_sigma_dcaxy_truth);
     gr1->SetName("gr_truthseed");
	gr1->SetMarkerStyle(25);
	gr1->SetMarkerColor(kMagenta);
	gr1->SetMarkerSize(2.0);
	gr1->SetTitle(";p_{T} (GeV/c); #sigma_{DCA_{T}} (#mum)");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	
     TGraphErrors *gr2 = new TGraphErrors(size_real,pt_real,sigma_dcaxy_real,err_pt_real,err_sigma_dcaxy_real);
     gr2->SetName("gr_realseed");
	gr2->SetMarkerStyle(34);
	gr2->SetMarkerColor(kRed);
	gr2->SetMarkerSize(2.0);
	gr2->SetTitle(";p_{T} (GeV/c); #sigma_{DCA_{T}} (#mum)");
	gr2->GetXaxis()->CenterTitle();
	gr2->GetYaxis()->CenterTitle();
	
//	mgDCAT->Add(gr1);
	mgDCAT->Add(gr2);
	c_dcaxy->cd();
//	c_dcaxy->SetLogy();
     mgDCAT->GetXaxis()->SetRangeUser(0.18,mgDCAT->GetXaxis()->GetXmax()); 
	   mgDCAT->GetYaxis()->SetRangeUser(0.0,1.50*TMath::MaxElement(gr2->GetN(),gr2->GetY())); 
     mgDCAT->Draw("AP");
    // lDCAT->AddEntry(gr1,"Truth Seeding");
     lDCAT->AddEntry(gr2,"Realistic Seeding");
     lDCAT->Draw("same");
     draw_req_DCA(etamin,etamax,0.,mgDCAT->GetXaxis()->GetXmax());
     c_dcaxy->SaveAs(Form("Final_Results/%s/dca/dcaxy_resol_%1.1f_eta_%1.1f.png",particle.Data(),etamin,etamax));
 
   // Write the numbers in output file for comparisons
   outfile<<"ePIC"<<setw(20)<<epic.Data()<<setw(20)<<"EICRecon"<<setw(20)<<eicrecon.Data()<<endl;
   outfile<<"Etamin"<<setw(20)<<"Etamax"<<setw(20)<<"Pt (GeV/c) \t"<<setw(20)<<"DCAT Resol #mum (Real)"<<endl;
   for (Int_t i = 0; i<gr1->GetN(); ++i){
   double x,ytrue, yreal;
   gr2->GetPoint(i,x,yreal);  
   outfile<<etamin<<setw(20)<<etamax<<setw(20)<<x<<setw(20)<<yreal<<endl;
   } 
    outfile.close();
     
     fout->cd();
     mgDCAT->SetName(Form("dcaxy_resol_%1.1f_eta_%1.1f",etamin,etamax));
     mgDCAT->Write();
     fout->Close();
}

//From Yellow report from section 11.2.2
//===Fit Pointing Resolution
float FitPointingAngle(Double_t *x, Double_t *par)
{
  float func = sqrt((par[0]*par[0])/(x[0]*x[0])+par[1]*par[1]);
  return func;
}

void draw_req_DCA(double etamin, double etamax, double xmin=0., double xmax=0.)
{
  TF1 *PWGReq_DCA2D;
  if (etamin >= -3.5 && etamax <= -2.5) PWGReq_DCA2D = new TF1("PWGReq_DCA2D", "TMath::Sqrt((30./x)^2+40.0^2)",xmin,xmax);
  else if (etamin >= -2.5 && etamax <= -1.0) PWGReq_DCA2D = new TF1("PWGReq_DCA2D", "TMath::Sqrt((30./x)^2+20.0^2)",xmin,xmax);
  else if (etamin >= -1.0 && etamax <= 1.0) PWGReq_DCA2D = new TF1("PWGReq_DCA2D", "TMath::Sqrt((20./x)^2+5.0^2)",xmin,xmax);
  else if (etamin >= 1.0 && etamax <= 2.5) PWGReq_DCA2D = new TF1("PWGReq_DCA2D", "TMath::Sqrt((30./x)^2+20.0^2)",xmin,xmax);
  else if (etamin >= 2.5 && etamax <= 3.5) PWGReq_DCA2D = new TF1("PWGReq_DCA2D", "TMath::Sqrt((30./x)^2+40.0^2)",xmin,xmax);
  else return;
  PWGReq_DCA2D->SetLineStyle(7);
  PWGReq_DCA2D->SetLineWidth(3.0);
  PWGReq_DCA2D->SetLineColor(kBlue);
  PWGReq_DCA2D->Draw("same");
		
  TLegend *l= new TLegend(0.70,0.75,0.90,0.80);
  l->SetTextSize(0.03);
  l->SetBorderSize(0);
  l->AddEntry(PWGReq_DCA2D,"PWGReq","l");
  l->Draw("same");
}
