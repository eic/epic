// Code to compare the tracking performances: Truth seeding vs real seeding
// Shyam Kumar; shyam.kumar@ba.infn.it; shyam055119@gmail.com

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#define mpi 0.139  // 1.864 GeV/c^2

void draw_req_Mom(double etamin, double etamax, double xmin=0., double xmax=0.);
void doCompare_truth_real_widebins_mom(TString particle = "pi-",double etamin=-1.0, double etamax=1.0, double range =0.3, Bool_t drawreq=1, TString extra_legend = "") // name = p, pt for getting p or pt dependence fitted results
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
  
   const Int_t nfiles = 15;
   double mom[nfiles] = { 0.50, 0.75, 1.0, 1.25, 1.75, 2.0, 2.50, 3.0, 4.0, 5.0, 7.0, 8.5, 10.0, 12.5, 15.0};
   std::vector<double> momV_truth, momV_real, momresolV_truth, err_momresolV_truth, momresolV_real, err_momresolV_real;
   momV_truth.clear(); momV_real.clear(); momresolV_truth.clear(); err_momresolV_truth.clear(); momresolV_real.clear(); err_momresolV_real.clear();
   TString symbolname = "";
   if (particle == "pi-") symbolname = "#pi^{-}"; 
   else symbolname = particle; 
   ofstream outfile;
   outfile.open ("Mom_resol.txt",ios_base::app);  
   
   TF1 *f1=new TF1("f1","FitMomentumResolution",0.,30.0,2);
   f1->SetParLimits(0,0.,0.1);	
   f1->SetParLimits(1,0.,5.0);	
   
   TCanvas *c_mom = new TCanvas("cmom","cmom",1400,1000);
   c_mom->SetMargin(0.10, 0.05 ,0.1,0.05);
   c_mom->SetGridy();

 //Reading the root file
  TFile *fmom_truth[nfiles], *fmom_real[nfiles];
  TGraphErrors *gr_mom_truth, *gr_mom_real;
 
  TMultiGraph *mgMom; 
  TLegend *lmom; 
  mgMom = new TMultiGraph("mgMom",";p (GeV/c); #sigmap/p %");
  
  lmom = new TLegend(0.65,0.80,0.90,0.93);
  lmom->SetTextSize(0.03);
  lmom->SetBorderSize(0);
  lmom->SetHeader(extra_legend.Data(), "C");
  lmom->AddEntry((TObject*)0, Form("%s, %1.1f < #eta < %1.1f", symbolname.Data(), etamin, etamax), "C");
  
  TF1 *func_truth = new TF1("func_truth","gaus",-0.5,0.5);
  TF1 *func_real = new TF1("func_real","gaus",-0.5,0.5);
 
  for (int i =0; i<nfiles; ++i){
      
   TCanvas *cp = new TCanvas("cp","cp",1400,1000);
   cp->SetMargin(0.10, 0.05 ,0.1,0.07);

	 fmom_truth[i] = TFile::Open(Form("./truthseed/pi-/mom/Performances_mom_%1.1f_mom_resol_truth_%s.root",mom[i],particle.Data()));
	 fmom_real[i] = TFile::Open(Form("./realseed/pi-/mom/Performances_mom_%1.1f_mom_resol_realseed_%s.root",mom[i],particle.Data()));
	 
	 TH1D *hist_truth = (TH1D*) fmom_truth[i]->Get(Form("hist_mom_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	 hist_truth->Rebin(2);
	 hist_truth->SetName(Form("truthseed_hist_mom_%1.1f_%1.1f_eta_%1.1f_%s",mom[i],etamin,etamax,particle.Data()));
	 hist_truth->SetTitle(Form("Truth Seed (%s): Momentum = %1.1f && %1.1f<#eta<%1.1f;#Delta p/p; Entries(a.u.)",particle.Data(),mom[i],etamin,etamax));

	 double mu_truth = hist_truth->GetMean(); 
	 double sigma_truth = hist_truth->GetStdDev();
         hist_truth->GetXaxis()->SetRangeUser(-1.0*range,1.0*range);
	 func_truth->SetRange(mu_truth-2.0*sigma_truth,mu_truth+2.0*sigma_truth); // fit with in 2 sigma range
	 hist_truth->Fit(func_truth,"NR+");
	 mu_truth = func_truth->GetParameter(1); 
	 sigma_truth = func_truth->GetParameter(2);
	 func_truth->SetRange(mu_truth-2.0*sigma_truth,mu_truth+2.0*sigma_truth);
	 hist_truth->Fit(func_truth,"R+");
	 float truth_par2 = func_truth->GetParameter(2)*100;
	 float truth_par2_err = func_truth->GetParError(2)*100;
	 momV_truth.push_back(mom[i]);
	 momresolV_truth.push_back(truth_par2);
	 err_momresolV_truth.push_back(truth_par2_err);

	 TH1D *hist_real = (TH1D*) fmom_real[i]->Get(Form("hist_mom_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	 hist_real->Rebin(2);
	 hist_real->SetName(Form("realseed_hist_mom_%1.1f_%1.1f_eta_%1.1f_%s",mom[i],etamin,etamax,particle.Data()));
	 hist_real->SetTitle(Form("Realistic Seed (%s): Momentum = %1.1f && %1.1f<#eta<%1.1f;#Delta p/p; Entries(a.u.)",particle.Data(),mom[i],etamin,etamax));
	 
	 double mu_real = hist_real->GetMean(); 
	 double sigma_real = hist_real->GetStdDev();
         hist_real->GetXaxis()->SetRangeUser(-1.0*range,1.0*range);
         func_real->SetRange(mu_real-2.0*sigma_real,mu_real+2.0*sigma_real); // fit with in 2 sigma range
	 hist_real->Fit(func_real,"NR+");
	 mu_real = func_real->GetParameter(1); 
	 sigma_real = func_real->GetParameter(2);
	 func_real->SetRange(mu_real-2.0*sigma_real,mu_real+2.0*sigma_real);
	 hist_real->Fit(func_real,"R+");
	 float real_par2 = func_real->GetParameter(2)*100;
	 float real_par2_err = func_real->GetParError(2)*100;
	 momV_real.push_back(mom[i]);
	 momresolV_real.push_back(real_par2);
	 err_momresolV_real.push_back(real_par2_err);
	 cp->cd();
	 hist_truth->Draw();
	 cp->SaveAs(Form("Debug_Plots/truth/%s/mom/truth_mom_resol_mom%1.1f_%1.1f_eta_%1.1f.png",particle.Data(),mom[i],etamin,etamax));
	 cp->Clear();
	 cp->cd();
	 hist_real->Draw();
	 cp->SaveAs(Form("Debug_Plots/real/%s/mom/real_mom_resol_mom%1.1f_%1.1f_eta_%1.1f.png",particle.Data(),mom[i],etamin,etamax));
 } // all files
 	 
	const int size_truth = momV_truth.size();
	double p_truth[size_truth], err_p_truth[size_truth], sigma_p_truth[size_truth], err_sigma_p_truth[size_truth]; 
	
	for (int i=0; i<size_truth; i++){
	p_truth[i] = momV_truth.at(i);
	sigma_p_truth[i] = momresolV_truth.at(i);
	err_sigma_p_truth[i] = err_momresolV_truth.at(i);
	err_p_truth[i] = 0.;
	}
	
        const int size_real = momV_real.size();
	double p_real[size_real], err_p_real[size_real], sigma_p_real[size_real], err_sigma_p_real[size_real]; 
	
	for (int i=0; i<size_real; i++){
	p_real[i] = momV_real.at(i);
        sigma_p_real[i] = momresolV_real.at(i);
	err_sigma_p_real[i] = err_momresolV_real.at(i);
	err_p_real[i] = 0.;
	}

  TFile *fout = new TFile(Form("Final_Results/%s/mom/mom_resol_%1.1f_eta_%1.1f.root",particle.Data(),etamin,etamax),"recreate");
	TGraphErrors *gr1 = new TGraphErrors(size_truth,p_truth,sigma_p_truth,err_p_truth,err_sigma_p_truth);
        gr1->SetName("gr_truthseed");
	gr1->SetMarkerStyle(25);
	gr1->SetMarkerColor(kBlue);
	gr1->SetMarkerSize(2.0);
	gr1->SetTitle(";p (GeV/c);#sigmap/p");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	
        TGraphErrors *gr2 = new TGraphErrors(size_real,p_real,sigma_p_real,err_p_real,err_sigma_p_real);
        gr2->SetName("gr_realseed");
	gr2->SetMarkerStyle(34);
	gr2->SetMarkerColor(kRed);
	gr2->SetMarkerSize(2.0);
	gr2->SetTitle(";p (GeV/c);#sigmap/p");
	gr2->GetXaxis()->CenterTitle();
	gr2->GetYaxis()->CenterTitle();
	
	mgMom->Add(gr1);
	mgMom->Add(gr2);
	c_mom->cd();
	mgMom->GetXaxis()->SetRangeUser(0.40,20.2);
	mgMom->GetYaxis()->SetRangeUser(0.0,1.50*TMath::MaxElement(gr2->GetN(),gr2->GetY())); // 50% more of the maximum value on yaxis
	mgMom->Draw("AP");
	lmom->AddEntry(gr1,"Truth Seeding");
	lmom->AddEntry(gr2,"Realistic Seeding");
	lmom->Draw("same");
	draw_req_Mom(etamin,etamax,0.,mgMom->GetXaxis()->GetXmax());
	c_mom->SaveAs(Form("Final_Results/%s/mom/mom_resol_%1.1f_eta_%1.1f.png",particle.Data(),etamin,etamax));
	
	// Write the numbers in output file for comparisons
	outfile << extra_legend << endl;
  outfile<<"Etamin"<<setw(20)<<"Etamax"<<setw(20)<<"p (GeV/c) \t"<<setw(20)<<"Resol  #mum (Truth)"<<setw(20)<<"Resol #mum (Real)"<<endl;
  for (Int_t i = 0; i<gr1->GetN(); ++i){
  double x,ytrue, yreal;
  gr1->GetPoint(i,x,ytrue);    gr2->GetPoint(i,x,yreal);  
  outfile<<etamin<<setw(20)<<etamax<<setw(20)<<x<<setw(20)<<ytrue<<setw(20)<<yreal<<endl;
  }
  outfile.close();
	
	fout->cd();
	mgMom->SetName(Form("mom_resol_%1.1f_eta_%1.1f",etamin,etamax));
	mgMom->Write();
	fout->Close();
}

//===Fit Momentum Resolution
float FitMomentumResolution(Double_t *x, Double_t *par)
{
  float func = sqrt(par[0]*par[0]*x[0]*x[0]+par[1]*par[1]);
  return func;
}

//From Yellow report from section 11.2.2

void draw_req_Mom(double etamin, double etamax, double xmin=0., double xmax=0.)
{

   TF1 *dd4hep_p;
   if (etamin >= -3.5 && etamax <= -2.5) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.1*x)^2+2.0^2)",xmin,xmax);
   else if (etamin >= -2.5 && etamax <= -1.0) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.05*x)^2+1.0^2)",xmin,xmax);
   else if (etamin >= -1.0 && etamax <= 1.0) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.05*x)^2+0.5^2)",xmin,xmax);
   else if (etamin >= 1.0 && etamax <= 2.5) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.05*x)^2+1.0^2)",xmin,xmax);
   else if (etamin >= 2.5 && etamax <= 3.5) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.1*x)^2+2.0^2)",xmin,xmax);
   else return;
   dd4hep_p->SetLineStyle(7);
   dd4hep_p->SetLineColor(kMagenta);
   dd4hep_p->SetLineWidth(3.0);
   dd4hep_p->Draw("same");

  TLegend *l= new TLegend(0.70,0.75,0.90,0.80);
  l->SetTextSize(0.03);
  l->SetBorderSize(0);
  l->AddEntry(dd4hep_p,"PWGReq","l");
  l->Draw("same");
 }
