void TestOB(){
  TFile* eicFile = new TFile("reconTest.edm4eic.root");
  TTree* eicTree = (TTree*)eicFile->Get("events");
  eicTree->Draw("CentralCKFTrajectories.nMeasurements>>nM");
  auto h0 = (TH1I*)gPad->GetPrimitive("nM");
  cout << "Mean nMeasurements: " << h0->GetMean() << endl;
}
