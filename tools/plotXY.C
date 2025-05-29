void plotXY(string infile="debugVol.txt") {

   ifstream in;
   in.open(infile);

   Double_t x,y,z;
   
   TGraph* g1 = new TGraph();
   int np=0;

   while (1) {
      in >> x >> y >> z;
      if (!in.good()) break;
      if (np < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
           g1->SetPoint(np,x,y);
// g1->SetPoint(np,z,x);
      np++;
   }

   in.close();
TCanvas *c1 = new TCanvas("c1","hist",200,10,1400,1000);
  gStyle->SetOptStat(0);
g1->Draw("AP");
}

