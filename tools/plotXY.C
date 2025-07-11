void plotXY(string infile="debugVol.txt") {

   ifstream in;
   in.open(infile);

   Double_t x,y,z;
   
   TGraph* g1 = new TGraph();
TGraph* g2 = new TGraph();
   
   int np1=0, np2=0;
   
   Double_t minX=9999, maxX=-9999, minY=9999, maxY=-9999,minZ=9999,maxZ=-9999,sumX=0, sumY=0, sumZ=0;

   while (1) {
      in >> x >> y >> z;
      if (!in.good()) break;
      if (np1 < 5 && np2<5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
      if (x<minX) minX=x;
      if(x>maxX) maxX=x;
      if(y<minY) minY=y;
      if(y>maxY) maxY=y;
      if(z<minZ) minZ=z;
      if(z>maxZ) maxZ=z;
      sumX+=x; sumY+=y; sumZ+=z;
          g1->SetPoint(np1,x,z);     
          np1++;
 /*     
      if(y<0){
 
 g1->SetPoint(np1,z,x);
      np1++;
      }
      else {
g2->SetPoint(np2,z,x);
      np2++;
      
      }
      */
   }

   in.close();
   
   cout << "\tX\tY\tZ\n";
   cout << "Min\t" << minX << "\t" << minY << "\t" << minZ << endl;
   cout << "Max\t" << maxX << "\t" << maxY << "\t" << maxZ << endl;
cout << "Average\t" << sumX/Double_t(np1) << "\t" << sumY/Double_t(np1) << "\t" << sumZ/Double_t(np1)  << endl;   
TCanvas *c1 = new TCanvas("c1","hist",200,10,1400,1000);
  gStyle->SetOptStat(0);
//  g1->SetMarkerColor(3);
g1->Draw("AP");
//g2->SetMarkerColor(2);
//g2->Draw("P");
}

