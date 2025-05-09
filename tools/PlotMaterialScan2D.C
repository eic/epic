void PlotMaterialScan2D(string file = "g4mat_scan.csv") {
  TCanvas* c  = new TCanvas("c", "Material Scan", 0, 0, 700, 600);
  TGraph2D* g = new TGraph2D();
  g->SetTitle("Material Scan; eta; phi; X/X0");
  ifstream inputFile(file);
  string MaterialNames[300];

  int N_phi_points = 0;
  string xxx;
  float eta, phi, air;
  getline(inputFile, xxx);
  stringstream ss(xxx);
  string mat;
  int i = 0, Nmat = 0;
  while (getline(ss, mat, '\t')) {
    if (i > 3) {
      MaterialNames[i - 4] = mat;
      Nmat++;
    }
    i++;
  }
  cout << "Materials: ";
  for (i = 0; i < Nmat; i++)
    cout << MaterialNames[i] << ", ";
  cout << endl;
  inputFile >> xxx;
  int j     = 0;
  i         = 0;
  float max = 0.0, avcX0 = 0.0;
  float eta_at_max;
  float X0_by_material[300];
  float av_X0_by_material[300];
  for (i = 0; i < 300; ++i)
    av_X0_by_material[i] = 0.0;
  i = 0;
  do {
    inputFile >> eta >> phi >> air; // >> xxx;  // skip air value
                                    //          cout << "eta: " << eta << ", phi: " << phi;
    float cumX0 = 0.0;
    for (int k = 0; k < Nmat; ++k) {
      inputFile >> X0_by_material[k];
      cumX0 += X0_by_material[k];
    }
    inputFile >> xxx;
    avcX0 += cumX0;

    //          cout << eta << "\t" << phi << "\t" << cumX0 << endl;
    g->SetPoint(i++, eta, phi, cumX0);
    if (cumX0 > max) {
      max        = cumX0;
      eta_at_max = eta;
    }
  } while (xxx != "lambda");
  inputFile.close();
  g->SetNpy(360);
  g->SetNpx(268);
  g->Draw("COLZ");
  cout << "Maximum material thickness X/X0 = " << max << " at eta = " << eta_at_max << endl;
}
