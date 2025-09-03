void PlotMaterialScan2D(string file="g4mat_scan.csv"){
	TCanvas* c = new TCanvas("c", "Material Scan",0,0,700,600);
	c->SetRightMargin(0.21); // Z axis won't get cut off now
	c->SetLeftMargin(0.11);
	
	TGraph2D* g = new TGraph2D();
	g->SetTitle("Material Scan; #eta; #phi (#circ); X/X_{0}");
	
	gStyle->SetTitleFontSize(0.06);  //title size and font
    gStyle->SetTitleFont(42, "");     
	gStyle->SetPalette(kRainBow);     // slash out for default colours
	
	ifstream inputFile(file);
	string MaterialNames[300];

	int N_phi_points = 0;	
	string xxx;
	float eta, phi, air;
	getline(inputFile,xxx);
	stringstream ss(xxx);
	string mat;
	int i=0, Nmat=0;
	while(getline(ss,mat,'\t')){
		if(i>3){
			MaterialNames[i-4]=mat;
			Nmat++;
		}
		i++;
	}
	cout << "Materials: ";
	for(i=0; i<Nmat; i++)
		cout << MaterialNames[i] << ", ";
	cout << endl;
	inputFile >> xxx;
	int j=0;
	i=0;
	float max=0.0, avcX0=0.0;
	float eta_at_max;
	float X0_by_material[300];
	float av_X0_by_material[300];
	for(i=0; i<300; ++i) av_X0_by_material[i]=0.0;
	i=0;
	do{
		inputFile >> eta >> phi >> air; // >> xxx;  // skip air value
//		cout << "eta: " << eta << ", phi: " << phi;
		float cumX0=0.0;
		for (int k=0; k<Nmat; ++k){
			inputFile >> X0_by_material[k];
			cumX0 += X0_by_material[k];
		}
		inputFile >> xxx;
		avcX0 += cumX0;
		
//		cout << eta << "\t" << phi << "\t" << cumX0 << endl;
		g -> SetPoint(i++,eta,phi,cumX0);
		if(cumX0 > max){
			max = cumX0;
			eta_at_max = eta;
		}
	} while (xxx != "lambda");
	inputFile.close();
	g->SetNpy(360);
	g->SetNpx(268);
	g->Draw("COLZ");

    gPad->Update();  

    // Style axes 
    TH1 *h2 = g->GetHistogram();
    if (h2) {
        h2->GetXaxis()->SetTitleSize(0.06);
        h2->GetYaxis()->SetTitleSize(0.06);
        h2->GetZaxis()->SetTitleSize(0.06);

        h2->GetXaxis()->SetLabelSize(0.045);
        h2->GetYaxis()->SetLabelSize(0.05);
        h2->GetZaxis()->SetLabelSize(0.05);

        h2->GetXaxis()->SetTitleOffset(0.8);
        h2->GetYaxis()->SetTitleOffset(0.8);
       // h2->GetZaxis()->SetTitleOffset(1);
		h2->GetZaxis()->SetRangeUser(0.0, max*0.35); // Set you z-axis range
    }

    // Style the Z-axis color palette
    TPaletteAxis *pal = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
    if (pal) {
        pal->SetLabelSize(0.05);
        pal->SetTitleSize(0.06);
        pal->SetTitleOffset(1.4);
    }

    c->Modified();
    c->Update();
	cout << "Maximum material thickness X/X0 = " << max << " at eta = " << eta_at_max << endl;
}
