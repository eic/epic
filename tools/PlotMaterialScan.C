void PlotMaterialScan(string file="g4mat_scan.csv"){
	TCanvas* c = new TCanvas("c", "Material Scan",0,0,700,600);
	c->SetLeftMargin(0.15);  // increasing the canvas size in the given direction
	c->SetBottomMargin(0.15); 
	c->SetGrid();

	gStyle->SetTitleFontSize(0.06);   // increase main title font size
	gStyle->SetTitleFont(42, "");     // optional: change font

	TGraph* g1d = new TGraph();
	g1d->SetTitle("Material Scan; #eta; X/X_{0}"); // edit title name for L3 or L4
	g1d->GetXaxis()->SetTitleSize(0.06);       // title and label sizes for the axis
	g1d->GetYaxis()->SetTitleSize(0.06);
	g1d->GetXaxis()->SetLabelSize(0.05);
	g1d->GetYaxis()->SetLabelSize(0.05);
	//g1d->GetXaxis()->SetLimits(-1,1);   // set eta range
	g1d->GetXaxis()->SetTitleOffset(0.8); // Moves the axis title
	g1d->GetYaxis()->SetTitleOffset(1.29);

	TGraph* graphs_by_material[300];
	ifstream inputFile(file);
	string MaterialNames[300];
	int N_phi_points = 0;	
	float eta,phi, air,prev_eta;
	string xxx;
	getline(inputFile,xxx);
	stringstream ss(xxx);
	string mat;
	int i=0, Nmat=0;
	while(getline(ss,mat,'\t')){
		if(i>3){
			MaterialNames[i-4]=mat;
			graphs_by_material[i-4] = new TGraph();
			graphs_by_material[i-4]->SetName(mat.c_str());
			Nmat++;
		}
		i++;
	}
	cout << "Materials: ";
	for(i=0; i<Nmat; i++)
		cout << MaterialNames[i] << ", ";
	cout << endl;
	float avcX0=0.0;
	float X0_by_material[300];
	float av_X0_by_material[300];
	for(i=0; i<300; ++i) 
		av_X0_by_material[i]=0.0;
	i=0;
	inputFile >> xxx >> eta >> phi >> air; 
	do{
		if(phi==85 && N_phi_points !=0){   //change phi==0 for either L3 | 85 or L4 | 87 for stave material scan
			for(int k=0; k<Nmat; ++k){
				av_X0_by_material[k] /= (float)N_phi_points;
				graphs_by_material[k]->SetPoint(i,prev_eta,av_X0_by_material[k]);
			}	
			avcX0 /= (float)N_phi_points;
			g1d->SetPoint(i++,prev_eta,avcX0); 
			g1d->SetLineWidth(2);
			N_phi_points=0;		
		}

		float cumX0=0.0;
		for (int k=0; k<Nmat; ++k){
			inputFile >> X0_by_material[k];
			cumX0 += X0_by_material[k];
		}
		avcX0 += cumX0;
		N_phi_points++;
		for(int k=0; k<Nmat; ++k){ 
			av_X0_by_material[k] +=  X0_by_material[k];
		}
		prev_eta=eta;
		inputFile >> xxx >> eta >> phi >> air; 
	} 	
	while (xxx != "lambda");
		inputFile.close();
        g1d->Draw();

	// custom legend colors
	int darkGreen = TColor::GetColor(0,200,0);   
	int colors[] = {kBlue, kRed, darkGreen, kMagenta, kCyan, kOrange};
	int Ncolors = sizeof(colors)/sizeof(int);

	// change legend size and position for your liking
	TLegend* legend = new TLegend(0.38,0.60,0.65,0.85);

	legend->SetTextSize(0.035);
	legend->SetTextFont(42); 
		for(int k=0; k<Nmat; ++k){
			graphs_by_material[k]->SetLineColor(colors[k % Ncolors]); //change to default colors with k+2 in colors[]
			graphs_by_material[k]->SetLineWidth(2);
			graphs_by_material[k]->Draw("L");
			legend->AddEntry(graphs_by_material[k],MaterialNames[k].c_str(),"l");
	}
	legend->AddEntry("g1d","All","l");
	legend->Draw();
}
