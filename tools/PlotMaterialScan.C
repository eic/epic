void PlotMaterialScan(string file="g4mat_scan.csv"){
	TCanvas* c = new TCanvas("c", "Material Scan",0,0,700,600);
	TGraph* g1d = new TGraph();
	TGraph* graphs_by_material[300];
	g1d->SetTitle("Material Scan; eta; X/X0");
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
		if(phi==0 && N_phi_points !=0){
			for(int k=0; k<Nmat; ++k) {
				av_X0_by_material[k] /= (float)N_phi_points;
				graphs_by_material[k]->SetPoint(i,prev_eta,av_X0_by_material[k]);
			}	
			avcX0 /= (float)N_phi_points;
			g1d->SetPoint(i++,eta,avcX0);
			N_phi_points=0;		
		}
		float cumX0=0.0;
		for (int k=0; k<Nmat; ++k){
			inputFile >> X0_by_material[k];
			cumX0 += X0_by_material[k];
		}
		avcX0 += cumX0;
		N_phi_points++;
		for(int k=0; k<Nmat; ++k) 
			av_X0_by_material[k] +=  X0_by_material[k];
		prev_eta=eta;
		inputFile >> xxx >> eta >> phi >> air; 
	} while (xxx != "lambda");
	inputFile.close();
	
        g1d->Draw();
	TLegend* legend = new TLegend(0.1,0.8,0.4,0.95);
        for(int k=0; k<Nmat; ++k){
        	graphs_by_material[k]->SetLineColor(k+2);
	        graphs_by_material[k]->Draw("L");
	        string gr_name = graphs_by_material[k]->GetName();
	        legend->AddEntry(gr_name.c_str(),MaterialNames[k].c_str(),"l");
	}
	legend->AddEntry("g1d","All","l");
	legend->Draw();
}
