# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		PlotMaterialScan.py
#   	 	Author：		LongLI <l.li.8@bham.ac.uk> <long.l@cern.ch>
#   		Time：			2024.11.08
#   		Description：
#
#======================================================================


from ROOT import *
import math
import array
import numpy as np
import os
import sys
import pandas as pd
import argparse
from array import array

def prepare_data(filename):
    
    
    isfirst = True
    # calculate the dimension of the data array
    count_eta = 0
    count_phi = 0
    pre_eta, pre_phi = None, None
    eta = set()
    phi = set()    
    data_list = []
    matName = []
    
    f = open(filename)
    # calculate the dimension of the data array
    count_eta = 0
    count_phi = 0
    count_mat = 0
    pre_eta, pre_phi = None, None
    
    for eachline in f:
        if isfirst:
            items = eachline.strip().split('\t')
            count_mat = len(items) - 4 # material air excludes
            for idx in range(count_mat):
                matName.append(items[idx+4])
            isfirst = False

        else:
            items = eachline.strip().split('\t')
            if(count_mat == 0): exit(-1)

            if items[0] == 'X0':
                items = eachline.strip().split('\t')
                if(count_mat == 0): exit(-1)
                for idx in range(count_mat):
                    data_list.append(float(items[idx+4]))
                
                if pre_eta == None or pre_phi == None:
                    count_eta+=1
                    count_phi+=1
                    pre_eta = items[1]
                    pre_phi = items[2]
                    eta.add(float(items[1]))                      
                
                else:
                    # for eta
                    if items[1] != pre_eta:
                        count_eta +=1
                        count_phi = 0
                        pre_eta = items[1]
                        eta.add(float(items[1]))
                        
                    else:
                        pass
                    
                    # for phi
                    if items[2] != pre_phi:
                        count_phi +=1
                        pre_phi = items[2]
                        phi.add(float(items[2]))
                    
                    else:
                        pass

            else:
                continue

    eta = sorted(eta)
    phi = sorted(phi)
    
    print(matName)
    np_eta = np.array(eta)
    np_phi = np.array(phi)
    np_data = np.array(data_list).reshape((count_eta, count_phi, count_mat))
    
    return np_eta, np_phi, np_data, matName    





def PlotMaterialScanEta(eta, data, matName, output):
    
    cvs = TCanvas('cvs', 'Material Scan', 800, 600)
    mg = TMultiGraph()
    
    
    data_array = array('d')
    eta_array = array('d')
    
    data_mean_phi = np.mean(data, axis=1)
    
    for idx in range(np.shape(data_mean_phi)[-1]):
        data_array = data_mean_phi[:,idx]*100
        eta_array = eta
        graph = TGraph(len(data_array), eta_array, data_array)

        graph.SetTitle(matName[idx])
        graph.SetLineColor(idx+2)
        graph.SetLineWidth(2)
        
        mg.Add(graph, 'AL')
    
    data_array = np.sum(data_mean_phi, axis=-1)*100 # unit of %
    graph = TGraph(len(data_array), eta_array, data_array)
    graph.SetTitle('TOT')
    graph.SetLineColor(1)
    graph.SetLineWidth(2)
    mg.Add(graph, 'AL')
        
        
    mg.GetXaxis().SetTitle("#eta")
    mg.GetYaxis().SetTitle("X/X0 [%]")
    mg.GetYaxis().SetRangeUser(0, 2.5)
    mg.Draw('AL')
    
    lgd = cvs.BuildLegend(0.4, 0.65, 0.6, 0.85)

    if not os.path.exists(output):
        os.makedirs(output)
        
    save_file = os.path.join(output, f'MaterialScan_vs_eta.pdf')

    cvs.Print(save_file)

    return eta_array, data_array


def PlotMaterialScan2D(eta, phi, data, output):
    
    gStyle.SetOptStat(0000)
    
    cvs = TCanvas("cvs", "Material Scan", 800, 600)
    cvs.SetRightMargin(0.2)
    
    h2d = TH2D("h2d", ";#eta;#Phi;X/X0 [%]", len(eta), eta[0], eta[-1], len(phi), phi[0], phi[-1])
    
    
    
    for e in range(len(eta)):
        for p in range(len(phi)):
            h2d.SetBinContent(e, p, data[e][p][-1])
    
    
    h2d.GetZaxis().SetTitleOffset(1.3)
    h2d.GetZaxis().SetRangeUser(0, 0.5)
    h2d.Draw("COLZ")
    
    if not os.path.exists(output):
        os.makedirs(output)
        
    save_file = os.path.join(output, f'MaterialScan2D_Silicon.pdf')

    cvs.Print(save_file)
    
    
    data_sum = np.sum(data, axis=-1)*100
    
    for e in range(len(eta)):
        for p in range(len(phi)):
            h2d.SetBinContent(e, p, data_sum[e][p])
    
    h2d.GetZaxis().SetTitleOffset(1.1)
    h2d.GetZaxis().SetRangeUser(0, 2.5)
    h2d.Draw("COLZ")
    
    if not os.path.exists(output):
        os.makedirs(output)
        
    save_file = os.path.join(output, f'MaterialScan2D_TOT.pdf')

    cvs.Print(save_file)
    
def MaterialComparisonEta(data, output):
    
    cvs = TCanvas("cvs", "Material Scan", 800, 600)
    cvs.SetRightMargin(0.2)
    
    mg = TMultiGraph()
    
    for i, item in enumerate(data.items()):
        
        graph = TGraph(len(item[-1][-1]), item[-1][0], item[-1][-1])
        graph.SetLineColor(i+1)
        graph.SetLineWidth(2)
        mg.Add(graph, 'AL')
        
        
    mg.GetXaxis().SetTitle('#eta')
    mg.GetYaxis().SetTitle('X/X0 [%]')
    mg.GetYaxis().SetRangeUser(0, 2.5)
    mg.Draw('AL')
    
    if not os.path.exists(output):
        os.makedirs(output)
        
    save_file = os.path.join(output, f'MaterialComparison.pdf')

    cvs.Print(save_file) 
    
        


def main(args):
    
    data_dict = {}
    
    for cPath, _, files in os.walk(args.data_path):
        
        for item in files:
            np_eta, np_phi, np_data, matName =  prepare_data(os.path.join(cPath, item))
            eta, data = PlotMaterialScanEta(np_eta, np_data, matName, args.output)
            PlotMaterialScan2D(np_eta, np_phi, np_data, args.output)
            
            title = item.split('.')[0]
            data_dict[title] = [eta, data]           
                        
                        
                        
    MaterialComparisonEta(data_dict, args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plotting scripts for the material scan")
    parser.add_argument('--data_path', type=str, default='../csv/', help='path to the raw data file (.csv)')
    parser.add_argument('--output', type=str, default='./plots/', help='output dir')
    
    args = parser.parse_args()
    
    main(args)