#!/usr/bin/env python

import script_utils as script_utils
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh
from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np 
import matplotlib.pylab as plt
import prettyplotlib as ppl

def plot_EC1_EC2(bolo_name, data_dir, tree_name):

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")
    standard_cuts = standard_cuts.split("&&")
    standard_cuts.remove("KTH>0")
    standard_cuts.remove("KTH<1")
    standard_cuts.remove("abs(EC1-EC2)<1")
    standard_cuts = " && ".join(standard_cuts)
    # raw_input()

    fidcut = "&& 0.5*(EIB+EID)>1 && EIA<1 && EIC<1 && Q>0.6"
    heatcut = "&& 0.5*(EIB+EID)<0.7 && EIA<2 && EIC<2"   

    # h2Dfid = TH2F("h2Dfid", "h2Dfid", 2000, 0, 20., 2000, 0, 20.)    
    # h2Dheat = TH2F("h2Dheat", "h2Dheat", 2000, 0, 20., 2000, 0, 20.)

    # tree.Project("h2Dfid", "EC2:EC1", standard_cuts + fidcut)    
    # tree.Project("h2Dheat", "EC2:EC1", standard_cuts + heatcut)

    # PyRPl.process_TH2(h2Dfid, X_title = "EC1", Y_title = "EC2", marker_color = kBlack)    
    # PyRPl.process_TH2(h2Dheat, X_title = "EC1", Y_title = "EC2", color = kRed)

    # l = TLine(0,0,20.,20.)
    # l.SetLineWidth(2)

    # c1 = TCanvas("c1", "c1")
    # h2Dfid.Draw()
    # h2Dheat.Draw("same")
    # l.Draw("same")
    # c1.Print("./Figures/FID837_EC1_EC2_heatonly_vs_fid.eps")
    # raw_input()

    h2Dfidchi = TH2F("h2Dfidchi", "h2Dfidchi", 2000, 0, 20., 2000, -0.2, 0.2)    
    h2Dheatchi = TH2F("h2Dheatchi", "h2Dheatchi", 2000, 0, 20., 2000, -0.2, 0.2)

    tree.Project("h2Dfidchi", "XOC1:EC", standard_cuts + fidcut)    
    tree.Project("h2Dheatchi", "XOC1:EC", standard_cuts + heatcut)

    PyRPl.process_TH2(h2Dfidchi, X_title = "EC", Y_title = "XOC1", marker_color = kBlack)    
    PyRPl.process_TH2(h2Dheatchi, X_title = "EC", Y_title = "XOC1", color = kRed)

    l = TLine(0,0,20.,20.)
    l.SetLineWidth(2)

    c2 = TCanvas("c2", "c2")
    h2Dfidchi.Draw()
    h2Dheatchi.Draw("same")
    l.Draw("same")
    c2.Print("./Figures/FID837_EC_XOC1_heatonly_vs_fid.eps")
    raw_input()


plot_EC1_EC2("FID837", "../Fond_ERA_merged/", "data")