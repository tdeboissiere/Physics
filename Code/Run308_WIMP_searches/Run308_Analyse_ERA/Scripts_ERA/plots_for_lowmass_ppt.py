#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import PyROOTPlots as PyRPl
import Analysis_utilities as Ana_ut
import numpy as np

def get_KTH_plot(bolo_name, data_dir, tree_name):

    """
    Do the plot of KTH vs time
    Also estimate the mean KTH

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the tree directory (containing Jules n-tuple)
    tree_name (str) the tree name in the tree file

    Outputs:

    """
    
    #Get tree 
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)
    nEntries    = tree.GetEntries()

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")
    l_kth = TEventList("l_kth", "l_kth")
    tree.Draw(">>l_kth",standard_cuts)
    pop_len = l_kth.GetN()
    list_kth = []
    for k in range(pop_len):
        counter = l_kth.GetEntry(k)
        tree.GetEntry(counter)
        list_kth.append(tree.KTH)

    print np.average(np.array(list_kth))

    #Initialize the start and end time 
    tree.GetEntry(0)
    all_tree_start_time =str(tree.JOUR)
    tree.GetEntry(nEntries-1)
    all_tree_end_time   =str(tree.JOUR)
    
    l_KTH = TLine(float(all_tree_start_time), 1, float(all_tree_end_time), 1)
    l_KTH.SetLineColor(kRed)
    l_KTH.SetLineWidth(2)

    c_KTH=TCanvas("c_KTH", "c_KTH")

    tree.Draw("KTH:JOUR>>KTH_hist(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,2)")
    PyRPl.process_TH2(KTH_hist, Y_title = "KTH", X_title = "Days")
    KTH_hist.Draw()
    l_KTH.Draw("same")

    raw_input()

    c_KTH.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_KTH_over_time.png")


def get_FWHM_plots(bolo_name, data_dir, tree_name):

    """
    Do the plot of baselines vs time

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the tree directory (containing Jules n-tuple)
    tree_name (str) the tree name in the tree file

    Outputs:

    """
    
    #Get tree 
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)
    nEntries    = tree.GetEntries()

    #Initialize the start and end time 
    tree.GetEntry(0)
    all_tree_start_time =str(tree.JOUR)
    tree.GetEntry(nEntries-1)
    all_tree_end_time   =str(tree.JOUR)
    
    l_ion = TLine(float(all_tree_start_time), 1.5, float(all_tree_end_time), 1.5)
    l_heat = TLine(float(all_tree_start_time), 1, float(all_tree_end_time), 1)

    l_ion.SetLineColor(kRed)
    l_heat.SetLineColor(kRed)

    l_ion.SetLineWidth(2)
    l_heat.SetLineWidth(2)

    c_fwhm_fid=TCanvas("c_fwhm_fid", "c_fwhm_fid")
    c_fwhm_fid.Divide(1,2)

    c_fwhm_fid.cd(1)
    tree.Draw("FWIB:JOUR>>fwhist2(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,4)")
    PyRPl.process_TH2(fwhist2, Y_title = "FWIB")
    fwhist2.Draw()
    l_ion.Draw("same")
    c_fwhm_fid.cd(2)
    tree.Draw("FWID:JOUR>>fwhist4(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,4)")
    PyRPl.process_TH2(fwhist4, Y_title = "FWID", X_title = "Days")
    fwhist4.Draw()
    l_ion.Draw("same")
    
    c_fwhm_vet=TCanvas("c_fwhm_vet", "c_fwhm_vet")
    c_fwhm_vet.Divide(1,2)
    c_fwhm_vet.cd(1)
    tree.Draw("FWIA:JOUR>>fwhist1(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,4)")
    PyRPl.process_TH2(fwhist1, Y_title = "FWIA")
    fwhist1.Draw()
    l_ion.Draw("same")
    c_fwhm_vet.cd(2)
    tree.Draw("FWIC:JOUR>>fwhist3(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,4)")
    PyRPl.process_TH2(fwhist3, Y_title = "FWIC", X_title = "Days")
    fwhist3.Draw()
    l_ion.Draw("same")


    c_fwhm_heat=TCanvas("c_fwhm_heat", "c_fwhm_heat")
    c_fwhm_heat.Divide(1,2)

    c_fwhm_heat.cd(1)
    tree.Draw("OWC1:JOUR>>fwhist5(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,4)")
    PyRPl.process_TH2(fwhist5, Y_title = "OWC1")
    fwhist5.Draw()
    l_heat.Draw("same")
    c_fwhm_heat.cd(2)
    tree.Draw("OWC2:JOUR>>fwhist6(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,0,4)")
    PyRPl.process_TH2(fwhist6, Y_title = "OWC2", X_title = "Days")
    fwhist6.Draw()
    l_heat.Draw("same")

    raw_input()

    c_fwhm_fid.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_collectrode_baselines.png")
    c_fwhm_vet.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_veto_baselines.png")
    c_fwhm_heat.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_heat_baselines.png")

bolo_name = "FID837"
data_dir = "../Fond_ERA_merged/"
tree_name = "data"
# get_FWHM_plots(bolo_name, data_dir, tree_name)
get_KTH_plot(bolo_name, data_dir, tree_name)