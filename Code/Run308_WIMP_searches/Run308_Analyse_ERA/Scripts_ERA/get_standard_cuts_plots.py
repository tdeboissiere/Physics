#!/usr/bin/env python

from ROOT import TCanvas, TFile, TTree, TH2F
import script_utils as script_utils
import os
import PyROOTPlots as PyRPl

def get_standard_cuts(bolo_name, data_dir, tree_name):

    """
    Creates .png files for the plots of the FWHM and CHI2 over time (resp. energy)

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the tree directory (containing Jules n-tuple)
    tree_name (str) the tree name in the tree file

    Outputs:

    .png files with the corresponding plots
    """
    
    #Get tree 
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)
    nEntries    = tree.GetEntries()

    #Initialize the start and end time 
    tree.GetEntry(0)
    all_tree_start_time =str(1E6*tree.UT1+tree.UT2)
    all_tree_end_time   =str(1E6*tree.UT1+tree.UT2)
    
    #Find the start and end time of the whole data tree
    list_time = []
    for i in range(1, nEntries):
        tree.GetEntry(i)
        time= 1E6*tree.UT1+tree.UT2
        list_time.append(time)
    
    all_tree_start_time = float(min(list_time))
    all_tree_end_time   = float(max(list_time))


    fwhist1 = TH2F("fwhist1", "fwhist1", 10000,float(all_tree_start_time),float(all_tree_end_time),1000,0,4)
    fwhist2 = TH2F("fwhist2", "fwhist2", 10000,float(all_tree_start_time),float(all_tree_end_time),1000,0,4)
    fwhist3 = TH2F("fwhist3", "fwhist3", 10000,float(all_tree_start_time),float(all_tree_end_time),1000,0,4)
    fwhist4 = TH2F("fwhist4", "fwhist4", 10000,float(all_tree_start_time),float(all_tree_end_time),1000,0,4)
    fwhist5 = TH2F("fwhist5", "fwhist5", 10000,float(all_tree_start_time),float(all_tree_end_time),1000,0,4)
    fwhist6 = TH2F("fwhist6", "fwhist6", 10000,float(all_tree_start_time),float(all_tree_end_time),1000,0,4)

    c_fwhm_ion=TCanvas("c_fwhm_ion", "c_fwhm_ion")
    c_fwhm_ion.Divide(2,2)

    tree.Project("fwhist1", "FWIA:1E6*UT1+UT2")
    tree.Project("fwhist2", "FWIB:1E6*UT1+UT2")
    tree.Project("fwhist3", "FWIC:1E6*UT1+UT2")
    tree.Project("fwhist4", "FWID:1E6*UT1+UT2")

    PyRPl.process_TH2(fwhist1, X_title = "Time", Y_title = "FWIA")
    PyRPl.process_TH2(fwhist2, X_title = "Time", Y_title = "FWIB")
    PyRPl.process_TH2(fwhist3, X_title = "Time", Y_title = "FWIC")
    PyRPl.process_TH2(fwhist4, X_title = "Time", Y_title = "FWID")

    c_fwhm_ion.cd(1)
    fwhist1.Draw()
    c_fwhm_ion.cd(2)
    fwhist2.Draw()
    c_fwhm_ion.cd(3)
    fwhist3.Draw()
    c_fwhm_ion.cd(4)
    fwhist4.Draw()

    raw_input()
    c_fwhm_ion.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_ion_FWHM_over_time.png")


    c_fwhm_heat=TCanvas("c_fwhm_heat", "c_fwhm_heat")
    c_fwhm_heat.Divide(1,2)

    tree.Project("fwhist5", "FWC1:1E6*UT1+UT2")
    tree.Project("fwhist6", "FWC2:1E6*UT1+UT2")

    PyRPl.process_TH2(fwhist5, X_title = "Time", Y_title = "FWC1")
    PyRPl.process_TH2(fwhist6, X_title = "Time", Y_title = "FWC2")

    c_fwhm_heat.cd(1)
    fwhist5.Draw()
    c_fwhm_heat.cd(2)
    fwhist6.Draw()

    raw_input()
    c_fwhm_heat.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_heat_FWHM_over_time.png")



    bolo_TCut_line = "FWIA<2&&FWIB<2&&FWIC<2&&FWID<2&&FWC1<2&&FWC2<2&&SDEL>1"
    bolo_TCut_line += "&&FWIA>0&&FWIB>0&&FWIC>0&&FWID>0&&FWC1>0&&FWC2>0"


    chist1 = TH2F("chist1", "chist1", 200,0,8,200,0,8)
    chist2 = TH2F("chist2", "chist2", 200,0,8,200,0,8)
    chist3 = TH2F("chist3", "chist3", 200,0,8,200,0,8)
    chist4 = TH2F("chist4", "chist4", 200,0,8,200,0,8)
    chist5 = TH2F("chist5", "chist5", 200,0,8,200,0,8)
    chist6 = TH2F("chist6", "chist6", 200,0,8,200,0,8)

    tree.Project("chist1", "CHIA:0.5*(EC1+EC2)", bolo_TCut_line)
    tree.Project("chist2", "CHIB:0.5*(EC1+EC2)", bolo_TCut_line)
    tree.Project("chist3", "CHIC:0.5*(EC1+EC2)", bolo_TCut_line)
    tree.Project("chist4", "CHID:0.5*(EC1+EC2)", bolo_TCut_line)

    PyRPl.process_TH2(chist1, X_title = "Heat (keV)", Y_title = "CHIA")
    PyRPl.process_TH2(chist2, X_title = "Heat (keV)", Y_title = "CHIB")
    PyRPl.process_TH2(chist3, X_title = "Heat (keV)", Y_title = "CHIC")
    PyRPl.process_TH2(chist4, X_title = "Heat (keV)", Y_title = "CHID")

    c_chi2_ion_of_energ=TCanvas("c_chi2_ion_of_energ", "c_chi2_ion_of_energ")
    c_chi2_ion_of_energ.Divide(2,2)

    c_chi2_ion_of_energ.cd(1)
    chist1.Draw("col")
    c_chi2_ion_of_energ.cd(2)
    chist2.Draw("col")
    c_chi2_ion_of_energ.cd(3)
    chist3.Draw("col")
    c_chi2_ion_of_energ.cd(4)
    chist4.Draw("col")

    raw_input()
    c_chi2_ion_of_energ.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_ion_chi2_over_energ.png")


    c_chi2_heat_of_energ=TCanvas("c_chi2_heat_of_energ", "c_chi2_heat_of_energ")
    c_chi2_heat_of_energ.Divide(1,2)

    tree.Project("chist5", "CHIC1:0.5*(EC1+EC2)", bolo_TCut_line)
    tree.Project("chist6", "CHIC2:0.5*(EC1+EC2)", bolo_TCut_line)

    PyRPl.process_TH2(chist1, X_title = "Heat (keV)", Y_title = "CHIC1")
    PyRPl.process_TH2(chist2, X_title = "Heat (keV)", Y_title = "CHIC2")


    c_chi2_heat_of_energ.cd(1)
    chist5.Draw("col")
    c_chi2_heat_of_energ.cd(2)
    chist6.Draw("col")

    raw_input()
    c_chi2_heat_of_energ.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_heat_chi2_over_energ.png")




get_standard_cuts("FID837", "../Fond_ERA_merged/", "data")




