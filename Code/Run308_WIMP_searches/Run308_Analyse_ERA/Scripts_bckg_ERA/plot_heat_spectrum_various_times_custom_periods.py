#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import os
import numpy as np
import Analysis_utilities as Ana_ut

def launch_plots(bolo_name, data_dir, tree_name = "data"):


    #Load estimators
    estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
    estimator_file_name =bolo_name + "_estimators.txt" 
    d_estimator         ={}
    assert(os.path.isfile(estimator_path_name + estimator_file_name) )
    with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
        list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
        for line in list_estimator_lines:
            d_estimator[line[0]] = line[1]

    # Load start and end times
    tmin, tmax = Ana_ut.open_start_and_end_file(bolo_name)
    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_estimator["HEAT"] + ">0"]

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)
    print cut_line

    bin_X, min_X, max_X = 50, 0, 6
    ytitle = str(float(max_X-min_X)/float(bin_X))

    arr_time =  [0,65,110,150,200]
         
    cut1 = cut_line + "&& JOUR>" + str(arr_time[0]) + "&& JOUR<" + str(arr_time[1])
    cut2 = cut_line + "&& JOUR>" + str(arr_time[1]) + "&& JOUR<" + str(arr_time[2])
    cut3 = cut_line + "&& JOUR>" + str(arr_time[2]) + "&& JOUR<" + str(arr_time[3])
    cut4 = cut_line + "&& JOUR>" + str(arr_time[3]) + "&& JOUR<" + str(arr_time[4])

    h1= TH1F("h1", "h1", bin_X, min_X, max_X)
    h2= TH1F("h2", "h2", bin_X, min_X, max_X)
    h3= TH1F("h3", "h3", bin_X, min_X, max_X)
    h4= TH1F("h4", "h4", bin_X, min_X, max_X)

    hnorm1= TH1F("hnorm1", "hnorm1", bin_X, min_X, max_X)
    hnorm2= TH1F("hnorm2", "hnorm2", bin_X, min_X, max_X)
    hnorm3= TH1F("hnorm3", "hnorm3", bin_X, min_X, max_X)
    hnorm4= TH1F("hnorm4", "hnorm4", bin_X, min_X, max_X)


    tree.Project("h1","0.5*(EC1+EC2)", cut1)
    tree.Project("h2","0.5*(EC1+EC2)", cut2)
    tree.Project("h3","0.5*(EC1+EC2)", cut3)
    tree.Project("h4","0.5*(EC1+EC2)", cut4)

    tree.Project("hnorm1","0.5*(EC1+EC2)", cut1)
    tree.Project("hnorm2","0.5*(EC1+EC2)", cut2)
    tree.Project("hnorm3","0.5*(EC1+EC2)", cut3)
    tree.Project("hnorm4","0.5*(EC1+EC2)", cut4)


    h1.SetLineColor(kRed)
    h2.SetLineColor(kOrange-8)
    h3.SetLineColor(kViolet+1)
    h4.SetLineColor(kGreen-3)

    #Save to TFile

    #Normalise to counts/keV
    h1.Scale(1./hnorm1.GetBinWidth(10))
    h2.Scale(1./hnorm2.GetBinWidth(10))
    h3.Scale(1./hnorm3.GetBinWidth(10))
    h4.Scale(1./hnorm4.GetBinWidth(10))

    g = TFile("./ROOT_files/" + bolo_name + "_countsperkeV_heat_spectrum.root", "recreate")
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    g.Close()

    hnorm1.SetLineColor(kRed)
    hnorm2.SetLineColor(kOrange-8)
    hnorm3.SetLineColor(kViolet+1)
    hnorm4.SetLineColor(kGreen-3)

    hnorm1.Scale(1./hnorm1.Integral())
    hnorm2.Scale(1./hnorm2.Integral())
    hnorm3.Scale(1./hnorm3.Integral())
    hnorm4.Scale(1./hnorm4.Integral())


    list_hist = [h1, h2, h3, h4]
    list_max = [hist.GetMaximum() for hist in list_hist]
    max_max = max(list_max)

    h= TH1F("h", "h", bin_X, min_X, max_X)
    h.SetMaximum(1.2*max_max)
    h.SetStats(0)

    h.GetXaxis().SetTitle("Heat (keVee)")
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    # h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle("Counts/" + ytitle + "keV")
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    # h.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    gPad.SetLogy()
    h.SetMaximum(2E4)
    h.Draw()   
    for hist in list_hist:
        hist.Draw("same")
    leg = TLegend(0.46,0.51,0.65,0.87)
    leg.AddEntry("h1", "Period 1", "leg")
    leg.AddEntry("h2", "Period 2", "leg")
    leg.AddEntry("h3", "Period 3", "leg")
    leg.AddEntry("h4", "Period 4", "leg")
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(0)
    leg.Draw("same")
    raw_input()
    # cc.Print("../Analyse_FID837/Figures/FID837_heat_spec_over_time_lowmass.eps")


    #Do the plots
    cc2=TCanvas("cc2","cc2")
    gPad.SetLogy()
    list_hist = [hnorm1, hnorm2, hnorm3, hnorm4]
    list_max = [hist.GetMaximum() for hist in list_hist]
    max_max = max(list_max)
    h.SetMaximum(1.2*max_max)
    h.GetYaxis().SetTitle("Arbitrary counts")

    #Save to TFile
    f = TFile("./ROOT_files/" + bolo_name + "_normalised_heat_spectrum.root", "recreate")
    hnorm1.Write()
    hnorm2.Write()
    hnorm3.Write()
    hnorm4.Write()
    f.Close()

    h.Draw()   
    for hist in list_hist:
        hist.Draw("same")
    leg = TLegend(0.46,0.51,0.65,0.87)
    leg.AddEntry("hnorm1", "Period 1", "leg")
    leg.AddEntry("hnorm2", "Period 2", "leg")
    leg.AddEntry("hnorm3", "Period 3", "leg")
    leg.AddEntry("hnorm4", "Period 4", "leg")
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(0)
    leg.Draw("same")
    raw_input()
    # cc2.Print("../Analyse_FID837/Figures/FID837_normalised_heat_spec_over_time_lowmasss.eps")


bolo_name, data_dir = "FID837", "../Fond_ERA_merged/"
launch_plots(bolo_name, data_dir)