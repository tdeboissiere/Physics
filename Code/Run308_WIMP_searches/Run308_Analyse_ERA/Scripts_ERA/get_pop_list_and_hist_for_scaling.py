#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh

def get_pop_list_for_scaling(bolo_name, data_dir, tree_name):

    """Get the population files for scaling
    
    Detail:
        Select each event type (Gamma, Beta, Pb, heatonly)
        also include Surface 1 Surface 2 and Fiducial distinction
        when applicable

        Get a txt file with the EC1 of each event for each event type 
        Also build ROOT TH1 histograms of EC1 

        They will be used for scaling the bckg in sidebands and
        estimating the bckg.

    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """



    tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")

    #Create background hist directory
    pop_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/")

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    d_est = BDT_fh.open_estimator_file(bolo_name)

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name)
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])


    
    l_heatonly = TEventList("l_heatonly")
    
    l_GammaFid = TEventList("l_GammaFid")
    l_GammaS1  = TEventList("l_GammaS1")
    l_GammaS2  = TEventList("l_GammaS2")
    
    l_BetaS1   = TEventList("l_BetaS1")
    l_BetaS2   = TEventList("l_BetaS2")
    
    l_PbS1     = TEventList("l_PbS1")
    l_PbS2     = TEventList("l_PbS2")


    print "Standard cuts are:  " , standard_cuts


    ###############################
    # Heatonly
    ###############################


    #Heatonly
    tree.Draw(">>l_heatonly",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] )
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close()

    ##################################
    #    G A M M A   E V E N T S
    ##################################

    #Fiducial gammas
    tree.Draw(">>l_GammaFid",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" +  d_est["Q_FID"]  + ">0.7")
    pop_len = l_GammaFid.GetN()
    pop_file_name = bolo_name + "_FidGamma.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_GammaFid.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close()   

    #S1 gammas
    tree.Draw(">>l_GammaS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&&" + d_est["Q_S1"] + ">0.65")
    pop_len = l_GammaS1.GetN()
    pop_file_name = bolo_name + "_S1Gamma.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_GammaS1.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close()   

    #S2 gammas
    tree.Draw(">>l_GammaS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" + d_est["Q_S2"] + ">0.65")
    pop_len = l_GammaS2.GetN()
    pop_file_name = bolo_name + "_S2Gamma.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_GammaS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close()  

    ##################################
    #    B E T A   E V E N T S
    ##################################

    #S1 beta
    tree.Draw(">>l_BetaS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&&" +  d_est["Q_S1"] + "<0.65 && " + d_est["Q_S1"] + ">0.2")
    pop_len = l_BetaS1.GetN()
    pop_file_name = bolo_name + "_S1Beta.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_BetaS1.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close()   

    #S2 beta
    tree.Draw(">>l_BetaS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" +  d_est["Q_S2"] + "<0.65 && " + d_est["Q_S2"] + ">0.2")
    pop_len = l_BetaS2.GetN()
    pop_file_name = bolo_name + "_S2Beta.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_BetaS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close()  



    ##################################
    #    P b    E V E N T S
    ##################################

    # S1 Pb
    tree.Draw(">>l_PbS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0 &&  EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&&" +  d_est["Q_S1"] + "<0.15 &&" +  d_est["Q_S1"] + ">0.04")
    print 
    pop_len = l_PbS1.GetN()
    pop_file_name = bolo_name + "_S1Pb.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_PbS1.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    
    pop_file.close()    

    # S2 Pb
    tree.Draw(">>l_PbS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" +  d_est["Q_S2"] + "<0.15 &&" +  d_est["Q_S2"] + ">0.04")
    pop_len = l_PbS2.GetN()
    pop_file_name = bolo_name + "_S2Pb.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_PbS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")

    pop_file.close() 

    del l_heatonly

    del l_GammaFid
    del l_GammaS1
    del l_GammaS2 

    del l_BetaS1 
    del l_BetaS2 

    del l_PbS1 
    del l_PbS2 


def get_hist(bolo_name, evt_type, bin, min_X, max_X):

    file_path="../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
    file_name_Fid = bolo_name + "_" + evt_type + ".txt"    
    arr = np.loadtxt(file_path + file_name_Fid)

    cc = TCanvas("cc", "cc")
    # gPad.SetLogy()
    h1 = TH1F ("h1", "h1", bin, min_X, max_X)
    for i in range(arr.size):
        h1.Fill(arr[i])


    PyRPl.process_TH1(h1, X_title = "Heat (keVee)", Y_title = "counts/" + str(h1.GetBinWidth(2))[:4] + " keV" )
    h1.Draw("")

    fig_path = script_utils.create_directory("../Analyse_" + bolo_name + "/Figures/Fig_scaling/")
    cc.Print(fig_path + bolo_name + "_spectrum_" + evt_type + ".png")

    hist_path = script_utils.create_directory("../Analyse_" + bolo_name + "/ROOT_files/ROOT_scaling/")
    fhist= TFile(hist_path + bolo_name + "_spectrum_" + evt_type +".root", "recreate")
    h1.Write()
    fhist.Close()



bolo_name = "FID837"


data_dir  = "../Fond_ERA_merged/"
# get_pop_list_for_scaling( bolo_name, data_dir, "data")

get_hist(bolo_name, "S1Beta",  100, 0, 60)
get_hist(bolo_name, "S2Beta",  100, 0, 60)

get_hist(bolo_name, "S1Gamma",  150, 0, 15)
get_hist(bolo_name, "S2Gamma",  150, 0, 15)

get_hist(bolo_name, "S1Pb",  150, 0, 50)
get_hist(bolo_name, "S2Pb",  150, 0, 50)

get_hist(bolo_name, "FidGamma",  200, 0, 15)
get_hist(bolo_name, "heatonly",  200, 0, 15)

