#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh

def get_pop_list_for_scaling_full_info(bolo_name, data_dir, tree_name):

    """Get the population files for scaling with full info
    
    Detail:
        Select each event type (Gamma, Beta, Pb, heatonly)
        also include Surface 1 Surface 2 and Fiducial distinction
        when applicable

        Get a txt file with the full info of each event for each event type
        full info = EC1 EC2, EIA, EIB, EIC, EID 

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


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

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

    
    l_heatonly = TEventList("l_all")
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
    # All
    ###############################


    #All
    tree.Draw(">>l_all",standard_cuts +  "&&" + d_est["HEAT"] + ">0 &&" + d_est["HEAT"] + "<60" )
    pop_len = l_all.GetN()
    pop_file_name = bolo_name + "_all_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_all.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()


    ###############################
    # Heatonly
    ###############################


    #Heatonly
    tree.Draw(">>l_heatonly",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] )
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()

    ##################################
    #    G A M M A   E V E N T S
    ##################################

    #Fiducial gammas
    tree.Draw(">>l_GammaFid",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_FID"]  + ">0.7")
    pop_len = l_GammaFid.GetN()
    pop_file_name = bolo_name + "_FidGamma_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_GammaFid.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    #S1 gammas
    tree.Draw(">>l_GammaS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" + d_est["Q_S1"] + ">0.65")
    pop_len = l_GammaS1.GetN()
    pop_file_name = bolo_name + "_S1Gamma_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_GammaS1.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    #S2 gammas
    tree.Draw(">>l_GammaS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" + d_est["Q_S2"] + ">0.65")
    pop_len = l_GammaS2.GetN()
    pop_file_name = bolo_name + "_S2Gamma_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_GammaS2.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  

    ##################################
    #    B E T A   E V E N T S
    ##################################

    #S1 beta
    tree.Draw(">>l_BetaS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" +  d_est["Q_S1"] + "<0.65 && " + d_est["Q_S1"] + ">0.2")
    pop_len = l_BetaS1.GetN()
    pop_file_name = bolo_name + "_S1Beta_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_BetaS1.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    #S2 beta
    tree.Draw(">>l_BetaS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_S2"] + "<0.65 && " + d_est["Q_S2"] + ">0.2")
    pop_len = l_BetaS2.GetN()
    pop_file_name = bolo_name + "_S2Beta_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_BetaS2.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  



    ##################################
    #    P b    E V E N T S
    ##################################

    # S1 Pb
    tree.Draw(">>l_PbS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0 &&  EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" +  d_est["Q_S1"] + "<0.15 &&" +  d_est["Q_S1"] + ">0.04")
    print 
    pop_len = l_PbS1.GetN()
    pop_file_name = bolo_name + "_S1Pb_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_PbS1.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()    

    # S2 Pb
    tree.Draw(">>l_PbS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_S2"] + "<0.15 &&" +  d_est["Q_S2"] + ">0.04")
    pop_len = l_PbS2.GetN()
    pop_file_name = bolo_name + "_S2Pb_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_PbS2.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close() 

    del l_heatonly

    del l_GammaFid
    del l_GammaS1
    del l_GammaS2 

    del l_BetaS1 
    del l_BetaS2 

    del l_PbS1 
    del l_PbS2 




bolo_name = "FID837"
data_dir  = "../Fond_ERA_merged/"
get_pop_list_for_scaling_full_info( bolo_name, data_dir, "data")



