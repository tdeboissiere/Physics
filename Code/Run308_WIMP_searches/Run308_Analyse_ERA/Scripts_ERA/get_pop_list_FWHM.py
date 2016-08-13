#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import BDT_file_handler as BDT_fh
import Analysis_utilities as Ana_ut


def get_pop_list_FWHM(bolo_name, d_cut,  data_dir, tree_name = "t_merged"):

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


    file_tree   = TFile(data_dir+bolo_name+"_lowmass_fond.root")
    tree        = file_tree.Get(tree_name)

    #Create background hist directory
    pop_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/")

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    d_est = BDT_fh.open_estimator_file(bolo_name)

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name)
    for key in ["FWC1_ERA", "FWC2_ERA", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])

    l_all      = TEventList("l_all")
    l_heatonly = TEventList("l_heatonly")
    l_FidGamma = TEventList("l_FidGamma")

    standard_cuts = standard_cuts

    # heat_cut = str(d_cut["ECinf"]) + "<" + d_est["HEAT"] + "&&" + str(d_cut["ECsup"]) + ">" + d_est["HEAT"] 
    # ion_cut = str(d_cut["EIinf"]) + "<" + d_est["FID"] + "&&" + str(d_cut["EIsup"]) + ">" + d_est["FID"]
    # veto_cut = "EIA<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIA"] + "&&" + "EIC<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIC"]
    
    # all_cuts = "&&".join([standard_cuts, heat_cut, ion_cut, veto_cut])


    ##################################
    #    A L L   E V E N T S
    ##################################

    #Fiducial gammas
    tree.Draw(">>l_all",standard_cuts )
    pop_len = l_all.GetN()
    pop_file_name = bolo_name + "_all_FWHM.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_all.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.FWC1_ERA) + "," + str(tree.FWC2_ERA) + "," + str(tree.FWIA) + "," + str(tree.FWIB) + "," + str(tree.FWIC) + "," + str(tree.FWID) + "\n")
    pop_file.close()   

    ##################################
    #    H E A T O N L Y   E V E N T S
    ##################################

    #Fiducial gammas
    tree.Draw(">>l_heatonly",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"])
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly_FWHM.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.FWC1_ERA) + "," + str(tree.FWC2_ERA) + "," + str(tree.FWIA) + "," + str(tree.FWIB) + "," + str(tree.FWIC) + "," + str(tree.FWID) + "\n")
    pop_file.close()   

    ##################################
    #    G A M M A   E V E N T S
    ##################################

    #Fiducial gammas
    tree.Draw(">>l_FidGamma",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" +  d_est["Q_FID"]  + ">0.7")
    pop_len = l_FidGamma.GetN()
    pop_file_name = bolo_name + "_FidGamma_FWHM.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_FidGamma.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.FWC1_ERA) + "," + str(tree.FWC2_ERA) + "," + str(tree.FWIA) + "," + str(tree.FWIB) + "," + str(tree.FWIC) + "," + str(tree.FWID) + "\n")
    pop_file.close()   


    del l_all
    del l_heatonly
    del l_FidGamma


bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
#1 sigma on heat only band=  EFid>0.2
d_cut       = {"ECinf": -2, "ECsup": 15, "EIinf": -2, "EIsup": 15, "sigma_vet": 5}
data_dir = "../Fond_ERA_merged/"
get_pop_list_FWHM( bolo_name, d_cut, data_dir, "t_merged")



