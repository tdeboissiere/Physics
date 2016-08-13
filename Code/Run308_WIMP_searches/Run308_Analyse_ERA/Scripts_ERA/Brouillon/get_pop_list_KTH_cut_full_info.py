#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import BDT_file_handler as BDT_fh

def get_pop_list_for_scaling(bolo_name, d_cut,   data_dir, tree_name = "t_merged"):


    """ Skim true event tree so that it passes the standard cut + analysis cuts
    
    Detail:
        Detailed description

    Args:
        bolo_name     = (str) bolometer name
        d_cut         = (dict) indicates the cuts (inf/sup) on heat and ion 
        data_dir      = (str) data directory

    Returns:
        void

    Raises:
        AssertionError 
    """
    
    file_tree   = TFile(data_dir+bolo_name+"_lowmass_fond.root")
    tree        = file_tree.Get(tree_name)

    #Create background hist directory
    pop_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/")

    #Load the estimator
    d_est             = BDT_fh.open_estimator_file(bolo_name)
    d_std_true_events = BDT_fh.open_true_event_FWHM_file(bolo_name)

    #Best estimator for heat: coefficients
    coeff_EC1, coeff_EC2 = str(d_est["HEAT"][:5]), str(1- float(d_est["HEAT"][:5]))
    coeff_EIA, coeff_EIB = str(d_est["S1"][:5]), str(1-float(d_est["S1"][:5]))
    coeff_EIC, coeff_EID = str(d_est["S2"][:5]), str(1-float(d_est["S2"][:5]))
    
    sigma_IA = str(d_std_true_events["FWIA"])
    sigma_IC = str(d_std_true_events["FWIC"])
    sigma_IB = str(d_std_true_events["FWIB"])
    sigma_ID = str(d_std_true_events["FWID"])

    #Load standard cuts
    TCut_path_name = script_utils.create_directory("../Cut_files/")  
    TCut_file_name ="TCuts.txt" 
    file_TCut      ="" 
    #Add an exception if the file does not exist
    try:
        file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
    except IOError:  
        script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
        sys.exit()
    
    # Load the cut values. 
    list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
    standard_cuts =""
    # Add a boolean flag to check if the bolo has its cuts in the file
    is_bolo_in_file =False
    for line in list_file_TCut_lines:
        if bolo_name == line[0]:
            standard_cuts = line[1]
            is_bolo_in_file = True
    assert(is_bolo_in_file)

    
    l_all      = TEventList("l_all")
    l_heatonly = TEventList("l_heatonly")
    
    l_FidGamma = TEventList("l_FidGamma")
    l_S1Gamma  = TEventList("l_S1Gamma")
    l_S2Gamma  = TEventList("l_S2Gamma")
    
    l_S1Beta   = TEventList("l_S1Beta")
    l_S2Beta   = TEventList("l_S2Beta")
    
    l_S1Pb     = TEventList("l_S1Pb")
    l_S2Pb     = TEventList("l_S2Pb")


    string_EC = coeff_EC1 + "*EC1_ERA+" + coeff_EC2 + "*EC2_ERA"
    string_EI = coeff_EIB + "*EIB+" + coeff_EID + "*EID"

    standard_cuts = standard_cuts + "&&KTH<1&&KTH>0"
    heat_cut = str(d_cut["ECinf"]) + "<" + string_EC + "&&" + str(d_cut["ECsup"]) + ">" + string_EC + "&& abs(EC1_ERA-EC2_ERA)<1"
    ion_cut = str(d_cut["EIinf"]) + "<" + string_EI + "&&" + str(d_cut["EIsup"]) + ">" + string_EI
    veto_cut = "EIA<" + str(d_cut["sigma_vet"]) + "*" + sigma_IA + "&&" + "EIC<" + str(d_cut["sigma_vet"]) + "*" + sigma_IC
    
    # all_cuts = "&&".join([standard_cuts, heat_cut, ion_cut, veto_cut])
    all_cuts = "&&".join([standard_cuts, heat_cut, ion_cut, veto_cut])

    # print tree
    # print all_cuts.split("&&")
    # raw_input()

    ###############################
    # All
    ###############################
    tree.Draw(">>l_all", all_cuts )
    pop_len = l_all.GetN()
    pop_file_name = bolo_name + "_all_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_all.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()

    ###############################
    # Heatonly
    ###############################
    tree.Draw(">>l_heatonly",all_cuts + " && EIA<2.7*" + sigma_IA +" && EIB<2.7*" + sigma_IB +"&& EIC<2.7*" + sigma_IC +"&& EID<2.7*" + sigma_ID)
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()


    ##################################
    #    G A M M A   E V E N T S
    ##################################
    #Fiducial gammas
    tree.Draw(">>l_FidGamma",all_cuts + " && EIA<2.7*" + sigma_IA +" && EIB>2.7*" + sigma_IB +"&& EIC<2.7*" + sigma_IC +"&& EID>2.7*" + sigma_ID + " &&" +  d_est["Q_FID"]  + ">0.7")
    pop_len = l_FidGamma.GetN()
    pop_file_name = bolo_name + "_FidGamma_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_FidGamma.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   
    #S1 gammas
    tree.Draw(">>l_S1Gamma",all_cuts +   "  && EIA>2.7*" + sigma_IA +" && EIB>2.7*" + sigma_IB +"&& EIC<2.7*" + sigma_IC +"&& EID<2.7*" + sigma_ID + " &&" + d_est["Q_S1"] + ">0.65")
    pop_len = l_S1Gamma.GetN()
    pop_file_name = bolo_name + "_S1Gamma_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_S1Gamma.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   
    #S2 gammas
    tree.Draw(">>l_S2Gamma",all_cuts +  " && EIA<2.7*" + sigma_IA +" && EIB<2.7*" + sigma_IB +"&& EIC>2.7*" + sigma_IC +"&& EID>2.7*" + sigma_ID + " &&" + d_est["Q_S2"] + ">0.65")
    pop_len = l_S2Gamma.GetN()
    pop_file_name = bolo_name + "_S2Gamma_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_S2Gamma.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  

    ##################################
    #    B E T A   E V E N T S
    ##################################
    #S1 beta
    tree.Draw(">>l_S1Beta",all_cuts +  "  && EIA>2.7*" + sigma_IA +" && EIB>2.7*" + sigma_IB +"&& EIC<2.7*" + sigma_IC +"&& EID<2.7*" + sigma_ID + " &&" +  d_est["Q_S1"] + "<0.65 && " + d_est["Q_S1"] + ">0.2")
    pop_len = l_S1Beta.GetN()
    pop_file_name = bolo_name + "_S1Beta_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_S1Beta.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   
    #S2 beta
    tree.Draw(">>l_S2Beta",all_cuts +  " && EIA<2.7*" + sigma_IA +" && EIB<2.7*" + sigma_IB +"&& EIC>2.7*" + sigma_IC +"&& EID>2.7*" + sigma_ID + " &&" +  d_est["Q_S2"] + "<0.65 && " + d_est["Q_S2"] + ">0.2")
    pop_len = l_S2Beta.GetN()
    pop_file_name = bolo_name + "_S2Beta_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_S2Beta.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  

    ##################################
    #    P b    E V E N T S
    ##################################
    # S1 Pb
    tree.Draw(">>l_S1Pb",all_cuts +  " &&  EIA>2.7*" + sigma_IA +" && EIB>2.7*" + sigma_IB +"&& EIC<2.7*" + sigma_IC +"&& EID<2.7*" + sigma_ID + " &&" +  d_est["Q_S1"] + "<0.15 &&" +  d_est["Q_S1"] + ">0.04")
    print 
    pop_len = l_S1Pb.GetN()
    pop_file_name = bolo_name + "_S1Pb_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_S1Pb.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()    
    # S2 Pb
    tree.Draw(">>l_S2Pb",all_cuts + " && EIA<2.7*" + sigma_IA +" && EIB<2.7*" + sigma_IB +"&& EIC>2.7*" + sigma_IC +"&& EID>2.7*" + sigma_ID + " &&" +  d_est["Q_S2"] + "<0.15 &&" +  d_est["Q_S2"] + ">0.04")
    pop_len = l_S2Pb.GetN()
    pop_file_name = bolo_name + "_S2Pb_KTH_cut_full_info.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    for k in range(pop_len):
        counter = l_S2Pb.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close() 

    list_list = [l_heatonly, l_FidGamma, l_S1Gamma, l_S2Gamma, l_S1Beta, l_S2Beta, l_S1Pb, l_S2Pb]
    list_num = [l.GetN() for l in list_list]
    list_ev = ["heatonly", "FidGamma", "S1Gamma", "S2Gamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb"]
    for ev, num in zip(list_ev, list_num):
        print ev, num

    print "all known", sum(list_num)
    print "all", l_all.GetN()

    del l_all
    del l_heatonly

    del l_FidGamma
    del l_S1Gamma
    del l_S2Gamma 

    del l_S1Beta 
    del l_S2Beta 

    del l_S1Pb 
    del l_S2Pb 


def plot_for_scaling_check(bolo_name):

    """ Superimpose all + identified bckg to see what is missing
    
    Detail:

    Args:
        bolo_name = (str) bolometer name

    Returns:
        void

    Raises:
        void
    """


    pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

    #Load the estimator
    d_est             = BDT_fh.open_estimator_file(bolo_name)

    #Best estimator for heat: coefficients
    coeff_EC1, coeff_EC2 = float(d_est["HEAT"][:5]), 1 - float(d_est["HEAT"][:5])
    coeff_EIB, coeff_EID = float(d_est["FID"][:5]), 1-float(d_est["FID"][:5])

    #Open event files
    data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}

    arr_heatonly = np.loadtxt(pop_path + bolo_name + "_heatonly_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_all = np.loadtxt(pop_path + bolo_name + "_all_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_FidGamma = np.loadtxt(pop_path + bolo_name + "_FidGamma_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_S1Gamma = np.loadtxt(pop_path + bolo_name + "_S1Gamma_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_S2Gamma = np.loadtxt(pop_path + bolo_name + "_S2Gamma_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_S1Beta = np.loadtxt(pop_path + bolo_name + "_S1Beta_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_S2Beta = np.loadtxt(pop_path + bolo_name + "_S2Beta_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_S1Pb   = np.loadtxt(pop_path + bolo_name + "_S1Pb_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)
    arr_S2Pb   = np.loadtxt(pop_path + bolo_name + "_S2Pb_KTH_cut_full_info.txt", delimiter=",",  dtype=data_types)

    arr_EI_heatonly, arr_EI_all = coeff_EIB*arr_heatonly["EIB"] + coeff_EID*arr_heatonly["EID"], coeff_EIB*arr_all["EIB"] + coeff_EID*arr_all["EID"]
    arr_EC_heatonly, arr_EC_all = coeff_EC1*arr_heatonly["EC1"] + coeff_EC2*arr_heatonly["EC2"], coeff_EC1*arr_all["EC1"] + coeff_EC2*arr_all["EC2"]
    arr_EI_FidGamma, arr_EC_FidGamma = coeff_EIB*arr_FidGamma["EIB"] + coeff_EID*arr_FidGamma["EID"], coeff_EC1*arr_FidGamma["EC1"] + coeff_EC2*arr_FidGamma["EC2"]
    arr_EI_S1Gamma, arr_EI_S2Gamma = coeff_EIB*arr_S1Gamma["EIB"] + coeff_EID*arr_S1Gamma["EID"], coeff_EIB*arr_S2Gamma["EIB"] + coeff_EID*arr_S2Gamma["EID"]
    arr_EC_S1Gamma, arr_EC_S2Gamma = coeff_EC1*arr_S1Gamma["EC1"] + coeff_EC2*arr_S1Gamma["EC2"], coeff_EC1*arr_S2Gamma["EC1"] + coeff_EC2*arr_S2Gamma["EC2"]
    arr_EI_S1Beta, arr_EI_S2Beta = coeff_EIB*arr_S1Beta["EIB"] + coeff_EID*arr_S1Beta["EID"], coeff_EIB*arr_S2Beta["EIB"] + coeff_EID*arr_S2Beta["EID"]
    arr_EC_S1Beta, arr_EC_S2Beta = coeff_EC1*arr_S1Beta["EC1"] + coeff_EC2*arr_S1Beta["EC2"], coeff_EC1*arr_S2Beta["EC1"] + coeff_EC2*arr_S2Beta["EC2"]
    arr_EI_S1Pb, arr_EI_S2Pb     = coeff_EIB*arr_S1Pb["EIB"] + coeff_EID*arr_S1Pb["EID"], coeff_EIB*arr_S2Pb["EIB"] + coeff_EID*arr_S2Pb["EID"]
    arr_EC_S1Pb, arr_EC_S2Pb     = coeff_EC1*arr_S1Pb["EC1"] + coeff_EC2*arr_S1Pb["EC2"], coeff_EC1*arr_S2Pb["EC1"] + coeff_EC2*arr_S2Pb["EC2"]

    lS1Beta, lS2Beta, lS1Pb, lS2Pb = np.where(arr_EC_S1Beta<15), np.where(arr_EC_S2Beta<15), np.where(arr_EC_S1Pb<15), np.where(arr_EC_S2Pb<15)
    lS1Gamma, lS2Gamma, lFidGamma = np.where(arr_EC_S1Gamma<15), np.where(arr_EC_S2Gamma<15), np.where(arr_EC_FidGamma<15)
    lheatonly, lall = np.where(arr_EC_heatonly<15), np.where(arr_EC_all<15)

    arr_EI_heatonly, arr_EC_heatonly = arr_EI_heatonly[lheatonly], arr_EC_heatonly[lheatonly]
    arr_EI_all, arr_EC_all = arr_EI_all[lall], arr_EC_all[lall]
    arr_EI_FidGamma, arr_EC_FidGamma = arr_EI_FidGamma[lFidGamma], arr_EC_FidGamma[lFidGamma]
    arr_EI_S1Gamma, arr_EC_S1Gamma = arr_EI_S1Gamma[lS1Gamma], arr_EC_S1Gamma[lS1Gamma]
    arr_EI_S2Gamma, arr_EC_S2Gamma = arr_EI_S2Gamma[lS2Gamma], arr_EC_S2Gamma[lS2Gamma]
    arr_EI_S1Beta, arr_EC_S1Beta = arr_EI_S1Beta[lS1Beta], arr_EC_S1Beta[lS1Beta]
    arr_EI_S2Beta, arr_EC_S2Beta = arr_EI_S2Beta[lS2Beta], arr_EC_S2Beta[lS2Beta]
    arr_EI_S1Pb, arr_EC_S1Pb     = arr_EI_S1Pb[lS1Pb], arr_EC_S1Pb[lS1Pb]
    arr_EI_S2Pb, arr_EC_S2Pb     = arr_EI_S2Pb[lS2Pb], arr_EC_S2Pb[lS2Pb]

    arr_EI_all, arr_EC_all = np.array(arr_EI_all).astype(float), np.array(arr_EC_all).astype(float)
    arr_EI_heatonly, arr_EC_heatonly = np.array(arr_EI_heatonly).astype(float), np.array(arr_EC_heatonly).astype(float)
    arr_EI_FidGamma, arr_EC_FidGamma = np.array(arr_EI_FidGamma).astype(float), np.array(arr_EC_FidGamma).astype(float)
    arr_EI_S1Gamma, arr_EC_S1Gamma = np.array(arr_EI_S1Gamma).astype(float), np.array(arr_EC_S1Gamma).astype(float)
    arr_EI_S2Gamma, arr_EC_S2Gamma = np.array(arr_EI_S2Gamma).astype(float), np.array(arr_EC_S2Gamma).astype(float)    
    arr_EI_S1Beta, arr_EC_S1Beta = np.array(arr_EI_S1Beta).astype(float), np.array(arr_EC_S1Beta).astype(float)
    arr_EI_S2Beta, arr_EC_S2Beta = np.array(arr_EI_S2Beta).astype(float), np.array(arr_EC_S2Beta).astype(float)
    arr_EI_S1Pb, arr_EC_S1Pb     = np.array(arr_EI_S1Pb).astype(float), np.array(arr_EC_S1Pb).astype(float)
    arr_EI_S2Pb, arr_EC_S2Pb     = np.array(arr_EI_S2Pb).astype(float), np.array(arr_EC_S2Pb).astype(float)


    gr_heatonly   = TGraph(len(arr_EI_heatonly), arr_EC_heatonly, arr_EI_heatonly)
    gr_FidGamma, gr_all   = TGraph(len(arr_EI_FidGamma), arr_EC_FidGamma, arr_EI_FidGamma),  TGraph(len(arr_EI_all), arr_EC_all, arr_EI_all)
    gr_S1Gamma, gr_S2Gamma   = TGraph(len(arr_EI_S1Gamma), arr_EC_S1Gamma, arr_EI_S1Gamma),  TGraph(len(arr_EI_S2Gamma), arr_EC_S2Gamma, arr_EI_S2Gamma)
    gr_S1Beta, gr_S2Beta   = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta),  TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)
    gr_S1Pb, gr_S2Pb       = TGraph(len(arr_EI_S1Pb), arr_EC_S1Pb, arr_EI_S1Pb),  TGraph(len(arr_EI_S2Pb), arr_EC_S2Pb, arr_EI_S2Pb)

    PyRPl.process_TGraph(gr_all, X_title = "Heat", Y_title = "Ion", color=kRed, marker_style = 20, marker_size = 0.1)
    PyRPl.process_TGraph(gr_FidGamma, X_title = "Heat", Y_title = "Ion", color=kBlack), PyRPl.process_TGraph(gr_heatonly, X_title = "Heat", Y_title = "Ion", color=kBlack)
    PyRPl.process_TGraph(gr_S1Gamma, X_title = "Heat", Y_title = "Ion", color=kBlack), PyRPl.process_TGraph(gr_S2Gamma, X_title = "Heat", Y_title = "Ion", color=kBlack)
    PyRPl.process_TGraph(gr_S1Beta, X_title = "Heat", Y_title = "Ion", color=kBlack), PyRPl.process_TGraph(gr_S2Beta, X_title = "Heat", Y_title = "Ion", color=kBlack)
    PyRPl.process_TGraph(gr_S1Pb, X_title = "Heat", Y_title = "Ion", color=kRed), PyRPl.process_TGraph(gr_S2Pb, X_title = "Heat", Y_title = "Ion", color=kBlack)

    list_gr = [gr_all, gr_FidGamma, gr_S1Gamma, gr_S2Gamma, gr_S1Beta, gr_S2Beta, gr_S1Pb, gr_S2Pb, gr_heatonly]
    list_pts = [gr.GetN() for gr in list_gr[1:]]
    print gr_all.GetN(), sum(list_pts)
    h = TH2F("h", "h", 100, -5, 15, 100, -5, 15)
    PyRPl.process_TH2(h, X_title = "Heat", Y_title = "Ion")
    h.Draw()
    for gr in list_gr:
        gr.Draw("*same")

    raw_input()

    # arr_Q_S1Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
    # arr_Q_S2Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
    # arr_Q_S1Pb   = arr_EI_S1Pb/((1+8./3)*arr_EC_S1Pb - arr_EI_S1Pb*5.5/3)
    # arr_Q_S2Pb   = arr_EI_S2Pb/((1+8./3)*arr_EC_S2Pb - arr_EI_S2Pb*5.5/3)
    
    # gr_QS1Beta, gr_QS2Beta = TGraph(len(arr_Q_S1Beta), arr_EC_S1Beta, arr_Q_S1Beta),  TGraph(len(arr_Q_S2Beta), arr_EC_S2Beta, arr_Q_S2Beta)
    # gr_QS1Pb, gr_QS2Pb     = TGraph(len(arr_Q_S1Pb), arr_EC_S1Pb, arr_Q_S1Pb),  TGraph(len(arr_Q_S2Pb), arr_EC_S2Pb, arr_Q_S2Pb)


    # PyRPl.process_TGraph(gr_QS1Beta, X_title = "Heat", Y_title = "Q", color=kOrange-3), PyRPl.process_TGraph(gr_QS2Beta, X_title = "Heat", Y_title = "Q", color=kBlue)
    # PyRPl.process_TGraph(gr_QS1Pb, X_title = "Heat", Y_title = "Q", color=kRed), PyRPl.process_TGraph(gr_QS2Pb, X_title = "Heat", Y_title = "Q", color=kGreen+2)


bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
#1 sigma on heat only band=  EFid>0.2
d_cut       = {"ECinf": 0, "ECsup": 60, "EIinf": 0, "EIsup": 60, "sigma_vet": 1000}
data_dir = "../Fond_ERA_merged/"
get_pop_list_for_scaling( bolo_name, d_cut,  data_dir, "t_merged")
# plot_for_scaling_check(bolo_name)


