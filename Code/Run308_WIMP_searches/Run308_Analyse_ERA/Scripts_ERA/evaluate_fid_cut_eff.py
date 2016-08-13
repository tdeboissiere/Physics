#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Poisson_file_handler as Poisson_fh

def evaluate_fid_cut_eff(bolo_name, data_dir, tree_name):


    file_tree   = TFile(data_dir+bolo_name+"_lowmass_fond.root")
    tree        = file_tree.Get(tree_name)

    #Create background hist directory
    pop_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/")


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


    #Load estimators
    estimator_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Text_files/")  
    estimator_file_name =bolo_name + "_estimators.txt" 
    d_estimator         ={}
    assert(os.path.isfile(estimator_path_name + estimator_file_name) )
    with open(estimator_path_name + estimator_file_name, "r") as estimator_file: 
        list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
        for line in list_estimator_lines:
            d_estimator[line[0]] = line[1]

    #Load FWHM
    FWHM_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Text_files/")  
    FWHM_file_name =bolo_name + "_FWHM.txt" 
    d_FWHM         ={}
    assert(os.path.isfile(FWHM_path_name + FWHM_file_name) )
    with open(FWHM_path_name + FWHM_file_name, "r") as FWHM_file: 
        list_FWHM_lines= [elem.rstrip().split(",") for elem in FWHM_file.readlines()]
        for index, fwhm in enumerate(list_FWHM_lines[0]):
            d_FWHM[fwhm] = str(2.7/2.3548*float(list_FWHM_lines[1][index]))
    print d_FWHM

    # Load the cut values. 
    list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
    standard_cuts        =""
    # Add a boolean flag to check if the bolo has its cuts in the file
    is_bolo_in_file      =False
    for line in list_file_TCut_lines:
        if bolo_name == line[0]:
            standard_cuts = line[1]
            is_bolo_in_file = True
    assert(is_bolo_in_file)

    standard_cuts = standard_cuts + "&&KTH<1"

    #Load the estimator
    d_est             = Poisson_fh.open_estimator_file(bolo_name)

    #Best estimator for heat: coefficients
    coeff_EC1, coeff_EC2 = float(d_est["HEAT"][:5]), 1 - float(d_est["HEAT"][:5])
    coeff_EIB, coeff_EID = float(d_est["FID"][:5]), 1-float(d_est["FID"][:5])
    
    # l_GammaFid = TEventList("l_GammaFid")
    # pop_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Populations/Pop_for_fid_cut_eff/")
    # print "Standard cuts are:  " , standard_cuts


    # # #Fiducial gammas
    # tree.Draw(">>l_GammaFid",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<5 && EIC<5"  + "&& abs(EC1_ERA-EC2_ERA)<3 && EIB>8 && EIB<15 && EID>8 && EID<15&&" +  d_estimator["Q_FID"]  + ">0.7")
    # pop_file_name = bolo_name + "_FidGamma_for_fid_cut_eff.txt"
    # pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, "w")
    # pop_len = l_GammaFid.GetN()
    # print pop_len, tree.GetEntries()
    # raw_input()
    # for k in range(pop_len):
    #     counter = l_GammaFid.GetEntry(k)
    #     tree.GetEntry(counter)
    #     pop_file.write(str(tree.EC1_ERA) + "," + str(tree.EC2_ERA) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    # pop_file.close()   


    # #Open created text file as numpy array for convenience
    # data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}
    # arr_FidGamma = np.loadtxt(pop_path_name+ pop_file_name, delimiter=",",  dtype=data_types)

    # arr_EI_FidGamma = coeff_EIB*arr_FidGamma["EIB"] + coeff_EID*arr_FidGamma["EID"]
    # arr_EC_FidGamma = coeff_EC1*arr_FidGamma["EC1"] + coeff_EC2*arr_FidGamma["EC2"]

    # arr_Q_FidGamma = arr_EI_FidGamma/((1+8./3)*arr_EC_FidGamma - arr_EI_FidGamma*8./3)
    # gr_FidGamma= TGraph(len(arr_EI_FidGamma), np.array(arr_EC_FidGamma).astype(float), np.array(arr_EI_FidGamma).astype(float))
    # gr_QFidGamma = TGraph(len(arr_Q_FidGamma), np.array(arr_EC_FidGamma).astype(float), np.array(arr_Q_FidGamma).astype(float))

    # PyRPl.process_TGraph(gr_FidGamma, X_title = "Heat", Y_title = "Ion", color=kOrange-3)
    # PyRPl.process_TGraph(gr_QFidGamma, X_title = "Heat", Y_title = "Q", color=kOrange-3)

    # cc = TCanvas("cc", "cc")
    # cc.Divide(2,1)
    # cc.cd(1)
    # gr_FidGamma.Draw("A*")    
    # cc.cd(2)
    # gr_QFidGamma.Draw("A*")

    # raw_input()

    d_std_true_events = Poisson_fh.open_true_event_FWHM_file(bolo_name)
    sIA = str(d_std_true_events["FWIA"])
    sIB = str(d_std_true_events["FWIB"])
    sIC = str(d_std_true_events["FWIC"])
    sID = str(d_std_true_events["FWID"])
    sHeat_dif = str(TMath.Sqrt(d_std_true_events["FWC1_ERA"]**2 + d_std_true_events["FWC2_ERA"]**2))
    sFid_dif =str(TMath.Sqrt(d_std_true_events["FWIB"]**2 + d_std_true_events["FWIA"]**2))

    print sIA, sIB, sIC, sID, sHeat_dif, sFid_dif

    full_cuts =standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<5 && EIC<5"  + "&& abs(EC1_ERA-EC2_ERA)<3 && EIB>8 && EIB<15 && EID>8 && EID<15&&" +  d_estimator["Q_FID"]  + ">0.7"
    full_cuts_entries = float(tree.GetEntries(full_cuts))
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<5*" + sFid_dif)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<4*" + sFid_dif)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<3*" + sFid_dif)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<2*" + sFid_dif)/full_cuts_entries

    print

    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<4*" + sFid_dif + " && EIA<5*" + sIA + " && EIC<5*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<4*" + sFid_dif + " && EIA<4*" + sIA + " && EIC<4*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<4*" + sFid_dif + " && EIA<3*" + sIA + " && EIC<3*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<4*" + sFid_dif + " && EIA<2*" + sIA + " && EIC<2*" + sIC)/full_cuts_entries

    print

    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<3*" + sFid_dif + " && EIA<5*" + sIA + " && EIC<5*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<3*" + sFid_dif + " && EIA<4*" + sIA + " && EIC<4*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<3*" + sFid_dif + " && EIA<3*" + sIA + " && EIC<3*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<3*" + sFid_dif + " && EIA<2*" + sIA + " && EIC<2*" + sIC)/full_cuts_entries

    print

    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<2*" + sFid_dif + " && EIA<5*" + sIA + " && EIC<5*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<2*" + sFid_dif + " && EIA<4*" + sIA + " && EIC<4*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<2*" + sFid_dif + " && EIA<3*" + sIA + " && EIC<3*" + sIC)/full_cuts_entries
    print tree.GetEntries(full_cuts + "&&abs(EIB-EID)<2*" + sFid_dif + " && EIA<2*" + sIA + " && EIC<2*" + sIC)/full_cuts_entries

bolo_name = "FID837"
data_dir = "../Fond_ERA_merged/"
tree_name = "t_merged"
evaluate_fid_cut_eff(bolo_name, data_dir, tree_name)