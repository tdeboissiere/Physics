#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np

def get_population_list(bolo_name, data_dir, tree_name):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Create background hist directory
    pop_path_name = script_utils.create_directory("../Analyse_" + bolo_name + "/Populations/")


    #Load standard cuts
    TCut_path_name = script_utils.create_directory('../Cut_files/')  
    TCut_file_name ="TCuts.txt" 
    file_TCut      ="" 
    #Add an exception if the file does not exist
    try:
        file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
    except IOError:
        script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
        sys.exit()


    #Load estimators
    estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
    estimator_file_name =bolo_name + "_estimators.txt" 
    d_estimator         ={}
    assert(os.path.isfile(estimator_path_name + estimator_file_name) )
    with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
        list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
        for line in list_estimator_lines:
            d_estimator[line[0]] = line[1]

    #Load FWHM
    FWHM_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
    FWHM_file_name =bolo_name + "_FWHM.txt" 
    d_FWHM         ={}
    assert(os.path.isfile(FWHM_path_name + FWHM_file_name) )
    with open(FWHM_path_name + FWHM_file_name, 'r') as FWHM_file: 
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


    l_heatonly = TEventList("l_heatonly")
    l_BetaS1   = TEventList("l_BetaS1")
    l_BetaS2   = TEventList("l_BetaS2")
    
    l_PbS1     = TEventList("l_PbS1")
    l_PbS2     = TEventList("l_PbS2")

    print "Standard cuts are:  " , standard_cuts

    #Then create various .txt files for the population selection
 
    ##################################
    #
    #
    #    H E A T  E V E N T S
    #
    #
    ##################################

    #S1 beta
    tree.Draw(">>l_heatonly",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0  && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 ")
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly_for_comparison.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  

    del l_heatonly

    ##################################
    #
    #
    #    B E T A   E V E N T S
    #
    #
    ##################################

    #S1 beta
    tree.Draw(">>l_BetaS1",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0  && EIA>" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S1"] + "<0.75 && " + d_estimator["Q_S1"] + ">0.2")
    pop_len = l_BetaS1.GetN()
    pop_file_name = bolo_name + "_BetaS1_for_comparison.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_BetaS1.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    #S2 beta
    tree.Draw(">>l_BetaS2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC>" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S2"] + "<0.75 && " + d_estimator["Q_S2"] + ">0.2")
    pop_len = l_BetaS2.GetN()
    pop_file_name = bolo_name + "_BetaS2_for_comparison.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_BetaS2.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  


    del l_BetaS1 
    del l_BetaS2 


    ##################################
    #
    #
    #    P b  E V E N T S
    #
    #
    ##################################

    #S1 Pb
    tree.Draw(">>l_PbS1",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0  && EIA>" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S1"] + "<0.2 &&" +  d_estimator["Q_S1"] + ">0.04")
    pop_len = l_PbS1.GetN()
    pop_file_name = bolo_name + "_PbS1_for_comparison.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_PbS1.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    #S2 Pb
    tree.Draw(">>l_PbS2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC>" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S2"] + "<0.2 &&" +  d_estimator["Q_S2"] + ">0.04")
    pop_len = l_PbS2.GetN()
    pop_file_name = bolo_name + "_PbS2_for_comparison.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_PbS2.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()  


    del l_PbS1 
    del l_PbS2 



bolo_list=["FID837"]
data_dir  = "../Fond_ERA_merged/"

for bolo_name in bolo_list:
    get_population_list( bolo_name, data_dir, "data")

