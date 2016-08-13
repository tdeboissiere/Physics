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


    l_true_lowE                = TEventList("l_true_lowE")
    l_true_highE          = TEventList("l_true_highE")
    l_heatonly_above4    = TEventList("l_heatonly_above4")
    l_heatonly_below4    = TEventList("l_heatonly_below4")
    l_heatonly_ioncut1keV = TEventList("l_heatonly_ioncut1keV")
    
    print "Standard cuts are:  " , standard_cuts
    #Then create various .txt files for the population selection

    #All events up to rather low energies
    tree.Draw(">>l_true_lowE",standard_cuts +  "&&" + d_estimator["HEAT"] + "<30 &&" + d_estimator["HEAT"] + ">4 && 0.5*(EIA+EIB+EIC+EID)>4" + "&& abs(EC1-EC2)<2")
    pop_len = l_true_lowE.GetN()
    pop_file_name = bolo_name + "_true_lowE.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_true_lowE.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()

    #All events up to high energies
    tree.Draw(">>l_true_highE",standard_cuts +  "&&" + d_estimator["HEAT"] + "<200 &&" + d_estimator["HEAT"] + ">30 && 0.5*(EIA+EIB+EIC+EID)>4" + "&& abs(EC1-EC2)<2")
    pop_len = l_true_highE.GetN()
    pop_file_name = bolo_name + "_true_highE.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_true_highE.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()

    #All heatonly above 4
    tree.Draw(">>l_heatonly_above4",standard_cuts +  "&&" + d_estimator["HEAT"] + ">4 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2")
    pop_len = l_heatonly_above4.GetN()
    pop_file_name = bolo_name + "_heatonly_above4.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly_above4.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    #All heatonly below 4
    tree.Draw(">>l_heatonly_below4",standard_cuts +  "&&" + d_estimator["HEAT"] + "<4 &&"  +d_estimator["HEAT"] + "<4 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2")
    pop_len = l_heatonly_below4.GetN()
    pop_file_name = bolo_name + "_heatonly_below4.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly_below4.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) +"," + str(tree.CHIC1) + "," + str(tree.CHIC2)+ "\n")
    pop_file.close() 

    #All heatonly above 4
    tree.Draw(">>l_heatonly_ioncut1keV","RUN<=251100&&FWIA<1.4&&FWIB<1.1&&FWIC<1.1&&FWID<1&&FWC1<0.8&&FWC2<0.8&&CHIA<2.5&&CHIB<2.5&&CHIC<2.2&&CHID<2&&CHIC1<3&&CHIC2<3&&SDEL>1&&FWIA>0&&FWIB>0&&FWIC>0&&FWID>0&&FWC1>0&&FWC2>0 &&" + d_estimator["HEAT"] + "<15&&" + d_estimator["HEAT"] + ">1 && 0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2")
    pop_len = l_heatonly_ioncut1keV.GetN()
    pop_file_name = bolo_name + "_heatonly_ioncut1keV.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly_ioncut1keV.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()   

    del l_true_highE
    del l_true_lowE
    del l_heatonly_ioncut1keV
    del l_heatonly_below4
    del l_heatonly_above4

bolo_list=["FID837"]
duration_dict={}
#Load duration
for bolo_name in bolo_list:
    duration_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
    duration_file_name =bolo_name + "_polar_start_and_end_time.txt" 
    duration         = 0
    assert(os.path.isfile(duration_path_name + duration_file_name) )
    with open(duration_path_name + duration_file_name, 'r') as duration_file: 
        list_duration_lines                                   = duration_file.readlines()[0].rstrip().split(",") 
        duration                                              = float(list_duration_lines[3])/(24*3600)
        script_utils.print_utility(script_utils.COL("Duration (in days) = "  + str(duration) , "blue"))
        duration_dict[bolo_name] = duration


data_dir  = "../Fond_ERA_merged/"
for bolo_name in bolo_list:
    get_population_list( bolo_name, data_dir, "data")

