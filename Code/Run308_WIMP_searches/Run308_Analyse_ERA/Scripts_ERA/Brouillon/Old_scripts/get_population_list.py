#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *

def get_population_list(bolo_name, data_dir, tree_name):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Create population directory
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


    l_all = TEventList("l_all")
    l_all_highE= TEventList("l_all_highE")
    l_all_ion_cut = TEventList("l_all_ion_cut")

    l_heatonly        = TEventList("l_heatonly")
    l_heatonly_above2 = TEventList("l_heatonly_above2")


    l_GammaFid = TEventList("l_GammaFid")

    l_BetaS1 = TEventList("l_BetaS1")
    l_PbS1 = TEventList("l_PbS1")
    l_GammaS1 = TEventList("l_GammaS1")

    l_BetaS2 = TEventList("l_BetaS2")
    l_PbS2 = TEventList("l_PbS2")
    l_GammaS2 = TEventList("l_GammaS2")

    # #Get Fiducial gammas
    # tree.Draw(">>l_heatonly_above2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">2 && EIA<1 && EIB<1 && EIC<1 && EID<1 && abs(EC1-EC2)<2")
    # pop_len = l_heatonly_above2.GetN()
    # pop_file_name = bolo_name + "_4Vheatonly_above2.txt"
    # pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    # for k in range(pop_len):
    #     counter = l_heatonly_above2.GetEntry(k)
    #     tree.GetEntry(counter)
    #     pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    #     # est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
    #     # pop_file.write(str(est_Heat) + "\n")
    # pop_file.close()   

    # #Then create various .txt files for the population selection
    # tree.Draw(">>l_4Vheatonly_above20",standard_cuts + " && 0.5*(EIA+ EIB+ EIC+EID)<1  && abs(EC1-EC2)<2")
    # print standard_cuts + " && 0.5*(EIA+ EIB+ EIC+EID)<1  && abs(EC1-EC2)<2"
    # pop_len = l_4Vheatonly_above20.GetN()
    # pop_file_name = bolo_name + "_4Vheatonly_above20.txt"
    # pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    # for k in range(pop_len):
    #     counter = l_4Vheatonly_above20.GetEntry(k)
    #     tree.GetEntry(counter)
    #     pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    # pop_file.close()

    # #Then create various .txt files for the population selection
    # tree.Draw(">>l_4Vheatonly_above2",standard_cuts + "&& 0.5*(EIA+ EIB+ EIC+EID)<1  && abs(EC1-EC2)<2  && 0.5*(EC1+EC2)>2")
    # pop_len = l_4Vheatonly_above2.GetN()
    # pop_file_name = bolo_name + "_4Vheatonly_above2.txt"
    # pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    # for k in range(pop_len):
    #     counter = l_4Vheatonly_above2.GetEntry(k)
    #     tree.GetEntry(counter)
    #     pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    # pop_file.close()

    # #Then create various .txt files for the population selection
    # tree.Draw(">>l_4Vheatonly_above5",standard_cuts + "&& 0.5*(EIA+ EIB+ EIC+EID)<1  && abs(EC1-EC2)<2 && 0.5*(EC1+EC2)>5")
    # pop_len = l_4Vheatonly_above5.GetN()
    # pop_file_name = bolo_name + "_4Vheatonly_above5.txt"
    # pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    # for k in range(pop_len):
    #     counter = l_4Vheatonly_above5.GetEntry(k)
    #     tree.GetEntry(counter)
    #     pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    # pop_file.close()

    #Then create various .txt files for the population selection
    tree.Draw(">>l_all",standard_cuts +  "&&" + d_estimator["HEAT"] + "<40 &&" + d_estimator["HEAT"] + ">0 && abs(EC1-EC2)<2")
    pop_len = l_all.GetN()
    pop_file_name = bolo_name + "_all.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_all.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()

    #Then create various .txt files for the population selection
    tree.Draw(">>l_all_highE",standard_cuts +  "&&" + d_estimator["HEAT"] + "<400 &&" + d_estimator["HEAT"] + ">0 && abs(EC1-EC2)<2")
    pop_len = l_all_highE.GetN()
    pop_file_name = bolo_name + "_all_highE.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_all_highE.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()



    tree.Draw(">>l_heatonly",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2")
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()


    ##################################
    #
    #
    #    G A M M A   E V E N T S
    #
    #
    ##################################

    #Get Fiducial gammas
    tree.Draw(">>l_GammaFid",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_FID"]  + ">0.7")
    pop_len = l_GammaFid.GetN()
    pop_file_name = bolo_name + "_GammaFid.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_GammaFid.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close()   

    tree.Draw(">>l_GammaS1",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0  && EIA>" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" + d_estimator["Q_S1"] + ">0.75")
    # print standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA>1 && EIB>1 && EIC<1 && EID<1 && abs(EC1-EC2)<2 && " + d_estimator["Q_S1"] + ">0.75"
    pop_len = l_GammaS1.GetN()
    pop_file_name = bolo_name + "_GammaS1.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_GammaS1.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close()   

    tree.Draw(">>l_GammaS2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC>" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" + d_estimator["Q_S2"] + ">0.75")
    pop_len = l_GammaS2.GetN()
    pop_file_name = bolo_name + "_GammaS2.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_GammaS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close()  

    ##################################
    #
    #
    #    B E T A   E V E N T S
    #
    #
    ##################################

    #Get surface S1 beta, Pb and gamma events
    tree.Draw(">>l_BetaS1",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0  && EIA>" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S1"] + "<0.75 && " + d_estimator["Q_S1"] + ">0.2")
    # print standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA>1 && EIB>1 && EIC<1 && EID<1 && abs(EC1-EC2)<2 && " +  d_estimator["Q_S1"] + "<0.75 && " + d_estimator["Q_S1"] + ">0.2"
    pop_len = l_BetaS1.GetN()
    pop_file_name = bolo_name + "_BetaS1.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_BetaS1.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close()   


    tree.Draw(">>l_BetaS2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC>" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S2"] + "<0.75 && " + d_estimator["Q_S2"] + ">0.2")
    pop_len = l_BetaS2.GetN()
    pop_file_name = bolo_name + "_BetaS2.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_BetaS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close()  



    ##################################
    #
    #
    #    P b    E V E N T S
    #
    #
    ##################################

    # S 1     E V E N T S
    tree.Draw(">>l_PbS1",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 &&  EIA>" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S1"] + "<0.2 " )
    print 
    pop_len = l_PbS1.GetN()
    pop_file_name = bolo_name + "_PbS1.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_PbS1.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close()    

    # S 2     E V E N T S
    tree.Draw(">>l_PbS2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC>" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S2"] + "<0.2 " )
    pop_len = l_PbS2.GetN()
    pop_file_name = bolo_name + "_PbS2.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_PbS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close() 

    del l_all
    del l_all_highE
    del l_heatonly
    del l_GammaFid

    del l_BetaS1 
    del l_PbS1 
    del l_GammaS1

    del l_BetaS2 
    del l_GammaS2 
    del l_PbS2 

# bolo_name = "FID828"
# bolo_list=["FID825", "FID824","FID828", "FID827", "FID826", "FID823", "FID810", "FID837", "FID838", "FID839", "FID821", "FID841", "FID842", "FID844", "FID845"]
bolo_list = ["FID837"]
data_dir  = "../Fond_ERA_merged/"
for bolo_name in bolo_list:
    get_population_list( bolo_name, data_dir, "data")