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

    #Load duration
    duration_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
    duration_file_name =bolo_name + "_polar_start_and_end_time.txt" 
    duration         = 0
    assert(os.path.isfile(duration_path_name + duration_file_name) )
    with open(duration_path_name + duration_file_name, 'r') as duration_file: 
        list_duration_lines                                   = duration_file.readlines()[0].rstrip().split(",") 
        duration                                              = float(list_duration_lines[3])/(24*3600)
        script_utils.print_utility(script_utils.COL("Duration (in days) = "  + str(duration) , "blue"))

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


    l_heatonly_2D       = TEventList("l_heatonly_2D")
    l_heatonly        = TEventList("l_heatonly")


    l_GammaFid = TEventList("l_GammaFid")

    l_BetaS1 = TEventList("l_BetaS1")
    l_PbS1 = TEventList("l_PbS1")
    l_GammaS1 = TEventList("l_GammaS1")

    l_BetaS2 = TEventList("l_BetaS2")
    l_PbS2 = TEventList("l_PbS2")
    l_GammaS2 = TEventList("l_GammaS2")

    print standard_cuts

    #Then create various .txt files for the population selection
    tree.Draw(">>l_heatonly_2D",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2")
    print standard_cuts +  "&&" + d_estimator["HEAT"] + ">0  && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"]
    pop_len = l_heatonly_2D.GetN()
    pop_file_name = bolo_name + "_heatonly_2D.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly_2D.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
    pop_file.close()

    #Then create various .txt files for the population selection
    tree.Draw(">>l_heatonly",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2")
    pop_len = l_heatonly.GetN()
    pop_file_name = bolo_name + "_heatonly.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_heatonly.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
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
    tree.Draw(">>l_PbS1",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 &&  EIA>" + d_FWHM["FWIA"] +" && EIB>" + d_FWHM["FWIB"] +"&& EIC<" + d_FWHM["FWIC"] +"&& EID<" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S1"] + "<0.2 &&" +  d_estimator["Q_S1"] + ">0.04")
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
    tree.Draw(">>l_PbS2",standard_cuts +  "&&" + d_estimator["HEAT"] + ">0 && EIA<" + d_FWHM["FWIA"] +" && EIB<" + d_FWHM["FWIB"] +"&& EIC>" + d_FWHM["FWIC"] +"&& EID>" + d_FWHM["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_estimator["Q_S2"] + "<0.2 &&" +  d_estimator["Q_S2"] + ">0.04")
    pop_len = l_PbS2.GetN()
    pop_file_name = bolo_name + "_PbS2.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_PbS2.GetEntry(k)
        tree.GetEntry(counter)
        est_Heat = float(d_estimator["HEAT"][:5])*tree.EC1 + (1- float(d_estimator["HEAT"][:5]))*tree.EC2
        pop_file.write(str(est_Heat) + "\n")
    pop_file.close() 

    del l_heatonly
    del l_heatonly_2D
    del l_GammaFid

    del l_BetaS1 
    del l_PbS1 
    del l_GammaS1

    del l_BetaS2 
    del l_GammaS2 
    del l_PbS2 


def get_hist_surf(bolo_list, evt_type, duration):

    array = np.array([])
    for bolo_name in bolo_list:
        file_path="../Analyse_" + bolo_name + "/Populations/"
        file_name_S1 = bolo_name + "_" + evt_type + "S1.txt"
        file_name_S2 = bolo_name + "_" + evt_type + "S2.txt"
        
        arr_S1 = np.loadtxt(file_path + file_name_S1)
        arr_S2 = np.loadtxt(file_path + file_name_S2)

        array = np.hstack((array, arr_S1))
        array = np.hstack((array, arr_S2))

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    h1 = TH1F ("h1", "h1", 1000, 0, 400)
    for i in range(array.size):
        h1.Fill(array[i])

    h1.Scale(1/(duration*h1.GetBinWidth(10)))
    h1.GetXaxis().SetTitle("Heat in keV")
    h1.GetXaxis().CenterTitle(kTRUE)
    h1.GetXaxis().SetTitleSize(0.06)
    h1.GetXaxis().SetTitleOffset(0.8)

    h1.GetYaxis().SetTitle("Counts/d.keV")
    h1.GetYaxis().CenterTitle(kTRUE)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetYaxis().SetTitleOffset(0.8)
    h1.Draw("")

    if evt_type == "Gamma":
        cc.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_spectrum_Gamma_Surf.eps")
        fhist= TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_spectrum_Gamma_Surf_highmass.root", "recreate")
        h1.Write()
        fhist.Close()
        del h1
        del fhist
    else:
        cc.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name +"_spectrum_" + evt_type + ".eps")
        fhist= TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_spectrum_" + evt_type +"_highmass.root", "recreate")
        h1.Write()
        fhist.Close()
        del h1
        del fhist
def get_hist_fid(bolo_list, evt_type, duration):

    array = np.array([])
    for bolo_name in bolo_list:
        file_path="../Analyse_" + bolo_name + "/Populations/"
        file_name_Fid = bolo_name + "_" + evt_type + ".txt"
        
        arr_Fid = np.loadtxt(file_path + file_name_Fid)

        array = np.hstack((array, arr_Fid))

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    h1 = TH1F ("h1", "h1", 1000, 0, 400)
    for i in range(array.size):
        h1.Fill(array[i])
 
    h1.Scale(1/(duration*h1.GetBinWidth(10)))
    h1.GetXaxis().SetTitle("Heat in keV")
    h1.GetXaxis().CenterTitle(kTRUE)
    h1.GetXaxis().SetTitleSize(0.06)
    h1.GetXaxis().SetTitleOffset(0.8)

    h1.GetYaxis().SetTitle("Counts/d.keV")
    h1.GetYaxis().CenterTitle(kTRUE)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetYaxis().SetTitleOffset(0.8)
    h1.Draw("")

    cc.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_spectrum_" + evt_type + ".eps")
    fhist= TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_spectrum_" + evt_type +"_highmass.root", "recreate")
    h1.Write()
    fhist.Close()
    del h1
    del fhist


def get_hist2D_heatonly(bolo_name, evt_type):

# file_path     ="../Analyse_" + bolo_name + "/Populations/"
# file_name = bolo_name + "_" + evt_type + ".txt"
# arr_heat      = np.loadtxt(file_path + file_name, delimiter="," , usecols = (2,3))

    file_path     ="../Analyse_" + bolo_name + "/Populations/"
    file_name = bolo_name + "_" + evt_type + ".txt"
    arr_heat      = np.loadtxt(file_path + file_name, delimiter="," , usecols = (2,3))
    
    
    cc2 = TCanvas("cc2", "cc2")
    # gPad.SetLogy()
    h2 = TH2F ("h2", "h2", 1000, -2, 40, 1000, -2, 40)
    for i in range(arr_heat.shape[0]):
        h2.Fill(arr_heat[i][0],arr_heat[i][1])

    h2.GetXaxis().SetTitle("Heat EC1 (keV)")
    h2.GetXaxis().CenterTitle(kTRUE)
    h2.GetXaxis().SetTitleSize(0.06)
    h2.GetXaxis().SetTitleOffset(0.8)

    h2.GetYaxis().SetTitle("Heat EC2 (keV)")
    h2.GetYaxis().CenterTitle(kTRUE)
    h2.GetYaxis().SetTitleSize(0.06)
    h2.GetYaxis().SetTitleOffset(0.8)
    h2.Draw("")

    cc2.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_spectrum_" + evt_type + ".eps")
    fhist= TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_spectrum_" + evt_type +"_highmass.root", "recreate")
    h2.Write()
    fhist.Close()
    del h2
    del fhist


bolo_list=["FID837"]
duration = 70.54 #FID837
# duration = 22.9832767778 #FID824
# duration = 23.640195809#FID825
# duration = 31.796886 #FID841
# duration = 17.4513580069#FID823
# duration = 18.2513209491#FID838
# duration = 22.99351142131#FID827


# data_dir  = "../Fond_ERA_merged/"
# for bolo_name in bolo_list:
#     get_population_list( bolo_name, data_dir, "data")

get_hist_surf( bolo_list, "Beta", duration)
# get_hist_surf( bolo_list, "Gamma", duration)
get_hist_surf( bolo_list, "Pb", duration)
# get_hist_fid(bolo_list, "GammaFid", duration)
# get_hist_fid(bolo_list, "heatonly", duration)
# get_hist2D_heatonly(bolo_list[0], "heatonly_2D")

