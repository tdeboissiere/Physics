#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import os, sys
import PyROOTPlots as PyRPl

def get_plot(bolo_name, data_dir, tree_name):

    # Load start and end times and duration
    path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
    file_name= bolo_name + "_polar_start_and_end_time.txt"  
    file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")

    #Create lists to hold the contents of the .txt file + the list for each polar duration
    list_tmin, list_tmax, list_string_polar, list_duration = [], [], [], []

    #Read the file lines and fill the lists
    list_lines_start_and_end = [elem.rstrip().split(",") for elem in file_start_and_end.readlines()]
    
    for k in range(len(list_lines_start_and_end)):
        list_string_polar.append(list_lines_start_and_end[k][0])
        list_tmin.append(float(list_lines_start_and_end[k][1]))
        list_tmax.append(float(list_lines_start_and_end[k][2]))
        list_duration.append(float(list_lines_start_and_end[k][3]))
    file_start_and_end.close()

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
    
    #Load estimators
    estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
    estimator_file_name =bolo_name + "_estimators.txt" 
    d_estimator         ={}
    assert(os.path.isfile(estimator_path_name + estimator_file_name) )
    with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
        list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
        for line in list_estimator_lines:
            d_estimator[line[0]] = line[1]

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

    file_tree   = TFile(data_dir+bolo_name+"_lowmass_fond.root")
    tree        = file_tree.Get(tree_name)

    list_cuts = [standard_cuts]
    cut_line='&&'.join(list_cuts)

    list_ion_chan = ["IA", "IB", "IC", "ID"]
    list_heat_chan = ["C1_ERA", "C2_ERA"]

    list_hist_ion_time = [TH2F("FW" + chan + "_time", "FW" + chan + "_time", 100000, list_tmin[0], list_tmax[0], 100, 0, 2) for chan in list_ion_chan]
    list_hist_heat_time = [TH2F("FW" + chan + "_time", "FW" + chan + "_time", 100000, list_tmin[0], list_tmax[0], 100, 0, 2) for chan in list_heat_chan]

    list_hist_chi_ion = [TH2F("CHI" + chan + "_heat", "CHI" + chan + "_heat", 500, 0, 5, 500, 0, 5) for chan in list_ion_chan]
    list_hist_chi_heat = [TH2F("CHI" + chan + "_heat", "CHI" + chan + "_heat", 500, 0, 5, 500, 0, 5) for chan in list_heat_chan]

    for chan in list_ion_chan:
        tree.Project("FW" + chan + "_time", "FW" + chan + ":DateSec")

    for chan in list_heat_chan:
        tree.Project("FW" + chan + "_time", "FW" + chan + ":DateSec")

    for chan in list_ion_chan:
        tree.Project("CHI" + chan + "_heat", "CH" + chan + ":EC1_ERA")

    for chan in list_heat_chan:
        tree.Project("CHI" + chan + "_heat", "CHI" + chan + ":EC1_ERA")

    import PyROOTPlots as PyRPl
    for index, chan in enumerate(list_ion_chan):
        PyRPl.process_TH2(list_hist_ion_time[index], X_title = "Time", Y_title = "FW" + chan, marker_style = 2, marker_size = 0.5)
        PyRPl.process_TH2(list_hist_chi_ion[index], X_title = "ERA heat (keVee)", Y_title = "CH" + chan)

    for index, chan in enumerate(list_heat_chan):
        PyRPl.process_TH2(list_hist_heat_time[index], X_title = "Time", Y_title = "FW" + chan, marker_style = 2, marker_size = 0.5)
        PyRPl.process_TH2(list_hist_chi_heat[index], X_title = "ERA heat (keVee)", Y_title = "CHI" + chan)



    #Do the plots
    cc1=TCanvas("cc1","cc1")
    cc1.Divide(2,2)
    for i in range(4):
        cc1.cd(i+1)
        list_hist_ion_time[i].Draw()

    cc2=TCanvas("cc2","cc2")
    cc2.Divide(2,2)
    for i in range(4):
        cc2.cd(i+1)
        list_hist_chi_ion[i].Draw("colz")

    #Do the plots
    cc3=TCanvas("cc3","cc3")
    cc3.Divide(1,2)
    for i in range(2):
        cc3.cd(i+1)
        list_hist_heat_time[i].Draw()

    cc4=TCanvas("cc4","cc4")
    cc4.Divide(1,2)
    for i in range(2):
        cc4.cd(i+1)
        list_hist_chi_heat[i].Draw("colz")

    raw_input()
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc1.Print(figure_path_name + bolo_name + "_FWION_time.png")
    cc2.Print(figure_path_name + bolo_name + "_CHION_Heat.png")
    cc3.Print(figure_path_name + bolo_name + "_FWHeat_time.png")
    cc4.Print(figure_path_name + bolo_name + "_CHIHeat_Heat.png")
    
bolo_name, data_dir, tree_name = "FID837", "../Fond_ERA_merged/", "t_merged"
get_plot(bolo_name, data_dir, tree_name)