#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *

def get_1D_plot(bolo_name, data_dir, tree_name, list_cuts1, list_cuts2):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line1='&&'.join(list_cuts1)
    cut_line2='&&'.join(list_cuts2)
    
    bin, min_X, max_X, channel_X = 200, 0.2, 2, "FWC1"
    bin, min_X, max_X, channel_X = 200, -3, 3, "EID"

    h1= TH1F("h_lowE", "h_lowE", bin, min_X, max_X)
    tree.Project("h_lowE",channel_X, cut_line1)

    h2= TH1F("h_highE", "h_highE", bin, min_X, max_X)
    tree.Project("h_highE",channel_X, cut_line2)
    
    h1.SetStats(0)
    h1.Scale(1/float(h1.Integral()))
    h2.Scale(1/float(h2.Integral()))

    h1.GetXaxis().SetTitle(channel_X)
    h1.GetXaxis().CenterTitle(kTRUE)
    h1.GetXaxis().SetTitleSize(0.06)
    h1.GetXaxis().SetTitleOffset(0.8)

    h1.GetYaxis().SetTitle("Counts")
    h1.GetYaxis().CenterTitle(kTRUE)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    gPad.SetLogy()
    # cc.SetLogy()
    h1.Draw()   
    h2.SetLineColor(kRed)
    h2.Draw("same")
    raw_input()
    # # Define path for the .txt file. Create the directory if it does not exist, then open file
    # figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    # cc.Print(figure_path_name + bolo_name + "_" + fig_title + ".eps")

def launch_plots(bolo_name, data_dir, tree_name = "data"):

    
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

    list_cuts1 = [standard_cuts]  + [d_estimator["HEAT"] + "<2"] + ["0.5*(abs(EIA)+abs(EIB)+abs(EIC)+abs(EID))<1"]+ ["abs(EC1-EC2)<2"]
    list_cuts2 = [standard_cuts]  + [d_estimator["HEAT"] + ">2"] + ["0.5*(abs(EIA)+abs(EIB)+abs(EIC)+abs(EID))<1"]+ ["abs(EC1-EC2)<2"]

    get_1D_plot(bolo_name, data_dir, tree_name, list_cuts1, list_cuts2)

    return

bolo_list= ["FID837"]
data_dir  = "../Fond_ERA_merged/"
for bolo_name in bolo_list:
    launch_plots( bolo_name, data_dir)