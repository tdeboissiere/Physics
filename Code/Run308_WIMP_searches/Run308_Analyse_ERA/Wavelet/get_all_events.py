#!/usr/bin/env python

import script_utils as script_utils
import os
from ROOT import *

def get_population_list(bolo_name):

    """Get all events for wavelet ana
    
    Detail:
        Get the list of all events
        We will compute their wavelet stuff

    Args:
        bolo_name    = (str) bolo name 

    Returns:
        void

    Raises:
        void
    """

    data_dir = "../../Fond_ERA_merged/"
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get("t_merged")

    #Create background hist directory
    pop_path_name = script_utils.create_directory("./Populations/" + bolo_name + "/")

    #Load estimators
    estimator_path_name = script_utils.create_directory('../../Analyse_' + bolo_name + "/Text_files/")  
    estimator_file_name =bolo_name + "_estimators.txt" 
    d_estimator         ={}
    assert(os.path.isfile(estimator_path_name + estimator_file_name) )
    with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
        list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
        for line in list_estimator_lines:
            d_estimator[line[0]] = line[1]
    
    l_all_event = TEventList("l_all_event")

    cuts = "FWIA<2&&FWIB<2&&FWIC<2&&FWID<2&&FWC1_ERA<2&&FWC2_ERA<2"
    cuts = cuts + "&&CHIA<6&&CHIB<6&&CHIC<6&&CHID<6&&CHIC1_ERA<6&&CHIC2_ERA<6 &&SDEL>1"
    cuts = cuts +  "&&" + d_estimator["HEAT"] + ">0 &&" +  d_estimator["HEAT"] + "<15 "


    #Fill event list
    tree.Draw(">>l_all_event",cuts )
    pop_len = l_all_event.GetN()
    pop_file_name = bolo_name + "_all_events.txt"
    pop_file = script_utils.open_text_file(pop_path_name, pop_file_name, 'w')
    for k in range(pop_len):
        counter = l_all_event.GetEntry(k)
        tree.GetEntry(counter)
        pop_file.write(str(tree.RUN) + "," + str(tree.SN) + "\n")
    pop_file.close()   

    del l_all_event 


bolo_list=["FID837"]
for bolo_name in bolo_list:
    get_population_list( bolo_name)


