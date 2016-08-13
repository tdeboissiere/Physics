#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import os
import numpy as np
import PyROOTPlots as PyRPl

#The conclusion of this tests:
#Hard to say if its really exponential (may be several populations, single burst, ...)

def launch_plots(bolo_name, data_dir, tree_name = "t_merged"):


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

    list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_estimator["HEAT"] + ">0"]+ ["abs(EC1-EC2)<2"]

    file_tree   = TFile(data_dir+bolo_name+"_lowmass_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)
    print cut_line



    arr_time = np.linspace(1406280024.25,1420977368.25,7)
         

    #This part to get the time delays once 

    list_time1 = []
    list_time2 = []
    list_time3 = []
    list_time4 = []
    list_time5 = []
    list_time6 = []

    nEntries = tree.GetEntries()

    # for i in range(nEntries):
    #     tree.GetEntry(i)
    #     time = tree.DateSec

    #     bool_fwhm = tree.FWIA<2 and tree.FWIB<2 and tree.FWIC<2 and tree.FWID<2 and tree.FWC1<2 and tree.FWC2<2
    #     bool_chi2 = tree.CHIA<6 and tree.CHIB<6 and tree.CHIC<6 and tree.CHID<6 and tree.CHIC1_ERA<2 and tree.CHIC2_ERA<2 
    #     bool_sdel = tree.SDEL>1
    #     bool_heat = tree.EIA<1 and tree.EIB<1 and tree.EIC<1 and tree.EID<1 and 0.5*(tree.EC1_ERA+tree.EC2_ERA)>2 and abs(tree.EC1_ERA-tree.EC2_ERA)<2
    #     bool_period = tree.FWIA>0 and tree.FWIB>0 and tree.FWIC>0 and tree.FWID>0 and tree.FWC1_ERA>0 and tree.FWC2_ERA>0
    #     bool_tot = bool_fwhm and bool_chi2 and bool_sdel and bool_heat and bool_period

    #     if bool_tot:
    #         if arr_time[0]<time<arr_time[1]:
    #             list_time1.append(time)

    #         elif arr_time[1]<time<arr_time[2]:
    #             list_time2.append(time)

    #         elif arr_time[2]<time<arr_time[3]:
    #             list_time3.append(time)

    #         elif arr_time[3]<time<arr_time[4]:
    #             list_time4.append(time)

    #         elif arr_time[4]<time<arr_time[5]:
    #             list_time5.append(time)

    #         elif arr_time[5]<time<arr_time[6]:
    #             list_time6.append(time)

    # list_time1 = sorted(list_time1)
    # list_time2 = sorted(list_time2)
    # list_time3 = sorted(list_time3)
    # list_time4 = sorted(list_time4)
    # list_time5 = sorted(list_time5)
    # list_time6 = sorted(list_time6)

    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period1.txt", list_time1)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period2.txt", list_time2)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period3.txt", list_time3)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period4.txt", list_time4)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period5.txt", list_time5)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period6.txt", list_time6)


    # list_timediff1 = np.array([list_time1[i]-list_time1[i-1] for i in range(1,len(list_time1))])
    # list_timediff2 = np.array([list_time2[i]-list_time2[i-1] for i in range(1,len(list_time2))])
    # list_timediff3 = np.array([list_time3[i]-list_time3[i-1] for i in range(1,len(list_time3))])
    # list_timediff4 = np.array([list_time4[i]-list_time4[i-1] for i in range(1,len(list_time4))])
    # list_timediff5 = np.array([list_time5[i]-list_time5[i-1] for i in range(1,len(list_time5))])
    # list_timediff6 = np.array([list_time6[i]-list_time6[i-1] for i in range(1,len(list_time6))])

    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period1.txt", list_timediff1)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period2.txt", list_timediff2)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period3.txt", list_timediff3)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period4.txt", list_timediff4)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period5.txt", list_timediff5)
    # np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period6.txt", list_timediff6)

    # #Quick code block to check the time distri of heat only events
    # arr3 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period3.txt")
    # arr4 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period4.txt")
    # arr5 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period5.txt")
    # arr6 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_times_period6.txt")

    # arr = np.hstack((arr3,arr4))    
    # arr = np.hstack((arr,arr5))    
    # arr = np.hstack((arr,arr6))    

    # # # list_timediff_short = [arr[i]-arr[i-1] for i in range(1,len(arr)) if arr[i] <1.415606E9] #test on one period
    # # # list_timediff_short = [arr[i]-arr[i-1] for i in range(1,len(arr)) if arr[i] <1.41716E9] #test on another period
    # # list_timediff_short = [arr[i]-arr[i-1] for i in range(1,len(arr)) if 1.4148<arr[i] <1.4156E9] #test on another period
    # list_timediff_short = [arr[i]-arr[i-1] for i in range(1,len(arr)) if 1.4122<arr[i] <1.4132E9] #test on another period
    # # list_timediff_short = [arr[i]-arr[i-1] for i in range(1,len(arr)) if 1.4103E9<arr[i] <1.4132E9] #test on another period
    # # list_timediff_short = [arr[i]-arr[i-1] for i in range(1,len(arr)) if 1.415E9 <arr[i] <1.417E9] #test on another period
    # h = TH1F("h", "h",1000, 1.4122E9, 1.4132E9)
    # ytitle = str(float(-1.4122E9 + 1.4132E9)/float(1000*24*3600))[:4]

    # for elem in arr:
    #     h.Fill(elem)
    # PyRPl.process_TH1(h, X_title = "Time", Y_title = "Counts/" + ytitle + "days")
    # h.Draw()
    # raw_input()
    # c1.Print("../Analyse_FID837/Figures/FID837_heat_rate_short_time.png")

    # hdiff = TH1F("hdiff", "hdiff",200, 0, 1000)
    # for elem in list_timediff_short:
    #     hdiff.Fill(elem)
    # PyRPl.process_TH1(hdiff, X_title = "#Delta t between events", Y_title = "Counts")
    # hdiff.Draw()
    # hdiff.Fit("expo")
    # raw_input()
    # # c1.Print("../Analyse_FID837/Figures/FID837_heat_delta_t_events.png")

    h1= TH1F("h1", "h1", bin_X, min_X, max_X)
    h2= TH1F("h2", "h2", bin_X, min_X, max_X)
    h3= TH1F("h3", "h3", bin_X, min_X, max_X)
    h4= TH1F("h4", "h4", bin_X, min_X, max_X)
    h5= TH1F("h5", "h5", bin_X, min_X, max_X)
    h6= TH1F("h6", "h6", bin_X, min_X, max_X)

    hnorm1= TH1F("hnorm1", "hnorm1", bin_X, min_X, max_X)
    hnorm2= TH1F("hnorm2", "hnorm2", bin_X, min_X, max_X)
    hnorm3= TH1F("hnorm3", "hnorm3", bin_X, min_X, max_X)
    hnorm4= TH1F("hnorm4", "hnorm4", bin_X, min_X, max_X)
    hnorm5= TH1F("hnorm5", "hnorm5", bin_X, min_X, max_X)
    hnorm6= TH1F("hnorm6", "hnorm6", bin_X, min_X, max_X)

    list_timediff1 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period1.txt")
    list_timediff2 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period2.txt")
    list_timediff3 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period3.txt")
    list_timediff4 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period4.txt")
    list_timediff5 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period5.txt")
    list_timediff6 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_heatonly_delta_t_period6.txt")

    for elem in list_timediff1:
        h1.Fill(elem)
        hnorm1.Fill(elem)
    for elem in list_timediff2:
        h2.Fill(elem)
        hnorm2.Fill(elem)
    for elem in list_timediff3:
        h3.Fill(elem)
        hnorm3.Fill(elem)
    for elem in list_timediff4:
        h4.Fill(elem)
        hnorm4.Fill(elem)
    for elem in list_timediff5:
        h5.Fill(elem)
        hnorm5.Fill(elem)
    for elem in list_timediff6:
        h6.Fill(elem)
        hnorm6.Fill(elem)

    h1.SetLineColor(kRed)
    h2.SetLineColor(kOrange-8)
    h3.SetLineColor(kGray)
    h4.SetLineColor(kGreen-3)
    h5.SetLineColor(kBlue-7)
    h6.SetLineColor(kBlack)

    hnorm1.SetLineColor(kRed)
    hnorm2.SetLineColor(kOrange-8)
    hnorm3.SetLineColor(kGray)
    hnorm4.SetLineColor(kGreen-3)
    hnorm5.SetLineColor(kBlue-7)
    hnorm6.SetLineColor(kBlack)

    hnorm1.Scale(1./hnorm1.Integral())
    hnorm2.Scale(1./hnorm2.Integral())
    hnorm3.Scale(1./hnorm3.Integral())
    hnorm4.Scale(1./hnorm4.Integral())
    hnorm5.Scale(1./hnorm5.Integral())
    hnorm6.Scale(1./hnorm6.Integral())

    list_hist = [h1, h2, h3, h4, h5, h6]
    list_max = [hist.GetMaximum() for hist in list_hist]
    max_max = max(list_max)

    h= TH1F("h", "h", bin_X, min_X, max_X)
    h.SetMaximum(1.2*max_max)
    h.SetStats(0)

    h.GetXaxis().SetTitle("#delta t between events")
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle("Counts")
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)

    cc=TCanvas("cc","cc")
    cc.Divide(2,3)
    for i,hist in enumerate(list_hist):
        cc.cd(i+1)
        gPad.SetLogy()
        PyRPl.process_TH1(list_hist[i], X_title = "#delta t between events period " + str(i+1), Y_title = "counts")
        list_hist[i].Draw("same")
    raw_input()
    cc.Print("../Analyse_FID837/Figures/FID837_heat_delta_t_over_time.png")

    #Do the plots
    cc2=TCanvas("cc2","cc2")
    list_hist = [hnorm1, hnorm2, hnorm3, hnorm4, hnorm5, hnorm6]
    list_max = [hist.GetMaximum() for hist in list_hist]
    max_max = max(list_max)
    h.SetMaximum(1.2*max_max)
    h.GetYaxis().SetTitle("Normalised counts")

    cc1=TCanvas("cc1","cc1")
    cc1.Divide(2,3)
    for i,hist in enumerate(list_hist):
        cc1.cd(i+1)
        gPad.SetLogy()
        PyRPl.process_TH1(list_hist[i], X_title = "#delta t between events period " + str(i+1), Y_title = "normalised counts")
        list_hist[i].Draw("same")


    raw_input()
    cc2.Print("../Analyse_FID837/Figures/FID837_normalised_heat_delta_t_over_time.png")


bin_X, min_X, max_X = 100, 0, 5000

bolo_name, data_dir = "FID837", "../Fond_ERA_merged/"
launch_plots(bolo_name, data_dir)