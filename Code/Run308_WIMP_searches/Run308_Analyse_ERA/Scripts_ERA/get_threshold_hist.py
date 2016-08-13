#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import numpy as np
import sys
import Analysis_utilities as Ana_ut

def get_threshold_hist(bolo_name, data_dir, tree_name ):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A ROOT file bolo_name + "_threshold.ROOT"
    it contains the threshold histograms/graphs
    """
    
    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)

    #Define path for txt files. Create the directory if it does not exist
    path_name       = script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    polar_file_name = bolo_name + "_polar_start_and_end_time_prelim.txt"

    print path_name + polar_file_name
    start_and_end_array = np.loadtxt(path_name + polar_file_name ,delimiter=",", usecols=(1,2))
    script_utils.print_utility(script_utils.COL('Opening ' + path_name+polar_file_name + ' for reading', 'blue'))    

    #Define absolute min and max of time from the _polar_start_and_end_time_prelim file. 
    tmin = np.amin(start_and_end_array)
    tmax = np.amax(start_and_end_array)

    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts_forduration.txt")
    # bolo_TCut_line = "APAT>1&&KTH>0&&KTH<1"

    bolo_TCut_line_heat = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    # Heat rate histogram
    njour = int((tmax-tmin)/(86400.))
    hheat      = TH1F("hheat","hheat",njour,tmin,tmax)
    print hheat.GetBinWidth(20)/(86400)
    tree.Draw("1E6*UT1+UT2>>hheat", bolo_TCut_line + "&& EC>1.5 && 0.5*(EIB+EID)<0")

    # 2 D threshold histogram
    hist_thresh      = TH2F("hist_thresh","hist_thresh",10000,tmin,tmax,1000,0,20)
    tree.Draw("KTH:1E6*UT1+UT2>>hist_thresh", bolo_TCut_line)
    
    hprojY    = TH1D("hprojY","hprojY",1000,0,20)
    hprojY    = hist_thresh.ProjectionY()
    hprojY.SetName("hprojY")
    
    hprojX    = TH1D("hprojX","hprojX",1000,0,20)
    hprojX    = hist_thresh.ProjectionX()
    hprojX.SetName("hprojX")
    
    htime     = TH2F("htime","htime",10000,tmin,tmax,100,0,2000)
    tree.Draw("0.5*(EC1+EC2):1E6*UT1+UT2>>htime")

    htime_FWHM     = TH2F("htime_FWHM","htime_FWHM",10000,tmin,tmax,100,0,2000)
    tree.Draw("0.5*(EC1+EC2):1E6*UT1+UT2>>htime_FWHM", bolo_TCut_line)
    
    hprojtime = TH1D("hprojtime","hprojtime",10000,tmin,tmax)
    hprojtime = htime.ProjectionX()
    hprojtime.SetName("hprojtime")    

    hprojtime_FWHM = TH1D("hprojtime_FWHM","hprojtime_FWHM",10000,tmin,tmax)
    hprojtime_FWHM = htime_FWHM.ProjectionX()
    hprojtime_FWHM.SetName("hprojtime_FWHM")   

    #Fill time and threshold lists used to build the 1D threshold histhograms
    list_time, list_kth = [], []
    for ix in range(1,10000) : 
        temp_thresh=0
        for iy in range(1, 1000): 
            c=hprojY.GetBinCenter(iy)
            if ( hist_thresh.GetBinContent(ix,iy)!=0) :
                temp_thresh=c
        
        list_time.append(hprojX.GetBinCenter(ix))
        list_kth.append(temp_thresh)

    # ROOT says there is an issue with increasing value for the bins represented by time
    # I disagree.
    # Solution: use a graph to build the 1D threshold histhogram
    gr = TGraph(np.size(list_time), np.array(list_time), np.array(list_kth))
    hkth = TH1F("hkth", "hkth",np.size(list_time),tmin,tmax)
    for k in range(np.size(list_time)):
        hkth.SetBinContent(k,gr.Eval(list_time[k]))
    
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh.root"
    f_thresh       = script_utils.open_ROOT_file(ROOT_path_name, ROOT_file_name, "recreate")
    hkth.Write()
    hheat.Write()
    hprojtime.Write()
    hprojtime_FWHM.Write()
    gr.Write()
    f_thresh.Close()

def get_threshold_hist_for_plot(bolo_name, data_dir, tree_name ):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A ROOT file bolo_name + "_threshold.ROOT"
    it contains the threshold histograms/graphs
    """
    
    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)

    #Define path for txt files. Create the directory if it does not exist
    path_name       = script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    polar_file_name = bolo_name + "_polar_start_and_end_time_prelim.txt"

    print path_name + polar_file_name
    start_and_end_array = np.loadtxt(path_name + polar_file_name ,delimiter=",", usecols=(1,2))
    script_utils.print_utility(script_utils.COL('Opening ' + path_name+polar_file_name + ' for reading', 'blue'))    

    #Define absolute min and max of time from the _polar_start_and_end_time_prelim file. 
    tmin = np.amin(start_and_end_array)
    tmax = np.amax(start_and_end_array)

    tmin = 0
    tmax = 200

    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts_forduration.txt")
    # bolo_TCut_line = "APAT>1&&KTH>0&&KTH<1"

    # 2 D threshold histogram
    hist_thresh      = TH2F("hist_thresh","hist_thresh",10000,tmin,tmax,1000,0,20)
    tree.Draw("KTH:JOUR>>hist_thresh", bolo_TCut_line)
    
    hprojY    = TH1D("hprojY","hprojY",1000,0,20)
    hprojY    = hist_thresh.ProjectionY()
    hprojY.SetName("hprojY")
    
    hprojX    = TH1D("hprojX","hprojX",1000,0,20)
    hprojX    = hist_thresh.ProjectionX()
    hprojX.SetName("hprojX")
    
    htime     = TH2F("htime","htime",10000,tmin,tmax,100,0,2000)
    tree.Draw("0.5*(EC1+EC2):JOUR>>htime")

    htime_FWHM     = TH2F("htime_FWHM","htime_FWHM",10000,tmin,tmax,100,0,2000)
    tree.Draw("0.5*(EC1+EC2):JOUR>>htime_FWHM", bolo_TCut_line)
    
    hprojtime = TH1D("hprojtime","hprojtime",10000,tmin,tmax)
    hprojtime = htime.ProjectionX()
    hprojtime.SetName("hprojtime")    

    hprojtime_FWHM = TH1D("hprojtime_FWHM","hprojtime_FWHM",10000,tmin,tmax)
    hprojtime_FWHM = htime_FWHM.ProjectionX()
    hprojtime_FWHM.SetName("hprojtime_FWHM")   

    #Fill time and threshold lists used to build the 1D threshold histhograms
    list_time, list_kth = [], []
    for ix in range(1,10000) : 
        temp_thresh=0
        for iy in range(1, 1000): 
            c=hprojY.GetBinCenter(iy)
            if ( hist_thresh.GetBinContent(ix,iy)!=0) :
                temp_thresh=c
        
        list_time.append(hprojX.GetBinCenter(ix))
        list_kth.append(temp_thresh)

    # ROOT says there is an issue with increasing value for the bins represented by time
    # I disagree.
    # Solution: use a graph to build the 1D threshold histhogram
    gr = TGraph(np.size(list_time), np.array(list_time), np.array(list_kth))
    hkth = TH1F("hkth", "hkth",np.size(list_time),tmin,tmax)
    for k in range(np.size(list_time)):
        hkth.SetBinContent(k,gr.Eval(list_time[k]))
    
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh_day.root"
    f_thresh       = script_utils.open_ROOT_file(ROOT_path_name, ROOT_file_name, "recreate")
    hkth.Write()
    hprojtime.Write()
    hprojtime_FWHM.Write()
    gr.Write()
    f_thresh.Close()