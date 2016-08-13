#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import os, sys
import PyROOTPlots as PyRPl
import Analysis_utilities as Ana_ut
from ROOT import *
import numpy as np
import BDT_file_handler as BDT_fh

def test(bolo_name, data_dir, tree_name):

    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)   

    tree.GetEntry(0)

    tmax = 1E6*tree.UT1 + tree.UT2

    for i in range(1, tree.GetEntries()):
        tree.GetEntry(i)
        time = 1E6*tree.UT1 + tree.UT2
        if time<tmax:
            print i, "Fail", time, tmax
            raw_input()
        else:
            tmax = float(time)

def get_FW_over_time():


    """

    Detail:

    Arguments:

    Outputs:

    """

    pass

def plot_FWHM_for_presentation(bolo_name, data_dir, tree_name):


    """

    Detail:

    Arguments:

    Outputs:

    """

    # #Load the data
    # file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    # tree       = file_tree.Get(tree_name)

    # #Define path 
    # path_name       = script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    # polar_file_name = bolo_name + "_polar_start_and_end_time_prelim.txt"

    # print path_name + polar_file_name
    # start_and_end_array = np.loadtxt(path_name + polar_file_name ,delimiter=",", usecols=(1,2))
    # script_utils.print_utility(script_utils.COL('Opening ' + path_name+polar_file_name + ' for reading', 'blue'))    

    # #Define absolute min and max of time from the _polar_start_and_end_time_prelim file. 
    # tmin = np.amin(start_and_end_array)
    # tmax = np.amax(start_and_end_array)

    # #Load the prelim cuts
    # bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts_forduration.txt")

    # hFID = TH2F("hFID", "hFID", 100000, 0,200, 100, 0, 2)    
    # hHEAT = TH2F("hHEAT", "hHEAT", 100000, 0,200, 100, 0, 2)    

    # tree.Project("hFID", "0.432*FWIB+(1-0.432)*FWID:JOUR" , bolo_TCut_line)
    # tree.Project("hHEAT", "0.57*OWC1+(1-0.57)*OWC2:JOUR" , bolo_TCut_line)

    # hFID.Draw()

    # raw_input()
    # # Get KTH histo
    # KTH_path = "../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_thresh.root"
    # hkth, f_thresh = PyRPl.open_ROOT_object(KTH_path, "hkth")

    # list_time= []

    # for i in range(1, hkth.GetNbinsX() +1):
    #     list_time.append(hkth.GetBinLowEdge(i))

    # #Load the data
    # file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    # tree       = file_tree.Get(tree_name)


    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)

    #Define path 
    path_name       = script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    polar_file_name = bolo_name + "_polar_start_and_end_time_prelim.txt"

    print path_name + polar_file_name
    start_and_end_array = np.loadtxt(path_name + polar_file_name ,delimiter=",", usecols=(1,2))
    script_utils.print_utility(script_utils.COL('Opening ' + path_name+polar_file_name + ' for reading', 'blue'))    

    #Define absolute min and max of time from the _polar_start_and_end_time_prelim file. 
    tmin = np.amin(start_and_end_array)
    tmax = np.amax(start_and_end_array)

    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts_forduration.txt")
    print bolo_TCut_line
    
    # evt_list = TEventList("evt_list")
    # tree.Draw(">>evt_list", bolo_TCut_line )
    # pop_len = evt_list.GetN()
    # tree.GetEntry(0)
    # tmax = tree.JOUR
    # list_time = []
    # list_FWIB = []
    # for k in range(pop_len):
    #     counter = evt_list.GetEntry(k)
    #     tree.GetEntry(counter)
    #     list_time.append(tree.JOUR)
    #     list_FWIB.append(tree.FWIB)

    # print len(list_time)
    # gr = TGraph(len(list_time), np.array(list_time).astype(float), np.array(list_FWIB).astype(float))
    # gr.Draw("A*")
    # raw_input()

    # 2 D threshold histogram
    binning = 2000
    days = 200
    hist_thresh      = TH2F("hist_thresh","hist_thresh",binning,0,days,1000,0,2)
    hist_FWIB      = TH2F("hist_FWIB","hist_FWIB",binning,0,days,1000,0,2)
    hist_FWID      = TH2F("hist_FWID","hist_FWID",binning,0,days,1000,0,2)
    hist_OWC1      = TH2F("hist_OWC1","hist_OWC1",binning,0,days,1000,0,2)
    hist_OWC2      = TH2F("hist_OWC2","hist_OWC2",binning,0,days,1000,0,2)

    tree.Project( "hist_thresh", "KTH:JOUR" , bolo_TCut_line)
    tree.Project( "hist_FWIB",  "FWIB:JOUR",bolo_TCut_line)
    tree.Project( "hist_FWID",  "FWID:JOUR",bolo_TCut_line)
    tree.Project( "hist_OWC1",  "OWC1:JOUR",bolo_TCut_line)
    tree.Project( "hist_OWC2",  "OWC2:JOUR",bolo_TCut_line)
    
    # hprojX_KTH    = TH1D("hprojX_KTH","hprojX_KTH",1000,0,20)
    # hprojX_FWIB    = TH1D("hprojX_FWIB","hprojX_FWIB",1000,0,20)
    # hprojX_FWID    = TH1D("hprojX_FWID","hprojX_FWID",1000,0,20)
    # hprojX_OWC1    = TH1D("hprojX_OWC1","hprojX_OWC1",1000,0,20)
    # hprojX_OWC2    = TH1D("hprojX_OWC2","hprojX_OWC2",1000,0,20)
    # hprojY_KTH    = TH1D("hprojY_KTH","hprojY_KTH",1000,0,20)
    # hprojX_FWIB    = TH1D("hprojX_FWIB","hprojX_FWIB",1000,0,20)
    # hprojX_FWID    = TH1D("hprojX_FWID","hprojX_FWID",1000,0,20)
    # hprojX_OWC1    = TH1D("hprojX_OWC1","hprojX_OWC1",1000,0,20)
    # hprojX_OWC2    = TH1D("hprojX_OWC2","hprojX_OWC2",1000,0,20)
    # hprojY_KTH    = TH1D("hprojY_KTH","hprojY_KTH",1000,0,20)

    hprojX_KTH    = hist_thresh.ProjectionX()
    # hprojX_KTH.Draw()
    # print hprojX_KTH.GetBinCenter(10)
    # raw_input()
    hprojX_FWIB = hist_FWIB.ProjectionX()
    hprojX_FWID = hist_FWID.ProjectionX()
    hprojX_OWC1 = hist_OWC1.ProjectionX()
    hprojX_OWC2 = hist_OWC2.ProjectionX()

    hprojY_KTH    = hist_thresh.ProjectionY()
    hprojY_FWIB = hist_FWIB.ProjectionY()
    hprojY_FWID = hist_FWID.ProjectionY()
    hprojY_OWC1 = hist_OWC1.ProjectionY()
    hprojY_OWC2 = hist_OWC2.ProjectionY()

    # hist_FWIB.Draw()

    print hprojY_FWIB.GetNbinsX()
    print hprojY_FWID.GetNbinsX()
    print hprojY_OWC1.GetNbinsX()
    print hprojY_FWIB.GetNbinsX()
    print hprojY_KTH.GetNbinsX()

    # raw_input()

    #Fill time and threshold lists used to build the 1D threshold histhograms
    list_time, list_fwid, list_fwib, list_owc1, list_owc2 = [], [], [], [], []
    for ix in range(1,binning) : 
        temp_fwib, temp_fwid, temp_owc1, temp_owc2=0,0,0,0
        for iy in range(1, hprojY_KTH.GetNbinsX()): 
            fwib = hprojY_FWIB.GetBinCenter(iy)
            fwid = hprojY_FWID.GetBinCenter(iy)
            owc1 = hprojY_OWC1.GetBinCenter(iy)
            owc2 = hprojY_OWC2.GetBinCenter(iy)

            if ( hist_FWIB.GetBinContent(ix,iy)!=0) :
                temp_fwib=fwib

            if ( hist_FWID.GetBinContent(ix,iy)!=0) :
                temp_fwid=fwid

            if ( hist_OWC1.GetBinContent(ix,iy)!=0) :
                temp_owc1=owc1

            if ( hist_OWC2.GetBinContent(ix,iy)!=0) :
                temp_owc2=owc2

        list_time.append(hprojX_KTH.GetBinCenter(ix))
        list_fwib.append(temp_fwib)
        list_fwid.append(temp_fwid)
        list_owc1.append(temp_owc1)
        list_owc2.append(temp_owc2)

    # ROOT says there is an issue with increasing value for the bins represented by time
    # I disagree.
    # Solution: use a graph to build the 1D threshold histhogram

    gr_fwib = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwib))
    hfwib = TH1F("hfwib", "hfwib",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hfwib.SetBinContent(k,gr_fwib.Eval(list_time[k-1]))

    hfwib.Draw()
    hist_FWIB.Draw("same")
    raw_input()

    gr_fwid = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwid))
    hfwid = TH1F("hfwid", "hfwid",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hfwid.SetBinContent(k,gr_fwid.Eval(list_time[k-1]))

    gr_owc1 = TGraph(np.size(list_time), np.array(list_time), np.array(list_owc1))
    howc1 = TH1F("howc1", "howc1",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        howc1.SetBinContent(k,gr_owc1.Eval(list_time[k-1]))

    gr_owc2 = TGraph(np.size(list_time), np.array(list_time), np.array(list_owc2))
    howc2 = TH1F("howc2", "howc2",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        howc2.SetBinContent(k,gr_owc2.Eval(list_time[k-1]))


    def compute_resolution( FWHM1,  FWHM2) :
        """
        Compute the resolution (i.e FWHM not sigma) of the best combination of two channels
        """
        s1=float(FWHM1)/2.3548
        s2=float(FWHM2)/2.3548

        w1=(s2*s2)/(s1*s1+s2*s2)

        res=np.sqrt(pow(w1*s1,2)+pow((1-w1)*s2,2))
        return 2.3548*res

    hfid = TH1F("hfid", "hfid",np.size(list_time),min(list_time), max(list_time))
    hheat = TH1F("hheat", "hheat",np.size(list_time),min(list_time), max(list_time))

    for i in range(1, np.size(list_time) +1):
        if (hfwib.GetBinContent(i)!=0and hfwid.GetBinContent(i)!=0 and howc1.GetBinContent(i)!=0 and howc2.GetBinContent(i)!=0 ):
            hfid.SetBinContent(i, compute_resolution(hfwib.GetBinContent(i), hfwid.GetBinContent(i)))
            hheat.SetBinContent(i, compute_resolution(howc1.GetBinContent(i), howc2.GetBinContent(i)))


    PyRPl.process_TH1(hheat, color = kRed, Y_title = "Heat FWHM (keV)")
    PyRPl.process_TH1(hfid, color = kBlack, Y_title = "Fiducial ionisation FWHM (keV)", X_title = "Days")

    cc = TCanvas("cc", "cc")
    cc.Divide(1,2)
    cc.cd(1)
    hheat.Draw()

    cc.cd(2)
    hfid.Draw()

    raw_input()

    cc.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_combined_resolution_over_time.eps")

bolo_name, data_dir, tree_name = "FID837", "../Fond_ERA_merged/", "data"
plot_FWHM_for_presentation(bolo_name, data_dir, tree_name)
# test(bolo_name, data_dir, tree_name)