#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh



def get_fid_gamma_plot(bolo_name, data_dir, tree_name):

    """ Fid Gamma plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    # d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])


    print d_std 
    raw_input()

    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"

    h1D = TH1F("h1D", "h1D", 150, 0, 15)
    tree.Project("h1D", "EC", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut)   
    PyRPl.process_TH1(h1D, X_title = "Heat energy (keV)", Y_title = "Counts/" + str(h1D.GetBinWidth(2))[:4] + " keV") 
    h1D.Draw("")
    gPad.SetLogy()


def get_surf_bckg_plot(bolo_name, evt_type, bin, min_X, max_X):

    """ Surf plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        evt_type  = (str) event type
        bin, min_X, max_X = (int, float, float) TH1F parameters

    Returns:
        void

    Raises:
        void
    """

    arr_calib = []

    if evt_type == "Beta" or evt_type == "Pb":

        arr_S1Beta_ref = np.loadtxt("../Text_files/S1Beta_heatremoved.txt", delimiter=",").astype(float)
        arr_S2Beta_ref = np.loadtxt("../Text_files/S2Beta_heatremoved.txt", delimiter=",").astype(float)

        arr_S1Beta_ref = 0.5*(arr_S1Beta_ref[:,2] + arr_S1Beta_ref[:,3])
        arr_S2Beta_ref = 0.5*(arr_S2Beta_ref[:,2] + arr_S2Beta_ref[:,3])

        arr_S1Pb_ref = np.loadtxt("../Text_files/S1Pb_heatremoved.txt", delimiter=",").astype(float)
        arr_S2Pb_ref = np.loadtxt("../Text_files/S2Pb_heatremoved.txt", delimiter=",").astype(float)

        arr_S1Pb_ref = 0.5*(arr_S1Pb_ref[:,2] + arr_S1Pb_ref[:,3])
        arr_S2Pb_ref = 0.5*(arr_S2Pb_ref[:,2] + arr_S2Pb_ref[:,3])

        if evt_type =="Beta":
            arr_calib = np.hstack((arr_S1Beta_ref, arr_S2Beta_ref))
        elif evt_type =="Pb":
            arr_calib = np.hstack((arr_S1Pb_ref, arr_S2Pb_ref))

    file_path="../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
    fS1 = bolo_name + "_S1" + evt_type + ".txt"    
    fS2 = bolo_name + "_S2" + evt_type + ".txt"    
    arrS1 = np.loadtxt(file_path + fS1)
    arrS2 = np.loadtxt(file_path + fS2)

    arr = np.hstack((arrS1, arrS2))

    cc = TCanvas("cc", "cc")
    h1 = TH1F ("h1", "h1", bin, min_X, max_X)
    hcalib = TH1F ("hcalib", "hcalib", bin, min_X, max_X)

    PyRPl.fill_TH1(h1, arr)
    PyRPl.process_TH1(h1, X_title = "Heat (keVee)", Y_title = "Counts/" + str(h1.GetBinWidth(2))[:4] + " keV" )
    h1.Draw("")
    
    if evt_type =="Beta" or evt_type =="Pb":
        PyRPl.fill_TH1(hcalib, arr_calib)
        hcalib.Scale(h1.Integral()/hcalib.Integral())
        hcalib.SetLineColor(kRed)
        h1.SetMaximum(1.2*max( h1.GetMaximum(), hcalib.GetMaximum() ) )
        hcalib.Draw("same")

    fig_path = script_utils.create_directory("../Analyse_" + bolo_name + "/Figures/Fig_scaling/")
    cc.Print(fig_path + bolo_name + "_spectrum_" + evt_type + ".png")


def get_heatonly_bckg_plot(bolo_name, evt_type, bin, min_X, max_X):

    """ Surf plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        evt_type  = (str) event type
        bin, min_X, max_X = (int, float, float) TH1F parameters

    Returns:
        void

    Raises:
        void
    """


    file_path="../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
    f = bolo_name + "_" + evt_type + ".txt"  
    arr = np.loadtxt(file_path + f)

    cc = TCanvas("cc", "cc")
    h1 = TH1F ("h1", "h1", bin, min_X, max_X)
    PyRPl.fill_TH1(h1, arr)

    gPad.SetLogy()

    PyRPl.process_TH1(h1, X_title = "Heat (keVee)", Y_title = "Counts/" + str(h1.GetBinWidth(2))[:4] + " keV" )
    h1.Draw("")

    fig_path = script_utils.create_directory("../Analyse_" + bolo_name + "/Figures/Fig_scaling/")
    cc.Print(fig_path + bolo_name + "_spectrum_" + evt_type + ".png")

bolo_name = "FID837"
data_dir  = "../Fond_ERA_merged/"

get_surf_bckg_plot(bolo_name, "Beta",  200, 2.5, 60)
get_surf_bckg_plot(bolo_name, "Gamma",  200, 2, 15)
get_surf_bckg_plot(bolo_name, "Pb",  100, 7, 50)
get_heatonly_bckg_plot(bolo_name, "heatonly",  150, 0.5, 15)


# get_hist(bolo_name, "FidGamma",  200, 0, 15)
# get_hist(bolo_name, "heatonly",  200, 0, 15)