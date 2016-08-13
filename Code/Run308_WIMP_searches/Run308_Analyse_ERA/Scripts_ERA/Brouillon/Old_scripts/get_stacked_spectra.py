#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np

def get_stacked_hist_surf(bolo_list, evt_type):

    array = np.array([])
    for bolo_name in bolo_list:
        file_path="../Analyse_" + bolo_name + "/Populations/"
        file_name_S1 = bolo_name + "_4V" + evt_type + "S1.txt"
        file_name_S2 = bolo_name + "_4V" + evt_type + "S2.txt"
        
        arr_S1 = np.loadtxt(file_path + file_name_S1)
        arr_S2 = np.loadtxt(file_path + file_name_S2)

        array = np.hstack((array, arr_S1))
        array = np.hstack((array, arr_S2))

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    h1 = TH1F ("h1", "h1", 200, 0, 100)
    for i in range(array.size):
        h1.Fill(array[i])

    h1.Scale(1/(25.*len(bolo_list)*h1.GetBinWidth(10)))
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
        cc.Print("../Figures/Stacked_spectrum_Gamma_Surf.eps")
    else:
        cc.Print("../Figures/Stacked_spectrum_" + evt_type + ".eps")

def get_stacked_hist_fid(bolo_list, evt_type):

    array = np.array([])
    for bolo_name in bolo_list:
        file_path="../Analyse_" + bolo_name + "/Populations/"
        file_name_Fid = bolo_name + "_4V" + evt_type + ".txt"
        
        arr_Fid = np.loadtxt(file_path + file_name_Fid)

        array = np.hstack((array, arr_Fid))

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    h1 = TH1F ("h1", "h1", 200, 0, 100)
    for i in range(array.size):
        h1.Fill(array[i])

    h1.Scale(1/(25.*len(bolo_list)*h1.GetBinWidth(10)))
    h1.GetXaxis().SetTitle("Heat in keV")
    h1.GetXaxis().CenterTitle(kTRUE)
    h1.GetXaxis().SetTitleSize(0.06)
    h1.GetXaxis().SetTitleOffset(0.8)

    h1.GetYaxis().SetTitle("Counts/d.keV")
    h1.GetYaxis().CenterTitle(kTRUE)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetYaxis().SetTitleOffset(0.8)
    h1.Draw("")

    cc.Print("../Figures/Stacked_spectrum_" + evt_type + ".eps")

# bolo_name = "FID828"
bolo_list=["FID825", "FID824","FID828", "FID827", "FID826", "FID823", "FID837", "FID838", "FID821", "FID841", "FID842", "FID844", "FID845"]
get_stacked_hist_surf( bolo_list, "Beta")
get_stacked_hist_surf( bolo_list, "Gamma")
get_stacked_hist_surf( bolo_list, "Pb")
get_stacked_hist_fid(bolo_list, "GammaFid")