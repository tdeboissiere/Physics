#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 
import numpy as np





def plot_bckg_and_neutron_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

    """Plot the output of the BDT classification + data 

    
    Detail:
        Only plot the backgrounds
        They are scaled to their expected contribution

    Args:
        bolo_name           = (str) bolometer name
        mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
        analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
        bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
        simu_number         = (str) data tree simulated number
        exposure            = (float) exposure of the simulated data

    Returns:
        void

    Raises:
        void
    """       

    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type)


    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"

    #############################################
    #SGet all histograms
    ##########################################

    fWIMP     = TFile(gen_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root", "read")
    fneutron     = TFile(gen_path +  "real_neutron" + "_mass_" + str(mass) + "_tree.root", "read")
    
    tWIMP     = fWIMP.Get("tout") 
    tneutron     = fneutron.Get("tout") 
    

    ER = "(1+8./3)*0.5*(EC1+EC2)-0.33*(1.5*EIA+4*EIB+1.5*EIC+4*EID)"

    #Get the histograms for weight computation
    bin_X, min_X, max_X = 400, -1,1
    hneutron = TH1F("hneutron", "hneutron", bin_X, min_X, max_X)
    hWIMP    = TH1F("hWIMP", "hWIMP", bin_X, min_X, max_X)

    tWIMP.Project("hWIMP", "NN")
    tneutron.Project("hneutron", "NN", "weight")

    hneutron.Scale(1/hneutron.Integral())
    hWIMP.Scale(1/hWIMP.Integral())

    hWIMP.SetLineColor(kRed)
    hneutron.Draw()
    hWIMP.Draw("same")

    leg = TLegend(0.14,0.50,0.33,0.87)
    leg.AddEntry(hneutron.GetName(), "Weighted Neutron " + str(mass) + " GeV", "leg")
    leg.AddEntry(hWIMP.GetName(),"WIMP " + str(mass) + " GeV","f")
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(0)
    leg.Draw("same")

    raw_input()

    c1.Print("./Figures/" +bolo_name + "_BDT_bckg_and_neutron_hist_" +analysis_type + "_mass_" + str(mass)   + "_baseline_time.png")


bolo_name           = "FID837"
mass                = 5
# simu_number         = 1
analysis_type       = "ana_1.5_0_5"
bin_X, min_X, max_X = 100, -0.8, 1.2
exposure = 65.
d_cut       = {"ECinf": 1.5, "ECsup": 15, "EIinf": 0*0.2, "EIsup": 15, "sigma_vet": 5}

# plot_bckg_BDT_output(bolo_name, mass, analysis_type, bin_X, min_X, max_X)
# for mass in [5, 6, 7, 10, 25]:
#     for simu_number in range(1):
#         plot_bckg_and_data_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, simu_number, exposure)

for mass in [5, 6, 7, 10, 25]:
# for mass in [10, 25]:
    plot_bckg_and_neutron_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X,  exposure)
