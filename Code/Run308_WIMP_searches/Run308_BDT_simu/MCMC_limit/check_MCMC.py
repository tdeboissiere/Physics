#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 



def check_MCMC_code(bolo_name, mass, analysis_type, ERA_name, simu_number, FWHM_type):

    d_scaling = BDT_fh.open_scaling_file(bolo_name, analysis_type, FWHM_type, ERA_name)

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application" + ERA_name +"/"
    list_all_evt_class =["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly", "WIMP", "data"]
    list_bckg_evt_class =["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly"]

    #############################################
    #Store all the hist data in a dictionnary
    #Also rescale the data
    ##########################################
    d_hist={}
    d_hist_file_name={}
    for file_name in list_all_evt_class:
        if file_name == "data":
            d_hist[file_name], d_hist_file_name[file_name] = PyRPl.open_ROOT_object(gen_path +  file_name + "_mass_" + str(mass) + "_hist.root", "MVA_BDT" + str(simu_number))
        else:
            d_hist[file_name], d_hist_file_name[file_name] = PyRPl.open_ROOT_object(gen_path +  file_name + "_mass_" + str(mass) + "_hist.root", "MVA_BDT")


    sum_data = d_hist["data"].Integral()
    script_utils.print_utility("Expected heat:" + str(sum_data*float(d_scaling["prop_heatonly"])))

    #Rescale WIMP hist to 1
    d_hist["WIMP"].Scale(sum_data/float(d_hist["WIMP"].Integral()))

    #Rescale  bckg histos to their expected value
    for evt_class in list_bckg_evt_class:
        d_hist[evt_class].Scale(float(d_scaling["prop_" + evt_class])*sum_data/float(d_hist[evt_class].Integral()))

    PyRPl.process_TH1(d_hist["S1Pb"], use_fill_bool = True, color = kOrange-7)
    PyRPl.process_TH1(d_hist["S2Pb"], use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(d_hist["S1Beta"], use_fill_bool = True, color = kGreen+3)
    PyRPl.process_TH1(d_hist["S2Beta"], use_fill_bool = True, color = kGreen-3)
    PyRPl.process_TH1(d_hist["S1Gamma"], use_fill_bool = True, color = kBlue-7)
    PyRPl.process_TH1(d_hist["S2Gamma"], use_fill_bool = True, color = kBlue)
    PyRPl.process_TH1(d_hist["FidGamma"], use_fill_bool = True, color = kAzure+10)
    PyRPl.process_TH1(d_hist["heatonly"], use_fill_bool = True, color = kRed)
    PyRPl.process_TH1(d_hist["WIMP"], use_fill_bool = True, color = kGray)


    hs=THStack("hs", "hs")

    for evt_class in list_bckg_evt_class:
        hs.Add(d_hist[evt_class])

    hs.Add(d_hist["WIMP"])

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    h=TH1F("h","h", 100,-2,2)
    PyRPl.process_TH1(h, X_title="MVA ouput", Y_title = "Events", min_Y = 0.1, max_Y = 3E4)
    h.GetXaxis().SetTitle("MVA response")
    h.GetYaxis().SetTitle("Events")

    h.Draw()
    hs.Draw("same")
    d_hist["data"].Draw("sameE1")

    raw_input()

    hbckg=TH1F("hbckg","hbckg", 100,-2,2)
    for evt_class in list_bckg_evt_class:
        for i in range(1,101):
            hbckg.SetBinContent(i, d_hist[evt_class].GetBinContent(i) + hbckg.GetBinContent(i))

    # class likelihood_hist:
    #     def __call__( self, x, par ):
    #         likelihood =0
    #         for i in range(1,101):
    #             N_expected =x[1]*hbckg.GetBinContent(i)+x[0]*d_hist["WIMP"].GetBinContent(i)
    #             N_obs      = d_hist["data"].GetBinContent(i)
    #             if N_expected>0:
    #                 likelihood +=-N_expected+N_obs*math.log(N_expected)
    #         return -(likelihood + par[0])

    # f2 = TF2("r", likelihood_hist(), 0,200, 1000, 3000, 1)
    # f2.SetParameter(0,0)
    # f2.SetContour(99)
    # f2.Draw("cont4z")
    # x,y=c_double(0.), c_double(0.)
    # f2.GetMinimumXY(x,y)
    # print x.value, y.value

    raw_input()
    

bolo_name = "FID837"
mass=6
analysis_type = "standard_resolution_no_ion_cut"
FWHM_type = "standard_resolution"
ERA_name = "_ERA"
# ERA_nmae = ""
simu_number = 0

check_MCMC_code(bolo_name, mass, analysis_type, ERA_name, simu_number, FWHM_type)
