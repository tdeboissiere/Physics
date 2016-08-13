#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 
import numpy as np


def plot_bckg_BDT_output(bolo_name, mass, analysis_type, bin_X, min_X, max_X):

    """Plot the output of the BDT classification

    
    Detail:
        Only plot the backgrounds
        They are scaled to their expected contribution

    Args:
        bolo_name           = (str) bolometer name
        mass                = (int) WIMP mass
        analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
        bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram

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

    fS1Pb     = TFile(gen_path +  "S1Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Pb     = TFile(gen_path +  "S2Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Beta   = TFile(gen_path +  "S1Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Beta   = TFile(gen_path +  "S2Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Gamma  = TFile(gen_path +  "S1Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Gamma  = TFile(gen_path +  "S2Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fFidGamma = TFile(gen_path +  "FidGamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fheat     = TFile(gen_path +  "heatonly" + "_mass_" + str(mass) + "_tree.root", "read")
    fWIMP     = TFile(gen_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root", "read")
    
    tS1Pb     = fS1Pb.Get("tout")
    tS2Pb     = fS2Pb.Get("tout") 
    tS1Beta   = fS1Beta.Get("tout") 
    tS2Beta   = fS2Beta.Get("tout") 
    tS1Gamma  = fS1Gamma.Get("tout") 
    tS2Gamma  = fS2Gamma.Get("tout") 
    tFidGamma = fFidGamma.Get("tout") 
    theat     = fheat.Get("tout")
    tWIMP     = fWIMP.Get("tout") 
    
    hS1Pb     = TH1F("hS1Pb", "hS1Pb", bin_X, min_X, max_X)
    hS2Pb     = TH1F("hS2Pb", "hS2Pb", bin_X, min_X, max_X)
    hS1Beta   = TH1F("hS1Beta", "hS1Beta", bin_X, min_X, max_X)
    hS2Beta   = TH1F("hS2Beta", "hS2Beta", bin_X, min_X, max_X)
    hS1Gamma  = TH1F("hS1Gamma", "hS1Gamma", bin_X, min_X, max_X)
    hS2Gamma  = TH1F("hS2Gamma", "hS2Gamma", bin_X, min_X, max_X)
    hFidGamma = TH1F("hFidGamma", "hFidGamma", bin_X, min_X, max_X)
    hheat     = TH1F("hheat", "hheat", bin_X, min_X, max_X)
    hWIMP     = TH1F("hWIMP", "hWIMP", bin_X, min_X, max_X)

    tS1Pb.Project("hS1Pb", "NN")
    tS2Pb.Project("hS2Pb", "NN")
    tS1Beta.Project("hS1Beta", "NN")
    tS2Beta.Project("hS2Beta", "NN")
    tS1Gamma.Project("hS1Gamma", "NN")
    tS2Gamma.Project("hS2Gamma", "NN")
    tFidGamma.Project("hFidGamma", "NN")
    theat.Project("hheat", "NN")
    tWIMP.Project("hWIMP", "NN")

    #Rescale WIMP  to 1
    hWIMP.Scale(1/float(hWIMP.Integral()))

    #Rescale  bckg histos to their expected value
    hheat.Scale(float(d_scaling["prop_heatonly"])*1./float(hheat.Integral()))
    hFidGamma.Scale(float(d_scaling["prop_FidGamma"])*1./float(hFidGamma.Integral()))
    hS1Gamma.Scale(float(d_scaling["prop_S1Gamma"])*1./float(hS1Gamma.Integral()))
    hS2Gamma.Scale(float(d_scaling["prop_S2Gamma"])*1./float(hS2Gamma.Integral()))
    hS1Beta.Scale(float(d_scaling["prop_S1Beta"])*1./float(hS1Beta.Integral()))
    hS2Beta.Scale(float(d_scaling["prop_S2Beta"])*1./float(hS2Beta.Integral()))
    hS1Pb.Scale(float(d_scaling["prop_S1Pb"])*1./float(hS1Pb.Integral()))
    hS2Pb.Scale(float(d_scaling["prop_S2Pb"])*1./float(hS2Pb.Integral()))

    PyRPl.process_TH1(hS1Pb, use_fill_bool = True, color = kOrange-7)
    PyRPl.process_TH1(hS2Pb, use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(hS1Beta, use_fill_bool = True, color = kGreen+3)
    PyRPl.process_TH1(hS2Beta, use_fill_bool = True, color = kGreen-3)
    PyRPl.process_TH1(hS1Gamma, use_fill_bool = True, color = kBlue-7)
    PyRPl.process_TH1(hS2Gamma, use_fill_bool = True, color = kBlue)
    PyRPl.process_TH1(hFidGamma, use_fill_bool = True, color = kAzure+10)
    PyRPl.process_TH1(hheat, use_fill_bool = True, color = kRed)
    PyRPl.process_TH1(hWIMP, use_fill_bool = True, color = kGray)

    list_hist =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat, hWIMP]

    hs=THStack("hs", "hs")

    for hist in list_hist:
        hs.Add(hist)

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    h=TH1F("h","h", bin_X, min_X, max_X)
    PyRPl.process_TH1(h, X_title="MVA ouput", Y_title = "Events", min_Y = 1E-5, max_Y = 0.5)
    h.GetXaxis().SetTitle("MVA response")
    h.GetYaxis().SetTitle("Events")

    h.Draw()
    hs.Draw("same")

    raw_input()

    cc.Print("./Figures/" +bolo_name + "_BDT_bckg_hist_" +analysis_type + "_mass_" + str(mass) + ".png")


def plot_bckg_and_data_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, simu_number, exposure):

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

    fS1Pb     = TFile(gen_path +  "S1Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Pb     = TFile(gen_path +  "S2Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Beta   = TFile(gen_path +  "S1Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Beta   = TFile(gen_path +  "S2Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Gamma  = TFile(gen_path +  "S1Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Gamma  = TFile(gen_path +  "S2Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fFidGamma = TFile(gen_path +  "FidGamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fheat     = TFile(gen_path +  "heatonly" + "_mass_" + str(mass) + "_tree.root", "read")
    fWIMP     = TFile(gen_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root", "read")
    fdata     = TFile(gen_path +  "data" + "_mass_" + str(mass) + "_tree.root", "read")
    
    tS1Pb     = fS1Pb.Get("tout")
    tS2Pb     = fS2Pb.Get("tout") 
    tS1Beta   = fS1Beta.Get("tout") 
    tS2Beta   = fS2Beta.Get("tout") 
    tS1Gamma  = fS1Gamma.Get("tout") 
    tS2Gamma  = fS2Gamma.Get("tout") 
    tFidGamma = fFidGamma.Get("tout") 
    theat     = fheat.Get("tout")
    tWIMP     = fWIMP.Get("tout") 
    tdata     = fdata.Get("tout" + str(simu_number)) 
    
    hS1Pb     = TH1F("hS1Pb", "hS1Pb", bin_X, min_X, max_X)
    hS2Pb     = TH1F("hS2Pb", "hS2Pb", bin_X, min_X, max_X)
    hS1Beta   = TH1F("hS1Beta", "hS1Beta", bin_X, min_X, max_X)
    hS2Beta   = TH1F("hS2Beta", "hS2Beta", bin_X, min_X, max_X)
    hS1Gamma  = TH1F("hS1Gamma", "hS1Gamma", bin_X, min_X, max_X)
    hS2Gamma  = TH1F("hS2Gamma", "hS2Gamma", bin_X, min_X, max_X)
    hFidGamma = TH1F("hFidGamma", "hFidGamma", bin_X, min_X, max_X)
    hheat     = TH1F("hheat", "hheat", bin_X, min_X, max_X)
    hWIMP     = TH1F("hWIMP", "hWIMP", bin_X, min_X, max_X)
    hdata     = TH1F("hdata", "hdata", bin_X, min_X, max_X)

    tS1Pb.Project("hS1Pb", "NN")
    tS2Pb.Project("hS2Pb", "NN")
    tS1Beta.Project("hS1Beta", "NN")
    tS2Beta.Project("hS2Beta", "NN")
    tS1Gamma.Project("hS1Gamma", "NN")
    tS2Gamma.Project("hS2Gamma", "NN")
    tFidGamma.Project("hFidGamma", "NN")
    theat.Project("hheat", "NN")
    tWIMP.Project("hWIMP", "NN")
    tdata.Project("hdata", "NN")

    sum_data = float(hdata.Integral())

    #Rescale WIMP  to hdata
    hWIMP.Scale(sum_data/float(hWIMP.Integral()))

    #Rescale  bckg histos to their expected value
    hheat.Scale(float(d_scaling["prop_heatonly"])*float(d_scaling["exp_per_day"])*exposure/float(hheat.Integral()))
    hFidGamma.Scale(float(d_scaling["prop_FidGamma"])*float(d_scaling["exp_per_day"])*exposure/float(hFidGamma.Integral()))
    hS1Gamma.Scale(float(d_scaling["prop_S1Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Gamma.Integral()))
    hS2Gamma.Scale(float(d_scaling["prop_S2Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Gamma.Integral()))
    hS1Beta.Scale(float(d_scaling["prop_S1Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Beta.Integral()))
    hS2Beta.Scale(float(d_scaling["prop_S2Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Beta.Integral()))
    hS1Pb.Scale(float(d_scaling["prop_S1Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Pb.Integral()))
    hS2Pb.Scale(float(d_scaling["prop_S2Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Pb.Integral()))


    print hdata.Integral()
    print sum([hS1Pb.Integral(), hS2Pb.Integral(), hS1Beta.Integral(), hS2Beta.Integral(), hS1Gamma.Integral(), hS2Gamma.Integral(), hFidGamma.Integral(), hheat.Integral()])

    # #Rescale  bckg histos to their expected value
    # hheat.Scale(0.915494830893*sum_data/float(hheat.Integral()))
    # hFidGamma.Scale(0.077775346672*sum_data/float(hFidGamma.Integral()))
    # hS1Gamma.Scale(0.000129377041922*sum_data/float(hS1Gamma.Integral()))
    # hS2Gamma.Scale(2.52615823815e-05*sum_data/float(hS2Gamma.Integral()))
    # hS1Beta.Scale(0.00306448788106*sum_data/float(hS1Beta.Integral()))
    # hS2Beta.Scale(0.0033224289936*sum_data/float(hS2Beta.Integral()))
    # hS1Pb.Scale(0.00014262881331*sum_data/float(hS1Pb.Integral()))
    # hS2Pb.Scale(4.56381229001e-05*sum_data/float(hS2Pb.Integral()))


    PyRPl.process_TH1(hS1Pb, use_fill_bool = True, color = kOrange-7)
    PyRPl.process_TH1(hS2Pb, use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(hS1Beta, use_fill_bool = True, color = kGreen+3)
    PyRPl.process_TH1(hS2Beta, use_fill_bool = True, color = kGreen-3)
    PyRPl.process_TH1(hS1Gamma, use_fill_bool = True, color = kBlue-7)
    PyRPl.process_TH1(hS2Gamma, use_fill_bool = True, color = kBlue)
    PyRPl.process_TH1(hFidGamma, use_fill_bool = True, color = kAzure+10)
    PyRPl.process_TH1(hheat, use_fill_bool = True, color = kRed)
    PyRPl.process_TH1(hWIMP, use_fill_bool = True, color = kGray)

    list_hist =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat, hWIMP]
    list_hist_bckg =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat]

    hsum=TH1F("hsum","hsum", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hsum.SetBinContent(i, sumcontent)

    hdatasimu = TH1F("hdatasimu","hdatasimu", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        hdatasimu.SetBinContent(i, np.random.poisson(hsum.GetBinContent(i)))

    hs=THStack("hs", "hs")

    for hist in list_hist:
        hs.Add(hist)


    cc = TCanvas("cc", "cc")
    cc.Divide(1,2)
    h1=TH1F("h1","h1", bin_X, min_X, max_X)
    PyRPl.process_TH1(h1, X_title="MVA ouput Poisson simulated data hist", Y_title = "Events", min_Y = 1E-1, max_Y = 1000)
    h2=TH1F("h2","h2", bin_X, min_X, max_X)
    PyRPl.process_TH1(h2, X_title="MVA ouput true simulated data hist", Y_title = "Events", min_Y = 1E-1, max_Y = 1000)
 
    cc.cd(1)
    gPad.SetLogy()
    h1.Draw()
    hs.Draw("same")
    hdatasimu.Draw("sameE1")
    cc.cd(2)
    gPad.SetLogy()
    h2.Draw()
    hs.Draw("same")
    hdata.Draw("sameE1")
    # raw_input()

    cc.Print("./Figures/" +bolo_name + "_BDT_bckg_and_data_hist_" +analysis_type + "_mass_" + str(mass) + "_simu" + str(simu_number) + ".png")




def plot_bckg_and_real_data_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

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

    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")


    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"

    #############################################
    #SGet all histograms
    ##########################################

    fS1Pb     = TFile(gen_path +  "S1Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Pb     = TFile(gen_path +  "S2Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Beta   = TFile(gen_path +  "S1Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Beta   = TFile(gen_path +  "S2Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Gamma  = TFile(gen_path +  "S1Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Gamma  = TFile(gen_path +  "S2Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fFidGamma = TFile(gen_path +  "FidGamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fheat     = TFile(gen_path +  "heatonly" + "_mass_" + str(mass) + "_tree.root", "read")
    fWIMP     = TFile(gen_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root", "read")
    fdata     = TFile(gen_path +  "real_data" + "_mass_" + str(mass) + "_tree.root", "read")
    
    tS1Pb     = fS1Pb.Get("tout")
    tS2Pb     = fS2Pb.Get("tout") 
    tS1Beta   = fS1Beta.Get("tout") 
    tS2Beta   = fS2Beta.Get("tout") 
    tS1Gamma  = fS1Gamma.Get("tout") 
    tS2Gamma  = fS2Gamma.Get("tout") 
    tFidGamma = fFidGamma.Get("tout") 
    theat     = fheat.Get("tout")
    tWIMP     = fWIMP.Get("tout") 
    tdata     = fdata.Get("tout") 
    
    hS1Pb     = TH1F("hS1Pb", "hS1Pb", bin_X, min_X, max_X)
    hS2Pb     = TH1F("hS2Pb", "hS2Pb", bin_X, min_X, max_X)
    hS1Beta   = TH1F("hS1Beta", "hS1Beta", bin_X, min_X, max_X)
    hS2Beta   = TH1F("hS2Beta", "hS2Beta", bin_X, min_X, max_X)
    hS1Gamma  = TH1F("hS1Gamma", "hS1Gamma", bin_X, min_X, max_X)
    hS2Gamma  = TH1F("hS2Gamma", "hS2Gamma", bin_X, min_X, max_X)
    hFidGamma = TH1F("hFidGamma", "hFidGamma", bin_X, min_X, max_X)
    hheat     = TH1F("hheat", "hheat", bin_X, min_X, max_X)
    hWIMP     = TH1F("hWIMP", "hWIMP", bin_X, min_X, max_X)
    hdata     = TH1F("hdata", "hdata", bin_X, min_X, max_X)

    tS1Pb.Project("hS1Pb", "NN")
    tS2Pb.Project("hS2Pb", "NN")
    tS1Beta.Project("hS1Beta", "NN")
    tS2Beta.Project("hS2Beta", "NN")
    tS1Gamma.Project("hS1Gamma", "NN")
    tS2Gamma.Project("hS2Gamma", "NN")
    tFidGamma.Project("hFidGamma", "NN")
    theat.Project("hheat", "NN")
    tWIMP.Project("hWIMP", "NN")
    tdata.Project("hdata", "NN")

    sum_data = float(hdata.Integral())

    #Rescale WIMP  to hdata
    hWIMP.Scale(sum_data/float(hWIMP.Integral()))

    # print d_scaling

    hheat.Scale(float(d_scaling["prop_heatonly"])*float(d_scaling["exp_per_day"])*exposure/float(hheat.Integral()))
    hFidGamma.Scale(float(d_scaling["prop_FidGamma"])*float(d_scaling["exp_per_day"])*exposure/float(hFidGamma.Integral()))
    hS1Gamma.Scale(float(d_scaling["prop_S1Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Gamma.Integral()))
    hS2Gamma.Scale(float(d_scaling["prop_S2Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Gamma.Integral()))
    hS1Beta.Scale(float(d_scaling["prop_S1Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Beta.Integral()))
    hS2Beta.Scale(float(d_scaling["prop_S2Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Beta.Integral()))
    hS1Pb.Scale(float(d_scaling["prop_S1Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Pb.Integral()))
    hS2Pb.Scale(float(d_scaling["prop_S2Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Pb.Integral()))


    print hdata.Integral()
    print sum([hS1Pb.Integral(), hS2Pb.Integral(), hS1Beta.Integral(), hS2Beta.Integral(), hS1Gamma.Integral(), hS2Gamma.Integral(), hFidGamma.Integral(), hheat.Integral()])

    PyRPl.process_TH1(hS1Pb, h_title = "hS1Pb", use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(hS2Pb, h_title = "hS2Pb", use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(hS1Beta, h_title = "hS1Beta", use_fill_bool = True, color = kGreen+2)
    PyRPl.process_TH1(hS2Beta, h_title = "hS2Beta", use_fill_bool = True, color = kGreen-3)
    PyRPl.process_TH1(hS1Gamma, h_title = "hS1Gamma", use_fill_bool = True, color = kBlue-7)
    PyRPl.process_TH1(hS2Gamma, h_title = "hS2Gamma", use_fill_bool = True, color = kBlue)
    PyRPl.process_TH1(hFidGamma, h_title = "hFidGamma", use_fill_bool = True, color = kAzure+10)
    PyRPl.process_TH1(hheat, h_title = "hheat", use_fill_bool = True, color = kRed)
    PyRPl.process_TH1(hWIMP, h_title = "hWIMP", use_fill_bool = True, color = kGray)

    hS1Pb.Add(hS2Pb)
    hS1Beta.Add(hS2Beta)
    hFidGamma.Add(hS1Gamma)
    hFidGamma.Add(hS2Gamma)

    list_hist =[hS1Pb, hS1Beta, hFidGamma, hheat, hWIMP]

    hs=THStack("hs", "hs")

    for hist in list_hist:
        hs.Add(hist)

    cc = TCanvas("cc", "cc")
    h1=TH1F("h1","h1", bin_X, min_X, max_X)
    PyRPl.process_TH1(h1, X_title="BDT ouput", min_Y = 1E-1, max_Y = 20000)
 
    # Compute p-value
    hsum_bckg = TH1F("hsum_bckg", "hsum_bckg", bin_X,min_X,max_X)
    for i in range(1,bin_X+1):
        bin_content=0
        for histo in [hS1Pb, hS1Beta, hFidGamma, hheat]:
            bin_content+=histo.GetBinContent(i)
        hsum_bckg.SetBinContent(i, bin_content)

    # hsum_bckg.Draw()
    # hdata.Draw("same")
    # raw_input()
    # print hdata.Chi2Test(hsum_bckg,"UWP")
    print hdata.KolmogorovTest(hsum_bckg,"M")

    # gPad.SetLogy()
    # h1.Draw()
    # hs.Draw("same")
    # hdata.SetMarkerStyle(1)
    # hdata.Draw("sameE1")

    # leg = TLegend(0.14,0.50,0.33,0.87)
    # leg.AddEntry(hdata.GetName(), "Data", "leg")
    # leg.AddEntry(hS1Pb.GetName(),"Lead" ,"f")
    # leg.AddEntry(hS1Beta.GetName(),"Beta", "f")
    # leg.AddEntry(hFidGamma.GetName(),"Gamma", "f")
    # leg.AddEntry(hheat.GetName(),"Heat-only", "f")
    # leg.AddEntry(hWIMP.GetName(),"WIMP " + str(mass) + " GeV","f")
    # leg.SetFillColor(kWhite)
    # leg.SetBorderSize(0)
    # leg.Draw("same")

    # raw_input()

    # out_dir = script_utils.create_directory("./Figures/" +bolo_name + "/" + analysis_type + "/")
    # cc.Print(out_dir + bolo_name + "_BDT_" +analysis_type + "_mass_" + str(mass)   + ".png")

bolo_name           = "FID837"
mass                = 5
# simu_number         = 1
analysis_type       = "ana_1.5_0_5"
bin_X, min_X, max_X = 100, -2, 1.2
exposure = 65.
d_cut       = {"ECinf": 1.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}

# plot_bckg_BDT_output(bolo_name, mass, analysis_type, bin_X, min_X, max_X)
# for mass in [5, 6, 7, 10, 25]:
#     for simu_number in range(1):
#         plot_bckg_and_data_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, simu_number, exposure)

# for mass in [3, 4, 5, 6, 7, 10, 25]:
for mass in [5, 6, 7, 10, 25]:
    plot_bckg_and_real_data_BDT_output(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X,  exposure)
