#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 
import numpy as np
import Poisson_90CL as PoissonCL

def plot_BDT_sensi(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

    """Plot BDT sensitivity nicely

    
    Detail:
        void

    Args:
        bolo_name           = (str) bolometer name
        mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
        analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
        bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
        exposure            = (float) exposure of the simulated data

    Returns:
        void

    Raises:
        void
    """       

    d_scaling = BDT_fh.open_scaling_file(bolo_name, d_cut, "")
    d_bckg_cut_eff = BDT_fh.open_bckg_cuteff_file(bolo_name, analysis_type, "")

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


    #Rescale WIMP  to hdata
    hWIMP.Scale(1./float(hWIMP.Integral()))

    #Rescale  bckg histos to their expected value
    hheat.Scale(float(d_scaling["rate_heatonly"])*exposure*float(d_bckg_cut_eff["heatonly"])/float(hheat.Integral()))
    hFidGamma.Scale(float(d_scaling["rate_FidGamma"])*exposure*float(d_bckg_cut_eff["FidGamma"])/float(hFidGamma.Integral()))
    hS1Gamma.Scale(float(d_scaling["rate_S1Gamma"])*exposure*float(d_bckg_cut_eff["S1Gamma"])/float(hS1Gamma.Integral()))
    hS2Gamma.Scale(float(d_scaling["rate_S2Gamma"])*exposure*float(d_bckg_cut_eff["S2Gamma"])/float(hS2Gamma.Integral()))
    hS1Beta.Scale(float(d_scaling["rate_S1Beta"])*exposure*float(d_bckg_cut_eff["S1Beta"])/float(hS1Beta.Integral()))
    hS2Beta.Scale(float(d_scaling["rate_S2Beta"])*exposure*float(d_bckg_cut_eff["S2Beta"])/float(hS2Beta.Integral()))
    hS1Pb.Scale(float(d_scaling["rate_S1Pb"])*exposure*float(d_bckg_cut_eff["S1Pb"])/float(hS1Pb.Integral()))
    hS2Pb.Scale(float(d_scaling["rate_S2Pb"])*exposure*float(d_bckg_cut_eff["S2Pb"])/float(hS2Pb.Integral()))

    list_hist_bckg =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat]

    hsum=TH1F("hsum","hsum", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hsum.SetBinContent(i, sumcontent)


    class limit_ratio:
       def __call__( self, x, par ):

            bin_number_sig = hWIMP.FindBin(x[0])
            bin_number_bckg = hsum.FindBin(x[0])
            #Signal eff for this cut value
            eff_sig = float(hWIMP.Integral(bin_number_sig, bin_X))
            #expected number of bckg events for this cut value
            exp_bckg = hsum.Integral(bin_number_bckg, bin_X)
            #Compute expected 90 CL poisson number of events
            vec_proba = [TMath.PoissonI(i, exp_bckg) for i in range(2000)]    
            lim_Poisson_bckg = np.sum(np.array([PoissonCL.compute_90CL_limit(i)*vec_proba[i] for i in range(2000)]))

            if eff_sig<=0:
                return 1E10
            else:
                return lim_Poisson_bckg/eff_sig + par[0]

    func_sensi = TF1("func_sensi", limit_ratio(), -1,10,1)
    func_sensi.SetParameter(0,0)
    PyRPl.process_TF1(func_sensi, color = kOrange+7, line_width = 3 )



    h = TH1F("h", "h",100, -1,1.8)
    h.SetMaximum(1E5)    
    h.SetMinimum(0.1)
    PyRPl.process_TH1(h, X_title = "BDT cut", Y_title = "Ratio R (a.u.)")

    cc = TCanvas("cc", "cc")
    h.Draw()
    gPad.SetLogy()
    func_sensi.Draw("same")

    l1 = TLine(-1,2.3,1.8,2.3)
    l1.SetLineColor(kBlack)
    l1.SetLineWidth(2)
    l1.SetLineStyle(7)
    l1.Draw("same")

    min_of_fun = func_sensi.GetMinimum()
    l2 = TLine(-1,min_of_fun,1.8,min_of_fun)
    l2.SetLineColor(kRed)
    l2.SetLineWidth(2)
    l2.SetLineStyle(7)
    l2.Draw("same")

    minX_of_fun = func_sensi.GetMinimumX()
    l3 = TLine(minX_of_fun,0.1,minX_of_fun,1E5)
    l3.SetLineColor(kOrange+8)
    l3.SetLineWidth(2)
    l3.SetLineStyle(7)
    l3.Draw("same")

    #Some legend
    tpoisson = TText(0.3,1E-10,"truc")
    tpoisson.SetTextSize(0.04)
    tpoisson.SetTextColor(kBlack)
    tpoisson.DrawText(-0.9,1,"Bckg free sensitivity")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kRed)
    tBDT.DrawText(-0.9,11,"Best BDT sensitivity")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kGray+2)
    tBDT.SetTextAngle(80)
    tBDT.DrawText(0.72,29.68,"Lose signal efficiency")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kGray+2)
    tBDT.DrawText(-0.9,2000,"Background dominated")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kOrange+8)
    tBDT.DrawText(0.4,0.5,"Best BDT cut value = 0.37")

    # raw_input()
    fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/")
    cc.Print(fig_dir + bolo_name + "_mass_" + str(mass) + "_sensitivity.eps")

bolo_name           = "FID837"
list_mass = [5,6,7,10,25]
list_mass = [5]
analysis_type       = "ana_1.5_0_5"
bin_X, min_X, max_X = 2000, -1, 1
exposure = 66
d_cut       = {"ECinf": 1.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}

for mass in list_mass:
    plot_BDT_sensi(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure)