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


def optimize_BDT_cut_billard(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

    """Find the cut efficiency for background and signal

    
    Detail:
        Get the BDT cut efficiency for signal
        And the bckg expected number of events after cut

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

    # hsum.Scale(1./float(hsum.Integral()))

    class Signal_eff:
       def __call__( self, x, par ):

            bin_number = hWIMP.FindBin(x[0])
            integ = float(hWIMP.Integral(bin_number, max_X))
            return par[0] + integ

    class Bckg_exp:
       def __call__( self, x, par ):

            bin_number = hsum.FindBin(x[0])
            integ = float(hsum.Integral(bin_number, max_X))
            return par[0] + integ

    class limit_ratio:
       def __call__( self, x, par ):

            bin_number_sig = hWIMP.FindBin(x[0])
            bin_number_bckg = hsum.FindBin(x[0])
            #Signal eff for this cut value
            eff_sig = float(hWIMP.Integral(bin_number_sig, max_X))
            #expected number of bckg events for this cut value
            exp_bckg = int(hsum.Integral(bin_number_bckg, max_X))
            #Compute expected 90 CL poisson number of events
            vec_proba = [TMath.PoissonI(i, exp_bckg) for i in range(400)]    
            # raw_input()        
            # lim_Poisson_bckg = PoissonCL.compute_90CL_limit(exp_bckg)
            # print lim_Poisson_bckg, np.sum(np.array([PoissonCL.compute_90CL_limit(i)*vec_proba[i] for i in range(100)]))
            # raw_input()
            lim_Poisson_bckg = np.sum(np.array([PoissonCL.compute_90CL_limit(i)*vec_proba[i] for i in range(100)]))

            if eff_sig<=0:
                return 1E10
            else:
                return lim_Poisson_bckg/eff_sig + par[0]

    fbckg_exp = TF1("fbckg_exp", Bckg_exp(),-1,1,1)
    fbckg_exp.SetParameter(0,0)
    # print fbckg_exp.Eval(0.2)
    # print fbckg_exp.Eval(0.5)

    fopt = TF1("fopt", limit_ratio(), 0.1, 0.5, 1)
    fopt.SetParameter(0,0)
    # fopt.SetNpx(500)

    fopt.Draw()
    raw_input()

    fsig_eff = TF1("fsig_eff", Signal_eff(), -1, 1, 1)
    fsig_eff.SetParameter(0,0)
    fsig_eff.SetNpx(500)
    
    print "Cut_val:", fopt.GetMinimumX(), "For mass:", mass, "Expect:", fbckg_exp(fopt.GetMinimumX()), "bckg events for Signal eff:", fsig_eff.Eval(fopt.GetMinimumX()), "%"

    return [fopt.GetMinimumX(),fsig_eff.Eval(fopt.GetMinimumX())]

def optimize_BDT_cut(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

    """Find the cut efficiency for background and signal

    
    Detail:
        Get the BDT cut efficiency for signal
        And the bckg expected number of events after cut

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

    # #Rescale  bckg histos to their expected value
    # hheat.Scale(float(d_scaling["rate_heatonly"])*exposure*float(d_bckg_cut_eff["heatonly"])/float(hheat.Integral()))
    # hFidGamma.Scale(float(d_scaling["rate_FidGamma"])*exposure*float(d_bckg_cut_eff["FidGamma"])/float(hFidGamma.Integral()))
    # hS1Gamma.Scale(float(d_scaling["rate_S1Gamma"])*exposure*float(d_bckg_cut_eff["S1Gamma"])/float(hS1Gamma.Integral()))
    # hS2Gamma.Scale(float(d_scaling["rate_S2Gamma"])*exposure*float(d_bckg_cut_eff["S2Gamma"])/float(hS2Gamma.Integral()))
    # hS1Beta.Scale(float(d_scaling["rate_S1Beta"])*exposure*float(d_bckg_cut_eff["S1Beta"])/float(hS1Beta.Integral()))
    # hS2Beta.Scale(float(d_scaling["rate_S2Beta"])*exposure*float(d_bckg_cut_eff["S2Beta"])/float(hS2Beta.Integral()))
    # hS1Pb.Scale(float(d_scaling["rate_S1Pb"])*exposure*float(d_bckg_cut_eff["S1Pb"])/float(hS1Pb.Integral()))
    # hS2Pb.Scale(float(d_scaling["rate_S2Pb"])*exposure*float(d_bckg_cut_eff["S2Pb"])/float(hS2Pb.Integral()))

    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")

    hheat.Scale(float(d_scaling["prop_heatonly"])*float(d_scaling["exp_per_day"])*exposure/float(hheat.Integral()))
    hFidGamma.Scale(float(d_scaling["prop_FidGamma"])*float(d_scaling["exp_per_day"])*exposure/float(hFidGamma.Integral()))
    hS1Gamma.Scale(float(d_scaling["prop_S1Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Gamma.Integral()))
    hS2Gamma.Scale(float(d_scaling["prop_S2Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Gamma.Integral()))
    hS1Beta.Scale(float(d_scaling["prop_S1Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Beta.Integral()))
    hS2Beta.Scale(float(d_scaling["prop_S2Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Beta.Integral()))
    hS1Pb.Scale(float(d_scaling["prop_S1Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Pb.Integral()))
    hS2Pb.Scale(float(d_scaling["prop_S2Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Pb.Integral()))


    list_hist_bckg =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat]

    hsum=TH1F("hsum","hsum", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hsum.SetBinContent(i, sumcontent)

    # hsum.Scale(1./float(hsum.Integral()))

    class Signal_eff:
       def __call__( self, x, par ):

            bin_number = hWIMP.FindBin(x[0])
            integ = float(hWIMP.Integral(bin_number, max_X))
            return par[0] + integ

    class Bckg_exp:
       def __call__( self, x, par ):

            bin_number = hsum.FindBin(x[0])
            integ = float(hWIMP.Integral(bin_number, max_X))
            return par[0] + integ

    class limit_ratio:
       def __call__( self, x, par ):

            bin_number_sig = hWIMP.FindBin(x[0])
            bin_number_bckg = hsum.FindBin(x[0])
            #Signal eff for this cut value
            eff_sig = float(hWIMP.Integral(bin_number_sig, max_X))
            #expected number of bckg events for this cut value
            exp_bckg = int(hsum.Integral(bin_number_bckg, max_X))
            #Compute expected 90 CL poisson number of events
            lim_Poisson_bckg = PoissonCL.compute_90CL_limit(exp_bckg)
            if eff_sig<=0:
                return 1E10
            else:
                return lim_Poisson_bckg/eff_sig + par[0]


    fopt = TF1("fopt", limit_ratio(), -1, 1, 1)
    fopt.SetParameter(0,0)
    fopt.SetNpx(100)

    fopt.Draw()
    raw_input()

    fsig_eff = TF1("fsig_eff", Signal_eff(), -1, 1, 1)
    fsig_eff.SetParameter(0,0)
    fsig_eff.SetNpx(500)

    fbckg_exp = TF1("fbckg_exp", Bckg_exp(),-1,1,1)
    fbckg_exp.SetParameter(0,0)

    print "For mass:", mass, "Expect:", fbckg_exp(fopt.GetMinimumX()), "bckg events for Signal eff:", fsig_eff.Eval(fopt.GetMinimumX()), "%"

    return [fopt.GetMinimumX(),fsig_eff.Eval(fopt.GetMinimumX())]


def save_optimize(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):
    """Find the cut efficiency for background and signal

    
    Detail:
        Get the BDT cut efficiency for signal
        And the bckg expected number of events after cut

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

    d_BDT_cut={}
    for m in list_mass:
        d_BDT_cut[m] = optimize_BDT_cut(bolo_name, m, d_cut, analysis_type, bin_X, min_X, max_X, exposure)


    # with open("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + "_" + str(exposure) + ".txt", "w") as fc:
    with open("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + ".txt", "w") as fc:
        for key in sorted(d_BDT_cut.keys()):
            fc.write(",".join( [str(key), str( d_BDT_cut[key][0] ), str( d_BDT_cut[key][1] )]) + "\n")

bolo_name           = "FID837"
# mass                = 6
list_mass = [3, 4, 5,6,7,10,25]
# list_mass = [5,6,7,10,25]
analysis_type       = "ana_0.5_0.3_5"
bin_X, min_X, max_X = 200, -1, 1
exposure = 66
d_cut       = {"ECinf": 0.5, "ECsup": 15, "EIinf": 0.3, "EIsup": 15, "sigma_vet": 5}

# optimize_BDT_cut(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure)
# optimize_BDT_cut_billard(bolo_name, 5, d_cut, analysis_type, bin_X, min_X, max_X, exposure)
save_optimize(bolo_name, list_mass,  d_cut, analysis_type, bin_X, min_X, max_X, exposure)
