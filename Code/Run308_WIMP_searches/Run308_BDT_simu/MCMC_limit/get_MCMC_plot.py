#!/usr/bin/env python
# -*- coding: utf-8 -*-


import emcee
import triangle
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from ROOT import *
import script_utils as script_utils
import os
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_utils as BDT_ut
import BDT_file_handler as BDT_fh

def log_factorial(k):
    if k==0:
        return 0
    else:
        return 0.5*math.log(2*math.pi*k) + k*math.log(k) - k + math.log(1+1./(12*k) + 1/(288.*k**2) -139./(51840*k**3)-571./(2488320*k**4) + 163879./(209018880*k**5))

def add_array(array_list):

    #Add all arrays in list element wise

    # There should be more than 1 array
    assert(len(array_list) > 1)
    added_array = array_list[0] +  array_list[1]
    for arr in array_list[2:]:
        added_array = added_array +  arr
    return added_array


    
def get_MCMC_limit(mass,bolo_name, bin_X, min_X, max_X, analysis_type, simu_number, FWHM_type, ERA_name):

    d_scaling = BDT_fh.open_scaling_file(bolo_name, analysis_type, FWHM_type, ERA_name)

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application" + ERA_name +"/"

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
    hdata     = TH1F("hdata", "hdata", bin_X, min_X, max_X)
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
    tdata.Project("hdata", "NN")

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

    #Sum the surface histograms
    hS1Pb.Add(hS1Pb, hS2Pb)
    hS1Beta.Add(hS1Beta, hS2Beta)
    hS1Gamma.Add(hS1Gamma, hS2Gamma)

    hbckg = TH1F ("hbckg", "hbckg", bin_X, min_X, max_X)
    hbckg.Add(hbckg, hS1Pb)
    hbckg.Add(hbckg, hS1Beta)
    hbckg.Add(hbckg, hS1Gamma)
    hbckg.Add(hbckg, hFidGamma)
    hbckg.Add(hbckg, hheat)
    hbckg.SetLineColor(kBlue)

    arr_data = np.array([hdata.GetBinContent(i) for i in range(1,bin_X+1)])
    arr_bckg = np.array([hbckg.GetBinContent(i) for i in range(1,bin_X+1)])
    arr_log_data = np.array([log_factorial(hdata.GetBinContent(i)) for i in range(1,bin_X+1)])
    arr_WIMP = np.array([hWIMP.GetBinContent(i) for i in range(1, bin_X+1)])

    #new method to remove log() evaluating to zero
    arr_exp = arr_bckg + arr_WIMP
    list_zero_index = np.where(arr_exp==0)

    arr_data = np.delete(arr_data, list_zero_index)
    arr_WIMP = np.delete(arr_WIMP, list_zero_index)
    arr_log_data = np.delete(arr_log_data, list_zero_index)
    arr_bckg = np.delete(arr_bckg, list_zero_index)

    def lnprior(theta):
        bckg_norm, WIMP_norm = theta
        if 1000 < bckg_norm < 6000  and 1E-2< WIMP_norm < 1000:
            return 0.0
        return -np.inf


    def lnlike(theta, data, log_data, bckg, WIMP):
        bckg_norm, WIMP_norm = theta
        exp = bckg_norm*bckg +WIMP_norm*WIMP
        return np.sum(-exp + data*np.log(exp) - log_data)

    def lnprob(theta, data, log_data, bckg, WIMP):
        # print "roger3"
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, data, log_data, bckg, WIMP)


    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [4800, 1], bounds = ((1000, 6000), (1E-2, 1000)), args=(arr_data, arr_log_data, arr_bckg, arr_WIMP))
    m_ml= result["x"]
    print m_ml

    ndim, nwalkers = 2, 400
    pos = [[4800, 1] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(arr_data, arr_log_data, arr_bckg, arr_WIMP))

    # Clear and run the production chain.
    script_utils.print_utility("  Running MCMC  ")
    sampler.run_mcmc(pos, 500, rstate0=np.random.get_state())
    script_utils.print_utility("  Done.  ")

    pl.clf()
    fig, axes = pl.subplots(2, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].set_ylabel("Bckg events")

    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].set_ylabel("WIMP events")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("step number")

    fig.tight_layout(h_pad=0.0)
    fig.savefig("MCMC_samples.png")

    # Make the triangle plot.
    burnin = 200
    print sampler.chain.shape
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

    fig = triangle.corner(samples[0:samples.shape[0]:10,:], labels=["Bckg events", "WIMP events"])
    fig.savefig("MCMC_distri.png")
    pl.close("all")
    # print np.percentile(samples[0:samples.shape[0]:10,1], 90)
    return np.percentile(samples[0:samples.shape[0]:10,1], 90)



bolo_name = "FID837"
bin_X, min_X, max_X = 100, -2, 2
analysis_type = "standard_resolution_no_ion_cut"
ERA_name = ""
FWHM_type = "standard_resolution"


mass_list=[25]
list_result=[]

for mass in mass_list: 
    list_result=[]
    for simu_number in range(1):
        print simu_number
        result = get_MCMC_limit(str(mass), bolo_name, bin_X, min_X, max_X, analysis_type, simu_number, FWHM_type, ERA_name)