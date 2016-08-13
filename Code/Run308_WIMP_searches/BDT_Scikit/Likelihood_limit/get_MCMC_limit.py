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
import os,sys, math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_utils as BDT_ut
import BDT_file_handler as BDT_fh
sys.path.append("../")
import data_preparation as dp
import xgboost as xgb


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

def get_projected_MCMC_limit(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):

    """
    Detail:
        Plot results

    Args:
        d_test (dict)                         = dict with test data
        d_eval (dict)                         = dict with eval data 
        d_event_dir (dict)                    = dict to get the proper directory of each event class 
        WIMP_mass (str)                       = WIMP mass
        bolo_name (str)                       = bolometer name 
        analysis_type (str)                   = type of analysis (which box cut)
        MVA_tag (str)                         = indicates which scaling file to use
        exposure (float)                      = exposure in days 
        bin_X, min_X, max_X (int float float) = TH1F parameters

    Returns:
        void

    Raises:
        void     
    """

    try:
        kwargs["weight_dir"]
    except KeyError:
        sys.exit()

    #Get scaling dict to set the weights
    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

    key_heat = ""
    for key in d_test.keys():
        if "heat" in key:
            key_heat = key

    # Get classifier
    model_dir = script_utils.create_directory("../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
    modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
    if kwargs.has_key("classifier_name"):
        modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
    bst = xgb.Booster({'nthread':16}, model_file = modelfile)

    #Get predictions on test sample
    d_pred = {}
    d_hist = {}
    d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
                     "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, key_heat: kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}

    for event_type in d_test.keys():
        d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values) )
        d_hist[event_type] = TH1F("h" + event_type+ WIMP_mass, "h" + event_type+ WIMP_mass, bin_X, min_X, max_X)
        PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
        PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
        if "WIMP" not in event_type:
            d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
            print "Event type:", event_type, "\tExpected #:", float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure

    #get predictions on data
    hdata = TH1F("hdata" + WIMP_mass, "hdata" + WIMP_mass, bin_X, min_X, max_X)
    PyRPl.fill_TH1(hdata, bst.predict( xgb.DMatrix(d_eval["realdata"].iloc[:,:].values) ) )
    d_hist["WIMP_mass_" + WIMP_mass].Scale(1./d_hist["WIMP_mass_" + WIMP_mass].Integral())

    d_hist["S1Pb"].Add(d_hist["S2Pb"])
    d_hist["S1Beta"].Add(d_hist["S2Beta"])
    d_hist["FidGamma"].Add(d_hist["S1Gamma"])
    d_hist["FidGamma"].Add(d_hist["S2Gamma"])

    list_hist =[d_hist["S1Pb"], d_hist["S1Beta"], d_hist["FidGamma"], d_hist[key_heat], d_hist["WIMP_mass_" + WIMP_mass]]
    hs=THStack("hs" + WIMP_mass, "hs" + WIMP_mass)
    for hist in list_hist:
        hs.Add(hist)

    list_hist =[d_hist["S1Pb"], d_hist["S1Beta"], d_hist["FidGamma"], d_hist[key_heat], d_hist["WIMP_mass_" + WIMP_mass]]
    hsum_bckg=TH1F("hsum_bckg"+WIMP_mass,"hsum_bckg"+WIMP_mass, bin_X, min_X, max_X)
    for i in range(1, bin_X+1): 
        hsum_bckg.SetBinContent(i, sum([h.GetBinContent(i) for h in list_hist[:-1]]))

    list_lim = []

    for niter in range(5):

        print WIMP_mass, niter
        hobs_sim = TH1F("hobs_sim"+str(niter)+WIMP_mass, "hobs_sim"+str(niter)+WIMP_mass, bin_X, min_X, max_X)
        for i in range(1, bin_X+1): 
            hobs_sim.SetBinContent(i, np.random.poisson(hsum_bckg.GetBinContent(i)))

        arr_data = np.array([hobs_sim.GetBinContent(i) for i in range(1,bin_X+1)])
        arr_bckg = np.array([hsum_bckg.GetBinContent(i) for i in range(1,bin_X+1)])
        arr_log_data = np.array([log_factorial(hobs_sim.GetBinContent(i)) for i in range(1,bin_X+1)])
        arr_WIMP = np.array([d_hist["WIMP_mass_" + WIMP_mass].GetBinContent(i) for i in range(1, bin_X+1)])

        #new method to remove log() evaluating to zero
        arr_exp = arr_bckg + arr_WIMP
        list_zero_index = np.where(arr_exp==0)

        arr_data = np.delete(arr_data, list_zero_index)
        arr_WIMP = np.delete(arr_WIMP, list_zero_index)
        arr_log_data = np.delete(arr_log_data, list_zero_index)
        arr_bckg = np.delete(arr_bckg, list_zero_index)

        def lnprior(theta):
            bckg_norm, WIMP_norm = theta
            if 0.9 < bckg_norm < 1.1  and 1E-2< WIMP_norm < 1000:
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
        result = op.minimize(nll, [1, 1], bounds = ((1E-2,100), (1E-2, 1000)), args=(arr_data, arr_log_data, arr_bckg, arr_WIMP))
        m_ml= result["x"]
        print m_ml

        ndim, nwalkers = 2, 100
        pos = [[1, 1] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(arr_data, arr_log_data, arr_bckg, arr_WIMP))

        # Clear and run the production chain.
        script_utils.print_utility("  Running MCMC  ")
        sampler.run_mcmc(pos, 500, rstate0=np.random.get_state())
        script_utils.print_utility("  Done.  ")

        pl.clf()
        fig, axes = pl.subplots(2, 1, sharex=True, figsize=(8, 9))
        axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
        axes[0].yaxis.set_major_locator(MaxNLocator(5))
        axes[0].set_ylabel("$Bckg$")

        axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
        axes[1].yaxis.set_major_locator(MaxNLocator(5))
        axes[1].set_ylabel("$WIMP$")
        axes[1].set_yscale("log")
        axes[1].set_xlabel("step number")

        fig.tight_layout()
        fig.savefig("line-time.png")

        # Make the triangle plot.
        burnin = 200
        print sampler.chain.shape
        samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

        fig = triangle.corner(samples[0:samples.shape[0]:10,:], labels=["$Bckg$", "$WIMP$"])
        fig.savefig("line-triangle.png")
        pl.close("all")
        print np.percentile(samples[0:samples.shape[0]:10,1], 90)
        # print np.percentile(samples[0:samples.shape[0]:100,1], 90)
        list_lim.append(np.percentile(samples[0:samples.shape[0]:10,1], 90))
        raw_input()

    return list_lim


bolo_name = "FID837"
data_dir = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
analysis_type = "ana_0.5_0_5"
exposure = 65.
bin_X, min_X, max_X = 500, -20, 20
MVA_tag = ""
list_mass = ["4", "5", "6", "7","10", "25"]
list_mass = ["25"]
# list_mass = ["10"]
list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test"]
d_list_variables = {}
d_list_variables["3"] =  ["EC1","EC2", "EIA", "EIC", "EFID"]
d_list_variables["4"] =  ["EC1","EC2", "EIA", "EIC", "EFID"]
d_list_variables["5"] = ["EC1","EC2", "EIA", "EIC", "EFID"]
d_list_variables["6"] = ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
d_list_variables["7"] =  ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
d_list_variables["10"] =  ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
d_list_variables["25"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"] 

d_lim = {}

for mass in list_mass:
    d_list_variables[mass] = ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]

#Loop over masses
for WIMP_mass in list_mass:
    d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                    "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

    script_utils.print_utility("Processing WIMP mass: " + str(WIMP_mass) + " GeV")

    #Load data
    d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[WIMP_mass], datasplit=1)
    d_eval  = dp.get_eval_array(bolo_name, analysis_type, "", d_list_variables[WIMP_mass])

    d_lim[WIMP_mass] = get_projected_MCMC_limit(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, weight_dir = "With_weights")
    print WIMP_mass, d_lim[WIMP_mass]

# with open("./Text_files/" + bolo_name + "_" + analysis_type + "_NWIMP_at_limit_MCMC.txt", "w") as f:
#     for WIMP_mass in list_mass:
#         f.write(WIMP_mass + "," + ",".join([str(elem) for elem in d_lim[WIMP_mass]]) + "\n")