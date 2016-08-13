#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math, sys
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 
sys.path.append("../")
import data_preparation as dp
import xgboost as xgb
import numpy as np
import matplotlib.pylab as plt

def get_projected_likelihood_limit(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):

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


    list_min_X = []

    for niter in range(200):

        print WIMP_mass, niter
        hobs_sim = TH1F("hobs_sim"+str(niter)+WIMP_mass, "hobs_sim"+str(niter)+WIMP_mass, bin_X, min_X, max_X)
        for i in range(1, bin_X+1): 
            hobs_sim.SetBinContent(i, 1.01*np.random.poisson(hsum_bckg.GetBinContent(i)))

        class likelihood_hist:
            def __call__( self, x, par ):
                likelihood =0
                for i in range(1,bin_X+1):
                    N_expected =hsum_bckg.GetBinContent(i)+x[0]*d_hist["WIMP_mass_" + WIMP_mass].GetBinContent(i)
                    N_obs      = hobs_sim.GetBinContent(i)
                    if N_expected>0:
                        likelihood +=-N_expected+N_obs*math.log(N_expected)
                return -(likelihood + par[0])

        f1 = TF1(str(niter), likelihood_hist(), 0,2000, 1)
        f1.SetParameter(0,0)
        f1.Draw()
        raw_input()

        list_min_X.append(f1.GetMinimumX()*d_hist["WIMP_mass_" + WIMP_mass].Integral())

    list_min_X = np.array(list_min_X)
    p= np.percentile(np.array(list_min_X), 90)

    # l = np.where(list_min_X<p)
    # print len(list_min_X[l])/float(len(list_min_X))
    # plt.hist(list_min_X, bins=50)
    # plt.show()
    # raw_input()
    raw_input()

    return p

bolo_name = "FID837"
data_dir = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
analysis_type = "ana_0.5_0_5"
exposure = 65.
bin_X, min_X, max_X = 500, -20, 20
MVA_tag = ""
list_mass = [ "4", "5", "6", "7","10", "25"]
# list_mass = ["7"]
list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test"]
d_list_variables = {}
d_list_variables["3"] =  ["EC1","EC2", "EIA", "EIC", "EFID"]
d_list_variables["4"] =  ["EC1","EC2", "EIA", "EIC", "EFID"]
d_list_variables["5"] = ["EC1","EC2", "EIA", "EIC", "EFID"]
d_list_variables["6"] = ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
d_list_variables["7"] =  ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
d_list_variables["10"] =  ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
d_list_variables["25"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"] 

for mass in list_mass:
    d_list_variables[mass] = ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]

d_lim = {}

#Loop over masses
for WIMP_mass in list_mass:
    d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                    "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

    script_utils.print_utility("Processing WIMP mass: " + str(WIMP_mass) + " GeV")

    #Load data
    d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[WIMP_mass], datasplit=1)
    d_eval  = dp.get_eval_array(bolo_name, analysis_type, "", d_list_variables[WIMP_mass])

    d_lim[WIMP_mass] = get_projected_likelihood_limit(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, weight_dir = "With_weights")
    print WIMP_mass, d_lim[WIMP_mass]

# with open("./Text_files/" + bolo_name + "_" + analysis_type + "_NWIMP_at_limit.txt", "w") as f:
#     for WIMP_mass in list_mass:
#         f.write(WIMP_mass + "," + str(d_lim[WIMP_mass]) + "\n")

# 3,87.9727185458
# 4,104.439615032
# 5,58.6568101614
# 6,32.0067650672
# 7,19.1987715304
# 10,1.85864952439
# 25 0
