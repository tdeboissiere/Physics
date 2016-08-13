#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math,os,sys
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 
import numpy as np
import Poisson_90CL as PoissonCL
import xgboost as xgb
sys.path.append("../../")
import data_preparation as dp

def optimize_BDT_cut(bolo_name, WIMP_mass, d_cut, analysis_type, MVA_tag, bin_X, min_X, max_X, exposure, list_variables, d_bounds, **kwargs):

    """Find the cut efficiency for background and signal

    
    Detail:
        Get the BDT cut efficiency for signal
        And the bckg expected number of events after cut

    Args:
        bolo_name             = (str) bolometer name
        mass                  = (int) WIMP mass
        d_cut                 = (dict) analysis cut dict
        analysis_type         = (str) name of analysis (name indicates which ion cut, which resolution...)
        MVA_tag (str)         = indicates which scaling file to use
        bin_X, min_X, max_X   = (int, float, float) = settings for BDT histogram
        exposure              = (float) exposure of the simulated data
        list_variables (list) = list of variables to retain for BDT
        d_bounds (list)         = lower and upper bound of the sensitivity function

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

    # if WIMP_mass == 5:
    #     list_variables = ["EIA", "EIB", "EID", "EIC", "EC1", "EC2"]

    d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                         "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly"}
    WIMP_mass = str(WIMP_mass)
    d_event_dir["WIMP_mass_" + WIMP_mass] = "WIMP"

    d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit=1)

    # Get classifier
    model_dir = script_utils.create_directory("../../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
    modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
    if kwargs.has_key("classifier_name"):
        modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
    bst = xgb.Booster({'nthread':16}, model_file = modelfile)

    #Get predictions on test sample
    d_pred = {}
    d_hist = {}
    d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
                     "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, "heatonly":kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}
    for event_type in d_test.keys():
        d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values) )
        d_hist[event_type] = TH1F("h" + event_type + str(WIMP_mass), "h" + event_type + str(WIMP_mass), bin_X, min_X, max_X)
        PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
        PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
        if "WIMP" not in event_type:
            d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
    d_hist["WIMP_mass_" + WIMP_mass].Scale(1./d_hist["WIMP_mass_" + WIMP_mass].Integral())

    list_hist_bckg =[d_hist["S1Pb"], d_hist["S2Pb"], d_hist["S1Beta"], d_hist["S2Beta"], d_hist["S1Gamma"], d_hist["S2Gamma"], d_hist["FidGamma"], d_hist["heatonly"]]

    hsum_bckg=TH1F("hsum_bckg" + str(WIMP_mass),"hsum_bckg" + str(WIMP_mass), bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hsum_bckg.SetBinContent(i, sumcontent)

    # hsum_bckg.Draw()
    # raw_input()

    class Signal_eff:
       def __call__( self, x, par ):

            bin_number = d_hist["WIMP_mass_" + WIMP_mass].FindBin(x[0])
            integ = float(d_hist["WIMP_mass_" + WIMP_mass].Integral(bin_number, bin_X))
            return par[0] + integ

    class Bckg_exp:
       def __call__( self, x, par ):

            bin_number = hsum_bckg.FindBin(x[0])
            integ = float(hsum_bckg.Integral(bin_number, bin_X))
            return par[0] + integ

    class Sensitivity:
       def __call__( self, x, par ):

            bin_number_sig = d_hist["WIMP_mass_" + WIMP_mass].FindBin(x[0])
            bin_number_bckg = hsum_bckg.FindBin(x[0])
            eff_sig = float(d_hist["WIMP_mass_" + WIMP_mass].Integral(bin_number_sig, bin_X))
            exp_bckg = hsum_bckg.Integral(bin_number_bckg, bin_X)

            vec_proba = [TMath.PoissonI(i, exp_bckg) for i in range(200)]    
            lim_Poisson_bckg = np.sum(np.array([PoissonCL.compute_90CL_limit(i)*vec_proba[i] for i in range(200)]))

            if eff_sig<=0:
                return 1E10
            else:
                return lim_Poisson_bckg/eff_sig + par[0]


    h = TH1F("h", "h",100, d_bounds[WIMP_mass][0], d_bounds[WIMP_mass][1])
    PyRPl.process_TH1(h, X_title = "BDT cut", Y_title = "Sensitivity (a.u.)")
    h.SetMinimum(1)
    h.SetMaximum(1E6)
    h.Draw()

    fopt = TF1("fopt", Sensitivity(), d_bounds[WIMP_mass][0], d_bounds[WIMP_mass][1], 1)
    fopt.SetParameter(0,0)
    fopt.SetNpx(100)
    fopt.Draw("same")

    gPad.SetLogy()
    raw_input()

    fsig_eff = TF1("fsig_eff", Signal_eff(), 0,5, 1)
    fsig_eff.SetParameter(0,0)
    fsig_eff.SetNpx(500)

    fbckg_exp = TF1("fbckg_exp", Bckg_exp(),0,5,1)
    fbckg_exp.SetParameter(0,0)

    min_val_X = fopt.GetMinimumX(d_bounds[WIMP_mass][0], d_bounds[WIMP_mass][1])

    print "For mass:", WIMP_mass, "Expect:", fbckg_exp(min_val_X), "bckg events for Signal eff:", fsig_eff.Eval(min_val_X), "%"
    return [min_val_X,fsig_eff.Eval(min_val_X)]


def save_optimize(bolo_name, WIMP_mass, d_cut, analysis_type,MVA_tag, bin_X, min_X, max_X, exposure, d_list_variables, d_bounds):
    """Find the cut efficiency for background and signal

    
    Detail:
        Get the BDT cut efficiency for signal
        And the bckg expected number of events after cut

    Args:
        bolo_name           = (str) bolometer name
        mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
        analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
        MVA_tag (str)       = indicates which scaling file to use
        bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
        exposure            = (float) exposure of the simulated data
        d_list_variables (dict)      = list of variables to retain for BDT
        d_bounds (list)         = lower and upper bound of the sensitivity function

    Returns:
        void

    Raises:
        void
    """       

    d_BDT_cut={}
    for m in list_mass:
        d_BDT_cut[m] = optimize_BDT_cut(bolo_name, m, d_cut, analysis_type,MVA_tag, bin_X, min_X, max_X, exposure, d_list_variables[str(m)], d_bounds, weight_dir = "With_weights")


    # with open("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + "_" + str(exposure) + ".txt", "w") as fc:
    with open("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + ".txt", "w") as fc:
        for key in sorted(d_BDT_cut.keys()):
            fc.write(",".join( [str(key), str( d_BDT_cut[key][0] ), str( d_BDT_cut[key][1] )]) + "\n")

bolo_name           = "FID837"
list_mass           = [3,4,5,6,7,10,25]
# list_mass           = [5,6,7,10,25]
# list_mass         = [25]
analysis_type       = "ana_0.5_0_5"
bin_X, min_X, max_X = 2000, -20, 20
d_bounds_0_5 = {"3":[2,3], "4":[2.5,4], "5":[3,8], "6":[3,6], "7":[3,10], "10":[5,10], "25":[4,12]}
# d_bounds_1_5 = {"3":[0, 10], "4":[0, 10], "5":[0, 10], "6":[0, 10], "7":[0, 10], "10":[0, 10], "25":[0, 10]}
exposure            = 66
MVA_tag             = ""
list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test", "HR"]
d_list_variables = {}
d_list_variables["3"] =  ["EC1", "EC2", "EIA", "EIC", "EFID"]
d_list_variables["4"] =  ["EC1", "EC2", "EIA", "EIC", "EFID"]
d_list_variables["5"] = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
d_list_variables["6"] = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
d_list_variables["7"] =  ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
d_list_variables["10"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]
d_list_variables["25"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test", "prod"]
# for mass in list_mass:
#         d_list_variables[str(mass)] = ["EC1", "EC2", "EIA", "EIB","EIC", "EID"]    


d_cut               = {"ECinf": 0.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}


save_optimize(bolo_name, list_mass,  d_cut, analysis_type,MVA_tag, bin_X, min_X, max_X, exposure, d_list_variables, d_bounds_0_5)