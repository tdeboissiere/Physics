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


def get_BDT_efficiency_curves(d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):

    """
    Detail:
        Obtain BDT efficiency curves as a function of BDT cut

    Args:
        d_test (dict)                         = dict with test data
        d_event_dir (dict)                    = dict to get the proper directory of each event class 
        classifier_type (str)                 = type of classifier
        WIMP_mass (str)                       = WIMP mass
        bolo_name (str)                       = bolometer name 
        analysis_type (str)                   = type of analysis (which box cut)
        exposure (float)                      = exposure in days 
        bin_X, min_X, max_X (int float float) = TH1F parameters

    Returns:
        void

    Raises:
        void     
    """

    #Get scaling dict to set the weights
    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

    try:
        kwargs["weight_dir"]
    except KeyError:
        sys.exit()

    #Get heat _fraction
    heat_fraction = kwargs["classifier_name"][13:]

    key_heat = ""
    for key in d_test.keys():
        if "heat" in key:
            key_heat = key

    # Get classifier
    model_dir = script_utils.create_directory("../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")    
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
        d_hist[event_type] = TH1F("h" + event_type, "h" + event_type, bin_X, min_X, max_X)
        PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
        PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
        if "WIMP" not in event_type:
            d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
        else:
            d_hist["WIMP_mass_" + WIMP_mass].Scale(1./d_hist["WIMP_mass_" + WIMP_mass].Integral())

    list_hist_bckg =[d_hist["S1Pb"], d_hist["S2Pb"], d_hist["S1Beta"], d_hist["S2Beta"], d_hist["S1Gamma"], d_hist["S2Gamma"], d_hist["FidGamma"], d_hist[key_heat]]

    hsum_bckg=TH1F("hsum_bckg","hsum_bckg", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hsum_bckg.SetBinContent(i, sumcontent)

    ###########
    #For tests
    ###########
    # hsum_bckg.Smooth(500)
    # ssum_bckg = TSpline5(hsum_bckg)
    # ssum_bckg.SetLineColor(kRed)

    # d_hist["WIMP_mass_" + WIMP_mass].Smooth(500)
    # sWIMP = TSpline5(d_hist["WIMP_mass_" + WIMP_mass])
    # sWIMP.SetLineColor(kRed)

    class Signal_eff:
       def __call__( self, x, par ):

            bin_number = d_hist["WIMP_mass_" + WIMP_mass].FindBin(x[0])
            integ = float(d_hist["WIMP_mass_" + WIMP_mass].Integral(bin_number, bin_X))
            return par[0] + integ

    class Bckg_eff:
       def __call__( self, x, par ):

            bin_number = hsum_bckg.FindBin(x[0])
            integ = float(hsum_bckg.Integral(bin_number, in_X))/float(hsum_bckg.Integral())
            return par[0] + integ

    feff_bckg = TF1("eff_bckg_mass_" + WIMP_mass + "_heat_fraction" + heat_fraction , Bckg_eff(), -10,8, 1)
    feff_bckg.SetParameter(0,0)
    feff_bckg.SetNpx(500)

    feff_WIMP = TF1("eff_signal_mass_" + WIMP_mass + "_heat_fraction" + heat_fraction , Signal_eff(), -10,8, 1)
    feff_WIMP.SetParameter(0,0)
    feff_WIMP.SetNpx(500)

    eff_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
    feff = TFile(eff_dir + bolo_name + "_WIMP_mass_" + WIMP_mass + "_BDT_cut_eff.root", "update")
    feff_WIMP.Write()
    feff_bckg.Write()
    feff.Close()


def get_BDT_sensitivity_curve(d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):

    """
    Detail:
        Obtain BDT sensitivity curve as a function of BDT cut

    Args:
        d_test (dict)                         = dict with test data
        d_event_dir (dict)                    = dict to get the proper directory of each event class 
        classifier_type (str)                 = type of classifier
        WIMP_mass (str)                       = WIMP mass
        bolo_name (str)                       = bolometer name 
        analysis_type (str)                   = type of analysis (which box cut)
        exposure (float)                      = exposure in days 
        bin_X, min_X, max_X (int float float) = TH1F parameters

    Returns:
        void

    Raises:
        void     
    """

    #Get scaling dict to set the weights
    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

    try:
        kwargs["weight_dir"]
    except KeyError:
        sys.exit()

    #Get heat _fraction
    heat_fraction = kwargs["classifier_name"][13:]

    key_heat = ""
    for key in d_test.keys():
        if "heat" in key:
            key_heat = key

    # Get classifier
    model_dir = script_utils.create_directory("../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")    
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
        d_hist[event_type] = TH1F("h" + event_type, "h" + event_type, bin_X, min_X, max_X)
        PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
        PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
        if "WIMP" not in event_type:
            d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
        else:
            d_hist["WIMP_mass_" + WIMP_mass].Scale(1./d_hist["WIMP_mass_" + WIMP_mass].Integral())

    list_hist_bckg =[d_hist["S1Pb"], d_hist["S2Pb"], d_hist["S1Beta"], d_hist["S2Beta"], d_hist["S1Gamma"], d_hist["S2Gamma"], d_hist["FidGamma"], d_hist[key_heat]]

    hsum_bckg=TH1F("hsum_bckg","hsum_bckg", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hsum_bckg.SetBinContent(i, sumcontent)

    class Sensitivity:
       def __call__( self, x, par ):

            bin_number_sig = d_hist["WIMP_mass_" + WIMP_mass].FindBin(x[0])
            bin_number_bckg = hsum_bckg.FindBin(x[0])
            eff_sig = float(d_hist["WIMP_mass_" + WIMP_mass].Integral(bin_number_sig, bin_X))
            exp_bckg = int(hsum_bckg.Integral(bin_number_bckg, bin_X))
            lim_Poisson_bckg = PoissonCL.compute_90CL_limit(exp_bckg)

            if eff_sig<=0:
                return 1E10
            else:
                return lim_Poisson_bckg/eff_sig + par[0]

    class Sensitivity_billard:
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

    fsensi = TF1("sensivity_mass_" + WIMP_mass + "_heat_fraction" + heat_fraction, Sensitivity(), 0,7, 1)
    fsensi.SetParameter(0,0)
    fsensi_billard = TF1("sensivity_billard_mass_" + WIMP_mass + "_heat_fraction" + heat_fraction, Sensitivity_billard(), 0,7, 1)
    fsensi_billard.SetLineColor(kBlue)
    fsensi_billard.SetParameter(0,0)
    fsensi.SetNpx(1000)
    fsensi_billard.SetNpx(1000)

    # fsensi.Draw()
    # fsensi_billard.Draw("same")
    # gPad.SetLogy()
    # raw_input()

    sensi_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
    file_sensi = TFile(sensi_dir + bolo_name + "_WIMP_mass_" + WIMP_mass + "_BDT_sensitivity.root", "update")
    fsensi.Write()
    fsensi_billard.Write()
    file_sensi.Close()

    del fsensi 
    del fsensi_billard 
    del file_sensi

def get_plot_BDT_sensitivity_curve(WIMP_mass, bolo_name, analysis_type, bin_X, min_X, max_X, list_heat_fraction):

    """
    Detail:
        Obtain BDT sensitivity curve as a function of BDT cut

    Args:
        WIMP_mass (str)                       = WIMP mass
        bolo_name (str)                       = bolometer name 
        analysis_type (str)                   = type of analysis (which box cut)
        bin_X, min_X, max_X (int float float) = TH1F parameters
        list_heat_fraction (list)             = list of heat fraction

    Returns:
        void

    Raises:
        void     
    """

    sensi_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
    file_sensi = TFile(sensi_dir + bolo_name + "_WIMP_mass_" + WIMP_mass + "_BDT_sensitivity.root", "read")

    d_func= {}

    list_color = [kOrange-8, kGreen+2, kBlue-7, kRed, kBlack, kMagenta, kAzure+10, kGreen-3, kOrange-9]

    for index, heat_fraction in enumerate(list_heat_fraction):
        d_func[heat_fraction] = file_sensi.Get("sensivity_billard_mass_" + WIMP_mass + "_heat_fraction_" + heat_fraction)
        PyRPl.process_TF1(d_func[heat_fraction], color = list_color[index], line_width = 3 )

    h = TH1F("h", "h",100, 0,7)
    h.SetMaximum(1E4)    
    h.SetMinimum(1)
    PyRPl.process_TH1(h, X_title = "BDT cut", Y_title = "Sensitivity (a.u.)")

    leg = TLegend(0.36,0.52,0.71,0.89)
    l1 = TLine(0,2.3,7,2.3)
    l1.SetLineColor(kBlack)
    l1.SetLineWidth(2)
    l1.SetLineStyle(7)


    cc = TCanvas("cc", "cc")
    h.Draw()
    gPad.SetLogy()
    for heat_fraction in list_heat_fraction[::-1]:
        d_func[heat_fraction].Draw("same")
        leg.AddEntry(d_func[heat_fraction].GetName(), "<" + heat_fraction +" *max rate","l")

    leg.SetFillColor(kWhite)
    leg.SetBorderSize(0)
    leg.Draw("same")
    l1.Draw("same")

    #Some legend
    tpoisson = TText(0.3,1E-10,"truc")
    tpoisson.SetTextSize(0.05)
    tpoisson.SetTextColor(kBlack)
    tpoisson.DrawText(1,2.8,"Bckg free")

    fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/")
    cc.Print(fig_dir + bolo_name + "_mass_" + WIMP_mass + "_sensitivity_various_heat_rate_cut.png")

    raw_input()


def get_plot_BDT_sensitivity_curve_PhD(bolo_name, WIMP_mass, analysis_type, MVA_tag, bin_X, min_X, max_X, exposure, list_variables, d_bounds, **kwargs):

    """Plot BDT sensitivity nicely

    
    Detail:
        void

    Args:
        bolo_name             = (str) bolometer name
        mass                  = (int) WIMP mass
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
    model_dir = script_utils.create_directory("../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
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

            vec_proba = [TMath.PoissonI(i, exp_bckg) for i in range(2000)]    
            lim_Poisson_bckg = np.sum(np.array([PoissonCL.compute_90CL_limit(i)*vec_proba[i] for i in range(2000)]))

            if eff_sig<=0:
                return 1E10
            else:
                return lim_Poisson_bckg/eff_sig + par[0]

    func_sensi = TF1("func_sensi", Sensitivity(), d_bounds[WIMP_mass][0], d_bounds[WIMP_mass][1],1)
    func_sensi.SetParameter(0,0)
    PyRPl.process_TF1(func_sensi, color = kOrange+7, line_width = 3 )

    h = TH1F("h", "h",100, d_bounds[WIMP_mass][0], d_bounds[WIMP_mass][1])
    h.SetMaximum(1E4)    
    h.SetMinimum(0.1)
    PyRPl.process_TH1(h, X_title = "BDT cut", Y_title = "Ratio R (a.u.)")

    cc = TCanvas("cc", "cc")
    h.Draw()
    gPad.SetLogy()
    func_sensi.Draw("same")

    l1 = TLine(d_bounds[WIMP_mass][0],2.3,d_bounds[WIMP_mass][1],2.3)
    l1.SetLineColor(kBlack)
    l1.SetLineWidth(2)
    l1.SetLineStyle(7)
    l1.Draw("same")

    min_of_fun = func_sensi.GetMinimum()
    l2 = TLine(d_bounds[WIMP_mass][0],min_of_fun,d_bounds[WIMP_mass][1],min_of_fun)
    l2.SetLineColor(kRed)
    l2.SetLineWidth(2)
    l2.SetLineStyle(7)
    l2.Draw("same")

    minX_of_fun = func_sensi.GetMinimumX()
    l3 = TLine(minX_of_fun,0.1,minX_of_fun,1E4)
    l3.SetLineColor(kOrange+8)
    l3.SetLineWidth(2)
    l3.SetLineStyle(7)
    l3.Draw("same")

    #Some legend
    tpoisson = TText(0.3,1E-10,"truc")
    tpoisson.SetTextSize(0.04)
    tpoisson.SetTextColor(kBlack)
    tpoisson.DrawText(.5,1,"Bckg free sensitivity")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kRed)
    tBDT.DrawText(0.5,10,"Best BDT sensitivity")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kGray+2)
    tBDT.SetTextAngle(80)
    tBDT.DrawText(8.4,25,"Lose signal efficiency")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kGray+2)
    tBDT.DrawText(.5,270,"Background dominated")

    #Some legend
    tBDT = TText(0.3,1E-10,"truc")
    tBDT.SetTextSize(0.04)
    tBDT.SetTextColor(kOrange+8)
    tBDT.DrawText(6.3,0.71,"Best BDT cut value = " + str(minX_of_fun)[:4])

    raw_input()
    fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/")
    cc.Print(fig_dir + bolo_name + "_mass_" + str(WIMP_mass) + "_sensitivity.eps")

