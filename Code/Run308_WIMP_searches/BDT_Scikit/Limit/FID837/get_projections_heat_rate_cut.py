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

from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Poisson_90CL as Poisson90

def get_limit_graph(file_name, line_width, color):

    """Open the file to get the graph of the limit 
    for the given file name
    
    Detail:
        Return the limit graph of given file_name 
        Also modify line_width and color

    Args:
        file_name (str)    = the name of the limit file 
        line_width (int)   = the line width 
        color (ROOT color) = the color for the line
        
    Returns:
        void

    Raises:
        void
    """

    if "cdms" in file_name:
        arr = np.loadtxt(file_name, delimiter = ",") 
        n_points                  =int(arr[:,0].shape[0])
        gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-10*np.array(list(arr[:,1])))
        gr.SetLineColor(color)
        gr.SetLineWidth(line_width)
    else:
        arr = np.loadtxt(file_name, delimiter = ",") 
        n_points                  =int(arr[:,0].shape[0])
        gr = TGraph(n_points, np.array(list(arr[:,0])), np.array(list(arr[:,1])))
        gr.SetLineColor(color)
        gr.SetLineWidth(line_width)

    return gr

def get_point_limit(bolo_name, nobs, WIMP_mass,analysis_type, heat_fraction, exposure, detector_mass = 0.6):

    """See detail
    
    Detail:
        From a single observation, get the 90 CL limit

    Args:
        bolo_name     = (str) bolometer name
        arr           = (np array) array of 90 CL values on the number of WIMPs
        WIMP_mass     = (str) WIMP mass in GeV
        analysis      = (str) type of analysis cut
        heat_fraction = (str) fraction of heat only events kept
        exposure      = (float) exposure in days
        detector_mass = (float) fiducial mass
        
    Returns:
        d (dict) = the percentile dictionnary

    Raises:
        void
    """

    #Convert the obs counts in arr to a 90CL Poisson limit
    nlim = Poisson90.compute_90CL_limit(nobs)

    #Get the efficiency 
    root_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
    file_root = TFile(root_dir + bolo_name + "_sensi_eff_curves_heat_fraction" + heat_fraction + "_mass_" + str(WIMP_mass) + ".root", "read")

    #Get the efficiency
    fsensi = file_root.Get("sensitivity_expo_" + str(exposure))
    feff = file_root.Get("signal_eff_expo_" + str(exposure))
    cut_val = fsensi.GetMinimumX(2,10)
    efficiency = feff.Eval(cut_val)

    #For normalisation
    count_arr = np.loadtxt("./Text_files/WIMP_counts_for_1kgday_" + analysis_type + "_2D.txt", delimiter = ",")
    d_conv={}
    for mass, Nevents in count_arr:
        d_conv[str(int(mass))] = float(Nevents)

    return nlim*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass* efficiency )


def get_simulated_event_limit(bolo_name, list_mass, analysis_type, heat_fraction, exposure, detector_mass = 0.6):

    """Get expected sensitivity confidence bands
    
    Detail:
        Return filled area graph for the 68 and 95 CI on the 
        expected sensitivity given a choice of analysis

    Args:
        bolo_name     = (str) bolometer name
        arr           = (np array) array of 90 CL values on the number of WIMPs
        WIMP_mass     = (str) WIMP mass in GeV
        analysis      = (str) type of analysis cut
        heat_fraction = (str) fraction of heat only events kept
        exposure      = (float) exposure in days
        detector_mass = (float) fiducial mass

    Returns:
        gr_lim = average limit

    Raises:
        void
    """

    d_mass_limit={}

    txt_dir = script_utils.create_directory("./Text_files/Simulated_sensitivity/")

    for WIMP_mass in list_mass:
        fin = open(txt_dir + "/simulated_events_passing_cut_heat_fraction" + heat_fraction + "_mass_" + str(WIMP_mass) + ".txt", "r")
        stuff = fin.readlines()
        list_lines = [map(float,line.rstrip().split(",")) for line in stuff[1:]]
        list_event_pass_cut = []
        for line in list_lines:
            if line[1] == exposure:
                list_event_pass_cut = line[2:]
        d_mass_limit[int(WIMP_mass)] = np.mean(np.array([get_point_limit(bolo_name, nobs, WIMP_mass, analysis_type, heat_fraction, exposure, detector_mass) for nobs in list_event_pass_cut]))

    list_lim = [d_mass_limit[int(mass)] for mass in list_mass]

    gr_lim = TGraph(len(list_mass), np.array(list_mass).astype(float), np.array(list_lim).astype(float))

    return gr_lim


def plot_limit(bolo_name, list_mass, analysis_type, exposure, detector_mass = 0.6):

    """Plot all the desired limits
    
    Detail:
        Modify the code to add more limits if needed

    Args:
        bolo_name = (str) bolometer name
        list_mass = (list) list of WIMP masses
        analysis_type      = (str) type of analysis cut
        exposure      = (float) exposure in days
        detector_mass = (float) fiducial mass
        
    Returns:
        void
    Raises:
        void
    """

    d_graph = {}
    list_color = [kOrange-8, kGreen+2, kBlue-7, kRed, kBlack, kMagenta, kAzure+10, kGreen-3, kOrange-9]

    for index, heat_fraction in enumerate(["0.3","0.4","0.5","0.8","1"]):
        d_graph[heat_fraction] = get_simulated_event_limit(bolo_name, list_mass, analysis_type, "_" + heat_fraction, exposure, detector_mass = 0.6)
        d_graph[heat_fraction].SetName(heat_fraction)
        PyRPl.process_TGraph(d_graph[heat_fraction], color = list_color[index])

    gr_edw_poisson = get_limit_graph("./Text_files/edw3_ana_1.5_0_5_poisson.txt", 2, kBlack)
    gr_edw_low = get_limit_graph("./Text_files/Published_limits/edw_lowmass_2012.txt", 2, kRed)
    gr_edw_low.SetLineStyle(7)
    gr_cdms = get_limit_graph("./Text_files/Published_limits/cdms_limit.txt", 2, kBlue)

    h = TH1F("h", "", 100, 3,25)
    PyRPl.process_TH1(h, X_title = "Mass (GeV)", Y_title = "#sigma (pb)", X_title_size = .06, Y_title_size = .06, X_title_offset = .98, Y_title_offset = .95)


    gr_edw_low.SetName("gr_edw_low")
    gr_edw_poisson.SetName("gr_edw_poisson")
    gr_cdms.SetName("gr_cdms")

    cc = TCanvas("cc", "cc")
    gPad.SetLogy()
    gPad.SetLogx()
    h.SetMaximum(1E-1)
    h.SetMinimum(4E-8)
    h.Draw()

    gr_cdms.Draw("sameC")
    gr_edw_poisson.Draw("sameC")
    gr_edw_low.Draw("sameC")

    for index, heat_fraction in enumerate(["0.3","0.4","0.5","0.8","1"]):
        d_graph[heat_fraction].Draw("sameC")

    leg =TLegend(0.564,0.584,0.83,0.857)
    leg.AddEntry("gr_cdms", "SCDMS" , "l")
    leg.AddEntry("gr_edw_low", "EDW II" , "l")
    leg.AddEntry("gr_edw_poisson", "EDW III Poisson" , "l")
    for index, heat_fraction in enumerate(["0.3","0.4","0.5","0.8","1"]):
        leg.AddEntry( d_graph[heat_fraction].GetName(), heat_fraction , "l")

    leg.SetFillColor(kWhite)
    leg.SetLineColor(kWhite)
    leg.Draw()
    raw_input()

def get_sensi_eff_curves_various_exp(bolo_name, WIMP_mass, d_cut, analysis_type, MVA_tag, bin_X, min_X, max_X, list_variables, **kwargs):

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
        list_variables (list)      = list of variables to retain for BDT

    Returns:
        void

    Raises:
        void
    """       

    try:
        kwargs["weight_dir"]
    except KeyError:
        sys.exit()

    #Get heat _fraction
    heat_fraction = kwargs["classifier_name"][13:]

    #Get scaling dict to set the weights
    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)
    # print d_scaling

    d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                            "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", 
                            "heatonly_heat_fraction" + heat_fraction: "Heatonly", "WIMP_mass_" + str(WIMP_mass): "WIMP"}
    key_heat = "heatonly_heat_fraction" + heat_fraction

    #Load data
    d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), 1, list_variables, datasplit = 1)

    # Get classifier
    model_dir = script_utils.create_directory("../../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")    
    if kwargs.has_key("classifier_name"):
        modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
    bst = xgb.Booster({'nthread':16}, model_file = modelfile)

    #Get predictions on test sample
    d_pred = {}
    d_hist = {}
    d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
                     "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, key_heat: kRed, "WIMP_mass_" + str(WIMP_mass):kGray, "neutron":kMagenta}

    #ROOT out_dir 
    root_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
    file_root = TFile(root_dir + bolo_name + "_sensi_eff_curves_heat_fraction" + heat_fraction + "_mass_" + str(WIMP_mass) + ".root", "recreate")

    #Loop over possible exposure values
    # for exposure in [10, 50, 100, 500]:
    for exposure in [66]:
        script_utils.print_utility("Getting sensi + eff for exposure of " + str(exposure) + " mass of " + str(WIMP_mass))
        for event_type in d_test.keys():
            d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values) )
            d_hist[event_type] = TH1F("h" + event_type + str(exposure), "h" + event_type + str(exposure), bin_X, min_X, max_X)
            PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
            PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
            if "WIMP" not in event_type:
                d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
            else:
                d_hist["WIMP_mass_" + str(WIMP_mass)].Scale(8000./d_hist["WIMP_mass_" + str(WIMP_mass)].Integral())

        list_hist_bckg =[d_hist["S1Pb"], d_hist["S2Pb"], d_hist["S1Beta"], d_hist["S2Beta"], d_hist["S1Gamma"], d_hist["S2Gamma"], d_hist["FidGamma"], d_hist[key_heat]]

        hsum_bckg=TH1F("hsum_bckg" + str(exposure),"hsum_bckg" + str(exposure), bin_X, min_X, max_X)
        for i in range(1,bin_X+1):
            sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
            hsum_bckg.SetBinContent(i, sumcontent)

        # print hsum_bckg.Integral(hsum_bckg.FindBin(3.5), bin_X)
        # print d_hist["WIMP_mass_" + str(WIMP_mass)].Integral(d_hist["WIMP_mass_" + str(WIMP_mass)].FindBin(3.5), bin_X)/d_hist["WIMP_mass_" + str(WIMP_mass)].Integral()

        hs=THStack("hs", "hs")
        for hist in list_hist_bckg + [d_hist["WIMP_mass_" + str(WIMP_mass)]]:
            hs.Add(hist)

        # cc = TCanvas("cc", "cc")
        # h1=TH1F("h1","h1", bin_X, min_X, max_X)
        # PyRPl.process_TH1(h1, X_title="BDT ouput", min_Y = 1E-1, max_Y = 20000)
        
        # gPad.SetLogy()
        # h1.Draw()
        # hs.Draw("same")
        # raw_input()

        class Sensitivity:
           def __call__( self, x, par ):

                bin_number_sig = d_hist["WIMP_mass_" + str(WIMP_mass)].FindBin(x[0])
                bin_number_bckg = hsum_bckg.FindBin(x[0])
                eff_sig = float(d_hist["WIMP_mass_" + str(WIMP_mass)].Integral(bin_number_sig, bin_X))
                exp_bckg = hsum_bckg.Integral(bin_number_bckg, bin_X)

                vec_proba = [TMath.PoissonI(i, exp_bckg) for i in range(500)]    
                lim_Poisson_bckg = np.sum(np.array([PoissonCL.compute_90CL_limit(i)*vec_proba[i] for i in range(500)]))

                if eff_sig<=0:
                    return 1E10
                else:
                    return lim_Poisson_bckg/eff_sig + par[0]

        class Signal_eff:
           def __call__( self, x, par ):

                bin_number = d_hist["WIMP_mass_" + str(WIMP_mass)].FindBin(x[0])
                integ = float(d_hist["WIMP_mass_" + str(WIMP_mass)].Integral(bin_number, bin_X))/float(d_hist["WIMP_mass_" + str(WIMP_mass)].Integral())
                return par[0] + integ

        h = TH1F("h", "h",100, 0, 10)
        PyRPl.process_TH1(h, X_title = "BDT cut", Y_title = "Sensitivity (a.u.)")
        h.SetMinimum(1)
        h.SetMaximum(1E3)
        # h.Draw()

        fopt = TF1("sensitivity_expo_" + str(exposure), Sensitivity(), 0,10, 1)
        fopt.SetParameter(0,0)
        fopt.SetNpx(100)
        # fopt.Draw("same")

        fsig_eff = TF1("signal_eff_expo_" + str(exposure), Signal_eff(), 0,10, 1)
        fsig_eff.SetParameter(0,0)
        fsig_eff.SetNpx(500)

        min_X = fopt.GetMinimumX(2,10)
        print "signal eff", fsig_eff.Eval(min_X)
        print "bckg_exp", hsum_bckg.Integral(hsum_bckg.FindBin(min_X), bin_X)

        # fopt.Write()
        # fsig_eff.Write()

        # gPad.SetLogy()
        # raw_input()
        # del h 

    # file_root.Close()

def get_events_passing_cuts(bolo_name, WIMP_mass, d_cut, analysis_type, MVA_tag, bin_X, min_X, max_X, list_variables, **kwargs):

    """Simulate data and find the number of events that pass
    the BDT cut

    
    Detail:
        void

    Args:
        bolo_name           = (str) bolometer name
        mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
        analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
        MVA_tag (str)       = indicates which scaling file to use
        bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
        list_variables (list)      = list of variables to retain for BDT

    Returns:
        void

    Raises:
        void
    """       

    try:
        kwargs["weight_dir"]
    except KeyError:
        sys.exit()

    #Get heat _fraction
    heat_fraction = kwargs["classifier_name"][13:]

    #Get scaling dict to set the weights
    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

    d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                            "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", 
                            "heatonly_heat_fraction" + heat_fraction: "Heatonly", "WIMP_mass_" + str(WIMP_mass): "WIMP"}
    key_heat = "heatonly_heat_fraction" + heat_fraction

    #Load data
    d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), 1, list_variables, datasplit = 1)

    # Get classifier
    model_dir = script_utils.create_directory("../../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")    
    if kwargs.has_key("classifier_name"):
        modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
    bst = xgb.Booster({'nthread':16}, model_file = modelfile)

    #Get predictions on test sample
    d_pred = {}
    d_hist = {}
    d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
                     "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, key_heat: kRed, "WIMP_mass_" + str(WIMP_mass):kGray, "neutron":kMagenta}

    #ROOT out_dir 
    root_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
    file_root = TFile(root_dir + bolo_name + "_sensi_eff_curves_heat_fraction" + heat_fraction + "_mass_" + str(WIMP_mass) + ".root", "read")

    #Write events that pass cut to a file 
    txt_dir = script_utils.create_directory("./Text_files/Simulated_sensitivity/")
    with open(txt_dir + "/simulated_events_passing_cut_heat_fraction" + heat_fraction + "_mass_" + str(WIMP_mass) + ".txt", "w") as fout:

        fout.write("heat_fraction,exposure,num_events_passing_cut\n")

        #Loop over possible exposure values
        for exposure in [10, 50, 100, 500]:
            script_utils.print_utility("Getting events passing cut for exposure of " + str(exposure) + " mass of " + str(WIMP_mass))
            for event_type in d_test.keys():
                d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values) )
                d_hist[event_type] = TH1F("h" + event_type + str(exposure), "h" + event_type + str(exposure), bin_X, min_X, max_X)
                PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
                PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
                if "WIMP" not in event_type:
                    d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
                else:
                    d_hist["WIMP_mass_" + str(WIMP_mass)].Scale(1./d_hist["WIMP_mass_" + str(WIMP_mass)].Integral())

            list_hist_bckg =[d_hist["S1Pb"], d_hist["S2Pb"], d_hist["S1Beta"], d_hist["S2Beta"], d_hist["S1Gamma"], d_hist["S2Gamma"], d_hist["FidGamma"], d_hist[key_heat]]

            hsum_bckg=TH1F("hsum_bckg" + str(exposure),"hsum_bckg" + str(exposure), bin_X, min_X, max_X)
            for i in range(1,bin_X+1):
                sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
                hsum_bckg.SetBinContent(i, sumcontent)

            fsensi = file_root.Get("sensitivity_expo_" + str(exposure))
            cut_val = fsensi.GetMinimumX(2,10)

            #Run Poisson simulations
            list_event_pass_cut=[]
            for nsimu in range(100):
                hdatasimu = TH1F("hdatasimu","hdatasimu", bin_X, min_X, max_X)
                for i in range(1,bin_X+1):
                    hdatasimu.SetBinContent(i, np.random.poisson(hsum_bckg.GetBinContent(i)))
                bin_cut = hdatasimu.FindBin(cut_val)
                num_entry_cut = int(hdatasimu.Integral(bin_cut, max_X))
                list_event_pass_cut.append(str(num_entry_cut))
                del hdatasimu
            fout.write(heat_fraction[1:] + "," + str(exposure) + "," + ",".join(list_event_pass_cut) + "\n")


bolo_name           = "FID837"
list_mass           = [5,6,7,10,25]
analysis_type       = "ana_1.5_0_5"
bin_X, min_X, max_X = 2000, -20, 20
exposure            = 66
MVA_tag             = ""
list_variables = ["EIA", "EIB", "EID", "EIC","EC1", "EC2",  "test"]

d_cut               = {"ECinf": 1.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}

# for WIMP_mass in list_mass:
for WIMP_mass in [5]:
    # for heat_fraction in ["0.3","0.4","0.5","0.8","1"]:
    for heat_fraction in ["0.2","0.3", "1"]:
        MVA_tag = "heat_fraction_" + heat_fraction
        get_sensi_eff_curves_various_exp(bolo_name, WIMP_mass, d_cut, analysis_type, MVA_tag, bin_X, min_X, max_X, list_variables,classifier_name=MVA_tag, weight_dir = "With_weights")
        # get_events_passing_cuts(bolo_name, WIMP_mass, d_cut, analysis_type, MVA_tag, bin_X, min_X, max_X, list_variables, classifier_name=MVA_tag, weight_dir = "With_weights")

# for exposure in [10]:
#     plot_limit(bolo_name, [5,6,7,10,25], analysis_type, exposure, detector_mass = 0.6)