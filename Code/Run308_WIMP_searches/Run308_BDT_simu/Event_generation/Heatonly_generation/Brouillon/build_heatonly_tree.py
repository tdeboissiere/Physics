#! /usr/bin/env python

from ROOT import *
from array import array
import BDT_file_handler as BDT_fh
import BDT_utils as BDT_ut
import PyROOTPlots as PyRPl
import numpy as np

def build_heatonly_tree(bolo_name, analysis_type, event_type, d_cut, num_event):

    """Build the training tree for heatonly events and the corresponding cut eff file
    
    Detail:
        Use a fitted 1D x3 exponential spectrum to generate the data
        Loop over events, select those which pass the cut and write to tree
        Also compute the efficiency of the event selection
        Write this efficiency to a .txt file

    Args:
        analysis_type = (str) name of analysis (name indicates which ion cut, which resolution...)
        event_type    = (str) class of event (S1Gamma S1Pb etc...)
        d_cut         = (dict) dictionnary which states which cut to use on heat/ion/veto
        num_event       = (int) number of events to simulate
    Returns:
        void

    Raises:
        void
    """

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"

    #Include NR computations directly
    util_path = "/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/"
    fEE_to_NR, file_EE_to_NR = PyRPl.open_ROOT_object(util_path +"conv_EE_to_NR.root", "conv")

    #Load the estimator, FWHM, for simulated events
    d_est             = BDT_fh.open_estimator_file(bolo_name)
    d_std_true_events = BDT_fh.open_true_event_FWHM_file(bolo_name,)
    # X+Y ~ N(0, sqrt(sx**2 + sy**2)) if X and Y follow iid gaussians of std dev sx and sy
    # Modify d_std_true_events to account for that since events are peaked from the spectrum of EC1
    for key in ["FWC2_ERA", "FWIA", "FWIB", "FWIC", "FWID"]:
         d_std_true_events[key] = TMath.Sqrt(d_std_true_events[key]**2 - d_std_true_events["FWC1_ERA"]**2)

    #Best estimator for heat: coefficients
    coeff_EC1         = float(d_est["HEAT"][:5])
    coeff_EC2         = 1 - coeff_EC1

    coeff_EIB         = float(d_est["FID"][:5])
    coeff_EID         = 1 - coeff_EIB

    # Sample the 1D Gamma bckg model
    Ana_path   = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"

    #Sample the 1D heat bckg model
    heatonly_path = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_heatonly_spectrum_extrapol.root"
    func_heatonly, file_heatonly = PyRPl.open_ROOT_object(heatonly_path, "heat_extra")
    # gRandom.SetSeed(0)

    #Create TFile + TTree to store the training data
    heatonly_out_path = gen_path + "BDT_" + bolo_name + "/" + analysis_type + "/Heatonly/ROOT_files/"
    file_tree = TFile(heatonly_out_path + bolo_name+"_heatonly_tree.root", "recreate")
    t_new = TTree("t_new", "t_new")

    # create 1 dimensional float arrays 
    EC1 = array("f", [0.])
    EC2 = array("f", [0.])
    EIA = array("f", [0.])
    EIB = array("f", [0.])
    EIC = array("f", [0.])
    EID = array("f", [0.])
    ENR = array("f", [0.])

    # create the branches and assign the fill-variables to them
    t_new.Branch("EC1", EC1, "EC1/F")
    t_new.Branch("EC2", EC2, "EC2/F")
    t_new.Branch("EIA", EIA, "EIA/F")
    t_new.Branch("EIB", EIB, "EIB/F")
    t_new.Branch("EIC", EIC, "EIC/F")
    t_new.Branch("EID", EID, "EID/F")
    t_new.Branch("ENR", ENR, "ENR/F")

    for i in range(num_event):

        evt_heat = func_heatonly.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))

        EIA[0] = np.random.normal(0, d_std_true_events["FWIA"])
        EIC[0] = np.random.normal(0, d_std_true_events["FWIC"])
        
        EIB[0] = np.random.normal(0, d_std_true_events["FWIB"])
        EID[0] = np.random.normal(0, d_std_true_events["FWID"])

        EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
        EC2[0]=evt_heat + np.random.normal(0, d_std_true_events["FWC2_ERA"])

        EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
        EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
        accept = 0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

        if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events) and np.random.binomial(1,0.9999*accept)==1):   
            ENR[0] = fEE_to_NR.Eval(EC)
            t_new.Fill()

    #Fill the bckg_cuteff file (the efficiency file for the given background)
    # See the BDT_file_handler code in /Desktop/Miscellaneous/Python/Useful_scripts/BDT_Analysis for details
    BDT_fh.fill_bckg_cuteff_file(bolo_name, analysis_type, t_new, "heatonly", num_event)

    t_new.Write()
    file_tree.Close()
  


