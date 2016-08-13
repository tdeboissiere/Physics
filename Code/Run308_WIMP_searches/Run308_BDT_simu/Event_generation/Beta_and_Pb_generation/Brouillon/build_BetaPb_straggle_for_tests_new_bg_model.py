#! /usr/bin/env python

from ROOT import *
from array import array
import BDT_file_handler as BDT_fh
import BDT_utils as BDT_ut
import PyROOTPlots as PyRPl
import numpy as np

def build_BetaPb_tree(bolo_name, event_type, num_event):

    """Build the training tree for Beta events and the corresponding cut eff file
    
    Detail:
        Use 1D fit to data as the background model 
        Loop over events, select those which pass the cut and write to tree
        Also compute the efficiency of the event selection
        Write this efficiency to a .txt file

    Args:
        bolo_name     = (str) bolometer name
        event_type    = (str) class of event (S1Gamma S1Pb etc...)
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

    #Include wavelet efficiency
    swave_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/" + bolo_name + "/"
    fwave_eff, file_wave = PyRPl.open_ROOT_object(swave_path + "PSA_cut_eff.root", "PSA_eff")

    #Load the estimator, FWHM, for simulated events
    d_est             = BDT_fh.open_estimator_file(bolo_name)
    d_std_true_events = BDT_fh.open_true_event_FWHM_file(bolo_name,)

    #Best estimator for heat: coefficients
    coeff_EC1         = float(d_est["HEAT"][:5])
    coeff_EC2         = 1 - coeff_EC1

    coeff_EIB         = float(d_est["FID"][:5])
    coeff_EID         = 1 - coeff_EIB



    #Get the function that give the standard deviation corresponding to the straggle
    Ana_path             = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
    strag_path = Ana_path +  "/ROOT_files/" + "FID808_straggle_std.root"
    func_stragS1Beta, file_stragS1Beta = PyRPl.open_ROOT_object(strag_path, "S1Beta_strag")
    func_stragS2Beta, file_stragS2Beta = PyRPl.open_ROOT_object(strag_path, "S2Beta_strag")
    func_stragS1Pb, file_stragS1Pb     = PyRPl.open_ROOT_object(strag_path, "S1Pb_strag")
    func_stragS2Pb, file_stragS2Pb     = PyRPl.open_ROOT_object(strag_path, "S2Pb_strag")

    # Sample the 1D Surf bckg model
    Ana_path                 = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
    func_S1Beta, file_S1Beta = PyRPl.open_ROOT_object(Ana_path+ "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Beta_from_XD_spectrum_extrapol.root", "S1Beta_extra")
    func_S2Beta, file_S2Beta = PyRPl.open_ROOT_object(Ana_path+ "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Beta_from_bolo_spectrum_extrapol.root", "S2Beta_extra")
    func_S1Pb, file_S1Pb     = PyRPl.open_ROOT_object(Ana_path+ "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Pb_from_bolo_spectrum_extrapol.root", "S1Pb_extra")
    func_S2Pb, file_S2Pb     = PyRPl.open_ROOT_object(Ana_path+ "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Pb_from_bolo_spectrum_extrapol.root", "S2Pb_extra")


    #Get the function that converts EI to EC for surface events
    Ana_path             = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
    convert_path = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_ion_heat_surface_relation.root"
    func_convS1Beta, file_convS1Beta = PyRPl.open_ROOT_object(convert_path, "S1Beta")
    func_convS2Beta, file_convS2Beta = PyRPl.open_ROOT_object(convert_path, "S2Beta")
    func_convS1Pb, file_convS1Pb     = PyRPl.open_ROOT_object(convert_path, "S1Pb")
    func_convS2Pb, file_convS2Pb     = PyRPl.open_ROOT_object(convert_path, "S2Pb")

    # gRandom.SetSeed(0)

    #Create TFile + TTree to store the training data
    out_path = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_" + event_type + "_tree_strag_test.root"
    file_tree = TFile( out_path, "recreate")
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


    #Fill tree differently (it depends on the event type)
    
    if event_type == "S1Beta":
        for i in range(num_event):

            evt_heat = func_S1Beta.GetRandom(0,30)
            strag = func_stragS1Beta.Eval(evt_heat)
            strag=0

            EIA[0]=func_convS1Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt(d_std_true_events["FWIA"]**2+strag**2))
            EIC[0]=np.random.normal(0, d_std_true_events["FWIC"])

            EIB[0]=func_convS1Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt(d_std_true_events["FWIB"]**2+strag**2))
            EID[0]=np.random.normal(0,d_std_true_events["FWID"])

            EC1[0]=evt_heat + np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, d_std_true_events["FWC2_ERA"])

            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.5) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (np.random.binomial(1,0.999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                t_new.Fill()

    if event_type == "S2Beta":
        for i in range(num_event):

            evt_heat = func_S2Beta.GetRandom(1,30)
            strag = func_stragS2Beta.Eval(evt_heat)
            strag=0

            EIA[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIA"])**2 - d_std_true_events["FWC1_ERA"]**2))
            EIC[0]=func_convS2Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIC"]+strag)**2 - d_std_true_events["FWC1_ERA"]**2))

            EIB[0]=np.random.normal(0,TMath.Sqrt((d_std_true_events["FWIB"])**2 - d_std_true_events["FWC1_ERA"]**2))
            EID[0]=func_convS2Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWID"]+strag)**2 - d_std_true_events["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWC2_ERA"])**2 - d_std_true_events["FWC1_ERA"]**2))
            
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.5) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (np.random.binomial(1,0.999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                t_new.Fill()


    if event_type == "S1Pb":
        for i in range(num_event):

            evt_heat = func_S1Pb.GetRandom(1,30)
            strag = func_stragS1Pb.Eval(evt_heat)
            strag=0

            EIA[0]=func_convS1Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIA"]+strag)**2 - d_std_true_events["FWC1_ERA"]**2))
            EIC[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIC"])**2 - d_std_true_events["FWC1_ERA"]**2))

            EIB[0]=func_convS1Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIB"]+strag)**2 - d_std_true_events["FWC1_ERA"]**2))
            EID[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events["FWID"])**2 - d_std_true_events["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWC2_ERA"])**2 - d_std_true_events["FWC1_ERA"]**2))

            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.5) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )
            
            if (np.random.binomial(1,0.999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                t_new.Fill()

    if event_type == "S2Pb":
        for i in range(num_event):

            evt_heat = func_S2Pb.GetRandom(1,30)
            strag = func_stragS2Pb.Eval(evt_heat)
            strag=0

            EIA[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIA"])**2 - d_std_true_events["FWC1_ERA"]**2))
            EIC[0]=func_convS2Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWIC"]+strag)**2 - d_std_true_events["FWC1_ERA"]**2))

            EIB[0]=np.random.normal(0,TMath.Sqrt((d_std_true_events["FWIB"])**2 - d_std_true_events["FWC1_ERA"]**2))
            EID[0]=func_convS2Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWID"]+strag)**2 - d_std_true_events["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events["FWC2_ERA"])**2 - d_std_true_events["FWC1_ERA"]**2))
            
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.5) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (np.random.binomial(1,0.999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                t_new.Fill()



    t_new.Write()
    file_tree.Close()


bolo_name = "FID837"
list_type = ["S1Beta", "S2Beta", "S1Pb", "S2Pb"]
list_numevt = [438, 587, 64, 93]
for index, event_type in enumerate(list_type):
    build_BetaPb_tree(bolo_name, event_type, list_numevt[index])
