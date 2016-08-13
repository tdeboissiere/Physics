#! /usr/bin/env python

from ROOT import *
from array import array
import BDT_file_handler as BDT_fh
import BDT_utils as BDT_ut
import PyROOTPlots as PyRPl
import numpy as np
import script_utils as script_utils

def build_true_events_tree(bolo_name, analysis_type, d_cut, nsimu, exposure):

    """Build the training tree for simulated data events 
    
    Detail:
        Reuse the codes written for bckg simulation. 
        Scale the data so as to respect the observed distribution
        Save to a tree
        Repeat nsimu times

    Args:
        bolo_name     = (str) bolometer name
        analysis_type = (str) name of analysis (name indicates which ion cut, which resolution...)
        d_cut         = (dict) dictionnary which states which cut to use on heat/ion/veto
        nsimu         = (int) the number of independent bckg realisations to simulate
        exposure      = (float) desired exposure

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


    #Load the estimator, FWHM, for simulated events. Different d_std for surfae to account for straggle
    d_est                   = BDT_fh.open_estimator_file(bolo_name)
    d_std_true_events       = BDT_fh.open_true_event_FWHM_file(bolo_name,)
    d_std_true_events_strag = BDT_fh.open_true_event_FWHM_file(bolo_name,)
    # X+Y ~ N(0, sqrt(sx**2 + sy**2)) if X and Y follow iid gaussians of std dev sx and sy
    # Modify d_std_true_events to account for that since events are peaked from the spectrum of EC1
    for key in ["FWC2_ERA", "FWIA", "FWIB", "FWIC", "FWID"]:
         d_std_true_events[key] = TMath.Sqrt(d_std_true_events[key]**2 - d_std_true_events["FWC1_ERA"]**2)

    #Best estimator for heat: coefficients
    coeff_EC1         = float(d_est["HEAT"][:5])
    coeff_EC2         = 1 - coeff_EC1

    coeff_EIB         = float(d_est["FID"][:5])
    coeff_EID         = 1 - coeff_EIB

    d_scaling = BDT_fh.open_scaling_file(bolo_name, d_cut)

    nFidGamma, nS1Gamma, nS2Gamma, nS1Beta, nS2Beta, nS1Pb, nS2Pb, nheatonly = BDT_ut.get_bckg_num_of_event(d_scaling, exposure)
    print nFidGamma, nS1Gamma, nS2Gamma, nS1Beta, nS2Beta, nS1Pb, nS2Pb, nheatonly
    raw_input()
    # #Get the function that give the standard deviation corresponding to the straggle
    # Ana_path             = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
    # strag_path = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_straggle_std.root"
    # func_stragS1Beta, file_stragS1Beta = PyRPl.open_ROOT_object(strag_path, "S1Beta_strag")
    # func_stragS2Beta, file_stragS2Beta = PyRPl.open_ROOT_object(strag_path, "S2Beta_strag")
    # func_stragS1Pb, file_stragS1Pb     = PyRPl.open_ROOT_object(strag_path, "S1Pb_strag")
    # func_stragS2Pb, file_stragS2Pb     = PyRPl.open_ROOT_object(strag_path, "S2Pb_strag")

    #Get the function that converts EI to EC for surface events
    Ana_path             = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
    convert_path = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_ion_heat_surface_relation.root"
    func_convS1Beta, file_convS1Beta = PyRPl.open_ROOT_object(convert_path, "S1Beta")
    func_convS2Beta, file_convS2Beta = PyRPl.open_ROOT_object(convert_path, "S2Beta")
    func_convS1Pb, file_convS1Pb     = PyRPl.open_ROOT_object(convert_path, "S1Pb")
    func_convS2Pb, file_convS2Pb     = PyRPl.open_ROOT_object(convert_path, "S2Pb")

    # Sample the 1D Gamma bckg model
    Ana_path                     = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name + "/ROOT_files/"
    func_FidGamma, file_FidGamma = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Gamma_spectrum_extrapol.root", "FidGamma_extra")
    func_S1Gamma, file_S1Gamma   = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Gamma_spectrum_extrapol.root", "S1Gamma_extra")
    func_S2Gamma, file_S2Gamma   = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Gamma_spectrum_extrapol.root", "S2Gamma_extra")
    func_S1Beta, file_S1Beta     = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Beta_from_bolo_spectrum_extrapol.root", "S1Beta_extra")
    func_S2Beta, file_S2Beta     = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Beta_from_bolo_spectrum_extrapol.root", "S2Beta_extra")
    func_S1Pb, file_S1Pb         = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Pb_from_bolo_spectrum_extrapol.root", "S1Pb_extra")
    func_S2Pb, file_S2Pb         = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_Pb_from_bolo_spectrum_extrapol.root", "S2Pb_extra")
    func_heatonly, file_heatonly = PyRPl.open_ROOT_object(Ana_path + bolo_name +  "_heatonly_spectrum_extrapol.root", "heat_extra")

    #Out path
    path_name = gen_path + "BDT_" + bolo_name + "/" + analysis_type + "/True_events/ROOT_files/"
    file_true=TFile(path_name + bolo_name + "_true_events_tree.root", "recreate")

    #Loop over number of times we want to simulate the experiment
    for i in range(nsimu):

        script_utils.print_utility("Starting simulation # " + str(i))
        cFidGamma, cS1Gamma, cS2Gamma, cS1Beta, cS2Beta, cS1Pb, cS2Pb, cheatonly =0,0,0,0,0,0,0,0

        # create 1 dimensional float arrays 
        EC1 = array("f", [0.])
        EC2 = array("f", [0.])
        EIA = array("f", [0.])
        EIB = array("f", [0.])
        EIC = array("f", [0.])
        EID = array("f", [0.])
        ENR = array("f", [0.])

        #New tree to be filled
        t_new = TTree("t_new" + str(i), "t_new" + str(i))

        # create the branches and assign the fill-variables to them
        t_new.Branch("EC1", EC1, "EC1/F")
        t_new.Branch("EC2", EC2, "EC2/F")
        t_new.Branch("EIA", EIA, "EIA/F")
        t_new.Branch("EIB", EIB, "EIB/F")
        t_new.Branch("EIC", EIC, "EIC/F")
        t_new.Branch("EID", EID, "EID/F")
        t_new.Branch("ENR", ENR, "ENR/F")

        #########################
        # Heatonly events
        #########################
        for k in range(nheatonly):

            evt_heat = func_heatonly.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))
            
            EIA[0]   = np.random.normal(0, d_std_true_events["FWIA"])
            EIC[0]   = np.random.normal(0, d_std_true_events["FWIC"])
            
            EIB[0]   = np.random.normal(0, d_std_true_events["FWIB"])
            EID[0]   = np.random.normal(0, d_std_true_events["FWID"])
            
            EC1[0]   = evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]   = evt_heat + np.random.normal(0, d_std_true_events["FWC2_ERA"])
            
            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            #No fwave eff, already taken into account in the fit
            accept = 0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events) and np.random.binomial(1,0.9999*accept)==1):   
                ENR[0] = fEE_to_NR.Eval(EC)
                cheatonly+=1
                t_new.Fill()

        #########################
        # Gamma events
        #########################
        Luke_factor = (1+8./3)/(1+5.5/3)

        for k in range(nFidGamma):

            evt_heat = func_FidGamma.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))

            EIA[0]=np.random.normal(0, d_std_true_events["FWIA"])
            EIC[0]=np.random.normal(0, d_std_true_events["FWIC"])

            EIB[0]=evt_heat + np.random.normal(0, d_std_true_events["FWIB"])
            EID[0]=evt_heat + np.random.normal(0, d_std_true_events["FWID"])

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, d_std_true_events["FWC2_ERA"])

            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events) and np.random.binomial(1,0.9999*accept)==1):   
                ENR[0] = fEE_to_NR.Eval(EC)
                cFidGamma+=1
                t_new.Fill()

        
        for k in range(nS1Gamma):

            evt_heat = func_S1Gamma.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))

            EIA[0]=Luke_factor*evt_heat + np.random.normal(0, d_std_true_events["FWIA"])
            EIC[0]=np.random.normal(0, d_std_true_events["FWIC"])

            EIB[0]=Luke_factor*evt_heat + np.random.normal(0, d_std_true_events["FWIB"])
            EID[0]=np.random.normal(0, d_std_true_events["FWID"])

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, d_std_true_events["FWC2_ERA"])

            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events) and np.random.binomial(1,0.9999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                cS1Gamma+=1
                t_new.Fill()

        for k in range(nS2Gamma):

            evt_heat = func_S2Gamma.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))

            EIA[0]=np.random.normal(0, d_std_true_events["FWIA"])
            EIC[0]=Luke_factor*evt_heat + np.random.normal(0, d_std_true_events["FWIC"])

            EIB[0]=np.random.normal(0, d_std_true_events["FWIB"])
            EID[0]=Luke_factor*evt_heat + np.random.normal(0, d_std_true_events["FWID"])

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, d_std_true_events["FWC2_ERA"])
            
            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events) and np.random.binomial(1,0.9999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                cS2Gamma+=1
                t_new.Fill()


        #####################
        # Beta
        #####################

        for k in range(nS1Beta):

            evt_heat = func_S1Beta.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))
            # strag = func_stragS1Beta.Eval(evt_heat)
            strag = 0.1*evt_heat

            EIA[0]=func_convS1Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIA"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EIC[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIC"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EIB[0]=func_convS1Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIB"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EID[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWID"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events_strag["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWC2_ERA"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events_strag["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events_strag) and np.random.binomial(1,0.9999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                cS1Beta+=1
                t_new.Fill()

        for k in range(nS2Beta):

            evt_heat = func_S2Beta.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))
            # strag = func_stragS2Beta.Eval(evt_heat)
            strag = 0.1*evt_heat

            EIA[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIA"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EIC[0]=func_convS2Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIC"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EIB[0]=np.random.normal(0,TMath.Sqrt((d_std_true_events_strag["FWIB"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EID[0]=func_convS2Beta.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWID"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events_strag["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWC2_ERA"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))
                
            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events_strag["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events_strag) and np.random.binomial(1,0.9999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                cS2Beta+=1
                t_new.Fill()

        #####################
        # Pb
        #####################

        for k in range(nS1Pb):

            evt_heat = func_S1Pb.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))
            # strag = func_stragS1Pb.Eval(evt_heat)
            strag = 0.05*evt_heat

            EIA[0]=func_convS1Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIA"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EIC[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIC"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EIB[0]=func_convS1Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIB"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EID[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWID"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events_strag["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWC2_ERA"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events_strag["FWHEAT"]) ) )
                
            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events_strag) and np.random.binomial(1,0.9999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                cS1Pb+=1
                t_new.Fill()

        for k in range(nS2Pb):

            evt_heat = func_S2Pb.GetRandom(float(d_cut["ECinf"]), float(d_cut["ECsup"]))
            # strag = func_stragS2Pb.Eval(evt_heat)
            strag = 0.05*evt_heat

            EIA[0]=np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIA"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EIC[0]=func_convS2Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWIC"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EIB[0]=np.random.normal(0,TMath.Sqrt((d_std_true_events_strag["FWIB"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))
            EID[0]=func_convS2Pb.Eval(evt_heat) + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWID"]+strag)**2 - d_std_true_events_strag["FWC1_ERA"]**2))

            EC1[0]=evt_heat #+ np.random.normal(0, d_std_true_events_strag["FWC1_ERA"])
            EC2[0]=evt_heat + np.random.normal(0, TMath.Sqrt((d_std_true_events_strag["FWC2_ERA"])**2 - d_std_true_events_strag["FWC1_ERA"]**2))
                
            EI = coeff_EIB*EIB[0] + coeff_EID*EID[0]
            EC = coeff_EC1*EC1[0] + coeff_EC2*EC2[0]
            accept = (fwave_eff.Eval(EC))*0.5*(1 + TMath.Erf( (EC-0.65) / (TMath.Sqrt(2)*d_std_true_events_strag["FWHEAT"]) ) )

            if (BDT_ut.bool_pass_cut(EC, EI, EIA[0], EIC[0], d_cut, d_std_true_events_strag) and np.random.binomial(1,0.9999*accept)==1):  
                ENR[0] = fEE_to_NR.Eval(EC)
                cS2Pb+=1
                t_new.Fill()

        t_new.Write()

        list_counts =[cFidGamma, cS1Gamma, cS2Gamma, cS1Beta, cS2Beta, cS1Pb, cS2Pb, cheatonly]
        print "total_counts", sum(list_counts)
        sum_counts = float(sum(list_counts))
        print cheatonly, cFidGamma, cS1Gamma, cS2Gamma, cS1Beta, cS2Beta, cS1Pb, cS2Pb
        print cheatonly/sum_counts, cFidGamma/sum_counts, cS1Gamma/sum_counts, cS2Gamma/sum_counts, cS1Beta/sum_counts, cS2Beta/sum_counts, cS1Pb/sum_counts, cS2Pb/sum_counts 

    file_true.Close()




        
