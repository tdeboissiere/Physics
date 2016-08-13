#! /usr/bin/env python
from ROOT import *
import BDT_file_handler as BDT_fh
import numpy as np
import sys 
sys.path.append("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/Event_generation/Heatonly_generation/")
import heatonly_KDE as heatonly_KDE

def adjust_MVA_scaling_file(bolo_name, d_cut, analysis_type,  exposure):

    """Change the value of the FWHM/ energycuts used in theget_MVA_scaling.C file
    
    Detail:
        Read the file once to save its content. 
        Then rewrite it with the appropriate modifications

    Args:
        bolo_name     = (str) bolometer name
        d_cut         = (dict) dictionnary which states which cut to use on heat/ion/veto
        analysis_type = (str) name of analysis (which FWHM, which cuts)
        exposure      = (int) exposure

    Returns:
        void

    Raises:
        void
    """

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"

    #Load the FWHM, scaling and estimator to update the EI and EC expression and also cut efficiency
    d_std_true_events = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    d_estimator = BDT_fh.open_estimator_file(bolo_name,"")
    d_scaling = BDT_fh.open_scaling_file(bolo_name, d_cut,"")
    d_cut_eff = BDT_fh.open_bckg_cuteff_file(bolo_name, analysis_type,"")


    nFidGamma = np.random.poisson(float(d_scaling["rate_FidGamma"])*exposure/d_cut_eff["FidGamma"])
    nS1Gamma  = np.random.poisson(float(d_scaling["rate_S1Gamma"])*exposure/d_cut_eff["S1Gamma"])
    nS2Gamma  = np.random.poisson(float(d_scaling["rate_S2Gamma"])*exposure/d_cut_eff["S2Gamma"])
    nS1Beta   = np.random.poisson(float(d_scaling["rate_S1Beta"])*exposure/d_cut_eff["S1Beta"])
    nS2Beta   = np.random.poisson(float(d_scaling["rate_S2Beta"])*exposure/d_cut_eff["S2Beta"])
    nS1Pb     = np.random.poisson(float(d_scaling["rate_S1Pb"])*exposure/d_cut_eff["S1Pb"])
    nS2Pb     = np.random.poisson(float(d_scaling["rate_S2Pb"])*exposure/d_cut_eff["S2Pb"])
    nheatonly = np.random.poisson(float(d_scaling["rate_heatonly"])*exposure/d_cut_eff["heatonly"])

    #Generate heatonly data
    # heatonly_KDE.generate_heatonly_from_KDE_scaling(bolo_name, analysis_type, nheatonly)
    # print "Heatonly KDE simulation complete"


    list_cpp_lines =[]

    with open(gen_path + "Event_generation/True_events_generation/get_MVA_scaling.C", "r") as f_cpp:
        list_cpp_lines = f_cpp.readlines()
   
    with open(gen_path + "Event_generation/True_events_generation/get_MVA_scaling.C", "w") as f_cpp:
        for line in list_cpp_lines:


            if "float s_heat=" in line:
                f_cpp.write("\tfloat s_heat=" + str(d_std_true_events["FWHEAT"]) + ";\n")
            elif "float s_EC1=" in line:
                f_cpp.write("\tfloat s_EC1=" + str(d_std_true_events["OWC1"]) + ";\n")
            elif "float s_EC2=" in line:
                f_cpp.write("\tfloat s_EC2=" + str(d_std_true_events["OWC2"]) + ";\n")

            elif "float s_IA=" in line:
                f_cpp.write("\tfloat s_IA=" + str(d_std_true_events["FWIA"]) + ";\n")
            elif "float s_IB=" in line:
                f_cpp.write("\tfloat s_IB=" + str(d_std_true_events["FWIB"]) + ";\n")
            elif "float s_IC=" in line:
                f_cpp.write("\tfloat s_IC=" + str(d_std_true_events["FWIC"]) + ";\n")
            elif "float s_ID=" in line:
                f_cpp.write("\tfloat s_ID=" + str(d_std_true_events["FWID"]) + ";\n")

            elif "float cut_vetA=" in line:
                f_cpp.write("\tfloat cut_vetA=" + str(d_cut["sigma_vet"]*d_std_true_events["FWIA"]) + ";\n")
            elif "float cut_vetC=" in line:
                f_cpp.write("\tfloat cut_vetC="  + str(d_cut["sigma_vet"]*d_std_true_events["FWIC"]) + ";\n")

            elif "float EC =" in line:
                f_cpp.write("\t\tfloat EC =" + str(d_estimator["HEAT"][:5]) + "*EC1+" +str(1-float(d_estimator["HEAT"][:5])) + "*EC2;\n")
            # elif "float EI =" in line:
            #     f_cpp.write("\t\t\tfloat EI =" + str(d_estimator["FID"][:5]) + "*EIB+" +str(1-float(d_estimator["FID"][:5]) )+ "*EID;\n")

            elif "float ECinf=" in line:
                f_cpp.write("\tfloat ECinf=" + str(d_cut["ECinf"]) + ";\n")
            elif "float ECsup=" in line:
                f_cpp.write("\tfloat ECsup=" + str(d_cut["ECsup"]) + ";\n")
            elif "float EIinf=" in line:
                f_cpp.write("\tfloat EIinf=" + str(d_cut["EIinf"]) + ";\n")
            elif "float EIsup=" in line:
                f_cpp.write("\tfloat EIsup=" + str(d_cut["EIsup"]) + ";\n")


            elif "TString bolo_name = TString" in line:
                f_cpp.write('\tTString bolo_name = TString("' + bolo_name + '");\n')
            elif "TString analysis_type = TString" in line:
                f_cpp.write('\tTString analysis_type = TString("' + analysis_type + '");\n')
            elif "float exposure =" in line:
                f_cpp.write("\tfloat exposure =" + str(exposure) + ";\n")

            elif "nFidGamma=" in line:
                f_cpp.write("\tnFidGamma=" +  str(nFidGamma) + ";\n")
            elif "nS1Gamma=" in line:
                f_cpp.write("\tnS1Gamma=" +  str(nS1Gamma)  + ";\n")
            elif "nS2Gamma=" in line:
                f_cpp.write("\tnS2Gamma=" +  str(nS2Gamma)  + ";\n")
            elif "nS1Beta=" in line:
                f_cpp.write("\tnS1Beta=" +  str(nS1Beta) + ";\n")
            elif "nS2Beta=" in line:
                f_cpp.write("\tnS2Beta=" + str(nS2Beta) + ";\n")
            elif "nS1Pb=" in line:
                f_cpp.write("\tnS1Pb=" +  str(nS1Pb)  +" ;\n")
            elif "nS2Pb=" in line:
                f_cpp.write("\tnS2Pb=" +  str(nS2Pb) + ";\n")
            elif "nheatonly=" in line:
                f_cpp.write("\tnheatonly=" +  str(nheatonly)  + ";\n")

            else:
                f_cpp.write(line)

