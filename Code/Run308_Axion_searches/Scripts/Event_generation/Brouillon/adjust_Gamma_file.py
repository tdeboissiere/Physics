#! /usr/bin/env python
from ROOT import *
import BDT_file_handler as BDT_fh

def adjust_Gamma_file(bolo_name):

    """Change the value of the FWHM/ energycuts used in the build_Gamma_tree.C file
    
    Detail:
        Read the file once to save its content. 
        Then rewrite it with the appropriate modifications

    Args:
        bolo_name              = (str) bolometer name
        d_cut                  = (dict) dictionnary which states which cut to use on heat/ion/veto
        analysis_type          = (str) name of analysis (which FWHM, which cuts)
        d_event_type_num_event = (dict) events for FidGamma, S1Gamma, S2Gamma 

    Returns:
        void

    Raises:
        void
    """

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Axion/Scripts/Event_generation/"

    list_cpp_lines =[]

    with open(gen_path + "build_Gamma_tree.C", "r") as f_cpp:
        list_cpp_lines = f_cpp.readlines()
   
    with open(gen_path + "build_Gamma_tree.C", "w") as f_cpp:
        for line in list_cpp_lines:

            if "TString bolo_name = TString" in line:
                f_cpp.write('\tTString bolo_name = TString("' + bolo_name + '");\n')
            else:
                f_cpp.write(line)

