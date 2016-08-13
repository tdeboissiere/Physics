#! /usr/bin/env python

import BDT_file_handler as BDT_fh


def adjust_MVA_files(bolo_name, analysis_type):

    """Adjust MVA files to bolo_name and analysis_type
    
    Detail:
        Scan the file to look for the fields to modify:
        scaling, bolo_name and anaysis_type

    Args:
        bolo_name = (str) bolometer type
        analysis_type = (str) type of analysis (which ion cut ...)

    Returns:
        void

    Raises:
        void
    """

    #Load the scaling
    d_scaling = BDT_fh.open_scaling_file(bolo_name, analysis_type)

    with open("./TMVAClassification.C", "r") as f_cpp:
        list_cpp_lines = f_cpp.readlines()
   
    with open("./TMVAClassification.C", "w") as f_cpp:
        print d_scaling
        for line in list_cpp_lines:

            if "factory->AddBackgroundTree(bckg_heat," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_heat," + d_scaling["prop_heatonly"]  + ");\n") 

            elif "factory->AddBackgroundTree(bckg_FidGamma," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_FidGamma," + d_scaling["prop_FidGamma"]  + ");\n")

            elif "factory->AddBackgroundTree(bckg_S1Gamma," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_S1Gamma," + d_scaling["prop_S1Gamma"]  + ");\n")

            elif "factory->AddBackgroundTree(bckg_S2Gamma," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_S2Gamma," + d_scaling["prop_S2Gamma"]  + ");\n")

            elif "factory->AddBackgroundTree(bckg_S1Beta," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_S1Beta," + d_scaling["prop_S1Beta"]  + ");\n") 

            elif "factory->AddBackgroundTree(bckg_S2Beta," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_S2Beta," + d_scaling["prop_S2Beta"]  + ");\n")

            elif "factory->AddBackgroundTree(bckg_S1Pb," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_S1Pb," + d_scaling["prop_S1Pb"]  + ");\n") 

            elif "factory->AddBackgroundTree(bckg_S2Pb," in line:
                f_cpp.write("\tfactory->AddBackgroundTree(bckg_S2Pb," + d_scaling["prop_S2Pb"]  + ");\n") 

            elif "TString bolo_name = TString" in line:
                f_cpp.write('\tTString bolo_name = TString("' + bolo_name  + '");\n') 

            elif "TString analysis_type = TString" in line:
                f_cpp.write('\tTString analysis_type = TString("' + analysis_type  + '");\n') 

            else:
                f_cpp.write(line)


    with open("./TMVAindividual.C", "r") as f_cppi:
        list_cpp_lines = f_cppi.readlines()
   
    with open("./TMVAindividual.C", "w") as f_cppi:
        print d_scaling
        for line in list_cpp_lines:

            if "TString bolo_name = TString" in line:
                f_cppi.write('\tTString bolo_name = TString("' + bolo_name  + '");\n') 

            elif "TString analysis_type = TString" in line:
                f_cppi.write('\tTString analysis_type = TString("' + analysis_type  + '");\n') 
                
            else:
                f_cppi.write(line)

    with open("./TMVAtruedata.C", "r") as f_cppd:
        list_cpp_lines = f_cppd.readlines()
   
    with open("./TMVAtruedata.C", "w") as f_cppd:
        print d_scaling
        for line in list_cpp_lines:

            if "TString bolo_name = TString" in line:
                f_cppd.write('\tTString bolo_name = TString("' + bolo_name  + '");\n') 

            elif "TString analysis_type = TString" in line:
                f_cppd.write('\tTString analysis_type = TString("' + analysis_type  + '");\n') 
                
            else:
                f_cppd.write(line)