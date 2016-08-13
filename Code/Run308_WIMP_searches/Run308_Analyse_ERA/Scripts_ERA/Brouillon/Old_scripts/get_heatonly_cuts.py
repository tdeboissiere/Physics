#!/usr/bin/env python

from ROOT import TCanvas, TFile, TTree
import script_utils as script_utils
import os
import get_plots as get_plots

def get_heatonly_cuts(bolo_name, data_dir, tree_name):

    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the tree directory (containing Jules n-tuple)
    tree_name (str) the tree name in the tree file

    Outputs:

    A .txt file "standard_cut_values".txt
    it contains the cut values on all the FWHM and CHI2 and is updated when a new one is added
    """
    
    #Do a ion heat fiducial plot
    fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ion", "Heat_ion" , "Heat", "Ion", d_estimator["HEAT"], 100,  -5,15, d_estimator["FID"], 100,  -5,15
    list_cuts=[standard_cuts] + ["abs(EC1-EC2)<1.5"] + ["0.5*(EIA+ EIB+ EIC+EID)<1 && SDEL>2"] + ["1E6*UT1+UT2>1409.2E6"] #+  ["EIA<1 && EIC<1"]
    get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)
    
    # Prevent TCanvas from closing
    raw_input(script_utils.print_utility(script_utils.COL("After determining the cuts, press any key. ", "blue")))

    # Manually fill in the CHI2 cut after looking at the data.
    str_cut_fwhm = raw_input('Enter (EIA,FWIB, FWIC, FWID, FWC1, FWC2) cut (a, b, c, d, e, f) format: ')

    # Create file if it does no exist
    cut_file_name='../Cut_files/standard_cut_values.txt'

    #2 cases: the file exists or does not.
    if not os.path.isfile(cut_file_name) :
        with open(cut_file_name, 'w') as cut_file: 
            #Include header
            cut_file.write("Bolo_name,FWIA,FWIB,FWIC,FWID,FWC1,FWC2,CHIA,CHIB,CHIC,CHID,CHIC1,CHIC2" + "\n")
            bolo_cut_line =bolo_name+','+ str_cut_fwhm + ',' + str_cut_chi2
            cut_file.write(bolo_cut_line+ '\n')
            script_utils.print_utility(script_utils.COL('Creating file with line ' + bolo_cut_line,'blue'))       

    else:
        # Read the whole file, save its contents.
        cut_lines=[]
        with open(cut_file_name, 'r') as cut_file:
            cut_lines=cut_file.readlines()

        #Rewrite it with the appropriate content
        with open(cut_file_name, 'w') as cut_file:
            #Include header
            cut_file.write("Bolo_name,FWIA,FWIB,FWIC,FWID,FWC1,FWC2,CHIA,CHIB,CHIC,CHID,CHIC1,CHIC2" + "\n")

            # Update if the line already exists. Else, add it to the file
            is_bolo_in_line=False
            for line in cut_lines[1:]:
                if bolo_name in line:
                    bolo_cut_line=bolo_name+','+ str_cut_fwhm + ',' + str_cut_chi2
                    cut_file.write(bolo_cut_line+'\n')
                    script_utils.print_utility(script_utils.COL('Updating line to ' + bolo_cut_line,'blue'))
                    is_bolo_in_line=True
                else:
                    cut_file.write(line)

            if not is_bolo_in_line:
                bolo_cut_line=bolo_name+','+ str_cut_fwhm + ',' + str_cut_chi2
                cut_file.write(bolo_cut_line+ '\n')
                script_utils.print_utility(script_utils.COL('Adding line ' + bolo_cut_line,'blue'))












