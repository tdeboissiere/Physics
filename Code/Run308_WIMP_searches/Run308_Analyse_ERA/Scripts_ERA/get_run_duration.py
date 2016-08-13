#!/usr/bin/env python

from ROOT import TTree, TFile, TH2F, TH1D, TGraph, TH1F
import script_utils as script_utils


def get_run_duration(bolo_name, data_dir, tree_name):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A .txt file bolo_name + "_run_start_and_end_time.txt
    it contains the duration for each polar
    """
    # Number of polarisations
    num_polar  = script_utils.get_howmany(bolo_name)
    
    # Open ROOT threshold files to be used for the determination of the run duration
    # Extract hkth, the threshold histogram, from it
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh.root"
    f_thresh = script_utils.open_ROOT_file(ROOT_path_name, ROOT_file_name, "read")
    hkth=f_thresh.Get("hprojtime_FWHM")

    # Define path for the file. Create the directory if it does not exist, then open file
    path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
    file_name= bolo_name + "_polar_start_and_end_time_prelim.txt"  
    file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")

    #Create lists to hold the contents of the .txt file + the list for each polar duration
    list_tmin, list_tmax, list_string_polar, list_duration = [], [], [], []

    #Read the file lines and fill the lists
    list_lines_start_and_end = [elem.rstrip().split(",") for elem in file_start_and_end.readlines()]
    
    for k in range(len(list_lines_start_and_end)):
        list_string_polar.append(list_lines_start_and_end[k][0])
        list_tmin.append(float(list_lines_start_and_end[k][1]))
        list_tmax.append(float(list_lines_start_and_end[k][2]))
    file_start_and_end.close()


    for k in range(num_polar):
        lowbin      = hkth.FindBin(list_tmin[k])
        upbin       = hkth.FindBin(list_tmax[k])
        
        running_bin = lowbin
        time        = 0
        # timebis     = 0
        while (running_bin+1 <=upbin) :
            if (hkth.GetBinContent(running_bin) !=0) :
                time    += hkth.GetBinCenter(running_bin+1)-hkth.GetBinCenter(running_bin)
                # timebis += hkth.GetXaxis().GetBinUpEdge(running_bin)-hkth.GetXaxis().GetBinLowEdge(running_bin)        
            running_bin+=1
        list_duration.append(time)

    outputfile = script_utils.open_text_file(path_name, bolo_name + "_polar_start_and_end_time.txt", "w")
    for k in range(num_polar):
        outputfile_line = str(list_string_polar[k]) + ","  + str(list_tmin[k])  + ","  +  str(list_tmax[k])  +  "," +  str(list_duration[k]) 
        outputfile.write(outputfile_line + '\n')
        #Print file_content
        script_utils.print_utility(outputfile_line) 


    outputfile.close()

def get_run_duration_KTH_cut(bolo_name, data_dir, tree_name):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A .txt file bolo_name + "_run_start_and_end_time.txt
    it contains the duration for each polar
    """
    # Number of polarisations
    num_polar  = script_utils.get_howmany(bolo_name)
    
    # Open ROOT threshold files to be used for the determination of the run duration
    # Extract hkth, the threshold histogram, from it
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh.root"
    f_thresh = script_utils.open_ROOT_file(ROOT_path_name, ROOT_file_name, "read")
    hkth=f_thresh.Get("hkth")

    # Define path for the file. Create the directory if it does not exist, then open file
    path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
    file_name= bolo_name + "_polar_start_and_end_time_prelim.txt"  
    file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")

    #Create lists to hold the contents of the .txt file + the list for each polar duration
    list_tmin, list_tmax, list_string_polar, list_duration = [], [], [], []

    #Read the file lines and fill the lists
    list_lines_start_and_end = [elem.rstrip().split(",") for elem in file_start_and_end.readlines()]
    
    for k in range(len(list_lines_start_and_end)):
        list_string_polar.append(list_lines_start_and_end[k][0])
        list_tmin.append(float(list_lines_start_and_end[k][1]))
        list_tmax.append(float(list_lines_start_and_end[k][2]))
    file_start_and_end.close()


    for k in range(num_polar):
        lowbin      = hkth.FindBin(list_tmin[k])
        upbin       = hkth.FindBin(list_tmax[k])
        
        running_bin = lowbin
        time        = 0
        # timebis     = 0
        while (running_bin+1 <=upbin) :
            # if (0<hkth.GetBinContent(running_bin) <1) :
            if (0<hkth.GetBinContent(running_bin)) :
                time    += hkth.GetBinCenter(running_bin+1)-hkth.GetBinCenter(running_bin)
                # timebis += hkth.GetXaxis().GetBinUpEdge(running_bin)-hkth.GetXaxis().GetBinLowEdge(running_bin)        
            running_bin+=1
        list_duration.append(time)

    print bolo_name, list_duration[0]/(86400.)

    outputfile = script_utils.open_text_file(path_name, bolo_name + "_duration_KTH_cut.txt", "w")
    for k in range(num_polar):
        outputfile_line = str(list_string_polar[k]) + ","  + str(list_tmin[k])  + ","  +  str(list_tmax[k])  +  "," +  str(list_duration[k]/(86400.)) 
        outputfile.write(outputfile_line + '\n')
        #Print file_content
        # script_utils.print_utility(outputfile_line) 

        
    outputfile.close()

def get_run_duration_KTH_heat_cut(bolo_name, data_dir, tree_name, heat_fraction):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file
    heat_fraction (float) heat_cut = heat_rate < heat_fraction*max(hheat)

    Outputs:

    A .txt file bolo_name + "_run_start_and_end_time.txt
    it contains the duration for each polar
    """
    # Number of polarisations
    num_polar  = script_utils.get_howmany(bolo_name)
    
    # Open ROOT threshold files to be used for the determination of the run duration
    # Extract hkth, the threshold histogram, from it
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh.root"
    f_thresh = script_utils.open_ROOT_file(ROOT_path_name, ROOT_file_name, "read")
    hkth=f_thresh.Get("hkth")
    hheat=f_thresh.Get("hheat")

    # Define path for the file. Create the directory if it does not exist, then open file
    path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
    file_name= bolo_name + "_polar_start_and_end_time_prelim.txt"  
    file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")

    #Create lists to hold the contents of the .txt file + the list for each polar duration
    list_tmin, list_tmax, list_string_polar, list_duration = [], [], [], []

    #Read the file lines and fill the lists
    list_lines_start_and_end = [elem.rstrip().split(",") for elem in file_start_and_end.readlines()]
    
    for k in range(len(list_lines_start_and_end)):
        list_string_polar.append(list_lines_start_and_end[k][0])
        list_tmin.append(float(list_lines_start_and_end[k][1]))
        list_tmax.append(float(list_lines_start_and_end[k][2]))
    file_start_and_end.close()

    heat_max = hheat.GetMaximum()

    for k in range(num_polar):
        lowbin      = hkth.FindBin(list_tmin[k])
        upbin       = hkth.FindBin(list_tmax[k])
        
        running_bin = lowbin
        time        = 0
        # timebis     = 0
        while (running_bin+1 <=upbin) :
            # if (0<hkth.GetBinContent(running_bin) <1) :
            if (0<hkth.GetBinContent(running_bin) and hheat.GetBinContent(hheat.FindBin(hkth.GetBinCenter(running_bin)))<heat_fraction*heat_max) :
                time    += hkth.GetBinCenter(running_bin+1)-hkth.GetBinCenter(running_bin)
                # timebis += hkth.GetXaxis().GetBinUpEdge(running_bin)-hkth.GetXaxis().GetBinLowEdge(running_bin)        
            running_bin+=1
        list_duration.append(time)

    print bolo_name, list_duration[0]/(86400.)

    outputfile = script_utils.open_text_file(path_name, bolo_name + "_duration_KTH_heat_cut_" + str(heat_fraction) + ".txt", "w")
    for k in range(num_polar):
        outputfile_line = str(list_string_polar[k]) + ","  + str(list_tmin[k])  + ","  +  str(list_tmax[k])  +  "," +  str(list_duration[k]/(86400.)) 
        outputfile.write(outputfile_line + '\n')
        #Print file_content
        # script_utils.print_utility(outputfile_line) 

        
    outputfile.close()