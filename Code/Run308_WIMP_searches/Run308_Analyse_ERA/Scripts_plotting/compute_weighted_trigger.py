#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import numpy as np
import sys
import matplotlib.pylab as plt
import PyROOTPlots as PyRPl

def get_threshold_hist(bolo_name, data_dir, tree_name ):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A ROOT file bolo_name + "_threshold.ROOT"
    it contains the threshold histograms/graphs
    """
    
    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)

    #Define path for txt files. Create the directory if it does not exist
    path_name       = script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    polar_file_name = bolo_name + "_polar_start_and_end_time_prelim.txt"

    print path_name + polar_file_name
    start_and_end_array = np.loadtxt(path_name + polar_file_name ,delimiter=",", usecols=(1,2))
    script_utils.print_utility(script_utils.COL('Opening ' + path_name+polar_file_name + ' for reading', 'blue'))    

    #Define absolute min and max of time from the _polar_start_and_end_time_prelim file. 
    tmin = np.amin(start_and_end_array)
    tmax = np.amax(start_and_end_array)

    # Define path for the .txt file. Create the directory if it does not exist, then open file
    TCut_path_name= script_utils.create_directory('../Cut_files/')  
    TCut_file_name="TCuts.txt" 
    file_TCut="" 
    #Add an exception if the file does not exist
    try:
        file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
    except IOError:
        script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
        sys.exit()

    # Load the cut values. 
    list_file_TCut_lines=[line.rstrip().split(",") for line in file_TCut.readlines()]
    bolo_TCut_line=""
    # Add a boolean flag to check if the bolo has its cuts in the file
    is_bolo_in_file=False
    for line in list_file_TCut_lines:
        if bolo_name == line[0]:
            bolo_TCut_line = line[1][:line[1].find("&&(CHIA")] #Do this to get only the FW cut
            is_bolo_in_file = True

    #Add an error message if the bolo is not found
    if not is_bolo_in_file:
        script_utils.print_utility(script_utils.COL("Bolo not in the cut file. Verify process", "fail"))
        sys.exit()


    # 2 D threshold histogram
    hist_thresh      = TH2F("hist_thresh","hist_thresh",10000,tmin,tmax,1000,0,20)
    tree.Draw("KTH:1E6*UT1+UT2>>hist_thresh")

    hist_OWC1      = TH2F("hist_OWC1","hist_OWC1",10000,tmin,tmax,1000,0,20)
    tree.Draw("OWC1:1E6*UT1+UT2>>hist_OWC1")

    hist_OWC2      = TH2F("hist_OWC2","hist_OWC2",10000,tmin,tmax,1000,0,20)
    tree.Draw("OWC2:1E6*UT1+UT2>>hist_OWC2")

    hprojY    = TH1D("hprojY","hprojY",1000,0,20)
    hprojY    = hist_thresh.ProjectionY()
    hprojY.SetName("hprojY")

    hprojY_OWC1    = TH1D("hprojY_OWC1","hprojY_OWC1",1000,0,20)
    hprojY_OWC1    = hist_OWC1.ProjectionY()
    hprojY_OWC1.SetName("hprojY_OWC1")

    hprojY_OWC2    = TH1D("hprojY_OWC2","hprojY_OWC2",1000,0,20)
    hprojY_OWC2    = hist_OWC2.ProjectionY()
    hprojY_OWC2.SetName("hprojY_OWC2")

    hprojX    = TH1D("hprojX","hprojX",1000,0,20)
    hprojX    = hist_thresh.ProjectionX()
    hprojX.SetName("hprojX")
    
    htime     = TH2F("htime","htime",10000,tmin,tmax,100,0,2000)
    tree.Draw("0.5*(EC1+EC2):1E6*UT1+UT2>>htime")

    htime_FWHM     = TH2F("htime_FWHM","htime_FWHM",10000,tmin,tmax,100,0,2000)
    tree.Draw("0.5*(EC1+EC2):1E6*UT1+UT2>>htime_FWHM", bolo_TCut_line)
    
    hprojtime = TH1D("hprojtime","hprojtime",10000,tmin,tmax)
    hprojtime = htime.ProjectionX()
    hprojtime.SetName("hprojtime")    

    hprojtime_FWHM = TH1D("hprojtime_FWHM","hprojtime_FWHM",10000,tmin,tmax)
    hprojtime_FWHM = htime_FWHM.ProjectionX()
    hprojtime_FWHM.SetName("hprojtime_FWHM")   

    #Fill time and threshold lists used to build the 1D threshold histhograms
    list_time, list_kth, list_OWC1, list_OWC2 = [], [], [], []
    for ix in range(1,10000) : 
        temp_thresh, temp_OWC1, temp_OWC2=0,0,0
        for iy in range(1, 1000): 
            c=hprojY.GetBinCenter(iy)
            if ( hist_thresh.GetBinContent(ix,iy)!=0) :
                temp_thresh=c

        for iy in range(1, 1000): 
            c=hprojY_OWC1.GetBinCenter(iy)
            if ( hist_OWC1.GetBinContent(ix,iy)!=0) :
                temp_OWC1=c

        for iy in range(1, 1000): 
            c=hprojY_OWC2.GetBinCenter(iy)
            if ( hist_OWC2.GetBinContent(ix,iy)!=0) :
                temp_OWC2=c

        list_time.append(hprojX.GetBinCenter(ix))
        list_kth.append(temp_thresh)
        list_OWC1.append(temp_OWC1)
        list_OWC2.append(temp_OWC2)

    # ROOT says there is an issue with increasing value for the bins represented by time
    # I disagree.
    # Solution: use a graph to build the 1D threshold histhogram
    gr = TGraph(np.size(list_time), np.array(list_time).astype(float), np.array(list_kth).astype(float))
    grOWC1 = TGraph(np.size(list_time), np.array(list_time).astype(float), np.array(list_OWC1).astype(float))
    grOWC2 = TGraph(np.size(list_time), np.array(list_time).astype(float), np.array(list_OWC2).astype(float))
    hkth = TH1F("hkth", "hkth",np.size(list_time),tmin,tmax)
    hOWC1 = TH1F("hOWC1", "hOWC1",np.size(list_time),tmin,tmax)
    hOWC2 = TH1F("hOWC2", "hOWC2",np.size(list_time),tmin,tmax)
    for k in range(np.size(list_time)):
        hkth.SetBinContent(k,gr.Eval(list_time[k]))
        hOWC1.SetBinContent(k,grOWC1.Eval(list_time[k]))
        hOWC2.SetBinContent(k,grOWC2.Eval(list_time[k]))
    
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh_for_trigger_eff.root"
    f_thresh       = script_utils.open_ROOT_file(ROOT_path_name, ROOT_file_name, "recreate")
    hkth.Write()
    hOWC1.Write()
    hOWC2.Write()
    f_thresh.Close()

def compute_resolution( FWHM1,  FWHM2) :
    """
    Compute the resolution (i.e FWHM not sigma) of the best combination of two channels
    """
    s1=float(FWHM1)/2.3548
    s2=float(FWHM2)/2.3548

    w1=(s2*s2)/(s1*s1+s2*s2)

    res=np.sqrt(pow(w1*s1,2)+pow((1-w1)*s2,2))
    return 2.3548*res


def get_KTH_and_baselines(bolo_name, data_dir, tree_name):
    """
    Get a txt file with time period, KTH, OWC1, OWC2

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:
        void
    """

    
    # Open ROOT threshold files to be used for the determination of the run duration
    # Extract hkth, the threshold histogram, from it
    ROOT_path_name = script_utils.create_directory('../Analyse_' + bolo_name + '/ROOT_files/')
    ROOT_file_name = bolo_name + "_thresh_for_trigger_eff.root"
    f = TFile(ROOT_path_name +  ROOT_file_name, "read")
    hkth=f.Get("hkth")
    hOWC1=f.Get("hOWC1")
    hOWC2=f.Get("hOWC2")

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

    lowbin      = hkth.FindBin(list_tmin[0])
    upbin       = hkth.FindBin(list_tmax[0])
        
    list_time_inf, list_time_sup, list_KTH, list_OWC1, list_OWC2, list_FWC_combined = [], [], [], [], [], []

    running_bin = lowbin
    while (running_bin+1 <=upbin) :
        if (0<hkth.GetBinContent(running_bin) <1) :
                list_time_inf.append(hkth.GetBinLowEdge(running_bin))
                list_time_sup.append(hkth.GetBinLowEdge(running_bin) + hkth.GetBinWidth(running_bin))
                list_KTH.append(hkth.GetBinContent(running_bin))
                list_OWC1.append(hOWC1.GetBinContent(running_bin))
                list_OWC2.append(hOWC2.GetBinContent(running_bin))
                list_FWC_combined.append(compute_resolution( hOWC1.GetBinContent(running_bin),  hOWC2.GetBinContent(running_bin)))
        running_bin+=1

    print np.average(list_KTH), np.average(list_OWC1), np.average(list_OWC2)
    print (np.sum(list_time_sup) - np.sum(list_time_inf))/(86400)

    with open("../Analyse_" + bolo_name + "/Text_files/" + bolo_name + "_kth_baseline_evol.txt", "w") as fevol:
        for tinf, tsup, kth, OWC1, OWC2, fwc in zip(list_time_inf, list_time_sup, list_KTH, list_OWC1, list_OWC2, list_FWC_combined):
            outputfile_line = str(tinf) + ","  + str(tsup)  + ","  +  str(kth)  +  "," +  str(OWC1) + "," + str(OWC2)+ "," + str(fwc) + "\n" 
            fevol.write(outputfile_line )


def compute_average_trigger_eff(bolo_name):
    """
    Compute weighted trigger eff and compare to average

    Detail:

    Arguments:
    bolo_name (str) the bolometer name

    Outputs:
        void
    """

    data_types = {"names": ("tinf", "tsup", "KTH", "OWC1", "OWC2", "FWC"), "formats": ("f", "f", "f", "f", "f", "f")}
    arr = np.loadtxt("../Analyse_" + bolo_name + "/Text_files/" + bolo_name + "_kth_baseline_evol.txt", delimiter=",",  dtype=data_types)

    arr_duration = arr["tsup"] - arr["tinf"]
    duration = float(np.sum(arr_duration))

    #EDW III keVee
    arr_xpoints = np.linspace(0,5,100)
    arr_eff = np.zeros(arr_xpoints.shape[0])
    for i in range(arr["tsup"].shape[0]):
        sigma = arr["FWC"][i]/2.3548
        kth = arr["KTH"][i]
        arr_eff = arr_eff +(arr_duration[i]/duration)*np.array([0.5*(1+TMath.Erf( (elem-kth) / (TMath.Sqrt(2)*sigma))) for elem in arr_xpoints]).astype(float)


    arr_eff = np.array(arr_eff).astype(float)
    arr_eff_avg = np.array([0.5*(1+TMath.Erf( (elem-0.65) / (TMath.Sqrt(2)*0.1648))) for elem in arr_xpoints]).astype(float)

    #Get a TGraph to save the efficiency as a function 
    grweighted_eff = TGraph(len(arr_xpoints), np.array(arr_xpoints).astype(float), arr_eff)

    class weighted_eff:
      def __call__( self, x, par ):
        if x[0]>5 :
            return grweighted_eff.Eval(4.99) + par[0]
        else: 
            return grweighted_eff.Eval(x[0]) + par[0]

    fweighted_eff = TF1("weighted_eff", weighted_eff(), 0, 100,1)
    fweighted_eff.SetParameter(0,0)
    fweighted_eff.SetNpx(1000)
    file_weight_eff = TFile("./ROOT_files/" + bolo_name + "_weighted_eff.root", "recreate")
    fweighted_eff.Write()
    file_weight_eff.Close()

    #Include NR computations directly
    util_path = "/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/"
    fEE_to_NR, file_EE_to_NR = PyRPl.open_ROOT_object(util_path +"conv_EE_to_NR.root", "conv")

    arr_NRxpoints = np.array([fEE_to_NR.Eval(elem) for elem in arr_xpoints]).astype(float)

    #EDW III keVNR
    #0.65 keVee == 1.625 keVNR
    #sigma_kevee @ 0.165 keVee  == sigma_kevNR @ 0.388 keVNR
    #This is "my" theoretical curve for NR
    arr_eff_NR = [0.5*(1+ TMath.Erf( (elem-1.625) / (TMath.Sqrt(2)*0.388))) for elem in arr_xpoints]

    plt.plot(arr_xpoints, arr_eff, 'r', label = "keVee Weighted", linewidth = 2)
    plt.plot(arr_xpoints, arr_eff_avg, 'b', label = "keVee Average", linewidth = 2)
    plt.plot(arr_NRxpoints, arr_eff, 'k', label = "keVNR Weighted", linewidth = 2)
    plt.plot(arr_NRxpoints, arr_eff_avg, 'g', label = "keVNR Average", linewidth = 2)
    # plt.plot(arr_xpoints, arr_eff_NR, 'c', label = "keVNR Average theoretical")
    plt.legend(loc = "lower right")
    plt.xlim([0,4])
    plt.savefig("../../FID837_trigger_eff.png")
    plt.show()
    raw_input()

    np.savetxt("/home/irfulx204/mnt/tmain/Desktop/FID837_eff_keVee_avg.txt", zip(arr_xpoints, arr_eff_avg))
    np.savetxt("/home/irfulx204/mnt/tmain/Desktop/FID837_eff_keVNR_avg.txt", zip(arr_NRxpoints, arr_eff_avg))

    np.savetxt("/home/irfulx204/mnt/tmain/Desktop/FID837_eff_keVee_weighted.txt", zip(arr_xpoints, arr_eff))
    np.savetxt("/home/irfulx204/mnt/tmain/Desktop/FID837_eff_keVNR_weighted.txt", zip(arr_NRxpoints, arr_eff))

bolo_name = "FID837"
data_dir  = "../Fond_ERA_merged/"
# get_threshold_hist( bolo_name, data_dir, tree_name = "data")
# get_KTH_and_baselines( bolo_name, data_dir, tree_name = "data")
compute_average_trigger_eff(bolo_name)