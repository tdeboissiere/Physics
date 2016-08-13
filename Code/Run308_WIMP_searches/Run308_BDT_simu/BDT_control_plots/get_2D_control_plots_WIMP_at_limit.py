from ROOT import *
import argparse
import math
import script_utils as script_utils
import os
import numpy as np
from root_numpy import root2array
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {"red": [], "green": [], "blue": []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict["red"].append([item, r1, r2])
            cdict["green"].append([item, g1, g2])
            cdict["blue"].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap("CustomMap", cdict)

def BDT_contol_plots(bolo_name, mass, analysis_type, simu_number, exposure):


    """Get 2D control plots (BDT value in EI EC plane)

    Detail:
    
    Args:
        bolo_name           = (str) bolometer name
        mass                = (int) WIMP mass
        analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
        simu_number         = (str) data tree simulated number
        exposure            = (float) exposure of the simulated data

    Returns:
        void

    Raises:
        void
    """       



    #############################################
    #SGet all histograms
    ##########################################



    fWIMP     = "./ROOT_files/" + bolo_name + "/" + analysis_type + "/" + bolo_name + "_WIMP_mass_" + str(mass) + "_tree_atlimit.root"
    frealdata     = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/FID837_ana_min4_min4_2_fond.root"

    file_WIMP = TFile(fWIMP, "read")
    file_realdata = TFile(frealdata, "read")

    tWIMP = file_WIMP.Get("t_new0")
    tdata = file_realdata.Get("data")

    print tWIMP.GetEntries("0.5*(EIB+EID)<-0.3")
    print tdata.GetEntries("0.5*(EIB+EID)<-0.3")

    # raw_input()

    arr_WIMP=root2array(fWIMP , "t_new0")
    arr_realdata=root2array(frealdata , "data" )


    list_arr= [arr_realdata] 
    list_title= [""]
    list_fig_name= ["Real_Data"]

    # print arr_data.shape
    print arr_realdata.shape


    for index, arr in enumerate(list_arr) :

        arr_heat=0.5*arr["EC1"]+ (1-0.5)*arr["EC2"]
        arr_ion_fid= 0.5*arr["EIB"]+ (1-0.5)*arr["EID"]

        arr_heat_WIMP=0.5*(arr_WIMP["EC1"]+arr_WIMP["EC2"])
        arr_ion_fid_WIMP= 0.5*(arr_WIMP["EIB"]+arr_WIMP["EID"])        

        plt.figure()
        plt.scatter(arr_heat, arr_ion_fid, color = "0.5", s=2)
        plt.scatter(arr_heat_WIMP, arr_ion_fid_WIMP, color = "red", s=2)

        plt.xlabel("Heat (keV)", fontsize = 20)
        plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
        plt.ylim([-1,0])
        plt.xlim([0.5,5])
        plt.grid(True)
        plt.savefig("./Figures/" + bolo_name + "/" + analysis_type + "/Limit_mass_" + str(mass) + "_" + list_fig_name[index]+   "_WIMP_at_limit.png")
        plt.close("all")


list_mass = [3, 4, 5, 6, 7, 10, 25]
# list_mass = [5, 6, 7, 10, 25]
# list_mass = [6, 25]
# list_mass = [5]
bolo_name = "FID837"
simu_number         = 2
analysis_type       = "ana_0.5_min4_5"
exposure = 65

for mass in list_mass:
    BDT_contol_plots(bolo_name, str(mass), analysis_type, simu_number, exposure)


