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


    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"

    #############################################
    #SGet all histograms
    ##########################################

    fS1Pb     = gen_path +  "S1Pb" + "_mass_" + str(mass) + "_tree.root"
    fS2Pb     = gen_path +  "S2Pb" + "_mass_" + str(mass) + "_tree.root"
    fS1Beta   = gen_path +  "S1Beta" + "_mass_" + str(mass) + "_tree.root"
    fS2Beta   = gen_path +  "S2Beta" + "_mass_" + str(mass) + "_tree.root"
    fS1Gamma  = gen_path +  "S1Gamma" + "_mass_" + str(mass) + "_tree.root"
    fS2Gamma  = gen_path +  "S2Gamma" + "_mass_" + str(mass) + "_tree.root"
    fFidGamma = gen_path +  "FidGamma" + "_mass_" + str(mass) + "_tree.root"
    fheat     = gen_path +  "heatonly" + "_mass_" + str(mass) + "_tree.root"
    fWIMP     = gen_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root"
    # fdata     = gen_path +  "data" + "_mass_" + str(mass) + "_tree.root"
    frealdata     = gen_path +  "real_data" + "_mass_" + str(mass) + "_tree.root"
    
    # tS1Pb     = fS1Pb.Get("tout")
    # tS2Pb     = fS2Pb.Get("tout") 
    # tS1Beta   = fS1Beta.Get("tout") 
    # tS2Beta   = fS2Beta.Get("tout") 
    # tS1Gamma  = fS1Gamma.Get("tout") 
    # tS2Gamma  = fS2Gamma.Get("tout") 
    # tFidGamma = fFidGamma.Get("tout") 
    # theat     = fheat.Get("tout")
    # tWIMP     = fWIMP.Get("tout") 
    # tdata     = fdata.Get("tout" + str(simu_number)) 

    arr_S1Pb=root2array(fS1Pb , "tout")
    arr_S2Pb=root2array(fS2Pb , "tout")
    arr_S1Beta=root2array(fS1Beta , "tout")
    arr_S2Beta=root2array(fS2Beta , "tout")
    arr_S1Gamma=root2array(fS1Gamma , "tout")
    arr_S2Gamma=root2array(fS2Gamma , "tout")
    arr_FidGamma=root2array(fFidGamma , "tout")
    arr_heat=root2array(fheat , "tout")
    arr_WIMP=root2array(fWIMP , "tout")
    # arr_data=root2array(fdata , "tout" + str(simu_number))
    arr_realdata=root2array(frealdata , "tout" )

    list_arr= [arr_S1Pb, arr_S2Pb, arr_S1Beta, arr_S2Beta, arr_S1Gamma, arr_S2Gamma, arr_FidGamma, arr_heat, arr_WIMP, arr_realdata] 
    list_bckg_arr= [arr_S1Pb, arr_S2Pb, arr_S1Beta, arr_S2Beta, arr_S1Gamma, arr_S2Gamma, arr_FidGamma, arr_heat]
    list_title= [ "S1 Pb", "S2 Pb", "S1 Beta", "S2 Beta", "S1 Gamma", "S2 Gamma", "Fiducial Gamma", "Heat", "WIMP", "Real Data"]
    list_fig_name= [ "S1_Pb", "S2_Pb", "S1_Beta", "S2_Beta", "S1_Gamma", "S2_Gamma", "Fiducial_Gamma", "Heat", "WIMP", "Real Data"]

    c = mcolors.ColorConverter().to_rgb
    rvb = make_colormap([c("red"), c("red"), 0.33, c("red"), c("green"), 0.86, c("green")])

    list_arr= [arr_realdata, arr_FidGamma] 
    list_title= ["", ""]
    list_fig_name= ["Real_Data", "FidGamma"]

    # print arr_data.shape
    print arr_realdata.shape


    for index, arr in enumerate(list_arr) :

        arr_heat=0.5*arr["EC1"]+ (1-0.5)*arr["EC2"]
        arr_ion_fid= 0.5*arr["EIB"]+ (1-0.5)*arr["EID"]

        arr_heat_WIMP=0.5*(arr_WIMP["EC1"]+arr_WIMP["EC2"])
        arr_ion_fid_WIMP= 0.5*(arr_WIMP["EIB"]+arr_WIMP["EID"])        

        plt.figure()
        plt.scatter(arr_heat_WIMP, arr_ion_fid_WIMP, color="0.6", s=1 )
        plt.scatter(arr_heat, arr_ion_fid, s=20, c=arr["NN"], cmap=rvb, vmin = -0.5, vmax = 1)

        # plt.colorbar()
        cbar = plt.colorbar()
        cbar.set_label("BDT output", labelpad = 15, fontsize= 18)


        plt.xlabel("Heat (keV)", fontsize = 20)
        plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
        plt.ylim([0,5])
        plt.xlim([0.5,5])
        plt.grid(True)
        plt.savefig("./Figures/mass_" + str(mass) + "_" + list_fig_name[index]+   ".png")
        plt.close("all")


    arr_heat_WIMP=0.5*(arr_WIMP["EC1"]+arr_WIMP["EC2"])
    arr_ion_fid_WIMP= 0.5*(arr_WIMP["EIB"]+arr_WIMP["EID"])        

    plt.figure()
    plt.scatter(arr_heat_WIMP, arr_ion_fid_WIMP, s=20, c=arr_WIMP["NN"], cmap=rvb, vmin = -0.5, vmax = 1)

    # plt.colorbar()
    cbar = plt.colorbar()
    cbar.set_label("BDT output", labelpad = 15, fontsize= 18)

    plt.xlabel("Heat (keV)", fontsize = 20)
    plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
    plt.ylim([0,5])
    plt.xlim([0.5,5])
    plt.grid(True)
    plt.savefig("./Figures/WIMP_mass_" + str(mass) + "_" + list_fig_name[index]+   ".png")
    plt.close("all")


list_mass = [3, 4, 5, 6, 7, 10, 25]
list_mass = [5, 6, 7, 10, 25]
# list_mass = [5]
bolo_name = "FID837"
simu_number         = 2
analysis_type       = "ana_1.5_0_5"
exposure = 65

for mass in list_mass:
    BDT_contol_plots(bolo_name, str(mass), analysis_type, simu_number, exposure)


