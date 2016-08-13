#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import PyROOTPlots as PyRPl
import matplotlib.pylab as plt 
import seaborn as sns
import numpy as np
"""

Unblinding trick:

compute expected number of events with PDF 2D for sigma = 1E-5

compute real events in sideband (sideband = x percent of WIMP signal, multiply observed events appropriately to have 100 per cent)

the sideband code is in test_WIMP_cut_region.py

generate events to superimpose with heat only

"""

def plot_heatonly_vs_WIMP(bolo_name, mass):
    """Ion/Heat 2D plot of heat + WIMP for given mass
    
    Detail:
        void

    Args:
        bolo_name = (str) bolo name 
        mass      = (str or int ) WIMP mass 
        
    Returns:
        void

    Raises:
        void
    """

    #Load standard cuts
    TCut_path_name = script_utils.create_directory('../Cut_files/')  
    TCut_file_name ="TCuts.txt" 
    file_TCut      ="" 
    #Add an exception if the file does not exist
    try:
        file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
    except IOError:
        script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
        sys.exit()

    # Load the cut values. 
    list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
    standard_cuts        =""
    # Add a boolean flag to check if the bolo has its cuts in the file
    is_bolo_in_file      =False
    for line in list_file_TCut_lines:
        if bolo_name == line[0]:
            standard_cuts = line[1]
            is_bolo_in_file = True
    assert(is_bolo_in_file)

    WIMP_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_min2_min2_5/WIMP/ROOT_files/"
    tWIMP, fWIMp = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_" + str(mass) + "_tree.root", "t_new0")

    ERA_path = "../Fond_ERA_merged/"
    tERA, fERA = PyRPl.open_ROOT_object(ERA_path + bolo_name + "_fond.root", "data")

    heat_cut = standard_cuts + "&& abs(EC1-EC2)<2 && 0.5*(EIB+EID)<0 && 0.5*(EIA+EIB+EIC+EID)<1"
    
    bin_X, min_X, max_X, bin_Y, min_Y, max_Y = 1000, 0.5, 3, 1000, -1, 2
    
    hWIMP= TH2F("hWIMP", "hWIMP", bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    hheat= TH2F("hheat", "hheat", bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    
    tERA.Project("hheat","0.5*(EIB+EID):0.5*(EC1+EC2)", heat_cut)
    tWIMP.Project("hWIMP","0.5*(EIB+EID):0.5*(EC1+EC2)")

    print mass, tWIMP.GetEntries("0.5*(EIB+EID)<0 && 0.5*(EC1+EC2)>2")

    lWIMP = TEventList("lWIMP")
    lheat = TEventList("lheat")

    tERA.Draw(">>lheat", standard_cuts + "&& abs(EC1-EC2)<2 && 0.5*(EIA+EIB+EIC+EID)<1 && 0.5*(EIB+EID)<0.8 && EC>0.5 && EC<15")
    tWIMP.Draw(">>lWIMP")

    # arr_WIMP, arr_heat = [], []
    # for i in range(lWIMP.GetN()) :
    #     c = lWIMP.GetEntry(i)
    #     tWIMP.GetEntry(c)
    #     arr_WIMP.append([0.5*(tWIMP.EC1 + tWIMP.EC2), 0.5*(tWIMP.EIB+tWIMP.EID)])

    # for i in range(lheat.GetN()) :
    #     c = lheat.GetEntry(i)
    #     tERA.GetEntry(c)
    #     arr_heat.append([0.5*(tERA.EC1 + tERA.EC2), 0.5*(tERA.EIB+tERA.EID)])


    # arr_WIMP = np.array(arr_WIMP)
    # arr_heat = np.array(arr_heat)

    # np.savetxt("./Text_files/WIMP_heat_ion_mass_" + str(mass) + ".txt",arr_WIMP, delimiter = ",")
    # np.savetxt("./Text_files/heatonly_heat_ion.txt",arr_heat, delimiter = ",")

    arr_WIMP = np.loadtxt("./Text_files/WIMP_heat_ion_mass_" + str(mass) + ".txt",delimiter = ",")
    arr_heat = np.loadtxt("./Text_files/heatonly_heat_ion.txt", delimiter = ",")

    arr_WIMP = arr_WIMP
    arr_heat = np.array(arr_heat)

    line = [[0.5,0], [3,0]]

    print "starting plot"
    sns.set_style("white", {'grid.color': '.5'})
    sns.kdeplot(arr_WIMP, cmap = sns.dark_palette(sns.xkcd_rgb["ocean"], as_cmap = True))
    plt.plot([0.5,3], [0,0], "k--")
    sns.regplot(arr_heat[:,0], arr_heat[:,1], color =sns.xkcd_rgb["orange red"], fit_reg=False, scatter_kws={"s": 5, "alpha":0.5})
    # plt.scatter(arr_heat[:,0], arr_heat[:,1], color = sns.xkcd_rgb["pale red"], s=2)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid()
    plt.xlabel("Combined Heat (keVee)", fontsize = 20)
    plt.ylabel("Fiducial ionisation (keVee)", fontsize = 20)
    plt.xlim([0.5, 3])
    plt.ylim([-0.5, 2])
    plt.tick_params(axis='both', which='major', labelsize=15) 
    plt.tight_layout()
    plt.savefig("FID837_heat_sideband.png")
    import sys
    sys.exit()

    hheat.SetLineColor(kRed)
    hWIMP.SetLineColor(kBlack)
    hheat.SetMarkerColor(kRed)
    hWIMP.SetMarkerColor(kBlack)
    # hWIMP.SetMarkerStyle(7)
    hheat.SetMarkerStyle(6)

    PyRPl.process_TH2(hheat, X_title = "Heat (keVee)", Y_title = "Ion Fid (keV)")
    PyRPl.process_TH2(hWIMP, X_title = "Heat (keVee)", Y_title = "Ion Fid (keV)")

    #Do the plots
    cc=TCanvas("cc","cc")
    hheat.Draw()
    hWIMP.Draw("same")

    raw_input()

    leg = TLegend(0.46,0.51,0.65,0.87)
    leg.AddEntry("hWIMP", "WIMP " + str(mass) + " GeV" , "leg")
    leg.AddEntry("hheat", "Heat only", "leg")
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(0)
    # leg.Draw("same")
    # raw_input()
    #Create fig directory
    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/WIMP/")
    cc.Print(fig_path + bolo_name + "_heat_ion_WIMP_" + str(mass) + "_vs_heatonly.png")

    #Define triple exp fit function
    class triple_exp:
        def __call__( self, x, par ):
            return par[0]*TMath.Exp(par[1]*x[0]) + par[2]*TMath.Exp(par[3]*x[0]) + par[4]*TMath.Exp(par[5]*x[0])

    f_triple_exp = TF1("heat_extra", triple_exp(), 1, 15, 6)
    f_triple_exp.SetParameters(1,-1,1,-1, 1, -1)
    f_triple_exp.SetParLimits(1,-10,0)
    f_triple_exp.SetParLimits(3,-10,0)
    f_triple_exp.SetParLimits(5,-10,0)

    f_triple_exp.FixParameter(4,0)

    heatspec= TH1F("heatspec", "heatspec", 100, 1, 4)
    tERA.Project("heatspec","0.5*(EC1+EC2)", heat_cut)
    PyRPl.process_TH1(heatspec, X_title = "Heat ERA (keVee)", Y_title = "Counts (keV)")
    heatspec.Draw()
    heatspec.Fit("heat_extra","","",1.3,4)
    f_triple_exp.Draw("same")
    raw_input()
    cc.Print(fig_path + bolo_name + "_heat_spec_fit.png")
    raw_input()

bolo_name = "FID837"
list_mass = [5,6,7,8,9,10,15,20,25,30]
list_mass = [10]

for mass in list_mass:
    plot_heatonly_vs_WIMP(bolo_name, mass)