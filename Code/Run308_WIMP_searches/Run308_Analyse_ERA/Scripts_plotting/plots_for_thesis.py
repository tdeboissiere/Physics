#!/usr/bin/env python

import script_utils as script_utils
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh
from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np 
import matplotlib.pylab as plt
import prettyplotlib as ppl

def plot_NTD_events(bolo_name, data_dir, tree_name):

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    txt_path = script_utils.create_directory("./Text_files/" + bolo_name + "/")
    figure_path = script_utils.create_directory("./Figures/" + bolo_name + "/")

    # cut_line = "FWIA<1.5&&FWIB<1.5&&FWIC<1.5&&FWID<1.5&&OWC1<1&&OWC2<1&&(CHIA-RCIA)<0.28&&(CHIB-RCIB)<0.22&&(CHIC-RCIC)<0.18&&(CHID-RCID)<0.23&&TBEF>0.5&&TAFT>0.5&&APAT>1&&MULT==1&&KTH>0&&KTH<1"
    # cut_line+="&& 0.5*(EIA+EIB+EIC+EID)<1"

    # l_heat     = TEventList("l_heat")
    # l_NTD     = TEventList("l_NTD")

    # # 0.063 => 2sigma cut
    # tree.Draw(">>l_NTD",cut_line + "&&XOC1>0.063||XOC2>0.063" )
    # tree.Draw(">>l_heat",cut_line + "&&XOC1<0.063&&XOC2<0.063")

    # arr_heat, arr_NTD = [], []

    # for i in range(l_heat.GetN()):
    #     counter = l_heat.GetEntry(i)
    #     tree.GetEntry(counter)
    #     arr_heat.append([tree.EC1, tree.EC2])

    # for i in range(l_NTD.GetN()):
    #     counter = l_NTD.GetEntry(i)
    #     tree.GetEntry(counter)
    #     arr_NTD.append([tree.EC1, tree.EC2])

    # np.savetxt(txt_path + bolo_name + "_heatonly_heat2D.txt", np.array(arr_heat))
    # np.savetxt(txt_path + bolo_name + "_NTD_heat2D.txt", np.array(arr_NTD))

    arr_heat = np.loadtxt(txt_path + bolo_name + "_heatonly_heat2D.txt")
    arr_NTD = np.loadtxt(txt_path + bolo_name + "_NTD_heat2D.txt")

    ######################
    # Matplotlib plot
    #####################
    #
    almost_black = '#262626'
    import seaborn as sns 

    sns.set_style("white", {"legend.frameon": True})

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(arr_NTD[:,0], arr_NTD[:,1], c=sns.xkcd_rgb["denim blue"], label = "Heat only events high $\chi^2$", edgecolor = almost_black,  linewidth=0.15, s=20)
    plt.scatter(arr_heat[:,0], arr_heat[:,1], c=sns.xkcd_rgb["pale red"], label = "Heat-only events low $\chi^2$", edgecolor = almost_black,  linewidth=0.15, s=20)

    
    plt.xlabel("EC1 (keVee)", fontsize = "20")
    plt.ylabel("EC2 (keVee)", fontsize = "20")
    plt.ylim([-2,10])
    plt.xlim([-2,10])

    plt.legend(loc =9, fontsize = 15)
    plt.savefig(figure_path + bolo_name + "_heatonly_NTD.png")
    print "ok"



def plot_color_Qplot(bolo_name, data_dir, tree_name):

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    fidcut = "&& 0.5*(EIB+EID)>1 && EIA<1 && EIC<1 && Q>0.6"
    vetcut = "&& (EIA>1 && EIB>1 && EID<1) || (EIC>1 && EID>1 && EIB<1)"
    heatcut = "&& 0.5*(EIA+EIB+EIC+EID)<1"

    txt_path = script_utils.create_directory("./Text_files/" + bolo_name + "/")
    figure_path = script_utils.create_directory("./Figures/" + bolo_name + "/")

    # print standard_cuts
    # standard_cuts = "FWIA<1.5&&FWIB<1.5&&FWIC<1.5&&FWID<1.5&&OWC1<1&&OWC2<1&&(CHIA-RCIA)<0.28&&(CHIB-RCIB)<0.22&&(CHIC-RCIC)<0.18&&(CHID-RCID)<0.23&&TBEF>0.5&&TAFT>0.5&&APAT>1&&MULT==1&&KTH>0&&KTH<1"

    # l_heat     = TEventList("l_heat")
    # l_vet     = TEventList("l_vet")
    # l_fid     = TEventList("l_fid")

    # tree.Draw(">>l_heat",standard_cuts + heatcut )
    # tree.Draw(">>l_fid",standard_cuts + fidcut )
    # tree.Draw(">>l_vet",standard_cuts + vetcut )

    # arr_heat, arr_fid, arr_vet = [], [], []

    # for i in range(l_heat.GetN()):
    #     counter = l_heat.GetEntry(i)
    #     tree.GetEntry(counter)
    #     arr_heat.append([tree.ER, tree.Q])

    # for i in range(l_vet.GetN()):
    #     counter = l_vet.GetEntry(i)
    #     tree.GetEntry(counter)
    #     arr_vet.append([tree.ER, tree.Q])

    # for i in range(l_fid.GetN()):
    #     counter = l_fid.GetEntry(i)
    #     tree.GetEntry(counter)
    #     arr_fid.append([tree.ER, tree.Q])

    # np.savetxt(txt_path + bolo_name + "_heatonly.txt", np.array(arr_heat))
    # np.savetxt(txt_path + bolo_name + "_fid.txt", np.array(arr_fid))
    # np.savetxt(txt_path + bolo_name + "_vet.txt", np.array(arr_vet))

    # raw_input()

    arr_heat = np.loadtxt(txt_path + bolo_name + "_heatonly.txt")
    arr_fid = np.loadtxt(txt_path + bolo_name + "_fid.txt")
    arr_vet = np.loadtxt(txt_path + bolo_name + "_vet.txt")

    ######################
    # Matplotlib plot
    #####################
    #
    almost_black = '#262626'
    import seaborn as sns 

    sns.set_style("white", {"legend.frameon": True})

    plt.scatter(arr_vet[:,0], arr_vet[:,1], c=sns.xkcd_rgb["medium green"], label = "Surface events", edgecolor = almost_black,  linewidth=0.15, s=20)
    plt.scatter(arr_fid[:,0], arr_fid[:,1], c=sns.xkcd_rgb["denim blue"], label = "Fiducial events", edgecolor = almost_black,  linewidth=0.15, s=20)
    plt.scatter(arr_heat[:,0], arr_heat[:,1], c = sns.xkcd_rgb["pale red"], label = "Heat-only events", edgecolor = almost_black,  linewidth=0.15, s=20)

    
    plt.xlabel("Recoil Energy (keV)", fontsize = "20")
    plt.ylabel("Ionisation Quenching", fontsize = "20")
    plt.ylim([-0.2,1.8])
    plt.xlim([0,100])

    plt.legend(loc ="best", fontsize = 15)
    # plt.show()
    plt.savefig(figure_path + bolo_name + "_Qplot_matplotlib.png")
    print "ok"
    # raw_input()


def plot_Pb_hist(bolo_name, data_dir, tree_name):

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    txt_path = script_utils.create_directory("./Text_files/" + bolo_name + "/")
    figure_path = script_utils.create_directory("./Figures/" + bolo_name + "/")

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    hPb= TH2F("hPb", "hPb", 100, 70, 100, 100, -0.2,1.2)
    hPbQ0= TH2F("hPbQ0", "hPbQ0", 100, 70, 100, 100, -0.2,1.2)

    tree.Project("hPb", "Q:ER", standard_cuts + " && Q>0.04 && Q<0.3")
    tree.Project("hPbQ0", "Q:ER", standard_cuts + " && Q<0.04")

    PyRPl.process_TH2(hPbQ0, X_title = "ER (keV)", Y_title = "Counts", color = kRed)
    PyRPl.process_TH2(hPb, X_title = "ER (keV)", Y_title = "Counts", color = kBlack)

    cc=TCanvas("cc","cc")
    hPb.Draw("")   
    hPbQ0.Draw("same")   

    raw_input()


    hPb= TH1F("hPb", "hPb", 20, 70, 100)
    hPbQ0= TH1F("hPbQ0", "hPbQ0", 20, 70, 100)

    tree.Project("hPb", "ER", standard_cuts + " && Q>0.04 && Q<0.3")
    tree.Project("hPbQ0", "ER", standard_cuts + " && Q<0.04")
    
    PyRPl.process_TH1(hPbQ0, X_title = "ER (keV)", Y_title = "Counts", color = kRed)
    PyRPl.process_TH1(hPb, X_title = "ER (keV)", Y_title = "Counts", color = kBlack)

    cc=TCanvas("cc","cc")
    hPb.Draw("")   
    hPbQ0.Draw("same")   

    print hPb.Integral(), hPbQ0.Integral()

    raw_input()
    figure_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    cc.Print(figure_path + bolo_name + "_heatonly_NTD.png")


def calibration_plots():

    """Get calibration plots
    
    Detail:

    Args:
        void

    Returns:
        void

    Raises:
        void
    """

    ##########################################################
    #First make a mock calibration ionisC versus ionisD plot
    ##########################################################
    x = np.random.uniform(low=0.0, high=2000, size=10000)
    y = [np.random.normal(0.25*elem,50) for elem in x]

    x_noise = np.random.uniform(low=0.0, high=2000, size=100)
    y_noise = np.random.uniform(low=0.0, high=2000, size=100)

    x_gaussian_noise = np.random.uniform(low=0.0, high=2000, size=4000)
    y_gaussian_noise = [np.random.normal(0.25*elem,400) for elem in x]

    x_g, y_g = [], []

    for elx, ely in zip(x_gaussian_noise, y_gaussian_noise):
        if ely < 0.25*elx :
            x_g.append(elx)
            y_g.append(ely)

    h =TH2F("h", "h", 100, 0, 2000, 100, 0,2000)
    h_noise =TH2F("h_noise", "h_noise", 100, 0, 2000, 100, 0,2000)
    h_gaussian_noise =TH2F("h_gaussian_noise", "h_gaussian_noise", 100, 0, 2000, 100, 0,2000)
    PyRPl.fill_TH2(h,x,y)
    PyRPl.fill_TH2(h_noise,x_noise,y_noise)
    PyRPl.fill_TH2(h_gaussian_noise,x_g,y_g)

    PyRPl.process_TH2(h, X_title = "Ionisation D (ADU)", Y_title = "Ionisation C (ADU)", X_title_offset=1.14, Y_title_offset = 1.14)

    cc = TCanvas("cc", "cc")

    h.Draw() 
    h.Fit("pol1")
    h_noise.Draw("same")
    h_gaussian_noise.Draw("same")

    figure_path = script_utils.create_directory("./Figures/")
    cc.Print(figure_path + "crosstalk_plot.png")

    ##########################################################
    #Then plot the 356 keV peak calibration
    ##########################################################

    h356, f356 = PyRPl.open_ROOT_object("./ROOT_files/FID845_ionisD_calibration.root", "h1")
    fitfunc = f356.Get("func1")

    class gaus_fit:
        def __call__( self, x, par ):
            return 5*fitfunc.Eval(x[0]/5.) + par[0]

    hcalib = TH1F("hcalib", "hcalib", 250, 1000,2250)
    for i in range(1,251):
        hcalib.SetBinContent(i, 5*h356.GetBinContent(i))

    fgaus = TF1("gaus_fit", gaus_fit(), 1655, 2000, 1)
    fgaus.SetParameter(0,0)

    PyRPl.process_TH1(hcalib, X_title = "Ionisation D (ADU)", Y_title = "Counts/ADU", X_title_offset=1.14, Y_title_offset = 1.14)
    hcalib.Draw()
    fgaus.Draw("same")
    figure_path = script_utils.create_directory("./Figures/")
    cc.Print(figure_path + "ionisD_calibration.png")

def get_heat_over_time(bolo_name, data_dir, tree_name, func_eff):

    """ Heat over time plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    # d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])


    heat_cut = "&&EI<1 && EC>1.5 && EC<15"
    heat_cut = standard_cuts +  heat_cut 

    l_heat   = TEventList("l_heat")
    tree.Draw(">>l_heat",heat_cut)
    arr = []
    arr_kth  =[]

    # tree.Draw("EC>>hist(150,0,15)", heat_cut )
    # raw_input()

    for i in range(l_heat.GetN()):
        counter = l_heat.GetEntry(i)
        tree.GetEntry(counter)
        arr.append(tree.JOUR)
        arr_kth.append(tree.KTH)

    arr = np.array(arr)
    bin, min_X, max_X = 150, 0, 15

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    trucx,trucy,o = plt.hist(arr, bins = bin, histtype = "step", color = "k")
    plt.xlabel("Temps (jours)", fontsize = 16, labelpad = 20)
    plt.ylabel("Counts", fontsize = 16, labelpad = 20)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.set_ylim([0,5])
    # plt.plot(np.linspace(0,100,1000), [func_eff.Eval(x) for x in np.linspace(0,100,1000)], "r--")
    ax2.set_ylabel("Online trigger (keVee)", fontsize = 16, labelpad = 20)
    plt.scatter(arr, np.array(arr_kth), color = "r", s=0.1)
    # plt.xlim([min_X, max_X])
    # ax.set_ylim([0.5, 1.2*max(trucx)])
    # ax.set_yscale("log")
    plt.tight_layout()
    # plt.show()


    # cc = TCanvas("cc", "cc")
    # h1D = TH1F("h1D", "h1D", 150, 0, 15)
    # heff = TH1F("heff", "h1D", 150, 0, 15)
    # for i in range(1,151):
    #     heff.SetBinContent(i, func_eff.Eval(heff.GetBinCenter(i)))
    # tree.Project("h1D", "EC", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + heat_cut)   
    # PyRPl.process_TH1(h1D, X_title = "Combined Heat EC (keVee)", Y_title = "Counts/" + str(h1D.GetBinWidth(2))[:4] + " keV") 
    # h1D.Divide(heff)
    # h1D.Draw("")
    # gPad.SetLogy()

    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    plt.savefig(fig_path + bolo_name + "_heat_over_time.png")


def get_fid_gamma_plot(bolo_name, data_dir, tree_name, func_eff):

    """ Fid Gamma plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    # d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])


    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"
    gamma_cut = standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut 
    gamma_cut = standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] +  " &&EIB>" + d_std["FWIB"] + "&& EID>" + d_std["FWID"] + "&& EC>0.5"

    l_gamma   = TEventList("l_gamma")
    tree.Draw(">>l_gamma",gamma_cut)
    arr = []

    # tree.Draw("EC>>hist(150,0,15)", gamma_cut )
    # raw_input()

    for i in range(l_gamma.GetN()):
        counter = l_gamma.GetEntry(i)
        tree.GetEntry(counter)
        arr.append(tree.EC)

    arr = np.array(arr)
    bin, min_X, max_X = 150, 0, 15
    l = np.where(np.logical_and(arr>min_X, arr<max_X))
    arr = arr[l]
    list_weights = [1./func_eff.Eval(x) for x in arr]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    trucx,trucy,o = plt.hist(arr, bins = bin, histtype = "step", weights=list_weights, color = "k")
    plt.xlabel("Combined Heat EC (keVee)", fontsize = 16, labelpad = 20)
    plt.ylabel("Counts/" + str((max_X-min_X)/float(bin))[:4] + " keV", fontsize = 16, labelpad = 20)
    ax = plt.gca()
    ax2 = ax.twinx()
    # ax2.set_ylim([0,10])
    plt.plot(np.linspace(0,100,1000), [func_eff.Eval(x) for x in np.linspace(0,100,1000)], "r--")
    ax2.set_ylabel("Online trigger efficiency", fontsize = 16, labelpad = 20)
    plt.xlim([min_X, max_X])
    ax.set_ylim([0.3, 1.2*max(trucx)])
    ax.set_yscale("log")
    plt.tight_layout()


    # cc = TCanvas("cc", "cc")
    # h1D = TH1F("h1D", "h1D", 150, 0, 15)
    # heff = TH1F("heff", "h1D", 150, 0, 15)
    # for i in range(1,151):
    #     heff.SetBinContent(i, func_eff.Eval(heff.GetBinCenter(i)))
    # tree.Project("h1D", "EC", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut)   
    # PyRPl.process_TH1(h1D, X_title = "Combined Heat EC (keVee)", Y_title = "Counts/" + str(h1D.GetBinWidth(2))[:4] + " keV") 
    # h1D.Divide(heff)
    # h1D.Draw("")
    # gPad.SetLogy()

    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    plt.savefig(fig_path + bolo_name + "_spectrum_FidGamma.png")



def get_fid_fit_gamma_plot(bolo_name, data_dir, tree_name, func_eff):

    """ Fid Gamma plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    # d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])


    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"
    gamma_cut = standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut 
    gamma_cut = standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] +  " &&EIB>" + d_std["FWIB"] + "&& EID>" + d_std["FWID"] + "&& EC>0.5"

    l_gamma   = TEventList("l_gamma")
    tree.Draw(">>l_gamma",gamma_cut)
    arr = []

    # tree.Draw("EC>>hist(150,0,15)", gamma_cut )
    # raw_input()

    for i in range(l_gamma.GetN()):
        counter = l_gamma.GetEntry(i)
        tree.GetEntry(counter)
        arr.append(tree.EC)

    arr = np.array(arr)
    bin, min_X, max_X = 150, 0, 15
    l = np.where(np.logical_and(arr>min_X, arr<max_X))
    arr = arr[l]
    list_weights = [1./func_eff.Eval(x) for x in arr]

    #Open fit function
    ffit, filefit = PyRPl.open_ROOT_object("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Gamma_spectrum.root", "FidGamma")
    arr_x = np.linspace(min_X, max_X, 2000)
    arr_y = [ffit.Eval(x) for x in arr_x]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    trucx,trucy,o = plt.hist(arr, bins = bin, histtype = "step", weights=list_weights, color = "k")
    plt.plot(arr_x,arr_y, "r-")
    plt.xlabel("Combined Heat EC (keVee)", fontsize = 16, labelpad = 20)
    plt.ylabel("Counts/" + str((max_X-min_X)/float(bin))[:4] + " keV", fontsize = 16, labelpad = 20)
    ax = plt.gca()
    ax2 = ax.twinx()
    # ax2.set_ylim([0,10])
    plt.plot(np.linspace(0,100,1000), [func_eff.Eval(x) for x in np.linspace(0,100,1000)], "r--")
    ax2.set_ylabel("Online trigger efficiency", fontsize = 16, labelpad = 20)
    plt.xlim([min_X, max_X])
    ax.set_ylim([0.3, 1.2*max(trucx)])
    ax.set_yscale("log")
    plt.tight_layout()


    # cc = TCanvas("cc", "cc")
    # h1D = TH1F("h1D", "h1D", 150, 0, 15)
    # heff = TH1F("heff", "h1D", 150, 0, 15)
    # for i in range(1,151):
    #     heff.SetBinContent(i, func_eff.Eval(heff.GetBinCenter(i)))
    # tree.Project("h1D", "EC", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut)   
    # PyRPl.process_TH1(h1D, X_title = "Combined Heat EC (keVee)", Y_title = "Counts/" + str(h1D.GetBinWidth(2))[:4] + " keV") 
    # h1D.Divide(heff)
    # h1D.Draw("")
    # gPad.SetLogy()

    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    plt.savefig(fig_path + bolo_name + "_spectrum_FidGamma_with_fit.png")

def get_surf_bckg_plot(bolo_name, evt_type, bin, min_X, max_X, func_eff):

    """ Surf plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        evt_type  = (str) event type
        bin, min_X, max_X = (int, float, float) TH1F parameters

    Returns:
        void

    Raises:
        void
    """

    arr_calib = []

    if evt_type == "Beta" or evt_type == "Pb":

        arr_S1Beta_ref = np.loadtxt("../Text_files/S1Beta_heatremoved.txt", delimiter=",").astype(float)
        arr_S2Beta_ref = np.loadtxt("../Text_files/S2Beta_heatremoved.txt", delimiter=",").astype(float)

        arr_S1Beta_ref = 0.5*(arr_S1Beta_ref[:,2] + arr_S1Beta_ref[:,3])
        arr_S2Beta_ref = 0.5*(arr_S2Beta_ref[:,2] + arr_S2Beta_ref[:,3])

        arr_S1Pb_ref = np.loadtxt("../Text_files/S1Pb_heatremoved.txt", delimiter=",").astype(float)
        arr_S2Pb_ref = np.loadtxt("../Text_files/S2Pb_heatremoved.txt", delimiter=",").astype(float)

        arr_S1Pb_ref = 0.5*(arr_S1Pb_ref[:,2] + arr_S1Pb_ref[:,3])
        arr_S2Pb_ref = 0.5*(arr_S2Pb_ref[:,2] + arr_S2Pb_ref[:,3])

        if evt_type =="Beta":
            arr_calib = np.hstack((arr_S1Beta_ref, arr_S2Beta_ref))
            l = np.where(np.logical_and(arr_calib>min_X, arr_calib<max_X))
            arr_calib = arr_calib[l]
        elif evt_type =="Pb":
            arr_calib = np.hstack((arr_S1Pb_ref, arr_S2Pb_ref))
            l = np.where(np.logical_and(arr_calib>min_X, arr_calib<max_X))
            arr_calib = arr_calib[l]

    file_path="../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
    fS1 = bolo_name + "_S1" + evt_type + ".txt"    
    fS2 = bolo_name + "_S2" + evt_type + ".txt"    
    arrS1 = np.loadtxt(file_path + fS1)
    arrS2 = np.loadtxt(file_path + fS2)

    arr = np.hstack((arrS1, arrS2))
    l = np.where(np.logical_and(arr>min_X, arr<max_X))
    arr = arr[l]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    trucx,trucy,o = plt.hist(arr, bins = bin, histtype = "step", color = "k")
    if evt_type == "Beta" or evt_type == "Pb":
        list_weights = (arr.shape[0]/float(arr_calib.shape[0]))*np.ones(arr_calib.shape[0])
        plt.hist(arr_calib, bins = bin, histtype = "step", weights = list_weights, color = "r")
    plt.xlabel("Combined Heat EC (keVee)", fontsize = 16, labelpad = 20)
    plt.ylabel("Counts/" + str((max_X-min_X)/float(bin))[:4] + " keV", fontsize = 16, labelpad = 20)
    ax = plt.gca()

    if evt_type == "Gamma":
        ax2 = ax.twinx()
        # ax2.set_ylim([0,10])
        plt.plot(np.linspace(0,100,1000), [func_eff.Eval(x) for x in np.linspace(0,100,1000)], "r--")
        ax2.set_ylabel("Online trigger efficiency", fontsize = 16, labelpad = 20)
    plt.xlim([min_X, max_X])
    ax.set_ylim([0, 1.2*max(trucx)])
    plt.tight_layout()
    # plt.show() 
    # raw_input()

    # cc = TCanvas("cc", "cc")
    # h1 = TH1F ("h1", "h1", bin, min_X, max_X)
    # hcalib = TH1F ("hcalib", "hcalib", bin, min_X, max_X)
    # heff = TH1F("heff", "h1D", bin, min_X, max_X)
    # for i in range(1,bin+1):
    #     heff.SetBinContent(i, func_eff.Eval(heff.GetBinCenter(i)))

    # PyRPl.fill_TH1(h1, arr)
    # PyRPl.process_TH1(h1, X_title = "Combined Heat EC (keVee)", Y_title = "Counts/" + str(h1.GetBinWidth(2))[:4] + " keV" )
    
    # if evt_type =="Beta" or evt_type =="Pb":
    #     PyRPl.fill_TH1(hcalib, arr_calib)
    #     hcalib.Scale(h1.Integral()/hcalib.Integral())
    #     hcalib.SetLineColor(kRed)
    #     h1.SetMaximum(1.2*max( h1.GetMaximum(), hcalib.GetMaximum() ) )

    # # hcalib.Draw("same")

    # axis =  TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),0,1,510,"+L")
    # # axis.SetLineColor(kRed)
    # # axis.SetLabelColor(kRed)
    # axis.Draw()

    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    plt.savefig(fig_path + bolo_name + "_spectrum_" + evt_type + ".png")


def get_heatonly_bckg_plot(bolo_name, evt_type, bin, min_X, max_X, func_eff):

    """ Surf plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        evt_type  = (str) event type
        bin, min_X, max_X = (int, float, float) TH1F parameters

    Returns:
        void

    Raises:
        void
    """


    file_path="../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
    f = bolo_name + "_" + evt_type + ".txt"  
    arr = np.loadtxt(file_path + f)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # machx, machy,o = plt.hist(arr,weights = 10*np.ones(arr.shape[0]), bins = bin, histtype = "step", color = "k")
    # plt.clf()
    trucx,trucy,o = plt.hist(arr, bins = bin, histtype = "step", color = "k")
    plt.xlabel("Combined Heat EC (keVee)", fontsize = 16, labelpad = 20)
    plt.ylabel("Counts/" + str((max_X-min_X)/float(bin))[:4] + " keV", fontsize = 16, labelpad = 20)
    ax = plt.gca()
    ax2 = ax.twinx()
    # ax2.set_ylim([0,10])
    plt.plot(np.linspace(0,100,1000), [func_eff.Eval(x) for x in np.linspace(0,100,1000)], "r--")
    ax2.set_ylabel("Online trigger efficiency", fontsize = 16, labelpad = 20)
    plt.xlim([min_X, max_X])
    ax.set_yscale("log")
    ax.set_ylim([0, 1.2*max(trucx)])
    plt.tight_layout()

    # cc = TCanvas("cc", "cc")
    # h1 = TH1F ("h1", "h1", bin, min_X, max_X)
    # PyRPl.fill_TH1(h1, arr)

    # heff = TH1F("heff", "h1D", bin, min_X, max_X)
    # for i in range(1,bin+1):
    #     heff.SetBinContent(i, func_eff.Eval(heff.GetBinCenter(i)))

    # gPad.SetLogy()

    # PyRPl.process_TH1(h1, X_title = "Combined Heat EC (keVee)", Y_title = "Counts/" + str(h1.GetBinWidth(2))[:4] + " keV" )
    # h1.Draw("")

    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    plt.savefig(fig_path + bolo_name + "_spectrum_" + evt_type + ".png")

def get_plot_HR(bolo_name):

    """ Heat over time plot
    
    Detail:


    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """


    file_path_heat = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_1.5_0_5/Heatonly/ROOT_files/"
    tree_heat, file_heat    = PyRPl.open_ROOT_object(file_path_heat+bolo_name+"_heatonly_tree.root", "t_new0")

    file_path_WIMP = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_1.5_0_5/WIMP/ROOT_files/"
    tree_WIMP,file_WIMP   = PyRPl.open_ROOT_object(file_path_WIMP+bolo_name+"_WIMP_mass_10_tree.root", "t_new0")

    HR_heat = TH1F("HR_heat", "HR_heat", 30, 0, 0.06)
    HR_WIMP = TH1F("HR_WIMP", "HR_WIMP", 30, 0, 0.06)

    tree_heat.Project("HR_heat", "HR")
    tree_WIMP.Project("HR_WIMP", "HR")

    PyRPl.process_TH1(HR_heat, X_title = "HR", Y_title = "Counts (a.u.)", color = kRed)
    PyRPl.process_TH1(HR_WIMP, X_title = "HR", Y_title = "Counts (a.u.)", color = kBlack)

    HR_heat.Scale(1./HR_heat.Integral())
    HR_WIMP.Scale(1./HR_WIMP.Integral())

    HR_heat.SetFillColor(kRed-10)
    HR_heat.SetFillStyle(3003)
    HR_WIMP.SetFillColor(kGray)

    cc = TCanvas("cc", "cc")
    HR_WIMP.Draw()
    HR_heat.Draw("same")

    leg = TLegend(0.55,0.57,0.84,0.88,"", "brNDC")
    leg.AddEntry(HR_heat.GetName(),"Heat-only" ,"f")
    leg.AddEntry(HR_WIMP.GetName(),"Non heat-only", "f")
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(0)
    leg.Draw("same")
    
    raw_input()

    fig_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    cc.Print(fig_path + bolo_name + "_HR.png")


def launch_plots(bolo_name, data_dir, tree_name = "data"):

    """Get general plots for the given detector
    
    Detail:
        Show Q plot, Ion/Heat plot, heatonly plot, trigger plots..

    Args:
        bolo_name = (str) bolometer type
        data_dir  = (str) the ROOT data tree directory
        tree_name = (str) the ROOT data tree name

    Returns:
        void

    Raises:
        void
    """

    
    # Load start and end times
    tmin, tmax = Ana_ut.open_start_and_end_file(bolo_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    d_est = BDT_fh.open_estimator_file(bolo_name, "")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name, "")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])

    #Load efficiency function 
    func_eff, file_eff = PyRPl.open_ROOT_object("./ROOT_files/" + bolo_name + "_weighted_eff.root", "weighted_eff")

    # plot_NTD_events(bolo_name, data_dir, tree_name)
    # plot_color_Qplot(bolo_name, data_dir, tree_name)
    # plot_Pb_hist(bolo_name, data_dir, tree_name)
    # calibration_plots()
    # get_surf_bckg_plot(bolo_name, "Beta",  200, 2.5, 60, func_eff)
    # get_surf_bckg_plot(bolo_name, "Gamma",  200, 2, 15, func_eff)
    # get_surf_bckg_plot(bolo_name, "Pb",  100, 7, 50, func_eff)
    # get_fid_gamma_plot(bolo_name, "../Fond_ERA_merged/", "data", func_eff)
    # get_fid_fit_gamma_plot(bolo_name, "../Fond_ERA_merged/", "data", func_eff)
    # get_heatonly_bckg_plot(bolo_name, "heatonly",  150, 0.5, 15, func_eff)
    # get_heat_over_time(bolo_name, data_dir, tree_name, func_eff)
    get_plot_HR(bolo_name)

bolo_list=["FID837"]
data_dir  = "../Fond_ERA_merged/"
for bolo_name in bolo_list:
    launch_plots( bolo_name, data_dir)


