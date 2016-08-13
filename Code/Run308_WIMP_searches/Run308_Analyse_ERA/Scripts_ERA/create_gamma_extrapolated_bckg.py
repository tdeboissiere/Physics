#!/usr/bin/env python

import script_utils as script_utils
import sys,os
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh

def compare_to_simu(bolo_name, data_dir, tree_name):

    """Compare the gamma data to the simulated data
    
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
    d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])

    #Get true data
    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"
    gamma_S1cut = "Q>0.7 && EIA>1 && EIB>1 && EIC<2 && EID<2"
    h1D = TH1F("h1D", "h1D", 150, 2, 15)
    # tree.Project("h1D", "EID", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut)   
    tree.Project("h1D", "EC1", standard_cuts +  "&&" + gamma_S1cut)   

    #Get simu data
    simu_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr_with_neutron/BDT_FID837/ana_0.5_0_100/Gamma/ROOT_files/"
    simu_file = TFile(simu_path + bolo_name + "_S1Gamma_tree.root", "read")
    simu_tree = simu_file.Get("t_new0")
    hsimu = TH1F("hsimu", "hsimu", 150, 2, 15)
    simu_tree.Project("hsimu", "EC1")
    hsimu.SetLineColor(kRed)

    h1D.Scale(1./h1D.Integral())
    hsimu.Scale(1./hsimu.Integral())

    h1D.SetMaximum(1.2*max(h1D.GetMaximum(), hsimu.GetMaximum()) )


    h1D.Draw()
    hsimu.Draw("same")
    gPad.SetLogy()

    raw_input()

def get_fid_cut_efficiency(bolo_name, data_dir, tree_name):

    """Get efficiency from silulated data
    
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

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])

    gamma_cut = "EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&& EIB>" + d_std["FWIB"] + "&& EID>" + d_std["FWID"] + "&& 0.5*(EC1+EC2)>0.5"

    #Get simu data
    simu_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_min2_min2_100/Gamma/ROOT_files/"
    simu_file = TFile(simu_path + bolo_name + "_FidGamma_tree.root", "read")
    simu_tree = simu_file.Get("t_new0")
    hsimunocut = TH1F("hsimunocut", "hsimunocut", 150, 0, 15)
    hsimucut = TH1F("hsimucut", "hsimucut", 150, 0, 15)
    simu_tree.Project("hsimunocut", "0.5*(EC1+EC2)")
    simu_tree.Project("hsimucut", "0.5*(EC1+EC2)", gamma_cut)
    hsimucut.SetLineColor(kRed)

    teff = TEfficiency(hsimucut, hsimunocut)
    teff.Draw()
    raw_input()
    
    hsimunocut.Draw()
    hsimucut.Draw("same")


    gPad.SetLogy()
    raw_input()


def prospective_analysis(bolo_name, data_dir, tree_name):

    """Check the 1.3 keV model for the relevant detector
    
    Detail:
        1.3 keV was fixed at 1/10 10.37 keV.
        Check if this is indeed the case

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
    d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(50*d_std[key])

    l_GammaFid = TEventList("l_GammaFid")

    print d_std 
    # raw_input()


    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"
    #First plot the 2D data
    h2D = TH2F("h2D", "h2D", 1000, 0.5, 5, 1000, 0, 5)
    h2Dcut = TH2F("h2Dcut", "h2Dcut", 1000, .5, 5, 1000, 0, 5)
    h2Dsurf = TH2F("h2Dsurf", "h2Dsurf", 1000, .5, 5, 1000, 0, 5)
    tree.Project("h2D", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"])
    tree.Project("h2Dcut", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts + "&& EIA<2.7/50*" + d_std["FWIA"] + "&& EIC <2.7/50*" + d_std["FWIC"] + "&& EIB>2.7/50*" + d_std["FWIB"] + "&& EID>2.7/50*" + d_std["FWID"] + "&& EC>0.5")
    tree.Project("h2Dsurf", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts + "&& (EIA>4/50*" + d_std["FWIA"] + "|| EIC >4/50*" + d_std["FWIC"]  + ") && EC>0.5")
    h2Dcut.SetMarkerColor(kRed)
    h2Dsurf.SetMarkerColor(kBlue)
    l = TLine(0,-.75,15,14.25)
    l2 = TLine(0,2,1.5,0)
    l3 = TLine(0,0.8,15,0.8)
    l4 = TLine(0,0,15,15)
    l4.SetLineColor(kRed)

    print tree.GetEntries(standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&EC>0.5 && 0.5*(EIB+EID)>0.9 && EC<5")
    print tree.GetEntries(standard_cuts  + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&& EIB>2.7/5*" + d_std["FWIB"] + "&& EID>2.7/5*" + d_std["FWID"] + "&& EC>0.5" + " && 0.5*(EIB+EID)>0.9 && EC<5")

    h2D.Draw()
    h2Dcut.Draw("same")
    h2Dsurf.Draw("same")
    l.Draw("same")
    l2.Draw("same")
    l3.Draw("same")
    l4.Draw("same")
    raw_input()
    # c1.Print("./Figures/" + bolo_name + "_biplot_1.3keV.png")
    # sys.exit()

    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"

    h1D = TH1F("h1D", "h1D", 150, 0, 15)
    tree.Project("h1D", "EC", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut)   
    PyRPl.process_TH1(h1D, X_title = "Heat energy (keV)", Y_title = "Counts") 
    h1D.Draw("")
    gPad.SetLogy()

    def reso(E):
        return TMath.Sqrt(0.1388**2+(1.69E-2*E)**2)

    def reso_0(E):
        return 1.69E-2*E

    def gaussian(x,A, mean, sigma):

        pi = TMath.Pi()
        fact = A/(sigma*TMath.Sqrt(2*pi))
        inner = (x-mean)/sigma

        return fact*TMath.Exp(-0.5*inner*inner)



    #Define a fit function
    class Gamma:
        def __call__( self, x, par ):
            flat_part = par[0] 
            peak_10_4 = par[1]*TMath.Gaus(x[0], par[9]*10.37, par[8])
            peak_9_66 = 0.097*par[1]*TMath.Gaus(x[0], par[9]*9.66, par[8])
            peak_8_98 = par[2]*TMath.Gaus(x[0], par[9]*8.98, par[8])
            peak_7_11 = par[3]*TMath.Gaus(x[0], par[9]*7.11, par[8])
            peak_6_54 = par[4]*TMath.Gaus(x[0], par[9]*6.54, par[8])
            peak_5_99 = par[5]*TMath.Gaus(x[0], par[9]*5.99, par[8])
            peak_5_46 = par[6]*TMath.Gaus(x[0], par[9]*5.46, par[8])
            peak_4_97 = par[7]*TMath.Gaus(x[0], par[9]*4.97, par[8])
            peak_1_3  = 0.114*par[1]*TMath.Gaus(x[0], par[9]*1.297, par[8]) +0.108*par[2]*TMath.Gaus(x[0], par[9]*1.0961, par[8])
            peak_1_3+= 0.1**2*par[1]*TMath.Gaus(x[0], par[9]*1.1930, par[8])
            list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
            return flat_part + sum(list_peaks)

    class Gamma_deconv:
        def __call__( self, x, par ):
            flat_part = par[8] 
            reso_10_4, reso_9_66, reso_8_98, reso_7_11, reso_6_54, reso_5_99 = reso_0(10.37), reso_0(9.66), reso_0(8.98), reso_0(7.11), reso_0(6.54), reso_0(5.99)
            reso_5_46, reso_4_97, reso_1_2977, reso_1_0961 = reso_0(5.46), reso_0(4.97), reso_0(1.2977), reso_0(1.0961)

            peak_10_4 = gaussian(x[0], par[1], par[0]*10.37, reso_10_4)
            peak_9_66 = gaussian(x[0], 0.097*par[1], par[0]*9.66, reso_9_66)
            peak_8_98 = gaussian(x[0], par[2], par[0]*8.98, reso_8_98)
            peak_7_11 = gaussian(x[0], par[3], par[0]*7.11, reso_7_11)
            peak_6_54 = gaussian(x[0], par[4], par[0]*6.54, reso_6_54)
            peak_5_99 = gaussian(x[0], par[5], par[0]*5.99, reso_5_99)
            peak_5_46 = gaussian(x[0], par[6], par[0]*5.46, reso_5_46)
            peak_4_97 = gaussian(x[0], par[7], par[0]*4.97, reso_4_97)
            peak_1_2977 = gaussian(x[0], 0.1175*par[1], par[0]*1.2977, reso_1_2977)
            peak_1_0961 = gaussian(x[0], 0.119*par[2], par[0]*1.0961, reso_1_0961)

            list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_2977, peak_1_0961]
            return flat_part + sum(list_peaks)

    class Gamma_reso:
        def __call__( self, x, par ):
            flat_part = par[8] 
            reso_10_4, reso_9_66, reso_8_98, reso_7_11, reso_6_54, reso_5_99 = reso(10.37), reso(9.66), reso(8.98), reso(7.11), reso(6.54), reso(5.99)
            reso_5_46, reso_4_97, reso_1_2977, reso_1_0961 = reso(5.46), reso(4.97), reso(1.2977), reso(1.0961)

            peak_10_4 = gaussian(x[0], par[1], par[0]*10.37, reso_10_4)
            peak_9_66 = gaussian(x[0], 0.097*par[1], par[0]*9.66, reso_9_66)
            peak_8_98 = gaussian(x[0], par[2], par[0]*8.98, reso_8_98)
            peak_7_11 = gaussian(x[0], par[3], par[0]*7.11, reso_7_11)
            peak_6_54 = gaussian(x[0], par[4], par[0]*6.54, reso_6_54)
            peak_5_99 = gaussian(x[0], par[5], par[0]*5.99, reso_5_99)
            peak_5_46 = gaussian(x[0], par[6], par[0]*5.46, reso_5_46)
            peak_4_97 = gaussian(x[0], par[7], par[0]*4.97, reso_4_97)
            peak_1_2977 = gaussian(x[0], 0.1175*par[1], par[0]*1.2977, reso_1_2977)
            peak_1_0961 = gaussian(x[0], 0.119*par[2], par[0]*1.0961, reso_1_0961)
            list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_2977, peak_1_0961]
            return flat_part + sum(list_peaks)



    ###################
    #  Fit 1
    ###################
    # Call FidGamma function to get derivative
    # fFidGamma = TF1( "FidGamma", Gamma(), 0, 15., 10 )
    # fFidGamma.SetNpx(500)
    # fFidGamma.SetParName(0,"p0")
    # fFidGamma.SetParName(1,"A_10.4")
    # fFidGamma.SetParName(2,"A_8.98")
    # fFidGamma.SetParName(3,"A_7.11")
    # fFidGamma.SetParName(4,"A_6.54")
    # fFidGamma.SetParName(5,"A_5.99")
    # fFidGamma.SetParName(6,"A_5.46")
    # fFidGamma.SetParName(7,"A_4.97")
    # fFidGamma.SetParName(8,"sigma")
    # fFidGamma.SetParameters(1,1,1,1,1,1,1,1,0.3,1)
    # for i in range(8):
    #     fFidGamma.SetParLimits(i, 0, 100)
    # cst = h1D.GetBinWidth(20)*h1D.Integral(h1D.FindBin(11.5),h1D.FindBin(14.5))/(14.5-11.5)
    # fFidGamma.FixParameter(0,cst)
    # h1D.Fit("FidGamma","LL","",9.7,15)


    ###################
    #  Fit better reso
    ###################
    fFidGamma = TF1( "FidGamma", Gamma_reso(), 0, 15., 9 )
    fFidGamma.SetNpx(500)
    fFidGamma.SetParName(0,"offset")
    fFidGamma.SetParName(1,"A_10.4")
    fFidGamma.SetParName(2,"A_8.98")
    fFidGamma.SetParName(3,"A_7.11")
    fFidGamma.SetParName(4,"A_6.54")
    fFidGamma.SetParName(5,"A_5.99")
    fFidGamma.SetParName(6,"A_5.46")
    fFidGamma.SetParName(7,"A_4.97")
    fFidGamma.SetParName(8,"flat_part")
    fFidGamma.SetParameters(1,1,1,1,1,1,1,1,1,1)
    for i in range(8):
        fFidGamma.SetParLimits(i, 0, 100)
    cst = h1D.GetBinWidth(20)*h1D.Integral(h1D.FindBin(11.5),h1D.FindBin(14.5))/(14.5-11.5)
    fFidGamma.FixParameter(8,cst)
    h1D.Fit("FidGamma","LL","",0.7,15)


    # fFidGamma_deconv = TF1( "FidGamma_deconv", Gamma_deconv(), 0, 15., 9 )
    # fFidGamma_deconv.SetNpx(500)
    # fFidGamma_deconv.SetParameter(0,fFidGamma.GetParameter(0))
    # fFidGamma_deconv.SetParameter(1,fFidGamma.GetParameter(1))
    # fFidGamma_deconv.SetParameter(2,fFidGamma.GetParameter(2))
    # fFidGamma_deconv.SetParameter(3,fFidGamma.GetParameter(3))
    # fFidGamma_deconv.SetParameter(4,fFidGamma.GetParameter(4))
    # fFidGamma_deconv.SetParameter(5,fFidGamma.GetParameter(5))
    # fFidGamma_deconv.SetParameter(6,fFidGamma.GetParameter(6))
    # fFidGamma_deconv.SetParameter(7,fFidGamma.GetParameter(7))
    # fFidGamma_deconv.SetParameter(8,fFidGamma.GetParameter(8))

    # ff = TFile("gamma.root", "recreate")
    # fFidGamma.SetNpx(1000)
    # fFidGamma_deconv.SetNpx(1000)
    # fFidGamma.Write()
    # fFidGamma_deconv.Write()
    # ff.Close()

    # class Gamma10keV:
    #     def __call__( self, x, par ):
    #         reso_10_4= reso(10.37)
    #         peak_10_4 = gaussian(x[0], par[1], par[0]*10.37, reso_10_4)
    #         return peak_10_4

    # class Gamma1keV:
    #     def __call__( self, x, par ):
    #         reso_1_2977 = reso(1.2977)
    #         peak_1_2977 = gaussian(x[0], 0.1175*par[1], par[0]*1.2977, reso_1_2977)
    #         return peak_1_2977

    # #Check that the ratio are ok
    # f10keV = TF1( "Gamma10keV", Gamma10keV(), 0, 15., 2 )
    # f1keV = TF1( "Gamma1keV", Gamma1keV(), 0, 15., 2 )

    # f10keV.SetParameters(fFidGamma.GetParameter(0), fFidGamma.GetParameter(1))
    # f1keV.SetParameters(fFidGamma.GetParameter(0), fFidGamma.GetParameter(1))

    # print f10keV.Integral(0,15)
    # print f1keV.Integral(0,15)
    # print f1keV.Integral(0,15)/f10keV.Integral(0,15)

    ###################
    #  Fit 10 keV
    ###################
    # f10keV = TF1( "Gamma10keV", Gamma10keV(), 0, 15., 2 )
    # f10keV.SetNpx(500)
    # f10keV.SetParName(0,"p0")
    # f10keV.SetParName(1,"A_10.4")
    # f10keV.SetParName(2,"sigma")
    # f10keV.SetParameters(1,0.3)
    # h1D.Fit("Gamma10keV", "LL","",9.7,10.5)


    raw_input("Print unzoomed")
    c1.Print("./Figures/" + bolo_name + "_1_3keV_fitted.png")
    raw_input("Print zoomed")
    c1.Print("./Figures/" + bolo_name + "_1_3keV_zoom.png")

def create_extrapolated_gamma(bolo_name, data_dir, tree_name):

    """Create extrapolated backg for gamma (fid and surf)
    
    Detail:
        Take the varying resolution into account
        Deconvolve by the base resolution sigma_0 = combined OWC1, OWC2 
        scale fid/surf as 600g/800g

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
    d_est = BDT_fh.open_estimator_file(bolo_name,"")

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
    output_dir = gen_path + "Analyse_" + bolo_name +"/ROOT_files/" + bolo_name 
    pop_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"


    #Get surface histograms
    arr_S1Gamma = np.loadtxt(pop_path + bolo_name + "_S1Gamma.txt", delimiter=",").astype(float)
    arr_S2Gamma = np.loadtxt(pop_path + bolo_name + "_S2Gamma.txt", delimiter=",").astype(float)
    hS1Gamma = TH1F("hS1Gamma", "hS1Gamma", 100, 0, 15)
    hS2Gamma = TH1F("hS2Gamma", "hS2Gamma", 100, 0, 15)
    PyRPl.fill_TH1(hS1Gamma, arr_S1Gamma)
    PyRPl.fill_TH1(hS2Gamma, arr_S2Gamma)


    gamma_cut = "0.5*(EIB+EID)>0.5*(EC1+EC2)-0.75 && 0.5*(EIB+EID)>-4./3.*(0.5*(EC1+EC2))+2 && 0.5*(EIB+EID)>0.8"

    hFidGamma = TH1F("hFidGamma", "hFidGamma", 150, 0, 15)
    tree.Project("hFidGamma", "EC", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&&" + gamma_cut)    
    hFidGamma.Draw()
    gPad.SetLogy()

    def reso(E):
        return TMath.Sqrt(d_std["FWHEAT"]**2+(1.69E-2*E)**2)

    def reso_0(E):
        return 1.69E-2*E

    def gaussian(x,A, mean, sigma):

        pi = TMath.Pi()
        fact = A/(sigma*TMath.Sqrt(2*pi))
        inner = (x-mean)/sigma

        return fact*TMath.Exp(-0.5*inner*inner)

    #Define a fit function
    class Gamma_reso:
        def __call__( self, x, par ):
            flat_part = par[8] 
            reso_10_4, reso_9_66, reso_8_98, reso_7_11, reso_6_54, reso_5_99 = reso(10.37), reso(9.66), reso(8.98), reso(7.11), reso(6.54), reso(5.99)
            reso_5_46, reso_4_97, reso_1_2977, reso_1_0961 = reso(5.46), reso(4.97), reso(1.2977), reso(1.0961)

            peak_10_4 = gaussian(x[0], par[1], par[0]*10.37, reso_10_4)
            peak_9_66 = gaussian(x[0], 0.097*par[1], par[0]*9.66, reso_9_66)
            peak_8_98 = gaussian(x[0], par[2], par[0]*8.98, reso_8_98)
            peak_7_11 = gaussian(x[0], par[3], par[0]*7.11, reso_7_11)
            peak_6_54 = gaussian(x[0], par[4], par[0]*6.54, reso_6_54)
            peak_5_99 = gaussian(x[0], par[5], par[0]*5.99, reso_5_99)
            peak_5_46 = gaussian(x[0], par[6], par[0]*5.46, reso_5_46)
            peak_4_97 = gaussian(x[0], par[7], par[0]*4.97, reso_4_97)
            peak_1_2977 = gaussian(x[0], 0.1175*par[1], par[0]*1.2977, reso_1_2977)
            peak_1_0961 = gaussian(x[0], 0.119*par[2], par[0]*1.0961, reso_1_0961)

            list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_2977, peak_1_0961]
            return flat_part + sum(list_peaks)

    #Define deconvoluted fit function
    class Gamma_deconv:
        def __call__( self, x, par ):
            flat_part = par[8] 
            reso_10_4, reso_9_66, reso_8_98, reso_7_11, reso_6_54, reso_5_99 = reso_0(10.37), reso_0(9.66), reso_0(8.98), reso_0(7.11), reso_0(6.54), reso_0(5.99)
            reso_5_46, reso_4_97, reso_1_2977, reso_1_0961 = reso_0(5.46), reso_0(4.97), reso_0(1.2977), reso_0(1.0961)

            peak_10_4 = gaussian(x[0], par[1], par[0]*10.37, reso_10_4)
            peak_9_66 = gaussian(x[0], 0.097*par[1], par[0]*9.66, reso_9_66)
            peak_8_98 = gaussian(x[0], par[2], par[0]*8.98, reso_8_98)
            peak_7_11 = gaussian(x[0], par[3], par[0]*7.11, reso_7_11)
            peak_6_54 = gaussian(x[0], par[4], par[0]*6.54, reso_6_54)
            peak_5_99 = gaussian(x[0], par[5], par[0]*5.99, reso_5_99)
            peak_5_46 = gaussian(x[0], par[6], par[0]*5.46, reso_5_46)
            peak_4_97 = gaussian(x[0], par[7], par[0]*4.97, reso_4_97)
            peak_1_2977 = gaussian(x[0], 0.1175*par[1], par[0]*1.2977, reso_1_2977)
            peak_1_0961 = gaussian(x[0], 0.119*par[2], par[0]*1.0961, reso_1_0961)

            list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_2977, peak_1_0961]
            return flat_part + sum(list_peaks)

    ###################
    #  Fit
    ###################
    fFidGamma = TF1( "FidGamma", Gamma_reso(), 0, 15., 9 )
    fFidGamma.SetNpx(500)
    fFidGamma.SetParName(0,"offset")
    fFidGamma.SetParName(1,"A_10.4")
    fFidGamma.SetParName(2,"A_8.98")
    fFidGamma.SetParName(3,"A_7.11")
    fFidGamma.SetParName(4,"A_6.54")
    fFidGamma.SetParName(5,"A_5.99")
    fFidGamma.SetParName(6,"A_5.46")
    fFidGamma.SetParName(7,"A_4.97")
    fFidGamma.SetParName(8,"flat_part")
    fFidGamma.SetParameters(1,1,1,1,1,1,1,1,1,1)
    for i in range(8):
        fFidGamma.SetParLimits(i, 0, 100)
    cst = hFidGamma.GetBinWidth(20)*hFidGamma.Integral(hFidGamma.FindBin(11.5),hFidGamma.FindBin(14.5))/(14.5-11.5)
    fFidGamma.FixParameter(8,cst)
    hFidGamma.Fit("FidGamma","LL","",0.7,15)

    #Deconvolve
    fFidGamma_deconv = TF1( "FidGamma_deconv", Gamma_deconv(), 0, 15., 9 )
    fFidGamma_deconv.SetNpx(500)
    fFidGamma_deconv.SetParameter(0,fFidGamma.GetParameter(0))
    fFidGamma_deconv.SetParameter(1,fFidGamma.GetParameter(1))
    fFidGamma_deconv.SetParameter(2,fFidGamma.GetParameter(2))
    fFidGamma_deconv.SetParameter(3,fFidGamma.GetParameter(3))
    fFidGamma_deconv.SetParameter(4,fFidGamma.GetParameter(4))
    fFidGamma_deconv.SetParameter(5,fFidGamma.GetParameter(5))
    fFidGamma_deconv.SetParameter(6,fFidGamma.GetParameter(6))
    fFidGamma_deconv.SetParameter(7,fFidGamma.GetParameter(7))
    fFidGamma_deconv.SetParameter(8,fFidGamma.GetParameter(8))

    #Surface Gamma function
    #Take 25% amplitude on peaks, subtract the flat part of FidGamma
    class Gamma_Surf:
        def __call__( self, x, par ):
            return 0.25*fFidGamma.Eval(x[0]*(1+8./3)/(1+5.5/3)) + par[0] -fFidGamma.GetParameter(8)

    fS1Gamma = TF1("S1Gamma", Gamma_Surf(), 0, 15, 1)
    fS2Gamma = TF1("S2Gamma", Gamma_Surf(), 0, 15, 1)
    cstS1 = hS1Gamma.GetBinWidth(20)*hS1Gamma.Integral(hS1Gamma.FindBin(11.5),hS1Gamma.FindBin(14.5))/(14.5-11.5)
    cstS2 = hS2Gamma.GetBinWidth(20)*hS2Gamma.Integral(hS2Gamma.FindBin(11.5),hS2Gamma.FindBin(14.5))/(14.5-11.5)
    fS1Gamma.FixParameter(0,cstS1)
    fS2Gamma.FixParameter(0,cstS2)

    #Store the non deconvoluted files
    out_file_Gamma =TFile(output_dir + "_Gamma_spectrum.root", "recreate")
    fFidGamma.SetNpx(1000)
    fS1Gamma.SetNpx(1000)
    fS2Gamma.SetNpx(1000)
    fFidGamma.Write()
    fS1Gamma.Write()
    fS2Gamma.Write()
    out_file_Gamma.Close()


    # Deconvolved Surface Gamma function
    # Take 25% amplitude on peaks, subtract the flat part of FidGamma
    class Gamma_Surf_deconv:
        def __call__( self, x, par ):
            return 0.25*fFidGamma_deconv.Eval(x[0]*(1+8./3)/(1+5.5/3)) + par[0] -fFidGamma_deconv.GetParameter(8)

    #Store the convoluted files
    fS1Gamma_deconv = TF1("S1Gamma_deconv", Gamma_Surf_deconv(), 0, 15, 1)
    fS2Gamma_deconv = TF1("S2Gamma_deconv", Gamma_Surf_deconv(), 0, 15, 1)
    fS1Gamma_deconv.FixParameter(0,cstS1)
    fS2Gamma_deconv.FixParameter(0,cstS2)

    #Store the deconvolved files
    out_file_Gamma_deconv =TFile(output_dir + "_Gamma_spectrum_deconv.root", "recreate")
    fFidGamma_deconv.SetNpx(1000)
    fS1Gamma_deconv.SetNpx(1000)
    fS2Gamma_deconv.SetNpx(1000)
    fFidGamma_deconv.Write()
    fS1Gamma_deconv.Write()
    fS2Gamma_deconv.Write()
    out_file_Gamma_deconv.Close()


bolo_name = "FID837"
data_dir  = "../Fond_ERA_merged/"
tree_name = "data"
prospective_analysis(bolo_name, data_dir, tree_name)
# compare_to_simu(bolo_name, data_dir, tree_name)
# get_fid_cut_efficiency(bolo_name, data_dir, tree_name)
# create_extrapolated_gamma(bolo_name, data_dir, tree_name)