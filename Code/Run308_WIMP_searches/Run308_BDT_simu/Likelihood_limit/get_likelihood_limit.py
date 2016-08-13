#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ROOT import *
import script_utils as script_utils
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 



def get_likelihood_limit(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):


    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"

    #############################################
    #SGet all histograms
    ##########################################

    fS1Pb     = TFile(gen_path +  "S1Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Pb     = TFile(gen_path +  "S2Pb" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Beta   = TFile(gen_path +  "S1Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Beta   = TFile(gen_path +  "S2Beta" + "_mass_" + str(mass) + "_tree.root", "read")
    fS1Gamma  = TFile(gen_path +  "S1Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fS2Gamma  = TFile(gen_path +  "S2Gamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fFidGamma = TFile(gen_path +  "FidGamma" + "_mass_" + str(mass) + "_tree.root", "read")
    fheat     = TFile(gen_path +  "heatonly" + "_mass_" + str(mass) + "_tree.root", "read")
    fWIMP     = TFile(gen_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root", "read")
    fdata     = TFile(gen_path +  "real_data" + "_mass_" + str(mass) + "_tree.root", "read")
    
    tS1Pb     = fS1Pb.Get("tout")
    tS2Pb     = fS2Pb.Get("tout") 
    tS1Beta   = fS1Beta.Get("tout") 
    tS2Beta   = fS2Beta.Get("tout") 
    tS1Gamma  = fS1Gamma.Get("tout") 
    tS2Gamma  = fS2Gamma.Get("tout") 
    tFidGamma = fFidGamma.Get("tout") 
    theat     = fheat.Get("tout")
    tWIMP     = fWIMP.Get("tout") 
    tdata     = fdata.Get("tout") 
    
    hS1Pb     = TH1F("hS1Pb", "hS1Pb", bin_X, min_X, max_X)
    hS2Pb     = TH1F("hS2Pb", "hS2Pb", bin_X, min_X, max_X)
    hS1Beta   = TH1F("hS1Beta", "hS1Beta", bin_X, min_X, max_X)
    hS2Beta   = TH1F("hS2Beta", "hS2Beta", bin_X, min_X, max_X)
    hS1Gamma  = TH1F("hS1Gamma", "hS1Gamma", bin_X, min_X, max_X)
    hS2Gamma  = TH1F("hS2Gamma", "hS2Gamma", bin_X, min_X, max_X)
    hFidGamma = TH1F("hFidGamma", "hFidGamma", bin_X, min_X, max_X)
    hheat     = TH1F("hheat", "hheat", bin_X, min_X, max_X)
    hWIMP     = TH1F("hWIMP", "hWIMP", bin_X, min_X, max_X)
    hdata     = TH1F("hdata", "hdata", bin_X, min_X, max_X)

    tS1Pb.Project("hS1Pb", "NN")
    tS2Pb.Project("hS2Pb", "NN")
    tS1Beta.Project("hS1Beta", "NN")
    tS2Beta.Project("hS2Beta", "NN")
    tS1Gamma.Project("hS1Gamma", "NN")
    tS2Gamma.Project("hS2Gamma", "NN")
    tFidGamma.Project("hFidGamma", "NN")
    theat.Project("hheat", "NN")
    tWIMP.Project("hWIMP", "NN")
    tdata.Project("hdata", "NN")


    #Rescale WIMP  to hdata
    hWIMP.Scale(1./float(hWIMP.Integral()))

    # print d_scaling

    hheat.Scale(float(d_scaling["prop_heatonly"])*float(d_scaling["exp_per_day"])*exposure/float(hheat.Integral()))
    hFidGamma.Scale(float(d_scaling["prop_FidGamma"])*float(d_scaling["exp_per_day"])*exposure/float(hFidGamma.Integral()))
    hS1Gamma.Scale(float(d_scaling["prop_S1Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Gamma.Integral()))
    hS2Gamma.Scale(float(d_scaling["prop_S2Gamma"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Gamma.Integral()))
    hS1Beta.Scale(float(d_scaling["prop_S1Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Beta.Integral()))
    hS2Beta.Scale(float(d_scaling["prop_S2Beta"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Beta.Integral()))
    hS1Pb.Scale(float(d_scaling["prop_S1Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS1Pb.Integral()))
    hS2Pb.Scale(float(d_scaling["prop_S2Pb"])*float(d_scaling["exp_per_day"])*exposure/float(hS2Pb.Integral()))

    PyRPl.process_TH1(hS1Pb, h_title = "hS1Pb", use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(hS2Pb, h_title = "hS2Pb", use_fill_bool = True, color = kOrange-8)
    PyRPl.process_TH1(hS1Beta, h_title = "hS1Beta", use_fill_bool = True, color = kGreen+2)
    PyRPl.process_TH1(hS2Beta, h_title = "hS2Beta", use_fill_bool = True, color = kGreen-3)
    PyRPl.process_TH1(hS1Gamma, h_title = "hS1Gamma", use_fill_bool = True, color = kBlue-7)
    PyRPl.process_TH1(hS2Gamma, h_title = "hS2Gamma", use_fill_bool = True, color = kBlue)
    PyRPl.process_TH1(hFidGamma, h_title = "hFidGamma", use_fill_bool = True, color = kAzure+10)
    PyRPl.process_TH1(hheat, h_title = "hheat", use_fill_bool = True, color = kRed)
    PyRPl.process_TH1(hWIMP, h_title = "hWIMP", use_fill_bool = True, color = kGray)



    list_hist_bckg =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat]

    hbckg=TH1F("hbckg","hbckg", bin_X, min_X, max_X)
    for i in range(1,bin_X+1):
        bckgcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
        hbckg.SetBinContent(i, bckgcontent)


    class likelihood_hist:
        def __call__( self, x, par ):
            likelihood =0
            #bin for 0 : 63 
            #bin for 0.1 66
            #bin for 0.2: 69
            for i in range(1,101):
                N_expected =1.1*hbckg.GetBinContent(i)+x[0]*hWIMP.GetBinContent(i)
                N_obs      = hdata.GetBinContent(i)
                if N_expected>0:
                    likelihood +=-N_expected+N_obs*math.log(N_expected)
            return -(likelihood + par[0])

    f1 = TF1("r", likelihood_hist(), 0,200, 1)
    f1.SetParameter(0,0)
    # f1.Draw("")
    min_x = f1.GetMinimumX()
    lim_x = f1.GetX(f1.Eval(min_x) + 1.64, min_x, min_x+100)

    print f1.GetMinimumX(), lim_x


    # f2 = TF2("r", likelihood_hist(), 0,200, 1000, 3000, 1)
    # f2.SetParameter(0,0)
    # f2.SetContour(99)
    # f2.Draw("cont4z")
    # x,y=c_double(0.), c_double(0.)
    # f2.GetMinimumXY(x,y)
    # print x.value, y.value

    # hS1Pb.Add(hS2Pb)
    # hS1Beta.Add(hS2Beta)
    # hFidGamma.Add(hS1Gamma)
    # hFidGamma.Add(hS2Gamma)
    # hWIMP.Scale(lim_x/hWIMP.Integral())
    # list_hist =[hS1Pb, hS1Beta, hFidGamma, hheat, hWIMP]

    # hs=THStack("hs", "hs")

    # for hist in list_hist:
    #     hs.Add(hist)

    # cc = TCanvas("cc", "cc")
    # h1=TH1F("h1","h1", bin_X, min_X, max_X)
    # PyRPl.process_TH1(h1, X_title="BDT ouput", min_Y = 1E-1, max_Y = 20000)
 

    # gPad.SetLogy()
    # h1.Draw()
    # hs.Draw("same")
    # hdata.SetMarkerStyle(1)
    # hdata.Draw("sameE1")

    # leg = TLegend(0.14,0.50,0.33,0.87)
    # leg.AddEntry(hdata.GetName(), "Data", "leg")
    # leg.AddEntry(hS1Pb.GetName(),"Lead" ,"f")
    # leg.AddEntry(hS1Beta.GetName(),"Beta", "f")
    # leg.AddEntry(hFidGamma.GetName(),"Gamma", "f")
    # leg.AddEntry(hheat.GetName(),"Heat-only", "f")
    # leg.AddEntry(hWIMP.GetName(),"WIMP " + str(mass) + " GeV","f")
    # leg.SetFillColor(kWhite)
    # leg.SetBorderSize(0)
    # leg.Draw("same")

    # raw_input()

    return lim_x

bolo_name           = "FID837"
mass                = 5
# simu_number         = 1
analysis_type       = "ana_0.5_0_5"
bin_X, min_X, max_X = 100, -2, 1.2
exposure = 65.
d_cut       = {"ECinf": 0.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}
# d_lim={}

# with open("./Text_files/" + bolo_name + "_" + analysis_type + "_NWIMP90.txt", "r") as f:
#     stuff = f.readlines()
#     stuff = [elem.rstrip().split(",") for elem in stuff]
#     for m, lim in stuff:
#         d_lim[int(m)] = float(lim) 

mass_list=[3,4,5,6,7,10,25]

with open("./Text_files/" + bolo_name + "_" + analysis_type + "_NWIMP90.txt", "w") as f:
    for mass in mass_list: 
        result = get_likelihood_limit(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure)
        f.write(str(mass)+","+str(result)+"\n")