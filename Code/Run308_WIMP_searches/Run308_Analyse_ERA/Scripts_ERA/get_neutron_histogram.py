#! /usr/bin/env python

import ROOT
import glob
import numpy as np
from ROOT import *
import PyROOTPlots as PyRPl

def get_neutron_simu_histogram(bolo_name, ECinf):

    """Build neutron simu histogram from MdJ data
    
    Detail:

    Args:
        bolo_name = (str) bolometer name
        ECinf = (float) min heat energy in analysis

    Returns:
        void

    Raises:
        void
    """

    datadir="../Neutron_simu"
    list_files=glob.glob(datadir+"/*-Run308.root")

    distris=[]
    names=[]
    for i,f in enumerate(list_files):
        ff=ROOT.TFile(f,"READ")
        names.append(f[len(datadir)+1:-5])
        h=ff.Get("htot")
        if i==0 : 
            nb=h.GetNbinsX()
            xvals=np.asarray([h.GetBinCenter(p+1) for p in range(nb)])
        else :
            if h.GetNbinsX()!=nb : print "Pbl",h.GetNbinsX(),f
            if h.GetBinCenter(1)!=xvals[0] : print "Pbl",h.GetBinCenter(1),f

        yvals=np.asarray([h.GetBinContent(p+1) for p in range(nb)])
        distris.append(yvals)

    vv=np.asarray(distris)
    fulldistri=np.asarray([np.sum(vv[:,p]) for p in range(len(xvals))])

    # print xvals 
    # print fulldistri

    conv_path    = "/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/conv_EE_to_NR.root"
    func_conv, fileconv = PyRPl.open_ROOT_object(conv_path, "conv")

    hneut = TH1F("hneut", "hneut",200, 0, 20)
    min_bin = hneut.FindBin(func_conv.Eval(ECinf))
    for i in range(200):
        if i>=min_bin:
            hneut.SetBinContent(i, fulldistri[i])
        else:
            hneut.SetBinContent(i,0)
    # hneut.Draw()
    # raw_input()
    f = TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_neutron_simu_hist.root", "recreate")
    hneut.Write()
    f.Close()

bolo_name = "FID837"
ECinf = 0.5
get_neutron_simu_histogram(bolo_name, ECinf)