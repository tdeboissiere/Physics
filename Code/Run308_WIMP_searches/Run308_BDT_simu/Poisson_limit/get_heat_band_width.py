from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import script_utils as script_utils
import math
from array import array


def get_heat_band_width_gauss(bolo_name):

    """Get width of heat band
    
    Detail:

    Args:
        bolo_name = (str) bolometer name

    Returns:
        void

    Raises:
        void
    """

    path_name = "../../Poisson_" + bolo_name + "/True_events/ROOT_files/"
    t, f = PyRPl.open_ROOT_object(path_name + bolo_name + "_500_days_true_events_tree.root", "t_new0")

    hist = TH1F("hist", "hist", 500, -3, 3)

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        EIB, EID = t.EIB, t.EID 
        EC1, EC2 = t.EC1, t.EC2
        EIA, EIC = t.EIA, t.EIC

        EI = 0.5*(EIB+EID)
        EC = 0.5*(EC1+EC2)
        ENR = t.ENR

        if (0.5*(EIA+EIB+EIC+EID)<2 and abs(EC1-EC2)<2 and ENR>4):
            dist = EI
            hist.Fill(dist)

    hist.Draw()
    hist.Fit("gaus", "", "", -3,3)
    # f = hist.GetFunction("gaus")
    # print f.GetParameter(0)
    # print f.GetParameter(1)
    # print f.GetParameter(2)

    raw_input()

def get_heat_band_width(bolo_name):

    """Get width of gamma band
    
    Detail:
        works in ENR plane

    Args:
        bolo_name = (str) bolometer name

    Returns:
        void

    Raises:
        void
    """

    path_name = "../../Poisson_" + bolo_name + "/True_events/ROOT_files/"
    t, f = PyRPl.open_ROOT_object(path_name + bolo_name + "_500_days_true_events_tree.root", "t_new0")
    list_trees = [f.Get("t_new" + str(i)) for i in range(10,12)]

    EC1 = array("f", [0.])
    EC2 = array("f", [0.])
    EIA = array("f", [0.])
    EIB = array("f", [0.])
    EIC = array("f", [0.])
    EID = array("f", [0.])
    ENR = array("f", [0.])

    #New tree to be filled
    t_new = TTree("t_new", "t_new")

    # create the branches and assign the fill-variables to them
    t_new.Branch("EC1", EC1, "EC1/F")
    t_new.Branch("EC2", EC2, "EC2/F")
    t_new.Branch("EIA", EIA, "EIA/F")
    t_new.Branch("EIB", EIB, "EIB/F")
    t_new.Branch("EIC", EIC, "EIC/F")
    t_new.Branch("EID", EID, "EID/F")
    t_new.Branch("ENR", ENR, "ENR/F")

    for tree in list_trees:
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)

            EC1[0], EC2[0], EIA[0], EIB[0] = tree.EC1, tree.EC2, tree.EIA, tree.EIB
            ENR[0], EIC[0], EID[0] = tree.ENR, tree.EIC, tree.EID


    tree.Draw("EIB:ENR>>hist(100,0,15,100,-2,5)", "0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>0")
    tree.Draw("EIB:ENR>>hist2(100,0,15,100,-2,5)", "0.5*(EIA+EIB+EIC+EID)<2 && abs(EC1-EC2)<2 && ENR>0")
    hist.SetMarkerColor(kRed)
    hist2.Draw()
    hist.Draw("same")
    raw_input()
    tot_entries = float(tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<2 && abs(EC1-EC2)<2 && ENR>4"))
    # print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<0.4")/tot_entries

    print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<2 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<0.512")/tot_entries
    print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<2 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<0.768")/tot_entries

    # print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<0.6")/tot_entries
    # print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<0.8")/tot_entries
    # print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<1")/tot_entries
    # print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<1.2")/tot_entries
    # print tree.GetEntries("0.5*(EIA+EIB+EIC+EID)<1 && abs(EC1-EC2)<2 && ENR>4 && 0.5*(EIB+EID)<1.4")/tot_entries

    #0.5ENR -1 : 95 %
    #0.5*ENR-1.2: 99%

bolo_name = "FID837"
# get_heat_band_width_gauss(bolo_name)
get_heat_band_width(bolo_name)