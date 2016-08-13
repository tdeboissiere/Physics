#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import sys
import numpy as np
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh

def compute_resolution( FWHM1,  FWHM2) :
    """
    Compute the resolution (i.e FWHM not sigma) of the best combination of two channels
    """
    s1=float(FWHM1)/2.3548
    s2=float(FWHM2)/2.3548

    w1=(s2*s2)/(s1*s1+s2*s2)

    res=np.sqrt(pow(w1*s1,2)+pow((1-w1)*s2,2))
    return 2.3548*res

def get_FWHM_hist(bolo_name, data_dir, tree_name):
    """
    Get the FWHM hist after cuts for BDT simulations

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A list : [mean_IA, mean_IB, mean_IC, mean_ID, mean_C1, mean_C2]
    Also prints the FWHM in a Figures repository
    """

    #Get tree 
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    # Define path for the .txt file. Create the directory if it does not exist, then open file
    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    hFWIA = TH1F("hFWIA", "hFWIA", 200, 0, 2)
    hFWIB = TH1F("hFWIB", "hFWIB", 200, 0, 2)
    hFWIC = TH1F("hFWIC", "hFWIC", 200, 0, 2)
    hFWID = TH1F("hFWID", "hFWID", 200, 0, 2)
    hOWC1 = TH1F("hOWC1", "hOWC1", 200, 0, 2)
    hOWC2 = TH1F("hOWC2", "hOWC2", 200, 0, 2)

    tree.Project("hFWIA", "FWIA", bolo_TCut_line)
    tree.Project("hFWIB", "FWIB", bolo_TCut_line)
    tree.Project("hFWIC", "FWIC", bolo_TCut_line)
    tree.Project("hFWID", "FWID", bolo_TCut_line)
    tree.Project("hOWC1", "OWC1", bolo_TCut_line)
    tree.Project("hOWC2", "OWC2", bolo_TCut_line)

    fout = TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_baseline_hist.root", "recreate")
    hFWIA.Write("FWIA")
    hFWIB.Write("FWIB")
    hFWIC.Write("FWIC")
    hFWID.Write("FWID")
    hOWC1.Write("OWC1")
    hOWC2.Write("OWC2")

def get_FWHM(bolo_name, data_dir, tree_name):
    """
    Return the list of the averaged FWHM in the 6 channels (IA,B,C,D then C1,C2)

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A list : [mean_IA, mean_IB, mean_IC, mean_ID, mean_C1, mean_C2]
    Also prints the FWHM in a Figures repository
    """

    #Get tree 
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    # Define path for the .txt file. Create the directory if it does not exist, then open file
    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    
    list_h=[TH1F() for i in range(6)]      #"FWHM","FWHM" ,100,0,5)
    list_channel = ["FWIA", "FWIB", "FWIC", "FWID", "OWC1", "OWC2"]
    for k, channel in enumerate(list_channel):
        htemp=TH1F("temp","temp" ,200, 0, 2)
        # print bolo_TCut_line
        tree.Project("temp",channel,bolo_TCut_line)
        list_h[k] = htemp.Clone(channel)
        list_h[k].SetStats(0)
        list_h[k].SetTitle(channel)
        list_h[k].GetXaxis().SetTitle(channel + " (keV)")
        list_h[k].GetXaxis().SetTitleSize(0.06)
        list_h[k].GetXaxis().SetTitleOffset(0.8)
        list_h[k].SetMarkerStyle(7)
        list_h[k].SetMarkerSize(0.7)
        list_h[k].GetXaxis().SetLabelSize(0.05)
        list_h[k].GetYaxis().SetLabelSize(0.05)
        del htemp

    #Do the plots
    c_ion=TCanvas("c_ion","c_ion")
    c_ion.Divide(2,2)
    for i in range(4):
        c_ion.cd(i+1)
        list_h[i].Draw()

    c_heat=TCanvas("c_heat", "c_heat")
    c_heat.Divide(1,2)
    for i in [0,1]:
        c_heat.cd(i+1)
        list_h[i+4].Draw()

    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory("../Analyse_" + bolo_name + "/Figures/")  
    c_ion.Print(figure_path_name + bolo_name + "_FWHM_ion.eps")
    c_heat.Print(figure_path_name + bolo_name + "_FWHM_heat.eps")

    list_mean= [ str(hist.GetMean())[:5] for hist in list_h]

    #Open polar file to access the bias voltages
    path_name=script_utils.create_directory("../Analyse_" + bolo_name + "/Text_files/")    
    file_polar = script_utils.open_text_file(path_name, bolo_name + "_all_polars_with_entries.txt", "r")
    list_lines_polar = [elem.rstrip().split(",") for elem in file_polar.readlines()]
    
    #Define and Fill in the polar cuts and an array of titles for the histograms
    VFID,VET = "",""

    #WARNING: assume only one polar: give error message if not the case.
    if len(list_lines_polar)>2:
        script_utils.print_utility(script_utils.COL("Warning, more than one polar", "fail"))
        sys.exit()

    for i in range(1,len(list_lines_polar)):
        VFID = list_lines_polar[i][3]
        VET  = list_lines_polar[i][4]



    FWHM_fid=str(compute_resolution(list_mean[1],list_mean[3]))[:5]
    FWHM_S1=str(compute_resolution(list_mean[0],list_mean[1]))[:5]
    FWHM_S2=str(compute_resolution(list_mean[2],list_mean[3]))[:5]
    FWHM_heat=str(compute_resolution(list_mean[3],list_mean[4]))[:5]
    FWHM_full=str(compute_resolution(FWHM_fid,FWHM_heat))[:5]

    # Create file directory if it does not exist and define the file name
    FWHM_file_path= script_utils.create_directory("../Analyse_" + bolo_name + "/Text_files/")
    FWHM_file_name= bolo_name + "_FWHM.txt"
   
    print VFID, VET

    with open(FWHM_file_path + FWHM_file_name, "w") as FWHM_file:
        FWHM_file.write("FWIA,FWIB,FWIC,FWID,OWC1,OWC2,VFID,VET,FWFID,FWS1,FWS2,FWHEAT,FWFULL\n")
        FWHM_line=",".join(list_mean) + "," + ",".join([VFID,VET]) + "," + ",".join([FWHM_fid, FWHM_S1, FWHM_S2, FWHM_heat, FWHM_full])
        FWHM_file.write(FWHM_line)
        script_utils.print_utility(script_utils.COL("Writing line " + FWHM_line,"blue"))




def get_FWHM_time_hist_for_BDT(bolo_name, data_dir, tree_name):

    """
    Return the list of the averaged FWHM in the 6 channels (IA,B,C,D then C1,C2)

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    void
    """

    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)

    l_events     = TEventList("l_events")
    tree.Draw(">>l_events" )
    pop_len = l_events.GetN()
    list_events = []
    for k in range(pop_len):
        counter = l_events.GetEntry(k)
        tree.GetEntry(counter)
        list_events.append(tree.JOUR)

    #Get min and max time, + safety margin
    tmin, tmax = min(list_events)-5, max(list_events)+5

    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts_forduration.txt")
    print bolo_TCut_line
    
    # 2 D threshold histogram
    binning = 2000
    hist_thresh   = TH2F("hist_thresh","hist_thresh",binning,tmin, tmax,1000,0,2)
    hist_FWIA     = TH2F("hist_FWIA","hist_FWIA",binning,tmin, tmax,1000,0,2)
    hist_FWIB     = TH2F("hist_FWIB","hist_FWIB",binning,tmin, tmax,1000,0,2)    
    hist_FWIC     = TH2F("hist_FWIC","hist_FWIC",binning,tmin, tmax,1000,0,2)
    hist_FWID     = TH2F("hist_FWID","hist_FWID",binning,tmin, tmax,1000,0,2)
    hist_OWC1     = TH2F("hist_OWC1","hist_OWC1",binning,tmin, tmax,1000,0,2)
    hist_OWC2     = TH2F("hist_OWC2","hist_OWC2",binning,tmin, tmax,1000,0,2)

    tree.Project( "hist_thresh", "KTH:JOUR" , bolo_TCut_line)
    tree.Project( "hist_FWIA",  "FWIA:JOUR",bolo_TCut_line)
    tree.Project( "hist_FWIB",  "FWIB:JOUR",bolo_TCut_line)    
    tree.Project( "hist_FWIC",  "FWIC:JOUR",bolo_TCut_line)
    tree.Project( "hist_FWID",  "FWID:JOUR",bolo_TCut_line)
    tree.Project( "hist_OWC1",  "OWC1:JOUR",bolo_TCut_line)
    tree.Project( "hist_OWC2",  "OWC2:JOUR",bolo_TCut_line)
    
    hprojX_KTH    = hist_thresh.ProjectionX()

    hprojY_KTH    = hist_thresh.ProjectionY()
    hprojY_FWIA = hist_FWIA.ProjectionY()
    hprojY_FWIB = hist_FWIB.ProjectionY()    
    hprojY_FWIC = hist_FWIC.ProjectionY()
    hprojY_FWID = hist_FWID.ProjectionY()
    hprojY_OWC1 = hist_OWC1.ProjectionY()
    hprojY_OWC2 = hist_OWC2.ProjectionY()


    print hprojY_FWIA.GetNbinsX()
    print hprojY_FWIB.GetNbinsX()    
    print hprojY_FWIC.GetNbinsX()
    print hprojY_FWID.GetNbinsX()
    print hprojY_OWC1.GetNbinsX()
    print hprojY_FWIB.GetNbinsX()
    print hprojY_KTH.GetNbinsX()


    #Fill time and threshold lists used to build the 1D threshold histhograms
    list_time, list_fwia, list_fwic, list_fwid, list_fwib, list_owc1, list_owc2, list_jour = [], [], [], [], [],[], [], []
    for ix in range(1,binning) : 
        temp_fwia, temp_fwic,temp_fwib, temp_fwid, temp_owc1, temp_owc2, temp_jour=0,0,0,0,0,0,0
        for iy in range(1, hprojY_KTH.GetNbinsX()): 
            fwia = hprojY_FWIA.GetBinCenter(iy)
            fwib = hprojY_FWIB.GetBinCenter(iy)            
            fwic = hprojY_FWIC.GetBinCenter(iy)
            fwid = hprojY_FWID.GetBinCenter(iy)
            owc1 = hprojY_OWC1.GetBinCenter(iy)
            owc2 = hprojY_OWC2.GetBinCenter(iy)

            if ( hist_FWIA.GetBinContent(ix,iy)!=0) :
                temp_fwia=fwia

            if ( hist_FWIB.GetBinContent(ix,iy)!=0) :
                temp_fwib=fwib

            if ( hist_FWIC.GetBinContent(ix,iy)!=0) :
                temp_fwic=fwic

            if ( hist_FWID.GetBinContent(ix,iy)!=0) :
                temp_fwid=fwid

            if ( hist_OWC1.GetBinContent(ix,iy)!=0) :
                temp_owc1=owc1

            if ( hist_OWC2.GetBinContent(ix,iy)!=0) :
                temp_owc2=owc2

        list_time.append(hprojX_KTH.GetBinCenter(ix))
        if (temp_fwia!=0 and temp_fwib!=0 and temp_fwic!=0 and temp_fwid!=0 and temp_owc1!=0 and temp_owc2!=0):
            list_fwia.append(temp_fwia)
            list_fwib.append(temp_fwib)        
            list_fwic.append(temp_fwic)
            list_fwid.append(temp_fwid)
            list_owc1.append(temp_owc1)
            list_owc2.append(temp_owc2)
            list_jour.append(1)
        else:
            list_fwia.append(0)
            list_fwib.append(0)        
            list_fwic.append(0)
            list_fwid.append(0)
            list_owc1.append(0)
            list_owc2.append(0)
            list_jour.append(0)            


    # ROOT says there is an issue with increasing value for the bins represented by time
    # I disagree.
    # Solution: use a graph to build the 1D threshold histhogram

    gr_fwia = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwia))
    hfwia = TH1F("hfwia", "hfwia",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hfwia.SetBinContent(k,gr_fwia.Eval(list_time[k-1]))

    gr_fwib = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwib))
    hfwib = TH1F("hfwib", "hfwib",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hfwib.SetBinContent(k,gr_fwib.Eval(list_time[k-1]))


    gr_fwic = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwic))
    hfwic = TH1F("hfwic", "hfwic",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hfwic.SetBinContent(k,gr_fwic.Eval(list_time[k-1]))

    gr_fwid = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwid))
    hfwid = TH1F("hfwid", "hfwid",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hfwid.SetBinContent(k,gr_fwid.Eval(list_time[k-1]))

    gr_owc1 = TGraph(np.size(list_time), np.array(list_time), np.array(list_owc1))
    howc1 = TH1F("howc1", "howc1",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        howc1.SetBinContent(k,gr_owc1.Eval(list_time[k-1]))

    gr_owc2 = TGraph(np.size(list_time), np.array(list_time), np.array(list_owc2))
    howc2 = TH1F("howc2", "howc2",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        howc2.SetBinContent(k,gr_owc2.Eval(list_time[k-1]))

    hist_jour = TH1F("hist_jour", "hist_jour", np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        # hist_jour.SetBinContent(k,gr_jour.Eval(list_time[k-1]))
        hist_jour.SetBinContent(k,list_jour[k-1])

    # hfwia.Draw()
    # hist_jour.SetLineColor(kRed)
    # hist_jour.Draw("same")
    # raw_input()

    root_path = "../Analyse_"+ bolo_name + "/ROOT_files/" + bolo_name + "_fwhm_time_hist_for_BDT.root"
    fout = TFile(root_path, "recreate")
    hist_jour.Write()
    hfwia.Write()
    hfwib.Write()
    hfwic.Write()
    hfwid.Write()
    howc1.Write()
    howc2.Write()
    fout.Close()


def get_FWHM_time_hist_for_heatonly_for_BDT(bolo_name, data_dir, tree_name):

    """
    Return the list of the averaged FWHM in the 2 heat channels
    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    void
    """

    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)

    #Define path 
    path_name       = script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    polar_file_name = bolo_name + "_polar_start_and_end_time_prelim.txt"

    print path_name + polar_file_name
    start_and_end_array = np.loadtxt(path_name + polar_file_name ,delimiter=",", usecols=(1,2))
    script_utils.print_utility(script_utils.COL('Opening ' + path_name+polar_file_name + ' for reading', 'blue'))    

    #Define absolute min and max of time from the _polar_start_and_end_time_prelim file. 
    tmin = np.amin(start_and_end_array)
    tmax = np.amax(start_and_end_array)

    bolo_TCut_line = Ana_ut.open_cut_file(bolo_name, "TCuts_forduration.txt")
    print bolo_TCut_line
    

    # 2 D threshold histogram
    binning = 20
    days = tmax + 10 #add 10 days just to have some "margin"
    hist_thresh   = TH2F("hist_thresh","hist_thresh",binning,tmin, tmax,1000,0,2)
    hist_OWC1     = TH2F("hist_OWC1","hist_OWC1",binning,tmin, tmax,1000,0,2)
    hist_OWC2     = TH2F("hist_OWC2","hist_OWC2",binning,tmin, tmax,1000,0,2)

    tree.Project( "hist_thresh", "KTH:1E6*UT1+UT2" , bolo_TCut_line)
    tree.Project( "hist_OWC1",  "OWC1:1E6*UT1+UT2",bolo_TCut_line)
    tree.Project( "hist_OWC2",  "OWC2:1E6*UT1+UT2",bolo_TCut_line)
    
    hprojX_KTH    = hist_thresh.ProjectionX()
    hprojY_KTH    = hist_thresh.ProjectionY()
    hprojY_OWC1 = hist_OWC1.ProjectionY()
    hprojY_OWC2 = hist_OWC2.ProjectionY()

    print hprojY_OWC1.GetNbinsX()
    print hprojY_OWC2.GetNbinsX()
    print hprojY_KTH.GetNbinsX()


    #Fill time and threshold lists used to build the 1D threshold histhograms
    list_time, list_owc1, list_owc2, list_jour = [], [], [], []
    for ix in range(1,binning) : 
        temp_owc1, temp_owc2, temp_jour=0,0,0
        for iy in range(1, hprojY_KTH.GetNbinsX()): 
            owc1 = hprojY_OWC1.GetBinCenter(iy)
            owc2 = hprojY_OWC2.GetBinCenter(iy)

            if ( hist_OWC1.GetBinContent(ix,iy)!=0) :
                temp_owc1=owc1

            if ( hist_OWC2.GetBinContent(ix,iy)!=0) :
                temp_owc2=owc2

        list_time.append(hprojX_KTH.GetBinCenter(ix))
        if (temp_owc1!=0 and temp_owc2!=0):
            list_owc1.append(temp_owc1)
            list_owc2.append(temp_owc2)
            list_jour.append(1)
        else:
            list_owc1.append(0)
            list_owc2.append(0)
            list_jour.append(0)            


    # ROOT says there is an issue with increasing value for the bins represented by time
    # I disagree.
    # Solution: use a graph to build the 1D threshold histhogram

    gr_owc1 = TGraph(np.size(list_time), np.array(list_time), np.array(list_owc1))
    howc1 = TH1F("howc1", "howc1",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        howc1.SetBinContent(k,gr_owc1.Eval(list_time[k-1]))

    gr_owc2 = TGraph(np.size(list_time), np.array(list_time), np.array(list_owc2))
    howc2 = TH1F("howc2", "howc2",np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        howc2.SetBinContent(k,gr_owc2.Eval(list_time[k-1]))

    hist_jour = TH1F("hist_jour", "hist_jour", np.size(list_time),min(list_time), max(list_time))
    for k in range(1, np.size(list_time)+1):
        hist_jour.SetBinContent(k,list_jour[k-1])

    root_path = "../Analyse_"+ bolo_name + "/ROOT_files/" + bolo_name + "_fwhm_time_hist_for_heatonly_for_BDT.root"
    fout = TFile(root_path, "recreate")
    hist_jour.Write()
    howc1.Write()
    howc2.Write()
    fout.Close()


def correct_FWHM(bolo_name, data_dir, tree_name):
    """
    Fit 10 keV to have correct estimate of baseline

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A list : [mean_IA, mean_IB, mean_IC, mean_ID, mean_C1, mean_C2]
    Also prints the FWHM in a Figures repository
    """

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    #Load standard cuts
    standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

    #Load estimators
    d_est = BDT_fh.open_estimator_file(bolo_name)

    #Load FWHM
    d_std = BDT_fh.open_true_event_FWHM_file(bolo_name)
    for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
        d_std[key] = str(2.7*d_std[key])

    l_FidGamma = TEventList("l_FidGamma")

    list_channel = ["C1", "C2", "IA", "IB", "IC", "ID"]
    list_hist = [TH1F("h" + channel, "h" + channel, 200,0,15) for channel in list_channel]
    hcomb = TH1F("hcomb", "hcomb", 200,0,15)
    list_FWHM = []

    print "Standard cuts are:  " , standard_cuts

    ##################################
    #    G A M M A   E V E N T S
    ##################################
    #Fiducial gammas
    tree.Draw(">>l_FidGamma",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&& abs(EC1-EC2)<2 &&" +  d_est["Q_FID"]  + ">0.7")
    pop_len = l_FidGamma.GetN()
    for k in range(pop_len):
        counter = l_FidGamma.GetEntry(k)
        tree.GetEntry(counter)
        list_val =[tree.EC1, tree.EC2, tree.EIB, tree.EID]
        list_hist[0].Fill(list_val[0])
        list_hist[1].Fill(list_val[1])
        list_hist[3].Fill(list_val[2])
        list_hist[5].Fill(list_val[3])
        hcomb.Fill(float(d_est["HEAT"][:5])*tree.EC1 + (1- float(d_est["HEAT"][:5]))*tree.EC2)

    #S1 gammas
    tree.Draw(">>l_GammaS1",standard_cuts +  "&&" + d_est["HEAT"] + ">0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&& abs(EC1-EC2)<2 &&" + d_est["Q_S1"] + ">0.65")
    pop_len = l_GammaS1.GetN()
    for k in range(pop_len):
        counter = l_GammaS1.GetEntry(k)
        tree.GetEntry(counter)
        list_val =[tree.EIA]
        list_hist[2].Fill(list_val[0])

    #S2 gammas
    tree.Draw(">>l_GammaS2",standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&& abs(EC1-EC2)<2 &&" + d_est["Q_S2"] + ">0.65")
    pop_len = l_GammaS2.GetN()
    for k in range(pop_len):
        counter = l_GammaS2.GetEntry(k)
        tree.GetEntry(counter)
        list_val =[tree.EIC]
        list_hist[4].Fill(list_val[0])


    class Gamma:
        def __call__( self, x, par ):
            peak_10_4 = par[0]*TMath.Gaus(x[0], par[1], par[2])
            return peak_10_4

    class Full_Gamma:
        def __call__( self, x, par ):
            peak_10_4 = par[0]*TMath.Gaus(x[0], par[1]*10.37, par[2])
            peak_9_66 = 0.1*par[0]*TMath.Gaus(x[0], par[1]*9.66, par[2])
            peak_8_98 = par[3]*TMath.Gaus(x[0], par[1]*8.98, par[2])
            return peak_10_4 + peak_8_98 + peak_9_66 +1

    # Call FidGamma function to get derivative
    fFidGamma = TF1( "FidGamma", Gamma(), 0, 14., 3 )
    fFidGamma.SetNpx(500)
    fFidGamma.SetParName(0,"A_10.4")
    fFidGamma.SetParName(1,"mu_10.4")
    fFidGamma.SetParName(2,"s_10.4")
    fFidGamma.SetParameters(1,10.37,0.3)
    fFidGamma.SetParLimits(2,0,1)
    fFidGamma.SetParLimits(1,10,12)
    # fFidGamma.FixParameter(1,10.37)

    # Call FidGamma function to get derivative
    fFidGamma_full = TF1( "FidGamma_full", Full_Gamma(), 0, 14., 4 )
    fFidGamma_full.SetNpx(500)
    fFidGamma_full.SetParName(0,"A_10.4")
    fFidGamma_full.SetParName(1,"calib_corr")
    fFidGamma_full.SetParName(2,"resolution")
    fFidGamma_full.SetParName(3,"A_8.98")
    fFidGamma_full.SetParameters(20,1,0.3,5)
    fFidGamma_full.SetParLimits(1,0.5,1.5)
    fFidGamma_full.SetParLimits(2,0.1,1)
    fFidGamma.FixParameter(1,1)


    c_ion=TCanvas("c_ion","c_ion")
    c_ion.Divide(2,2)
    for i in range(2,6):
        c_ion.cd(i-1)
        list_hist[i].Draw()
        list_hist[i].Fit("FidGamma_full", "LL", "", 8,12)
        print i, fFidGamma_full.GetParameter(1)*10.37
        list_FWHM.append(str(fFidGamma_full.GetParameter(2)*2.3548))

    c_heat=TCanvas("c_heat", "c_heat")
    c_heat.Divide(2,2)
    for i in range(2):
        c_heat.cd(i+1)
        list_hist[i].Draw()
        list_hist[i].Fit("FidGamma_full", "LL", "", 8,12)
        print i, fFidGamma_full.GetParameter(1)*10.37
        list_FWHM.append(str(fFidGamma_full.GetParameter(2)*2.3548))
    c_heat.cd(3)
    hcomb.Draw()
    hcomb.Fit("FidGamma_full", "LL", "", 8,12)
    res_comb = str(fFidGamma_full.GetParameter(2)*2.3548)

    # Create file directory if it does not exist and define the file name
    FWHM_file_path= script_utils.create_directory("../Analyse_" + bolo_name + "/Text_files/")
    FWHM_file_name= bolo_name + "_FWHM.txt"

    #First read the lines 
    list_lines = []
    with open(FWHM_file_path + FWHM_file_name, "r") as FWHM_file:
        list_lines = FWHM_file.readlines()

    print "old FWHM:", list_lines[1].rstrip().split(",")[:6]
    print "new FWHM:", list_FWHM
    print "res comb heat:", res_comb

    raw_input()
    sys.exit()

    list_FWHM = raw_input(" Tweak the fit until satisfied. Then Enter new FWHM (FWIA,FWIB,FWIC,FWID,OWC1,OWC2):   ")
    list_FWHM = list_FWHM.split(",")

    FWHM_fid=str(compute_resolution(list_FWHM[1],list_FWHM[3]))[:5]
    FWHM_S1=str(compute_resolution(list_FWHM[0],list_FWHM[1]))[:5]
    FWHM_S2=str(compute_resolution(list_FWHM[2],list_FWHM[3]))[:5]
    FWHM_heat=str(compute_resolution(list_FWHM[3],list_FWHM[4]))[:5]
    FWHM_full=str(compute_resolution(FWHM_fid,FWHM_heat))[:5]

    #Second update file
    with open(FWHM_file_path + FWHM_file_name, "w") as FWHM_file:
        FWHM_file.write(list_lines[0])
        new_line = list_lines[1].rstrip().split(",")
        new_line[0],new_line[1],new_line[2],new_line[3] = list_FWHM[0],list_FWHM[1],list_FWHM[2],list_FWHM[3]
        new_line[4],new_line[5],new_line[8],new_line[9], new_line[10] = list_FWHM[4],list_FWHM[5],FWHM_fid,FWHM_S1,FWHM_S2
        new_line[11],new_line[12] = FWHM_heat,FWHM_full
        new_line=",".join(new_line) 
        FWHM_file.write(new_line)
        script_utils.print_utility(script_utils.COL("Updating line to " + new_line,"blue"))


