#!/usr/bin/env python

from ROOT import *
import script_utils as script_utils
import PyROOTPlots as PyRPl


def get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)
    print cut_line

    h= TH1F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "_" + bolo_name,channel_X, cut_line)
    
    h.SetStats(0)

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    gPad.SetLogy()
    h.Draw()   
    print h.Integral()
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")

def get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)

    h= TH2F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    tree.Project(h_name + "_" + bolo_name,channel_Y+":"+channel_X, cut_line)
    
    h.SetStats(0)

    # WIMP_path = "../Scripts_bckg_ERA/ROOT_files/" + bolo_name + "/WIMP/"
    # tWIMP, fWIMp = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_25_tree.root", "t_new")

    # hWIMP= TH2F("hWIMP", "hWIMP", bin_X, min_X, max_X, bin_Y, min_Y, max_Y)    
    # tWIMP.Project("hWIMP","0.5*(EIB+EID):0.5*(EC1+EC2)")
    # hWIMP.SetMarkerColor(kRed)
    # hWIMP.SetMarkerStyle(1)

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    # h.SetMarkerStyle(20)
    h.Draw("")   
    if "vetocut" in fig_title:
        # hWIMP.Draw("same")
        h.Draw("same")

    raw_input()
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")

def get_2D_plot_col(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)

    h= TH2F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    tree.Project(h_name + "_" + bolo_name,channel_Y+":"+channel_X, cut_line)
    
    h.SetStats(0)

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    h.SetContour(50)
    h.Draw("cont4z")   
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")


def get_1D_plot_with_thresh(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)

    h= TH1F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "_" + bolo_name,channel_X, cut_line)
    
    # h.SetStats(0)

    PyRPl.process_TH1(h, X_title = channel_X_title, Y_title = channel_Y_title + " (counts/" + str(h.GetBinWidth(2))[:4] + " days)")

    #Get the threshold hist
    f_thresh = TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_thresh_day.root", "read")
    hkt= f_thresh.Get("hkth")
    rightmax = 3.5*hkt.GetMaximum()
    sep=506

    #Do the plots
    cc=TCanvas("cc","cc")
    h.Draw()   
    cc.Update()

    # gStyle.SetOptTitle(1)

    #Create new axis
    axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(),gPad.GetUymax(), 0,rightmax,sep,"+L")
    axis.SetTitle("Threshold (keV)")
    axis.SetTitleColor(kRed)
    axis.SetLineColor(kRed)
    axis.CenterTitle(kTRUE)

    #Scale hkt to the pad coordinates
    scale= gPad.GetUymax()/rightmax
    hkt.SetLineColor(kRed)
    hkt.SetMarkerColor(kRed)
    hkt.Scale(scale)

    #P option to draw the markers only
    hkt.Draw("Psame")
    axis.Draw("same")
    cc.Update()
    raw_input()
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_with_threshold_lowmass.eps")


def get_1D_plot_normalised(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts, norm):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)
    print cut_line

    h = TH1F("Heat only spectrum", "Heat only spectrum", bin_X, min_X, max_X)

    hplot= TH1F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "_" + bolo_name,channel_X, cut_line)
    
    h.SetStats(0)
    hplot.Scale(1/(norm*h.GetBinWidth(10)/(24*3600.)))

    script_utils.print_utility(script_utils.COL("NNNNNNNOOOOOOOORRRRRRRRRMMMMMMM: " + str(norm/(24*3600.)), "fail"))

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    cc.SetLogy()
    h.SetMaximum(200)
    h.SetMinimum(0.1)
    h.Draw()
    hplot.Draw("same")   
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")

def get_time_series_correlation_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, bin_time, min_time, max_time, list_cuts_X, list_cuts_Y):
    

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line_X='&&'.join(list_cuts_X)
    cut_line_Y='&&'.join(list_cuts_Y)

    #First  create the heat only plot
    h_heat= TH1F("heatonly", "heatonly", bin_time, min_time, max_time)
    tree.Project("heatonly","1E6*UT1+UT2", cut_line_X)

    # h_heat.Draw()
    # raw_input()

    #Then create the plot for KTH, FWHM etc. They will be averaged over the bin time period
    h_2D= TH2F("h_2D", "h_2D", bin_time, min_time, max_time, bin_Y, min_Y, max_Y)
    tree.Project("h_2D", channel_Y+":1E6*UT1+UT2", cut_line_Y)
    list_avg=[]
    for i in range(1, bin_time +1):
        counter=0
        sum_thresh=0
        for k in range(1, bin_Y+1):
            if h_2D.GetBinContent(i,k) != 0:
                counter+=1
                sum_thresh+=float(h_2D.GetYaxis().GetBinCenter(k))
        if counter!=0:
            list_avg.append(sum_thresh/float(counter))
        else:
            list_avg.append(0)



    # #First create a 2D time-channel_X /time-channel_Y histogram, then project it
    # h_X= TH2F("h_X", "h_X", bin_time, min_time, max_time, bin_X, min_X, max_X)
    # tree.Project("h_X", channel_X+":1E6*UT1+UT2", cut_line_X)

    # h_Y= TH2F("h_Y", "h_Y", bin_time, min_time, max_time, bin_Y, min_Y, max_Y)
    # tree.Project("h_Y", channel_Y+":1E6*UT1+UT2", cut_line_Y)
    
    # hprojX    = TH1D("hprojX","hprojX",bin_time, min_time, max_time)
    # hprojX    = h_X.ProjectionY()
    # hprojX.SetName("hprojX")

    # hprojY    = TH1D("hprojY","hprojY",bin_time, min_time, max_time)
    # hprojY    = h_Y.ProjectionY()
    # hprojY.SetName("hprojY")     


    h_corr = TH2F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    for k in range(1, bin_time +1):
        x_entry = h_heat.GetBinContent(k)
        y_entry = list_avg[k-1]
        h_corr.Fill(x_entry,y_entry)

    h_corr.SetStats(0)

    h_corr.GetXaxis().SetTitle(channel_X_title)
    h_corr.GetXaxis().CenterTitle(kTRUE)
    h_corr.GetXaxis().SetTitleSize(0.06)
    h_corr.GetXaxis().SetTitleOffset(0.8)

    h_corr.GetYaxis().SetTitle(channel_Y_title)
    h_corr.GetYaxis().CenterTitle(kTRUE)
    h_corr.GetYaxis().SetTitleSize(0.06)
    h_corr.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    h_corr.Draw()   
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")

def get_time_bolo_to_bolo(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, bin_time, min_time, max_time, list_cuts_X, list_cuts_Y):
    

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line_X='&&'.join(list_cuts_X)
    cut_line_Y='&&'.join(list_cuts_Y)

    #First  create the heat only plot
    h_heat= TH1F("heatonly", "heatonly", bin_time, min_time, max_time)
    tree.Project("heatonly","1E6*UT1+UT2", cut_line_X)

    # h_heat.Draw()
    # raw_input()

    #Then create the plot for KTH, FWHM etc. They will be averaged over the bin time period
    h_2D= TH2F("h_2D", "h_2D", bin_time, min_time, max_time, bin_Y, min_Y, max_Y)
    tree.Project("h_2D", channel_Y+":1E6*UT1+UT2", cut_line_Y)
    list_avg=[]
    for i in range(1, bin_time +1):
        counter=0
        sum_thresh=0
        for k in range(1, bin_Y+1):
            if h_2D.GetBinContent(i,k) != 0:
                counter+=1
                sum_thresh+=float(h_2D.GetYaxis().GetBinCenter(k))
        if counter!=0:
            list_avg.append(sum_thresh/float(counter))
        else:
            list_avg.append(0)



    # #First create a 2D time-channel_X /time-channel_Y histogram, then project it
    # h_X= TH2F("h_X", "h_X", bin_time, min_time, max_time, bin_X, min_X, max_X)
    # tree.Project("h_X", channel_X+":1E6*UT1+UT2", cut_line_X)

    # h_Y= TH2F("h_Y", "h_Y", bin_time, min_time, max_time, bin_Y, min_Y, max_Y)
    # tree.Project("h_Y", channel_Y+":1E6*UT1+UT2", cut_line_Y)
    
    # hprojX    = TH1D("hprojX","hprojX",bin_time, min_time, max_time)
    # hprojX    = h_X.ProjectionY()
    # hprojX.SetName("hprojX")

    # hprojY    = TH1D("hprojY","hprojY",bin_time, min_time, max_time)
    # hprojY    = h_Y.ProjectionY()
    # hprojY.SetName("hprojY")     


    h_corr = TH2F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    for k in range(1, bin_time +1):
        x_entry = h_heat.GetBinContent(k)
        y_entry = list_avg[k-1]
        h_corr.Fill(x_entry,y_entry)

    h_corr.SetStats(0)

    h_corr.GetXaxis().SetTitle(channel_X_title)
    h_corr.GetXaxis().CenterTitle(kTRUE)
    h_corr.GetXaxis().SetTitleSize(0.06)
    h_corr.GetXaxis().SetTitleOffset(0.8)

    h_corr.GetYaxis().SetTitle(channel_Y_title)
    h_corr.GetYaxis().CenterTitle(kTRUE)
    h_corr.GetYaxis().SetTitleSize(0.06)
    h_corr.GetYaxis().SetTitleOffset(0.8)

    #Do the plots
    cc=TCanvas("cc","cc")
    h_corr.Draw()   
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")





def get_1D_plot_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)
    print cut_line

    h= TH1F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "_" + bolo_name,channel_X, cut_line)
    
    h.SetStats(0)

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)

    # h.GetXaxis().SetTimeFormat("%d\/%m\/%y%F1969-12-31 23:00:00");  #(PARIS)
    h.GetXaxis().SetTimeDisplay(1)
    h.Draw()
    

    # raw_input()

    #Do the plots
    cc=TCanvas("cc","cc")
    # cc.SetLogy()
    h.Draw()   
    print h.Integral()
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")

def get_1D_plot_with_thresh_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)

    h = TH1F("Heatonly rate", "Heatonly rate", bin_X, min_X, max_X)

    hplot= TH1F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "_" + bolo_name,channel_X, cut_line)
    
    h.SetStats(0)

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)
    # h.GetXaxis().SetTimeDisplay(1)

    #Get the threshold hist
    f_thresh = TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_thresh.root", "read")
    hkt= f_thresh.Get("hkth")
    rightmax = 0.5*hkt.GetMaximum()
    sep=506

    #Do the plots
    cc=TCanvas("cc","cc")
    
    h.SetMinimum(1)
    h.SetMaximum(200)
    h.Draw()   
    # gPad.SetLogy()
    hplot.Draw("same")
    cc.Update()

    gStyle.SetOptTitle(0)

    #Create new axis
    axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(),gPad.GetUymax(), 0,rightmax,sep,"+L")
    axis.SetTitle("Threshold (keV)")
    axis.SetTitleColor(kRed)
    axis.SetLineColor(kRed)
    axis.CenterTitle(kTRUE)

    #Scale hkt to the pad coordinates
    scale= gPad.GetUymax()/rightmax
    hkt.SetLineColor(kRed)
    hkt.SetMarkerColor(kRed)
    hkt.Scale(scale)
    #P option to draw the markers only
    hkt.Draw("Psame")
    axis.Draw("same")
    cc.Update()
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/') 
    # raw_input() 
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_with_threshold_time_lowmass.png")


def get_2D_plot_col_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line='&&'.join(list_cuts)

    h= TH2F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X, bin_Y, min_Y, max_Y)
    tree.Project(h_name + "_" + bolo_name,channel_Y+":"+channel_X, cut_line)
    
    h.SetStats(0)

    h.GetXaxis().SetTitle(channel_X_title)
    h.GetXaxis().CenterTitle(kTRUE)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.8)

    h.GetYaxis().SetTitle(channel_Y_title)
    h.GetYaxis().CenterTitle(kTRUE)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.8)
    h.GetXaxis().SetTimeDisplay(1)

    #Do the plots
    cc=TCanvas("cc","cc")
    h.SetContour(99)
    h.Draw("")   
    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")

def compare_fid_distri_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts1, list_cuts2):


    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)

    cut_line1='&&'.join(list_cuts1)
    print cut_line1

    cut_line2='&&'.join(list_cuts2)
    print cut_line2

    h1= TH1F(h_name + "1_" + bolo_name, h_name + "1_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "1_" + bolo_name,channel_X, cut_line1)

    h2= TH1F(h_name + "2_" + bolo_name, h_name + "2_" + bolo_name, bin_X, min_X, max_X)
    tree.Project(h_name + "2_" + bolo_name,channel_X, cut_line2)
    
    h1.SetStats(0)

    h1.GetXaxis().SetTitle(channel_X_title)
    h1.GetXaxis().CenterTitle(kTRUE)
    h1.GetXaxis().SetTitleSize(0.06)
    h1.GetXaxis().SetTitleOffset(0.8)

    h1.GetYaxis().SetTitle(channel_Y_title)
    h1.GetYaxis().CenterTitle(kTRUE)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetYaxis().SetTitleOffset(0.8)

    h1.SetLineColor(kBlack)
    h2.SetLineColor(kRed)

    h1.SetMaximum(3*h1.GetMaximum())
    h2.Scale(h1.Integral()/float(h2.Integral()))
    #Do the plots
    cc=TCanvas("cc","cc")
    gPad.SetLogy()
    h1.Draw()   
    h2.Draw("same")

    # Define path for the .txt file. Create the directory if it does not exist, then open file
    figure_path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Figures/')  
    cc.Print(figure_path_name + bolo_name + "_" + fig_title + "_lowmass.png")