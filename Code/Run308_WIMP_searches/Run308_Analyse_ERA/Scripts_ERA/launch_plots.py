#!/usr/bin/env python

import script_utils as script_utils
import get_plots as get_plots
import sys,os
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh


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


	# #Plot CHIA as a function of heat
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "CHIA_vs_EC1", "CHIA_ERA_vs_EC1", "EC1 (keV)","CHIA", "EC1",200,  0,15, "CHIA", 200, 0, 5
	# list_cuts=[]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	# #Plot CHIB as a function of heat
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "CHIB_vs_EC1", "CHIB_vs_EC1", "EC1 (keV)","CHIB", "EC1",200,  0,15, "CHIB", 200, 0, 5
	# list_cuts=[]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	# #Plot CHIC as a function of heat
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "CHIC_vs_EC1", "CHIC_vs_EC1", "EC1 (keV)","CHIC", "EC1",200,  0,15, "CHIC", 200, 0, 5
	# list_cuts=[]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	# #Plot CHID as a function of heat
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "CHID_vs_EC1", "CHID_vs_EC1", "EC1 (keV)","CHID", "EC1",200,  0,15, "CHID", 200, 0, 5
	# list_cuts=[]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)



	# #Plot CHIC1 as a function of heat
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "XOC1_vs_EC1", "XOC1_vs_EC1", "EC1 (keV)","XOC1", "EC1",200,  0,15, "XOC1", 200, 0, 5
	# list_cuts=[]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	# #Plot CHIC2 as a function of heat
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "XOC2_vs_EC2", "XOC2_vs_EC2", "EC2 (keV)","XOC2", "EC2",200,  0,15, "XOC2", 200, 0, 5
	# list_cuts=[]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)


	# #Plot the FWHM versus time
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "OWC1_vs_time", "OWC1_vs_time", "Time","OWC1", "DateSec",100,  tmin, tmax, "OWC1", 200, 0, 5
	# list_cuts= []
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)


	# #Plot the chi2 of heat only
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "CHIC_heatonly", "CHIC_heatonly", "chi2", "Rate","CHIC", 300,  -2 ,5
	# list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_est["HEAT"] + ">0"]+ ["abs(EC1-EC2)<2"]
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "CHID_heatonly", "CHID_heatonly", "chi2", "Rate","CHID", 300,  -2 ,5
	# list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_est["HEAT"] + ">0"]+ ["abs(EC1-EC2)<2"]
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

 #    #Compare fiducial ionisation of noise events to that of E>4 heatonly events
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "ion_distri", "ion_distri", "IonFid (keV)", "Rate","0.5*(EIB+EID)", 200,-4,4
	# list_cuts1=[standard_cuts] + ["abs(EIB-EID)<2 && EIA<1 && EIC<1 "] +  [d_est["HEAT"] + "<2"]+ ["abs(EC1-EC2)<2"]
	# list_cuts2=[standard_cuts] + ["EIA<1 && EIC<1 && 0.5*(EIB+EID)<0.5*0.5*(EC1+EC2)"] + [d_est["HEAT"] + ">4"]+ ["abs(EC1-EC2)<2"]
	# get_plots.compare_fid_distri_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts1, list_cuts2)

 #    #Plot of the no heat only CHI2
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "CHIC_noheatonly", "CHIC_noheatonly", "chi2", "Rate","CHIA", 300,  -2 ,5
	# list_cuts=[standard_cuts] + ["0.5*(EIA+ EIB+ EIC+EID)>3 && 0.5*(EIA+ EIB+ EIC+EID)<30"] + ["abs(EC1-EC2)<2 && 0.5*(EC1+EC2)<30"] 
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "CHID_noheatonly", "CHID_noheatonly", "chi2", "Rate","CHIB", 300,  -2 ,5
	# list_cuts=[standard_cuts] + ["0.5*(EIA+ EIB+ EIC+EID)>3 && 0.5*(EIA+ EIB+ EIC+EID)<30"] + ["abs(EC1-EC2)<2 && 0.5*(EC1+EC2)<30"]
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

	# #Plot the time distribution of all events
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "time_distri_all", "time_distri_all", "Unix time", "Rate","DateSec", 300,  tmin, tmax
	# list_cuts=[standard_cuts] 
	# get_plots.get_1D_plot_with_thresh(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

 #    #Plot the time distribution of heat only events
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "time_distri_heatonly", "time_distri_heatonly", "Days", "Rate","JOUR", 300,  20,190
	# list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_est["HEAT"] + ">1.5"] #+ [d_est["HEAT"] + "<2"]
	# get_plots.get_1D_plot_with_thresh(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)
	# raw_input()
	# print tmin, tmax
	# bin_X = 600
	# ytitle = str((tmax-tmin)/(bin_X*24.*3600))[:4]
 #    # Plot the time distribution of heat only events
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "time_distri_heatonly", "time_distri_heatonly", "", "Rate (Events/" + ytitle+ " days)","DateSec", bin_X,  tmin, tmax
	# list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_est["HEAT"] + ">2"]+ ["abs(EC1-EC2)<2"]
	# get_plots.get_1D_plot_with_thresh_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

	# print tmin, tmax
 #    # Plot the time distribution of heat only events
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "time_distri_heatonly_from0keV", "time_distri_heatonly_from0keV", "", "Rate (Events/0.4 days)","DateSec", 300,  tmin, tmax
	# list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_est["HEAT"] + ">0"]+ ["abs(EC1-EC2)<2"]
	# get_plots.get_1D_plot_with_thresh_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)


 #    #Plot the 2D time energy distribution of heat only events
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "2D_energy_time_distri_heatonly", "2D_energy_time_distri_heatonly", "", "Energy","DateSec", 300,  tmin, tmax, d_est["HEAT"], 100, 1, 5
	# list_cuts=[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_est["HEAT"] + ">2"]+ ["abs(EC1-EC2)<2"]
	# get_plots.get_2D_plot_col_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

 
 #    #Do a S1 Q plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Qplot_S1", "Qplot_S1", "Heat", "Q", d_est["HEAT"], 200,  0,50, d_est["Q_S1"], 100,  0,1.2
	# # list_cuts=[standard_cuts]  + ["abs(EC1-EC2)<3 && EIC<0.1 && EID<0.1"] #+ ["EIA>0.50*.85*1.2 && EIB>0.50*.85*1.02 && EIC<2*0.85*1.01 && EID<2*0.85*0.91"]
	# list_cuts=[standard_cuts]  + [d_est["HEAT"] + ">0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&& abs(EC1-EC2)<2"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

 #    #Do a S2 Q plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Qplot_S2", "Qplot_S2", "Heat", "Q", d_est["HEAT"], 200,  0,50, d_est["Q_S2"], 100,  0,1.2
	# # list_cuts=[standard_cuts]   + ["abs(EC1-EC2)<3 && EIA<0.85*1.2 && EIB<0.85*1.02 && EIC>0.85*1.01 && EID>0.85*0.91"]
	# list_cuts=[standard_cuts]  + [d_est["HEAT"] + ">0  && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&& abs(EC1-EC2)<2"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

 #    #Do a S1 heat ion plot (multiply heat by (1+Vfid/3)/(1+5.5/3) = 1.29 to make sure Eh = Ei for gammas)
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ionS1", "Heat_ionS1" , "Heat", "IonS1", d_est["HEAT"], 1000,  -2,20.01, "0.5*(EIA+EIB+EIC+EID)", 1000,  -2,20.01
	# list_cuts=[standard_cuts]  + ["abs(EC1-EC2)<3 && EIC<1.5*0.85*1.01 && EIA<1.5*0.85*0.91"] #+ ["EIA>0.50*.85*1.2 && EIB>0.50*.85*1.02 && EIC<2*0.85*1.01 && EID<2*0.85*0.91"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

 #    #Do a S2 heat ion plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ionS2", "Heat_ionS2" , "Heat", "IonS2", d_est["HEAT"], 1000,  -2,20.01, "0.5*(EIA+EIB+EIC+EID)", 1000,  -2,20.01
	# list_cuts=[standard_cuts]   + ["abs(EC1-EC2)<3 && EIA<0.85*1.2 && EIB<0.85*1.02 && EIC>0.85*1.01 && EID>0.85*0.91"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)


	# #Plot the fiducial Gamma spectrum
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "Fid_Gamma_heat_1keV", "Fid_Gamma_heat_1keV", "Heat (keV)", "Counts",d_est["HEAT"], 100,  0, 5
	# list_cuts=[standard_cuts] + ["EIA<3 && EIC<3"] + [d_est["Q_FID"] + ">0.7"]
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

	# #Plot the fiducial Gamma spectrum
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "Fid_Gamma_heat_10keV", "Fid_Gamma_heat_10keV", "Heat (keV)", "Counts",d_est["HEAT"], 100,  8, 13
	# list_cuts=[standard_cuts] + ["EIA<3 && EIC<3"] + [d_est["Q_FID"] + ">0.7"]
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

	# #Plot the fiducial Gamma spectrum
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "Fid_Gamma_ion", "Fid_Gamma_ion", "Ion (keV)", "Counts",d_est["FID"], 200,  0, 5
	# list_cuts=[standard_cuts] + ["EIA<3 && EIC<3"] + [d_est["Q_FID"] + ">0.7"]
	# get_plots.get_1D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)


 #    #Do a fiducial Q plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Qplot_fid", "Qplot_fid", "EC", "Q", d_est["ER_FID"], 100,  0,100, d_est["Q_FID"], 100,  -1,2
	# list_cuts=[standard_cuts]  + [d_est["HEAT"] + ">0"] + ["EIA<1 && EIB>1 && EIC<1 && EID>1"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

 #    # Do a full Q plot_ERA
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Qplot_all", "Qplot_all", "ER", "Q", "((1+8./3)*0.5*(EC1+EC2)-0.33*(1.5*EIA+4*EIB+1.5*EIC+4*EID))", 200,  0,110, "0.5*(EIA+EIB+EIC+EID)/((1+8./3)*0.5*(EC1+EC2)-0.33*(1.5*EIA+4*EIB+1.5*EIC+4*EID))", 100,  -.1,1.2
	# list_cuts=[standard_cuts]  + [d_est["HEAT"] + ">0"]    #+ ["EIC<1 && EID<1"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	#Do a ion heat plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ion_all_ERA", "Heat_ion_all_ERA" , "Heat (keVee)", "Fiducial Ion (keVee)", d_est["HEAT"], 200,  1,10, d_est["FID"], 200,  1,10
	# list_cuts=[standard_cuts] + ["abs(EC1-EC2)<3"]#+ [d_est["HEAT"] + ">0"] 
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	#Do a ion heat plot
	fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ion_all", "Heat_ion_all" , "Heat (keVee)", "Fiducial Ion (keVee)", d_est["HEAT"], 200,  -1,10, "0.5*(EIB+EID)", 200,  -1,10
	list_cuts=[standard_cuts] + ["abs(EC1-EC2)<1"]
	get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	# #Do a ion heat plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ion_vetocut_all", "Heat_ion_vetocut_fid_all" , "Heat (keVee)", "Fiducial Ion (keVee)", d_est["HEAT"], 200,  1.5,15, "0.5*(EIB+EID)", 200,  -1,15
	# list_cuts=[standard_cuts] + ["abs(EC1-EC2)<1"]+ ["EIA<1 && EIC<1"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	# #Do a ion heat plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "Heat_ion_all_ANA", "Heat_ion_all_ANA" , "Heat (keVee)", "Fiducial Ion (keVee)", "0.5*(EC1+EC2)", 2000,  -2,20, "0.5*(EIB+EID)", 2000,  -2,20
	# list_cuts=[standard_cuts] + ["abs(EC1-EC2)<3"]#+ [d_est["HEAT"] + ">0"] 
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)


	# #Do a EC1 EC2 plot
	# fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y = "EC1_vs_EC2", "EC1_vs_EC2" , "EC1 (keV)", "EC2 (keV)", "EC1", 200,  -2,40, "EC2", 200,  -2, 40
	# list_cuts=[] #["0.5*(EIA+ EIB+ EIC+EID)>1"] # ["XOC1<6 && XOC2<6"] #["0.5*(EIA+ EIB+ EIC+EID)>1"] #["SDEL>4"]#["0.5*(EIA+ EIB+ EIC+EID)>1"]#[standard_cuts] #+ [d_est["HEAT"] + ">1"] #+ ["abs(EC1-EC2)<2"] + ['EIA<2.1 && EIC<1.9']  #+ ["0.5*(EIA+ EIB+ EIC+EID)<1 && SDEL>2"] + ["DateSec>1409.2E6"] #+  ["EIA<1 && EIC<1"]
	# get_plots.get_2D_plot(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, channel_Y, bin_Y, min_Y, max_Y, list_cuts)

	return


bolo_list=["FID837"]
data_dir  = "../Fond_ERA_merged/"
for bolo_name in bolo_list:
	launch_plots( bolo_name, data_dir)