#! /usr/bin/env python

from ROOT import *
import BDT_file_handler as BDT_fh
import PyROOTPlots as PyRPl
import os
import script_utils as script_utils
import numpy as np
import Analysis_utilities as Ana_ut
import matplotlib.pylab as plt


def get_scaling_heatonly(bolo_name, ECinf, ECsup):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	#Load tree and standard cuts
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")
	
	exp_heat = 2*tree.GetEntries(standard_cuts + "&& 0.5*(EIB+EID)<0 && EC>" + str(ECinf) + "&& EC<" + str(ECsup) )
	# tree.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)>>hist(100,0.5,15,100,-1,15)", standard_cuts)
	print exp_heat, tree.GetEntries(standard_cuts +"&& EIA<1.4 && EIC<1.4 && 0.5*(EIB+EID)<0.7&& EC>" + str(ECinf) + "&& EC<" + str(ECsup))

	return exp_heat

def get_scaling_FidGamma(bolo_name, bound, ECinf, ECsup):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	gamma_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name
	fFidGamma, file_Gamma = PyRPl.open_ROOT_object(gamma_path + "/ROOT_files/" + bolo_name + "_Gamma_spectrum.root", "FidGamma")

	#Load tree and standard cuts
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])

	gamma_cut = standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"] + "&& EIB>" + d_std["FWIB"] + "&& EID>" + d_std["FWID"] + "&& EC>0.5"
	########################
	#1) plot  (checks)
	########################
	# h2D = TH2F("h2D", "h2D", 1000, 0, 15, 1000, 0, 15)
	# h2Dcut = TH2F("h2Dcut", "h2Dcut", 1000, 0, 15, 1000, 0, 15)
	# tree.Project("h2D", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts + "&& EIA<" + d_std["FWIA"] + "&& EIC <" + d_std["FWIC"])
	# tree.Project("h2Dcut", "0.5*(EIB+EID):0.5*(EC1+EC2)", gamma_cut)
	# h2Dcut.SetMarkerColor(kRed)
	# h2D.Draw()
	# h2Dcut.Draw("same")

	########################
	#2) Compute
	########################

	ngamma_true = tree.GetEntries(gamma_cut + "&& EC>" + str(bound["low"]) + "&& 0.5*(EIB+EID)>" + str(bound["low"]) + "&& EC<" + str(bound["up"]) +"&& 0.5*(EIB+EID)<" + str(bound["up"]) )
	ngamma_func = fFidGamma.Integral(bound["low"],bound["up"])
	norm_gamma = ngamma_true/ngamma_func

	########################
	#3) Control
	########################

	# h1Dtrue = TH1F("h1Dtrue", "h1Dtrue", 150, 0, 15)
	# tree.Project("h1Dtrue", "EC", gamma_cut)
	# h1Dsimu = TH1F("h1Dsimu", "h1Dsimu", 150, 0, 15)
	# for i in range(1,151):
	# 	h1Dsimu.SetBinContent(i, fFidGamma.Eval(h1Dsimu.GetBinCenter(i)))
	# h1Dsimu.Scale(norm_gamma*fFidGamma.Integral(0,15)/h1Dsimu.Integral())
	# h1Dsimu.SetLineColor(kRed)
	# h1Dtrue.Draw() 
	# h1Dsimu.Draw("same")
	# raw_input()
	# print fFidGamma.Integral(3,bound["up"])*norm_gamma, tree.GetEntries(gamma_cut + "&& EC>3 && 0.5*(EIB+EID)>3 && EC<15 && 0.5*(EIB+EID)<15")

	# print fFidGamma.Integral(ECinf, ECsup)*norm_gamma, tree.GetEntries(gamma_cut + "&& EC>0.5 && 0.5*(EIB+EID)>0.5 && EC<15 && 0.5*(EIB+EID)<15")

	return 	fFidGamma.Integral(ECinf, ECsup)*norm_gamma

def get_scaling_SurfGamma(bolo_name, ECinf, ECsup):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	gamma_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name
	fFidGamma, file_Gamma = PyRPl.open_ROOT_object(gamma_path + "/ROOT_files/" + bolo_name + "_Gamma_spectrum.root", "FidGamma")
	fS1Gamma = file_Gamma.Get("S1Gamma")
	fS2Gamma = file_Gamma.Get("S2Gamma")

	#Load tree and standard cuts
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])
	#Load estimators
	d_est = BDT_fh.open_estimator_file(bolo_name, "")

	S1gamma_cut = standard_cuts +  "&& EC >0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&&" + d_est["Q_S1"] + ">0.65"
	S2gamma_cut = standard_cuts +  "&& EC >0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" + d_est["Q_S2"] + ">0.65"
	

	########################
	#1) plot  (checks)
	########################
	# h2D = TH2F("h2D", "h2D", 1000, 0, 15, 1000, 0, 15)
	# h2Dcut = TH2F("h2Dcut", "h2Dcut", 1000, 0, 15, 1000, 0, 15)
	# tree.Project("h2D", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts)
	# tree.Project("h2Dcut", "0.5*(EIB+EID):0.5*(EC1+EC2)", S2gamma_cut)
	# h2Dcut.SetMarkerColor(kRed)
	# h2D.Draw()
	# h2Dcut.Draw("same")
	# raw_input()

	########################
	#2) Compute
	########################

	nS1_true = tree.GetEntries(S1gamma_cut + "&& 0.5*(EIB+EID)>2 && 0.5*(EIB+EID)<6 && EC>5 && EC<10")
	nS2_true = tree.GetEntries(S2gamma_cut + "&& 0.5*(EIB+EID)>2 && 0.5*(EIB+EID)<6 && EC>5 && EC<10")
	nS1_func = fS1Gamma.Integral(5,10)
	nS2_func = fS2Gamma.Integral(5,10)
	norm_S1 = nS1_true/nS1_func
	norm_S2 = nS2_true/nS2_func

	########################
	#3) Control
	########################

	# h1Dtrue = TH1F("h1Dtrue", "h1Dtrue", 150, 0, 15)
	# tree.Project("h1Dtrue", "EC", S2gamma_cut)
	# h1Dsimu = TH1F("h1Dsimu", "h1Dsimu", 150, 0, 15)
	# for i in range(1,151):
	# 	h1Dsimu.SetBinContent(i, fS2Gamma.Eval(h1Dsimu.GetBinCenter(i)))
	# h1Dsimu.Scale(norm_S2*fS2Gamma.Integral(0,15)/h1Dsimu.Integral())
	# h1Dsimu.SetLineColor(kRed)
	# h1Dtrue.Draw() 
	# h1Dsimu.Draw("same")
	# raw_input()

	return 	fS1Gamma.Integral(ECinf, ECsup)*norm_S1, fS2Gamma.Integral(ECinf, ECsup)*norm_S2


def get_scaling_SurfBeta(bolo_name, ECinf, ECsup):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	Beta_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name
	file_Beta = TFile(Beta_path + "/ROOT_files/" + bolo_name + "_Beta_spectrum_extrapol.root", "read")
	fS1Beta = file_Beta.Get("S1Beta_extra")
	fS2Beta = file_Beta.Get("S2Beta_extra")

	#Load tree and standard cuts
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])
	#Load estimators
	d_est = BDT_fh.open_estimator_file(bolo_name, "")
	
	S1Beta_cut = standard_cuts +  "&&EC>0  && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&&" +  d_est["Q_S1"] + "<0.65 && " + d_est["Q_S1"] + ">0.2"
	S2Beta_cut = standard_cuts +  "&&EC>0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" +  d_est["Q_S2"] + "<0.65 && " + d_est["Q_S2"] + ">0.2"	

	########################
	#1) plot  (checks)
	########################
	# h2D = TH2F("h2D", "h2D", 1000, 0, 15, 1000, 0, 15)
	# h2Dcut = TH2F("h2Dcut", "h2Dcut", 1000, 0, 15, 1000, 0, 15)
	# tree.Project("h2D", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts)
	# tree.Project("h2Dcut", "0.5*(EIB+EID):0.5*(EC1+EC2)", S1Beta_cut)
	# h2Dcut.SetMarkerColor(kRed)
	# h2D.Draw()
	# h2Dcut.Draw("same")
	# raw_input()

	########################
	#2) Compute
	########################

	nS1_true = tree.GetEntries(S1Beta_cut + "&& 0.5*(EIB+EID)>1 && 0.5*(EIB+EID)<6 && EC>4 && EC<10")
	nS2_true = tree.GetEntries(S2Beta_cut + "&& 0.5*(EIB+EID)>1 && 0.5*(EIB+EID)<6 && EC>4 && EC<10")
	nS1_func = fS1Beta.Integral(4,10)
	nS2_func = fS2Beta.Integral(4,10)
	norm_S1 = nS1_true/nS1_func
	norm_S2 = nS2_true/nS2_func

	########################
	#3) Control
	########################

	# h1Dtrue = TH1F("h1Dtrue", "h1Dtrue", 150, 0, 40)
	# tree.Project("h1Dtrue", "EC", S1Beta_cut)
	# h1Dsimu = TH1F("h1Dsimu", "h1Dsimu", 150, 0, 40)
	# for i in range(1,151):
	# 	h1Dsimu.SetBinContent(i, fS1Beta.Eval(h1Dsimu.GetBinCenter(i)))
	# h1Dsimu.Scale(norm_S1*fS1Beta.Integral(0,15)/h1Dsimu.Integral())
	# h1Dsimu.SetLineColor(kRed)
	# h1Dtrue.Draw() 
	# h1Dsimu.Draw("same")
	# raw_input()

	# print fS2Beta.Integral(2,10)*norm_S2, tree.GetEntries(S2Beta_cut + "&& EC>2 && EC<10")

	return 	fS1Beta.Integral(ECinf, ECsup)*norm_S1, fS2Beta.Integral(ECinf, ECsup)*norm_S2

def get_scaling_SurfPb(bolo_name, ECinf, ECsup):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	Pb_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name
	file_Pb = TFile(Pb_path + "/ROOT_files/" + bolo_name + "_Pb_spectrum_extrapol.root", "read")
	fS1Pb = file_Pb.Get("S1Pb_extra")
	fS2Pb = file_Pb.Get("S2Pb_extra")

	#Load tree and standard cuts
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,"")
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])
	#Load estimators
	d_est = BDT_fh.open_estimator_file(bolo_name, "")
	
	S1Pb_cut = standard_cuts +  "&&" + d_est["HEAT"] + ">0 &&  EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + "&&" +  d_est["Q_S1"] + "<0.15 &&" +  d_est["Q_S1"] + ">0.04"
	S2Pb_cut = standard_cuts +  "&&" + d_est["HEAT"] + ">0 && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + "&&" +  d_est["Q_S2"] + "<0.15 &&" +  d_est["Q_S2"] + ">0.04"

	# #######################
	# #1) plot  (checks)
	# #######################
	# h2D = TH2F("h2D", "h2D", 1000, 0, 50, 1000, 0, 50)
	# h2Dcut = TH2F("h2Dcut", "h2Dcut", 1000, 0, 50, 1000, 0, 50)
	# tree.Project("h2D", "0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts)
	# tree.Project("h2Dcut", "0.5*(EIB+EID):0.5*(EC1+EC2)", S2Pb_cut)
	# h2Dcut.SetMarkerColor(kRed)
	# h2D.Draw()
	# h2Dcut.Draw("same")
	# raw_input()

	########################
	#2) Compute
	########################

	nS1_true = tree.GetEntries(S1Pb_cut + "&& 0.5*(EIB+EID)>1 && 0.5*(EIB+EID)<8 && EC>15 && EC<40")
	nS2_true = tree.GetEntries(S2Pb_cut + "&& 0.5*(EIB+EID)>1 && 0.5*(EIB+EID)<8 && EC>15 && EC<40")
	nS1_func = fS1Pb.Integral(15,40)
	nS2_func = fS2Pb.Integral(15,40)
	norm_S1 = nS1_true/nS1_func
	norm_S2 = nS2_true/nS2_func

	########################
	#3) Control
	########################

	# h1Dtrue = TH1F("h1Dtrue", "h1Dtrue", 30, 0, 40)
	# tree.Project("h1Dtrue", "EC", S1Pb_cut)
	# h1Dsimu = TH1F("h1Dsimu", "h1Dsimu", 30, 0, 40)
	# for i in range(1,151):
	# 	h1Dsimu.SetBinContent(i, fS1Pb.Eval(h1Dsimu.GetBinCenter(i)))
	# h1Dsimu.Scale(norm_S1*fS1Pb.Integral(0,15)/h1Dsimu.Integral())
	# h1Dsimu.SetLineColor(kRed)
	# h1Dtrue.Draw() 
	# h1Dsimu.Draw("same")
	# raw_input()

	print fS1Pb.Integral(2,10)*norm_S1, tree.GetEntries(S1Pb_cut + "&& EC>2 && EC<10")

	return 	fS1Pb.Integral(ECinf, ECsup)*norm_S1, fS2Pb.Integral(ECinf, ECsup)*norm_S2

def get_scaling():

	"""Get the scaling
    
	Detail:

	Args:
	Returns:
		void

	Raises:
		void
	"""

	bolo_name = "FID837"
	ECinf = .5
	ECsup = 15

	exp_heat = get_scaling_heatonly(bolo_name, ECinf, ECsup)
	exp_FidGamma = get_scaling_FidGamma(bolo_name, {"low":4,"up":15}, ECinf, ECsup)
	exp_S1Gamma, exp_S2Gamma = get_scaling_SurfGamma(bolo_name, ECinf, ECsup)
	exp_S1Beta, exp_S2Beta = get_scaling_SurfBeta(bolo_name, ECinf, ECsup)
	exp_S1Pb, exp_S2Pb = get_scaling_SurfPb(bolo_name, ECinf, ECsup)

	#Load duration
	duration = 0
	file_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
	duration_file_name =bolo_name + "_duration_KTH_cut.txt" 
	assert(os.path.isfile(file_path_name + duration_file_name) )
	with open(file_path_name + duration_file_name, 'r') as duration_file: 
	    list_duration_lines                                   = duration_file.readlines()[0].rstrip().split(",") 
	    duration                                              = float(list_duration_lines[3])
	    script_utils.print_utility(script_utils.COL("Duration (in days) = "  + str(duration) , "blue"))

	list_exp = [exp_heat, exp_FidGamma, exp_S1Gamma, exp_S2Gamma, exp_S1Beta, exp_S2Beta, exp_S1Pb, exp_S2Pb]
	list_pct = [float(elem)/float(sum(list_exp)) for elem in list_exp]
	list_rate = [float(elem)/float(duration) for elem in list_exp]
	for index, elem in enumerate(list_exp):
		print elem, list_pct[index], list_rate[index]

	with open(file_path_name + bolo_name + "_bckg_scaling_no_ion_cut_" + str(ECinf) + "_heat_thresh.txt", "w") as fs:
		fs.write("#no ion cut, 0.5*(EC1+EC2)>0, rate in days, " + str(ECinf)  + " to 15 keV in heat\n")
		fs.write("prop_heatonly,prop_FidGamma,prop_S1Gamma,prop_S2Gamma,prop_S1Beta,prop_S2Beta,prop_S1Pb,prop_S2Pb\n")
		fs.write(",".join([str(elem) for elem in list_pct]) + "\n")
		fs.write("rate_heatonly,rate_FidGamma,rate_S1Gamma,rate_S2Gamma,rate_S1Beta,rate_S2Beta,rate_S1Pb,rate_S2Pb\n")
		fs.write(",".join([str(elem) for elem in list_rate]) + "\n")

get_scaling()
