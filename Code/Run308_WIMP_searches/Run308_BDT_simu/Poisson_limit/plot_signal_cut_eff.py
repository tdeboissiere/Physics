

from ROOT import *
import script_utils as script_utils
import math
import PyROOTPlots as PyRPl
from ctypes import *
import BDT_file_handler as BDT_fh 
import numpy as np


def plot_ENR_spectra(bolo_name, mass, analysis_type):

	"""Plot ENR spectra of WIMP with different TCuts
	
	
	Detail:
		Plot the ENR spectra without cuts, with the pre sel cuts and with the BDT cut
		Do it as a function of the mass

	Args:
		bolo_name           = (str) bolometer name
		mass                = (int) WIMP mass
		analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
		bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram

	Returns:
		void

	Raises:
		void
	"""    

	#Load cut value on BDT output
	cut_val = 0     
	with open ("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + ".txt", "r") as fcut:
		stuff = [elem.rstrip().split(",") for elem in fcut.readlines()]
		for elem in stuff:
			mass_val = elem[0]
			if int(mass) ==int(mass_val):
				cut_val = float(elem[1])
	
	WIMP_gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/WIMP/ROOT_files/"
	Application_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"
	
	#Load the WIMP trees
	tnocut, fnocut   = PyRPl.open_ROOT_object(WIMP_gen_path + bolo_name + "_WIMP_mass_" + mass + "_tree.root", "t_newnocut1")
	tcut, fcut       = PyRPl.open_ROOT_object(WIMP_gen_path + bolo_name + "_WIMP_mass_" + mass + "_tree.root", "t_new1")
	tBDTcut, fBDTcut = PyRPl.open_ROOT_object(Application_path +  "WIMP" + "_mass_" + str(mass) + "_tree.root", "tout")

	tBDTcut.AddFriend(tcut)
	# tBDTcut.Draw("tout.NN:t_new1.ENR>>hist(100,2,8, 100,-2,2)")
	# raw_input()

	hWIMPnocut  = TH1F("hWIMPnocut", "hWIMPnocut", 100, 0, 15)
	hWIMPcut    = TH1F("hWIMPcut", "hWIMPcut", 100, 0, 15)
	hWIMPBDTcut = TH1F("hWIMPBDTcut", "hWIMPBDTcut", 100, 0, 15)

	tnocut.Project("hWIMPnocut", "ENR")
	tcut.Project("hWIMPcut", "ENR")
	tBDTcut.Project("hWIMPBDTcut", "t_new1.ENR", "tout.NN>" + str(cut_val))


	PyRPl.process_TH1(hWIMPnocut, X_title = "ENR (keV)", Y_title = "Efficiency", color = kBlack)
	PyRPl.process_TH1(hWIMPcut, X_title = "ENR (keV)", Y_title = "Efficiency", color = kRed)
	PyRPl.process_TH1(hWIMPBDTcut, X_title = "ENR (keV)", Y_title = "Efficiency", color = kBlue)

	#Set bin errors
	list_hist = [hWIMPnocut, hWIMPcut, hWIMPBDTcut]
	for h in list_hist:
		for i in range(1,101):
			h.SetBinError(i, TMath.Sqrt(h.GetBinContent(i)))

	# #Plot histograms as a check
	# hWIMPnocut.Draw("")	
	# hWIMPcut.Draw("same")	
	# hWIMPBDTcut.Draw("same")

	# raw_input()

	# plot the efficiency functions
	# hWIMPcut.Divide(hWIMPnocut)	
	# hWIMPBDTcut.Divide(hWIMPnocut)	
	# hWIMPBDTcut.SetMarkerStyle(1)
	# hWIMPcut.SetMarkerStyle(1)
	# hWIMPcut.SetMinimum(0)
	# hWIMPcut.SetMaximum(1.5)
	# hWIMPcut.Draw("E1")	
	# hWIMPBDTcut.Draw("E1same")
	
	# leg = TLegend(0.1934673,0.5813953,0.5037688,0.7043189,"","brNDC")
	# leg.AddEntry(hWIMPcut.GetName(), "Pre selection cut", "leg")
	# leg.AddEntry(hWIMPBDTcut.GetName(),"BDT cut" ,"leg")
	# leg.SetFillColor(kWhite)
	# leg.SetBorderSize(0)
	# leg.Draw("same")

	teff1 = TEfficiency(hWIMPcut, hWIMPnocut)
	teff2 = TEfficiency(hWIMPBDTcut, hWIMPnocut)


	teff1.SetName("teff1")
	teff2.SetName("teff2")

	teff1.SetLineColor(kRed)
	teff1.Draw()
	teff2.Draw("same")

	leg = TLegend(0.1934673,0.5813953,0.5037688,0.7043189,"","brNDC")
	leg.AddEntry(teff1.GetName(), "Pre selection cut", "leg")
	leg.AddEntry(teff2.GetName(),"Pre selection + BDT cut" ,"leg")
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)
	leg.Draw("same")
	
	raw_input()


	c1.Print("./Figures/" + bolo_name + "_efficiency_mass_" + str(mass) + "_baseline_time.eps")

bolo_name           = "FID837"
mass                = "25"
analysis_type       = "ana_0.5_0_5"
for mass in ["3", "4", "5", "6", "7", "10", "25"]:
	plot_ENR_spectra(bolo_name, mass, analysis_type)
