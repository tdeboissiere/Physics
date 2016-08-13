from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import script_utils as script_utils
import math
from array import array
import BDT_file_handler as BDT_fh


def get_events_passing_cuts(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

	"""Simulate data and find the number of events that pass
	the BDT cut

	
	Detail:
		void

	Args:
		bolo_name           = (str) bolometer name
		mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
		analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
		bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
		exposure            = (float) exposure of the simulated data

	Returns:
		void

	Raises:
		void
	"""       

	#Load cut value on BDT output
	cut_val = 0     
	# with open ("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + "_" + str(exposure) + ".txt", "r") as fcut:
	with open ("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + ".txt", "r") as fcut:
		stuff = [elem.rstrip().split(",") for elem in fcut.readlines()]
		for elem in stuff:
			mass_val = elem[0]
			if int(mass) ==int(mass_val):
				cut_val = float(elem[1])
	
	d_scaling = BDT_fh.open_scaling_file(bolo_name, d_cut, "")
	d_bckg_cut_eff = BDT_fh.open_bckg_cuteff_file(bolo_name, analysis_type, "")
	
	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"

	#############################################
	#Get all histograms
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
	
	tS1Pb     = fS1Pb.Get("tout")
	tS2Pb     = fS2Pb.Get("tout") 
	tS1Beta   = fS1Beta.Get("tout") 
	tS2Beta   = fS2Beta.Get("tout") 
	tS1Gamma  = fS1Gamma.Get("tout") 
	tS2Gamma  = fS2Gamma.Get("tout") 
	tFidGamma = fFidGamma.Get("tout") 
	theat     = fheat.Get("tout")
	tWIMP     = fWIMP.Get("tout") 
	
	hS1Pb     = TH1F("hS1Pb", "hS1Pb", bin_X, min_X, max_X)
	hS2Pb     = TH1F("hS2Pb", "hS2Pb", bin_X, min_X, max_X)
	hS1Beta   = TH1F("hS1Beta", "hS1Beta", bin_X, min_X, max_X)
	hS2Beta   = TH1F("hS2Beta", "hS2Beta", bin_X, min_X, max_X)
	hS1Gamma  = TH1F("hS1Gamma", "hS1Gamma", bin_X, min_X, max_X)
	hS2Gamma  = TH1F("hS2Gamma", "hS2Gamma", bin_X, min_X, max_X)
	hFidGamma = TH1F("hFidGamma", "hFidGamma", bin_X, min_X, max_X)
	hheat     = TH1F("hheat", "hheat", bin_X, min_X, max_X)
	hWIMP     = TH1F("hWIMP", "hWIMP", bin_X, min_X, max_X)
	
	tS1Pb.Project("hS1Pb", "NN")
	tS2Pb.Project("hS2Pb", "NN")
	tS1Beta.Project("hS1Beta", "NN")
	tS2Beta.Project("hS2Beta", "NN")
	tS1Gamma.Project("hS1Gamma", "NN")
	tS2Gamma.Project("hS2Gamma", "NN")
	tFidGamma.Project("hFidGamma", "NN")
	theat.Project("hheat", "NN")
	tWIMP.Project("hWIMP", "NN")
	
	
	#Rescale WIMP  to hdata
	hWIMP.Scale(1./float(hWIMP.Integral()))
	
	#Rescale  bckg histos to their expected value
	hheat.Scale(float(d_scaling["rate_heatonly"])*exposure*float(d_bckg_cut_eff["heatonly"])/float(hheat.Integral()))
	hFidGamma.Scale(float(d_scaling["rate_FidGamma"])*exposure*float(d_bckg_cut_eff["FidGamma"])/float(hFidGamma.Integral()))
	hS1Gamma.Scale(float(d_scaling["rate_S1Gamma"])*exposure*float(d_bckg_cut_eff["S1Gamma"])/float(hS1Gamma.Integral()))
	hS2Gamma.Scale(float(d_scaling["rate_S2Gamma"])*exposure*float(d_bckg_cut_eff["S2Gamma"])/float(hS2Gamma.Integral()))
	hS1Beta.Scale(float(d_scaling["rate_S1Beta"])*exposure*float(d_bckg_cut_eff["S1Beta"])/float(hS1Beta.Integral()))
	hS2Beta.Scale(float(d_scaling["rate_S2Beta"])*exposure*float(d_bckg_cut_eff["S2Beta"])/float(hS2Beta.Integral()))
	hS1Pb.Scale(float(d_scaling["rate_S1Pb"])*exposure*float(d_bckg_cut_eff["S1Pb"])/float(hS1Pb.Integral()))
	hS2Pb.Scale(float(d_scaling["rate_S2Pb"])*exposure*float(d_bckg_cut_eff["S2Pb"])/float(hS2Pb.Integral()))
	
	list_hist_bckg =[hS1Pb, hS2Pb, hS1Beta, hS2Beta, hS1Gamma, hS2Gamma, hFidGamma, hheat]
	
	hsum=TH1F("hsum","hsum", bin_X, min_X, max_X)
	for i in range(1,bin_X+1):
		sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
		hsum.SetBinContent(i, sumcontent)

	#Run Poisson simulations
	list_event_pass_cut=[]
	for nsimu in range(100):
		hdatasimu = TH1F("hdatasimu","hdatasimu", bin_X, min_X, max_X)
		for i in range(1,bin_X+1):
			hdatasimu.SetBinContent(i, np.random.poisson(hsum.GetBinContent(i)))
		bin_cut = hdatasimu.FindBin(cut_val)
		num_entry_cut = int(hdatasimu.Integral(bin_cut, max_X))
		list_event_pass_cut.append(str(num_entry_cut))
		del hdatasimu

	return list_event_pass_cut


def get_true_events_passing_cuts(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

	"""Simulate data and find the number of events that pass
	the BDT cut

	
	Detail:
		void

	Args:
		bolo_name           = (str) bolometer name
		mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
		analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
		bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
		exposure            = (float) exposure of the simulated data

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
	
	
	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type +"/Application/"

	#############################################
	#Get data
	##########################################
	
	frealdata     = TFile(gen_path +  "real_data_mass_" + str(mass) + "_tree.root", "read")
	trealdata     = frealdata.Get("tout")
	hrealdata     = TH1F("hrealdata", "hrealdata", bin_X, min_X, max_X)
	trealdata.Project("hrealdata", "NN")
	
	#Run Poisson simulations
	list_event_pass_cut=[]
	bin_cut = hrealdata.FindBin(cut_val)
	num_entry_cut = int(hrealdata.Integral(bin_cut, max_X))
	list_event_pass_cut.append(str(num_entry_cut))

	return list_event_pass_cut


bolo_name = "FID837"
list_mass = ["3", "4", "5", "6", "7", "10", "25"]
# list_mass = ["5", "6", "7", "10", "25"]
analysis_type = "ana_0.5_0.3_5"
bin_X, min_X, max_X = 200, -1, 1
exposure = 66
d_cut       = {"ECinf": 0.5, "ECsup": 15, "EIinf": 0.3, "EIsup": 15, "sigma_vet": 5}


# with open("./Text_files/events_passing_cuts_" + analysis_type + "_" + str(exposure) + ".txt", "w") as fev:
with open("./Text_files/events_passing_cuts_" + analysis_type +  ".txt", "w") as fev:
	for mass in list_mass:
		print mass
		list_ev = get_events_passing_cuts(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure)
		fev.write(mass + "," + ",".join(list_ev) + "\n")

with open("./Text_files/true_events_passing_cuts_" + analysis_type + ".txt", "w") as fev:
	for mass in list_mass:
		print mass
		list_ev = get_true_events_passing_cuts(bolo_name, mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure)
		fev.write(mass + "," + ",".join(list_ev) + "\n")