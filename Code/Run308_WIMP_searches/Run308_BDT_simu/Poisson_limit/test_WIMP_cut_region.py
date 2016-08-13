from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import script_utils as script_utils
import math
from array import array
import os, sys

def test_WIMP_cut_region(bolo_name, mass, exposure, analysis, nsimu):

	"""Get number of events passing Poisson cuts
	
	Detail:

	Args:
		bolo_name = (str) bolometer name
		mass      = (str) WIMP name
		exposure  = (str) exposure
		analysis  = (str) analysis type (which cut on fid and veto)
		nsimu     = (int) number of simulated data samples
	Returns:
		void

	Raises:
		void
	"""
	

	#Load 2D PDF
	fWIMP2D, f = PyRPl.open_ROOT_object("./ROOT_files/WIMP_PDF2D_" + analysis + ".root", "WIMP_" + mass + "_GeV")

	#Load cut value on PDF for 95% WIMP box
	cut_val = 0
	with open ("./Text_files/WIMP_PDF_95_cut_value_" + analysis + ".txt", "r") as fcut:
		stuff = [elem.rstrip().split(",") for elem in fcut.readlines()]
		for elem in stuff:
			mass_val = elem[0]
			if int(mass)==int(mass_val):
				cut_val = float(elem[1])
	

	wimp_path = "../../BDT_" + bolo_name + "/" + analysis + "/WIMP/ROOT_files/"
	filou = TFile(wimp_path + bolo_name + "_WIMP_mass_" + str(mass) + "_tree.root", "read")
	tree = filou.Get("t_new0")
	num_pass_cut =0

	for k in range(tree.GetEntries()):
		tree.GetEntry(k)
		if fWIMP2D.Eval(tree.ENR, 0.5*(tree.EIB+tree.EID))>cut_val:
			num_pass_cut+=1

	print float(num_pass_cut)/float(tree.GetEntries())

def test_WIMP_cut_region_on_true_data(bolo_name, mass, analysis):

	"""Get number of events passing Poisson cuts
	
	Detail:

	Args:
		bolo_name = (str) bolometer name
		mass      = (str) WIMP name
		exposure  = (str) exposure
		analysis  = (str) analysis type (which cut on fid and veto)
		nsimu     = (int) number of simulated data samples
	Returns:
		void

	Raises:
		void
	"""
	

	#Load 2D PDF
	fWIMP2D, f = PyRPl.open_ROOT_object("./ROOT_files/WIMP_PDF2D_" + analysis + ".root", "WIMP_" + mass + "_GeV")

	#Load cut value on PDF for 95% WIMP box
	cut_val_90, cut_val_99 = 0,0
	with open ("./Text_files/WIMP_PDF_90_and_99_cut_value_" + analysis + ".txt", "r") as fcut:
		stuff = [elem.rstrip().split(",") for elem in fcut.readlines()]
		for elem in stuff:
			mass_val = elem[0]
			if int(mass)==int(mass_val):
				cut_val_90 = float(elem[1])
				cut_val_99 = float(elem[2])
	

	data_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/"
	filou = TFile(data_path + bolo_name + "_" + analysis + "_fond.root", "read")
	tree = filou.Get("data")
	num_pass_cut =0

	hpass = TH2F("hpass", "hpass", 100, 0, 15, 100, 0, 15)

	# #T Check that the events are found where expected
	# arr1 = np.random.uniform(0,15,size=(200000,2))
	# for i in range(arr1.shape[0]):
	# 	PDF_val = fWIMP2D.Eval(arr1[i][0], arr1[i][1])
	# 	if (cut_val_99<PDF_val<cut_val_90):
	# 	# if (cut_val_99<PDF_val<cut_val_90):
	# 		num_pass_cut+=1
	# 		hpass.Fill(arr1[i][0], arr1[i][1])		

	# hpass.Draw()
	# raw_input()

	for k in range(tree.GetEntries()):
		tree.GetEntry(k)
		ER=(1+8./3)*0.5*(tree.EC1+tree.EC2)-0.33*(1.5*tree.EIA+4*tree.EIB+1.5*tree.EIC+4*tree.EID)
		PDF_val = fWIMP2D.Eval(ER, 0.5*(tree.EIB+tree.EID))
		if (cut_val_99<PDF_val<cut_val_90 and 0.5*(tree.EIB+tree.EID)>0.7):
		# if (cut_val_99<PDF_val<cut_val_90):
			num_pass_cut+=1
			hpass.Fill(0.5*(tree.EC1+tree.EC2), 0.5*(tree.EIB+tree.EID))

	print num_pass_cut
	hpass.Draw()
	raw_input()

bolo_name = "FID837"
list_mass = ["5", "6", "7", "10", "25"]
analysis = "ana_0.5_0_5"

for mass in list_mass:
	test_WIMP_cut_region_on_true_data(bolo_name, mass, analysis)
