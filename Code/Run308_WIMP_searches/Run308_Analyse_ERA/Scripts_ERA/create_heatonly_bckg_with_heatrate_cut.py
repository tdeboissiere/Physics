#! /usr/bin/env python

from ROOT import *
import Poisson_file_handler as Poisson_fh
import PyROOTPlots as PyRPl
import os,sys
import script_utils as script_utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from scipy import interpolate
import Analysis_utilities as Ana_ut
import BDT_file_handler as BDT_fh

def create_extrapolated_bckg(bolo_name):

	"""Create the extrapolated bckg for heat and gamma
    
	Detail:
		Use triple exp for heat and poly + gaussian peaks for gamma

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		Assertion Error if the initial histograms do not exist
	"""

	#Output dir
	output_heat_path = "../Analyse_" + bolo_name +"/ROOT_files/"

	##################
	#
	#    Heatonly 2D
	#
	##################
	
	hheat, file_hheat = PyRPl.open_ROOT_object("../Analyse_" + bolo_name +"/ROOT_files/" + bolo_name+"_thresh.root", "hheat")
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")	
	standard_cuts = standard_cuts + "&&KTH<1&&KTH>0&&0.5*(EIB+EID)<0"
	
	heat_max = hheat.GetMaximum()
	

	l_heatonly = TEventList("l_heatonly")
	tree.Draw(">>l_heatonly",standard_cuts)
	pop_len    = l_heatonly.GetN()

	try :
		os.remove(output_heat_path + bolo_name + "_heatonly_2D_with_heat_cut.root")
	except OSError:
		pass

	for fraction in [0.1,0.2,0.3,0.4,0.5,0.8,1]:
		print fraction
		file_heat2D = TFile(output_heat_path + bolo_name + "_heatonly_2D_with_heat_cut.root", "update")
		hheat2D = TH2F("heat2D_fraction_" + str(fraction), "heat2D_fraction_" + str(fraction), 200, -2, 15, 200, -2, 15)
		for k in range(pop_len):
			counter  = l_heatonly.GetEntry(k)
			tree.GetEntry(counter)
			time = 1E6*tree.UT1+ tree.UT2
			hr = hheat.GetBinContent(hheat.FindBin(time))
			if hr<= fraction*heat_max:
				hheat2D.Fill(tree.EC1, tree.EC2)

		hheat2D.Write()
		file_heat2D.Close() 
		del hheat2D 
		del file_heat2D
	


if __name__ == '__main__':

	bolo_name ="FID837"
	create_extrapolated_bckg(bolo_name)