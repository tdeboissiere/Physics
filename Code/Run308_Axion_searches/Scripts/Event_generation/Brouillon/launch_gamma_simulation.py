#! /usr/bin/env python

from ROOT import *
import script_utils as script_utils
import adjust_Gamma_file as adjust_Gamma
import subprocess as subp

def launch_analysis(bolo_name):

	"""Launch BDT simulations
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
		overwrite     = (dict) indicates if standard files need to be overwritten

	Returns:
		void

	Raises:
		AssertionError 
	"""


	#Get the correct root config options for compilation:
	cflags         = subp.Popen(["root-config", "--cflags"], stdout = subp.PIPE).communicate()[0].rstrip()
	libs           = subp.Popen(["root-config", "--libs"], stdout = subp.PIPE).communicate()[0].rstrip()

	#########
	#Simulate fiducial gammas
	#########
	script_utils.print_utility("Gamma")

	adjust_Gamma.adjust_Gamma_file(bolo_name, d_cut, analysis_type, d_event_type_num_event)

	list_subp_args = ["g++", "-o", "./essai.exe",  "./build_Gamma_tree.C", "-Wl,--no-as-needed"]
	list_subp_args.extend(cflags.split(" "))
	list_subp_args.extend(libs.split(" "))

	subp.call(list_subp_args)
	subp.call("../Event_generation/Gamma_generation/essai.exe")
	script_utils.print_utility("Done")


list_bolo_name = ["FID824", "FID825", "FID827", "FID837", "FID838", "FID839", "FID841", "FID842"]
for bolo_name in list_bolo_name :
	script_utils.print_utility("Processing bolo " + bolo_name)
	launch_analysis(bolo_name)
