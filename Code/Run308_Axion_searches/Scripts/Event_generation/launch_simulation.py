#! /usr/bin/env python

from ROOT import *
import script_utils as script_utils
import subprocess as subp

def launch_analysis():
	"""Launch FidGamma simulations
	"""

	#Get the correct root config options for compilation:
	cflags         = subp.Popen(["root-config", "--cflags"], stdout = subp.PIPE).communicate()[0].rstrip()
	libs           = subp.Popen(["root-config", "--libs"], stdout = subp.PIPE).communicate()[0].rstrip()

	# ###########################
	# #Simulate fiducial gammas
	# ###########################
	script_utils.print_utility("Gamma")

	list_subp_args = ["g++", "-o", "./essai.exe",  "./build_Gamma_tree.C", "-Wl,--no-as-needed"]
	list_subp_args.extend(cflags.split(" "))
	list_subp_args.extend(libs.split(" "))

	subp.call(list_subp_args)
	subp.call("./essai.exe")
	script_utils.print_utility("Done")

	###########################
	#Simulate signal
	###########################
	script_utils.print_utility("Signal")

	list_subp_args = ["g++", "-o", "./essai.exe",  "./build_signal_tree.C", "-Wl,--no-as-needed"]
	list_subp_args.extend(cflags.split(" "))
	list_subp_args.extend(libs.split(" "))

	subp.call(list_subp_args)
	subp.call("./essai.exe")
	script_utils.print_utility("Done")



if __name__ == "__main__":

	launch_analysis()
