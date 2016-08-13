#! /usr/bin/env python

from ROOT import *
import sys
import script_utils as script_utils
sys.path.append("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Scripts_ERA/") 
import create_extrapolated_bckg as create_bckg
import subprocess as subp
import BDT_utils as BDT_ut
import BDT_file_handler as BDT_fh
import PyROOTPlots as PyRPl

def launch_analysis(bolo_name, analysis_type, d_cut, d_overwrite, nsimu, exposure, d_event_type_num_event, list_mass):

	"""Launch BDT simulations
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		d_cut         = (dict) indicates the cuts (inf/sup) on heat and ion 
		nsimu         = (int) number of true event simulations
		overwrite     = (dict) indicates if standard files need to be overwritten
        exposure      = (float) desired exposure
        d_event_type_num_event = (dict) indicates how many events to simulate
        list_mass = (list) list of WIMP masses

	Returns:
		void

	Raises:
		AssertionError 
	"""


	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"

	#Create empty directories needed for chosen analysis
	BDT_ut.create_BDT_directories(bolo_name, analysis_type, "corr")

	#First check if the files for Heat + Gamma simulations have been created
	Ana_path   = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
	file_Gamma = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Gamma_spectrum_extrapol.root"
	file_heat = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_heatonly_spectrum_extrapol.root"
	try :
		assert( script_utils.bool_file_exist([file_heat, file_Gamma]) )
		if d_overwrite["Extra"]==1: create_bckg.create_extrapolated_bckg(bolo_name)
	except AssertionError:
		create_bckg.create_extrapolated_bckg(bolo_name)

	# Relevant imports
	sys.path.append("../Event_generation/Beta_and_Pb_generation/")
	sys.path.append("../Event_generation/Gamma_generation/")
	sys.path.append("../Event_generation/Heatonly_generation/")
	sys.path.append("../Event_generation/True_events_generation/")
	sys.path.append("../Event_generation/WIMP_generation/")
	sys.path.append("../MVA/")
	import adjust_BetaPb_file as adjust_BetaPb
	import adjust_Gamma_file as adjust_Gamma
	import adjust_heat_file as adjust_heat
	import heatonly_KDE as heatonly_KDE
	import adjust_true_events_file as adjust_true_events
	import adjust_WIMP_files as adjust_WIMP
	import adjust_MVA_scaling_file as adjust_MVA_scaling
	import adjust_MVA_files as adjust_MVA

	#Get the correct root config options for compilation:
	cflags         = subp.Popen(["root-config", "--cflags"], stdout = subp.PIPE).communicate()[0].rstrip()
	libs           = subp.Popen(["root-config", "--libs"], stdout = subp.PIPE).communicate()[0].rstrip()

	#Prepare the bckg cut file
	list_evt_type = ["heatonly", "FidGamma", "S1Gamma", "S2Gamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb"]
	file_eff = Ana_path + "Analyse_" + bolo_name + "/Text_files/" + bolo_name + "_bckg_cuteff_" + analysis_type + ".txt"
	try :
		assert( script_utils.bool_file_exist([file_eff]) )
	except AssertionError:
		with open(file_eff, "w") as feff:
			for event_type in list_evt_type:
				feff.write(event_type + "\n")


	# ##############
	# #Beta and Pb
	# ##############
	# script_utils.print_utility("Beta and Pb")

	# adjust_BetaPb.adjust_BetaPb_file(bolo_name, d_cut, analysis_type, d_event_type_num_event)
	# BetaPb_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Beta_and_Pb/ROOT_files/" + bolo_name + "_"
	# list_BetaPb_files = [BetaPb_path + event_type + "_tree.root" for event_type in d_event_type_num_event.keys() if ("Beta" in event_type or "Pb" in event_type) ]

	# list_subp_args = ["g++", "-o", "../Event_generation/Beta_and_Pb_generation/essai.exe",  "../Event_generation/Beta_and_Pb_generation/build_BetaPb_tree.C", "-Wl,--no-as-needed"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))

	# try :
	# 	assert( script_utils.bool_file_exist(list_BetaPb_files) )
	# 	if d_overwrite["BetaPb"] == 1:
	# 		subp.call(list_subp_args)
	# 		subp.call("../Event_generation/Beta_and_Pb_generation/essai.exe")

	# except AssertionError:
	# 	subp.call(list_subp_args)
	# 	subp.call("../Event_generation/Beta_and_Pb_generation/essai.exe")
	# script_utils.print_utility("Done")


	# #########
	# #Gamma
	# #########
	# script_utils.print_utility("Gamma")

	# adjust_Gamma.adjust_Gamma_file(bolo_name, d_cut, analysis_type, d_event_type_num_event)
	# Gamma_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Gamma/ROOT_files/" + bolo_name + "_"
	# list_Gamma_files = [Gamma_path + event_type + "_tree.root" for event_type in d_event_type_num_event.keys() if "Gamma" in event_type]

	# list_subp_args = ["g++", "-o", "../Event_generation/Gamma_generation/essai.exe",  "../Event_generation/Gamma_generation/build_Gamma_tree.C", "-Wl,--no-as-needed"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))

	# try :
	# 	assert( script_utils.bool_file_exist(list_Gamma_files) )
	# 	if d_overwrite["Gamma"] == 1:
	# 		subp.call(list_subp_args)
	# 		subp.call("../Event_generation/Gamma_generation/essai.exe")

	# except AssertionError:
	# 	subp.call(list_subp_args)
	# 	subp.call("../Event_generation/Gamma_generation/essai.exe")
	# script_utils.print_utility("Done")

	# sys.exit()

	# ###########
	# #Heatonly
	# ############
	# script_utils.print_utility("Heatonly")

	# adjust_heat.adjust_heat_file(bolo_name, d_cut, analysis_type, d_event_type_num_event["heatonly"])
	# heatonly_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Heatonly/ROOT_files/" + bolo_name + "_"
	# heatonly_txt_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Heatonly/Text_files/" + bolo_name + "_"
	# list_heatonly_files = [heatonly_path + "heatonly_tree.root", heatonly_txt_path + "heatonly_2D_time.txt"]	

	# list_subp_args = ["g++", "-o", "../Event_generation/Heatonly_generation/essai.exe",  "../Event_generation/Heatonly_generation/build_heatonly_tree_fromhist.C", "-Wl,--no-as-needed"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))

	# try :
	# 	assert( script_utils.bool_file_exist(list_heatonly_files) )
	# 	if d_overwrite["Heatonly"] == 1:
	# 		# heatonly_KDE.generate_heatonly_from_KDE(bolo_name, analysis_type, d_event_type_num_event["heatonly"])
	# 		subp.call(list_subp_args)
	# 		subp.call("../Event_generation/Heatonly_generation/essai.exe")

	# except AssertionError:
	# 	# heatonly_KDE.generate_heatonly_from_KDE(bolo_name, analysis_type, d_event_type_num_event["heatonly"])
	# 	subp.call(list_subp_args)
	# 	subp.call("../Event_generation/Heatonly_generation/essai.exe")
	# script_utils.print_utility("Done")


	# ############################
	# #Simulated data events
	# ############################
	# script_utils.print_utility("simulated data")

	# adjust_true_events.adjust_true_events_file(bolo_name, d_cut, analysis_type, nsimu, exposure)
	# true_events_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/True_events/ROOT_files/" + bolo_name + "_"
	# list_true_events_files = [true_events_path + "true_events_tree.root"]	

	# list_subp_args = ["g++", "-o", "../Event_generation/True_events_generation/essai.exe",  "../Event_generation/True_events_generation/build_true_events_tree.C", "-Wl,--no-as-needed"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))

	# try :
	# 	assert( script_utils.bool_file_exist(list_true_events_files) )
	# 	if d_overwrite["True"] == 1:
	# 		subp.call(list_subp_args)
	# 		subp.call("../Event_generation/True_events_generation/essai.exe")

	# except AssertionError:
	# 	subp.call(list_subp_args)
	# 	subp.call("../Event_generation/True_events_generation/essai.exe")
	# script_utils.print_utility("Done")


	# ############################
	# #MVA scaling
	# ############################
	# script_utils.print_utility("MVA scaling")

	# MVA_scaling_path = Ana_path + "Analyse_" + bolo_name + "/Text_files/" + bolo_name + "_MVA_scaling_" +analysis_type + ".txt"
	# list_MVA_scaling_file = [MVA_scaling_path]	

	# list_subp_args = ["g++", "-o", "../Event_generation/True_events_generation/essai.exe",  "../Event_generation/True_events_generation/get_MVA_scaling.C", "-Wl,--no-as-needed"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))

	# try :
	# 	assert( script_utils.bool_file_exist(list_MVA_scaling_file) )
	# 	if d_overwrite["MVA_scaling"] == 1:
	# 		adjust_MVA_scaling.adjust_MVA_scaling_file(bolo_name, d_cut, analysis_type, 10000)
	# 		subp.call(list_subp_args)
	# 		subp.call("../Event_generation/True_events_generation/essai.exe")

	# except AssertionError:
	# 	adjust_MVA_scaling.adjust_MVA_scaling_file(bolo_name, d_cut, analysis_type, 10000)
	# 	subp.call(list_subp_args)
	# 	subp.call("../Event_generation/True_events_generation/essai.exe")
	# script_utils.print_utility("Done")


	#################
	# WIMP events
	# ###############
	script_utils.print_utility("WIMP")

	adjust_WIMP.adjust_WIMP_files(bolo_name, d_cut, analysis_type, d_event_type_num_event)

	#Generate the event tree
	list_subp_args = ["g++", "-o", "../Event_generation/WIMP_generation/essai.exe",  "../Event_generation/WIMP_generation/build_WIMP_tree.C", "-Wl,--no-as-needed"]
	list_subp_args.extend(cflags.split(" "))
	list_subp_args.extend(libs.split(" "))

	if d_overwrite["WIMP"] == 1:
		subp.call(list_subp_args)
		subp.call("../Event_generation/WIMP_generation/essai.exe")

	sys.exit()

	# #Generate the no cut event tree
	# list_subp_args = ["g++", "-o", "../Event_generation/WIMP_generation/essai.exe",  "../Event_generation/WIMP_generation/build_WIMP_tree_nocut.C", "-Wl,--no-as-needed"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))

	# if d_overwrite["WIMP_nocut"] == 1:
	# 	subp.call(list_subp_args)
	# 	subp.call("../Event_generation/WIMP_generation/essai.exe")

	# script_utils.print_utility("Done")



	# #############
	# #MVA analysis
	# #############
	script_utils.print_utility("Launching MVA analysis")
	adjust_MVA.adjust_MVA_files(bolo_name, analysis_type, nsimu)

	list_subp_args = ["g++", "-o", "../MVA/essai.exe",  "../MVA/TMVAClassification.C"]
	list_subp_args.extend(cflags.split(" "))
	list_subp_args.extend(libs.split(" "))
	list_subp_args.append("-lTMVA")

	# Launch TMVAClassification
	subp.call(list_subp_args)
	subp.call("../MVA/essai.exe")

	#Launch TMVAindividual
	list_subp_args[3] = "../MVA/TMVAindividual.C"
	subp.call(list_subp_args)
	subp.call("../MVA/essai.exe")

	# #Launch TMVAtruedata
	# list_subp_args[3] = "../MVA/TMVAtruedata.C"
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")
	# script_utils.print_utility("Done")

	# #Launch TMVAtruedata
	# list_subp_args[3] = "../MVA/TMVAtruedata_xcheck.C"
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")
	# script_utils.print_utility("Done")

	#Launch TMVArealtruedata (= real non simulated data)
	list_subp_args[3] = "../MVA/TMVArealtruedata.C"
	subp.call(list_subp_args)
	subp.call("../MVA/essai.exe")
	script_utils.print_utility("Done")

	# #Launch TMVArealtruedata (= real non simulated data with strict veto cut)
	# list_subp_args[3] = "../MVA/TMVArealtruedata_xcheck.C"
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")
	# script_utils.print_utility("Done")

	# #Launch TMVArealtrueneutron (= real non simulated data with strict veto cut)
	# list_subp_args[3] = "../MVA/TMVArealtrueneutron.C"
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")
	# script_utils.print_utility("Done")

bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
analysis_type = "ana_min2_min2_5"
FWHM_type   = "standard_resolution"
d_cut       = {"ECinf": -2, "ECsup": 15, "EIinf": -2, "EIsup": 15, "sigma_vet": 5}
d_overwrite = {"Gamma": 1, "BetaPb": 1, "Heatonly":1, "WIMP":1, "WIMP_nocut":0, "True":0, "MVA_scaling":1, "Extra":0}


#Heat cut at min2 min2 keV  (enormous increase of 3GeV data)
d_event_type_num_event ={"S1Beta":600000, "S2Beta":600000, "S1Pb":300000, "S2Pb":250000, 
	                     "S1Gamma":1200000, "S2Gamma":1200000, "FidGamma":80000, "heatonly":200000,
	                      "3GeV":40, "4GeV":4800, "5GeV":14000, "6GeV":800000, "7GeV":400000, "10GeV":250000, "25GeV":120000}



# #Heat cut at 0.5 keV  (enormous increase of 3GeV data)
# d_event_type_num_event ={"S1Beta":600000, "S2Beta":600000, "S1Pb":300000, "S2Pb":250000, 
# 	                     "S1Gamma":1200000, "S2Gamma":1200000, "FidGamma":160000, "heatonly":800000,
# 	                      "3GeV":80000000, "4GeV":4800000, "5GeV":1400000, "6GeV":800000, "7GeV":400000, "10GeV":250000, "25GeV":120000}


# #Heat cut at 1.5 keV
# d_event_type_num_event ={"S1Beta":1200000, "S2Beta":2000000, "S1Pb":300000, "S2Pb":280000, 
# 	                     "S1Gamma":20000000, "S2Gamma":100000000, "FidGamma":100000, "heatonly":5000000,
# 	                     "3GeV":1, "4GeV":1, "5GeV":200000000, "6GeV":10000000, "7GeV":5000000, "10GeV":750000, "25GeV":150000} 

                 		 

list_mass = [3,4,5,6,7,10,25]
# list_mass = [5,6,7,10,25]

nsimu       = 10
exposure    = 66

launch_analysis(bolo_name, analysis_type, d_cut, d_overwrite, nsimu, exposure, d_event_type_num_event, list_mass)