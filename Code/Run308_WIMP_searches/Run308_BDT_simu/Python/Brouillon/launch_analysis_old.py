#! /usr/bin/env python

from ROOT import *
import sys
import script_utils as script_utils
sys.path.append("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Scripts_ERA/") 
import create_extrapolated_bckg as create_bckg
import subprocess as subp
import BDT_utils as BDT_ut
import BDT_file_handler as BDT_fh

def launch_analysis(bolo_name, analysis_type, d_cut, d_overwrite, nsimu, exposure):

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

	Returns:
		void

	Raises:
		AssertionError 
	"""


	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"

	#Create empty directories needed for chosen analysis
	BDT_ut.create_BDT_directories(bolo_name, analysis_type)

	#First check if the files for Heat + Gamma simulations have been created
	
	Ana_path   = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
	file_Gamma = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_Gamma_spectrum_extrapol.root"
	file_heat = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_heatonly_spectrum_extrapol.root"
	try :
		assert( script_utils.bool_file_exist([file_heat, file_Gamma]) )
		if d_overwrite["Extra"]==1: create_bckg.create_extrapolated_bckg(bolo_name)
	except AssertionError:
		create_bckg.create_extrapolated_bckg(bolo_name)

	# Then run the Event generation scripts



	#Beta and Pb
	script_utils.print_utility("Beta and Pb")
	sys.path.append("../Event_generation/Beta_and_Pb_generation/")
	import build_BetaPb_tree as build_BetaPb

	BetaPb_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Beta_and_Pb/ROOT_files/" + bolo_name + "_"
	d_event_type_num_event ={"S1Beta":50000, "S2Beta":50000, "S1Pb":10000, "S2Pb":30000}
	list_BetaPb_files = [BetaPb_path + event_type + "_tree.root" for event_type in d_event_type_num_event.keys()]
	try :
		assert( script_utils.bool_file_exist(list_BetaPb_files) )
		if d_overwrite["BetaPb"] == 1:
			for event_type, num_event in d_event_type_num_event.iteritems():
				build_BetaPb.build_BetaPb_tree(bolo_name, analysis_type, event_type, d_cut, num_event)
	except AssertionError:
		for event_type, num_event in d_event_type_num_event.iteritems():
			build_BetaPb.build_BetaPb_tree(bolo_name, analysis_type, event_type, d_cut, num_event)
	script_utils.print_utility("Done")




	#Gamma
	script_utils.print_utility("Gamma")
	sys.path.append("../Event_generation/Gamma_generation/")
	# import build_Gamma_tree as build_Gamma
	import adjust_Gamma_file as adjust_Gamma
	d_event_type_num_event ={"S1Gamma":1000000, "S2Gamma":25000000, "FidGamma":10000}

	# adjust_Gamma.adjust_Gamma_file(bolo_name, d_cut, analysis_type, d_event_type_num_event)

	Gamma_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Gamma/ROOT_files/" + bolo_name + "_"
	list_Gamma_files = [Gamma_path + event_type + "_tree.root" for event_type in d_event_type_num_event.keys()]
	# try :
	# 	assert( script_utils.bool_file_exist(list_Gamma_files) )
	# 	if d_overwrite["Gamma"] == 1:
	# 		for event_type, num_event in d_event_type_num_event.iteritems():
	# 			build_Gamma.build_Gamma_tree(bolo_name, analysis_type, event_type, d_cut, num_event)
	# except AssertionError:
	# 	for event_type, num_event in d_event_type_num_event.iteritems():
	# 		build_Gamma.build_Gamma_tree(bolo_name, analysis_type, event_type, d_cut, num_event)
	# script_utils.print_utility("Done")

	#Use this complicated pattern to solve a mysterious ROOT-based TypeError (linked to Processline)

	#Get the correct root config options:
	cflags         = subp.Popen(["root-config", "--cflags"], stdout = subp.PIPE).communicate()[0].rstrip()
	libs           = subp.Popen(["root-config", "--libs"], stdout = subp.PIPE).communicate()[0].rstrip()

	list_subp_args = ["g++", "-o", "../Event_generation/Gamma_generation/essai.exe",  "../Event_generation/Gamma_generation/build_Gamma_tree.C", "-Wl,--no-as-needed"]
	list_subp_args.extend(cflags.split(" "))
	list_subp_args.extend(libs.split(" "))

	try :
		assert( script_utils.bool_file_exist(list_Gamma_files) )
		if d_overwrite["Gamma"] == 1:
			subp.call(list_subp_args)
			subp.call("../Event_generation/Gamma_generation/essai.exe")

	except AssertionError:
		subp.call(list_subp_args)
		subp.call("../Event_generation/Gamma_generation/essai.exe")

	# #Fill bckg_cuteff file "manually". Needed when using the compiled version of the Gamma generation script
	# tFidGamma, fFidGamma = PyRPl.open_ROOT_object("", "t_new")
	# tS1Gamma, fS1Gamma = PyRPl.open_ROOT_object("", "t_new")
	# tS2Gamma, fS2Gamma = PyRPl.open_ROOT_object("", "t_new")


	# #Heatonly
	# script_utils.print_utility("Heatonly")
	# sys.path.append("../Event_generation/Heatonly_generation/")
	# import build_heatonly_tree as build_heatonly
	# heatonly_num_event = 10000

	# heatonly_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/Heatonly/ROOT_files/" + bolo_name + "_"
	# list_heatonly_files = [heatonly_path + "heatonly_tree.root"]	
	# try :
	# 	assert( script_utils.bool_file_exist(list_heatonly_files) )
	# 	if d_overwrite["Heatonly"] == 1:
	# 		build_heatonly.build_heatonly_tree(bolo_name, analysis_type, FWHM_type, d_cut, heatonly_num_event)
	# except AssertionError:
	# 	build_heatonly.build_heatonly_tree(bolo_name, analysis_type, FWHM_type, d_cut, heatonly_num_event)		
	# script_utils.print_utility("Done")

	# # Create MVA scaling file adapted to the analysis type 
	# # (We chose to create a file to read the coefficients from
	# # instead of recomputing it every time. This is to keep track of the scaling)
	# d_cut_eff = BDT_fh.open_bckg_cuteff_file(bolo_name, analysis_type)
	# d_scaling = BDT_fh.open_scaling_file(bolo_name, d_cut)
	# BDT_fh.fill_MVA_scaling_file(bolo_name, analysis_type, d_scaling, d_cut_eff)





	# #Simulated data events
	# script_utils.print_utility("simulated data")
	# sys.path.append("../Event_generation/True_events_generation/")
	# import build_true_events_tree as build_true_events

	# true_events_path = gen_path + "BDT_" + bolo_name + "/" +analysis_type + "/True_events/ROOT_files/" + bolo_name + "_"
	# list_true_events_files = [true_events_path + "true_events_tree.root"]		
	# try :
	# 	assert( script_utils.bool_file_exist(list_true_events_files) )
	# 	if d_overwrite ["True"] == 1:
	# 		build_true_events.build_true_events_tree(bolo_name, analysis_type, d_cut, nsimu, exposure)
	# except AssertionError:
	# 	build_true_events.build_true_events_tree(bolo_name, analysis_type, d_cut, nsimu, exposure)
	# script_utils.print_utility("Done")





	# # #WIMP events
	# script_utils.print_utility("WIMP")
	# #First adjust WIMP files to current analysis
	# sys.path.append("../Event_generation/WIMP_generation/")
	# import adjust_WIMP_files as adjust_WIMP
	# adjust_WIMP.adjust_WIMP_files(bolo_name, d_cut, analysis_type)

	# WIMP_path = gen_path + "BDT_" + bolo_name + "/"+ analysis_type + "/WIMP/ROOT_files/"
	# list_WIMP_files = [WIMP_path + "recoils_mass_" + str(mass) + "_GeV.txt" for mass in [5, 6, 7, 10, 25]]
	# list_WIMP_filesbis = [WIMP_path + "recoils_mass_" + str(mass) + "_tree.root" for mass in [5, 6, 7, 10, 25]]
	# list_WIMP_files.extend(list_WIMP_filesbis)

	# #Use this complicated pattern to solve a mysterious ROOT-based TypeError (linked to Processline)
	# try :
	# 	assert( script_utils.bool_file_exist(list_WIMP_files) )
	# 	if d_overwrite["WIMP"] == 1:
	# 		gROOT.ProcessLine(".L " + gen_path +"Event_generation/WIMP_generation/WIMP_recoil_generation.C")
	# 		try :
	# 			gROOT.ProcessLine(main())
	# 		except TypeError:
	# 			pass
	# 		gROOT.ProcessLine(".L " + gen_path +"Event_generation/WIMP_generation/build_WIMP_tree.C")
	# 		try :
	# 			gROOT.ProcessLine(main())
	# 		except TypeError:
	# 			pass
	# except AssertionError:
	# 	gROOT.ProcessLine(".L " + gen_path +"Event_generation/WIMP_generation/WIMP_recoil_generation.C")
	# 	try :
	# 		gROOT.ProcessLine(main())
	# 	except TypeError:
	# 		pass
	# 	gROOT.ProcessLine(".L " + gen_path +"Event_generation/WIMP_generation/build_WIMP_tree.C")
	# 	try :
	# 		gROOT.ProcessLine(main())
	# 	except TypeError:
	# 		pass
	# script_utils.print_utility("Done")



	# #MVA analysis
	# script_utils.print_utility("Launching MVA analysis")
	# sys.path.append("../MVA/")
	# # #First adjust MVA files to current analysis
	# import adjust_MVA_files as adjust_MVA
	# adjust_MVA.adjust_MVA_files(bolo_name, analysis_type, nsimu)

	# #Get the correct root config options:
	# cflags         = subp.Popen(["root-config", "--cflags"], stdout = subp.PIPE).communicate()[0].rstrip()
	# libs           = subp.Popen(["root-config", "--libs"], stdout = subp.PIPE).communicate()[0].rstrip()

	# list_subp_args = ["g++", "-o", "../MVA/essai.exe",  "../MVA/TMVAClassification.C"]
	# list_subp_args.extend(cflags.split(" "))
	# list_subp_args.extend(libs.split(" "))
	# list_subp_args.append("-lTMVA")

	# # Launch TMVAClassification
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")

	# #Launch TMVAindividual
	# list_subp_args[3] = "../MVA/TMVAindividual.C"
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")

	# #Launch TMVAtruedata
	# list_subp_args[3] = "../MVA/TMVAtruedata.C"
	# subp.call(list_subp_args)
	# subp.call("../MVA/essai.exe")
	# script_utils.print_utility("Done")


bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
analysis_type = "ana_1.5_0_4"
#1 sigma on heat only band=  EFid>0.2
FWHM_type   = "standard_resolution"
d_cut       = {"ECinf": 1.5, "ECsup": 15, "EIinf": 0*0.2, "EIsup": 15, "sigma_vet": 4}
d_overwrite = {"Gamma": 1, "BetaPb": 0, "True":0, "Heatonly":0, "WIMP":0, "Extra":0}
nsimu       = 10
exposure    = 70.

launch_analysis(bolo_name, analysis_type, d_cut, d_overwrite, nsimu, exposure)