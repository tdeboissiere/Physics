#! /usr/bin/env python

from ROOT import *
import BDT_file_handler as BDT_fh
import PyROOTPlots as PyRPl
import os
import script_utils as script_utils
import numpy as np
import Analysis_utilities as Ana_ut


def check_scaling(bolo_name):

	"""Code to make some checks before applying 
	    the  code for the scaling of the backgrounds
    
	Detail:
		Compute the expected number of background events
		which pass the standard quality cuts and have EC>1

		Strategy: compute # of events in a given energy band 
		Use normalised bckg model : 
		expected = #observed in band / integral of band

		We get the total number of expected events for the 
		relevant exposure

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		voidt
	"""

	#Generic path
	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name

	#Load background files
	# data_types = {"names": "EC", "formats": "f"}
	
	arr_heat     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_heatonly.txt", delimiter=",") 
	arr_FidGamma = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_FidGamma.txt", delimiter=",") 
	arr_S1Gamma  = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Gamma.txt", delimiter=",") 
	arr_S2Gamma  = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Gamma.txt", delimiter=",") 
	arr_S1Beta   = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Beta.txt", delimiter=",") 
	arr_S2Beta   = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Beta.txt", delimiter=",") 	
	arr_S1Pb     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Pb.txt", delimiter=",") 
	arr_S2Pb     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Pb.txt", delimiter=",") 
	
	#Load background models, compute integral
	fS1Beta, files_beta = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Beta_spectrum_extrapol.root", "S1Beta_extra")
	fS2Beta = files_beta.Get("S2Beta_extra")
	fS1Pb, files_Pb = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Pb_spectrum_extrapol.root", "S1Pb_extra")
	fS2Pb = files_Pb.Get("S2Pb_extra")

	fFidGamma, file_Gamma = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Gamma_spectrum_extrapol.root", "FidGamma_extra")
	fS1Gamma = file_Gamma.Get("S1Gamma_extra")
	fS2Gamma = file_Gamma.Get("S2Gamma_extra")

	fheat, file_heat = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_heatonly_spectrum_extrapol.root", "heat_extra")

	#The code below verifies that we are able to ~ accurately predict
	#the rate in a band, based on the scaling of the rate found in another band

	sum_S1Beta, sum_S2Beta   = float(fS1Beta.Integral(5,20)), float(fS1Beta.Integral(5,20))
	sum_S1Pb, sum_S2Pb       = float(fS1Pb.Integral(5,60)), float(fS2Pb.Integral(5,60))
	sum_S1Gamma, sum_S2Gamma = float(fS1Gamma.Integral(5,15)), float(fS2Gamma.Integral(5,15))
	sum_FidGamma, sum_heat   = float(fFidGamma.Integral(5,15)), float(fheat.Integral(1,15))
	Beta_inf, Beta_sup   = 5,6
	Pb_inf, Pb_sup       = 15,40
	Gamma_inf, Gamma_sup = 6,15
	heat_inf, heat_sup   = 3,5

	list_S1Beta_band   = [elem for elem in arr_S1Beta if Beta_inf<elem<Beta_sup]
	list_S2Beta_band   = [elem for elem in arr_S2Beta if Beta_inf<elem<Beta_sup]
	
	list_S1Pb_band     = [elem for elem in arr_S1Pb if Pb_inf<elem<Pb_sup]
	list_S2Pb_band     = [elem for elem in arr_S2Pb if Pb_inf<elem<Pb_sup]
	
	list_S1Gamma_band  = [elem for elem in arr_S1Gamma if Gamma_inf<elem<Gamma_sup]
	list_S2Gamma_band  = [elem for elem in arr_S2Gamma if Gamma_inf<elem<Gamma_sup]
	
	list_FidGamma_band = [elem for elem in arr_FidGamma if Gamma_inf<elem<Gamma_sup]
	list_heat_band     = [elem for elem in arr_heat if heat_inf<elem<heat_sup]
	
	# Count the real number of observed events in the band, compare it to the expected one
	len_S1Beta_check   = len([elem for elem in arr_S1Beta if 5<elem<20])
	len_S2Beta_check   = len([elem for elem in arr_S2Beta if 5<elem<20])
	
	len_S1Pb_check     = len([elem for elem in arr_S1Pb if 5<elem<60])
	len_S2Pb_check     = len([elem for elem in arr_S2Pb if 5<elem<60])
	
	len_S1Gamma_check  = len([elem for elem in arr_S1Gamma if 5<elem<15])
	len_S2Gamma_check  = len([elem for elem in arr_S2Gamma if 5<elem<15])
	
	len_FidGamma_check = len([elem for elem in arr_FidGamma if 5<elem<15])
	len_heat_check     = len([elem for elem in arr_heat if 1<elem<15])
	
	exp_S1Beta         = len(list_S1Beta_band)*sum_S1Beta/fS1Beta.Integral(Beta_inf,Beta_sup)
	exp_S2Beta         = len(list_S2Beta_band)*sum_S2Beta/fS2Beta.Integral(Beta_inf,Beta_sup)
	
	exp_S1Pb           = len(list_S1Pb_band)*sum_S1Pb/fS1Pb.Integral(Pb_inf,Pb_sup)
	exp_S2Pb           = len(list_S2Pb_band)*sum_S2Pb/fS2Pb.Integral(Pb_inf,Pb_sup)
	
	exp_S1Gamma        = len(list_S1Gamma_band)*sum_S1Gamma/fS1Gamma.Integral(Gamma_inf,Gamma_sup)
	exp_S2Gamma        = len(list_S2Gamma_band)*sum_S2Gamma/fS2Gamma.Integral(Gamma_inf,Gamma_sup)
	
	exp_FidGamma       = len(list_FidGamma_band)*sum_FidGamma/fFidGamma.Integral(Gamma_inf,Gamma_sup)
	exp_heat           = len(list_heat_band)*sum_heat/fheat.Integral(heat_inf,heat_sup)

	list_evt = ["S1Beta", "S2Beta", "S1Pb", "S2Pb", "S1Gamma", "S2Gamma", "FidGamma", "heatonly"]
	list_exp = [exp_S1Beta, exp_S2Beta, exp_S1Pb, exp_S2Pb, exp_S1Gamma, exp_S2Gamma, exp_FidGamma, exp_heat]
	list_check = [len_S1Beta_check, len_S2Beta_check, len_S1Pb_check, len_S2Pb_check, len_S1Gamma_check, len_S2Gamma_check, len_FidGamma_check, len_heat_check]
	for index, elem in enumerate(list_exp):
		print list_evt[index], elem, list_check[index]

	############################################################
	#Series of plots to show the models more or less fit the bg
	############################################################

	#Verif heat bckg
	class heat_bckg:
		def __call__( self, x, par ):
			return fheat.Eval(x[0])*par[0]
	funcheat = TF1("h", heat_bckg(), 1, 15, 1)
	funcheat.SetParameter(0,exp_heat/sum_heat)
	print funcheat.Integral(1,15)

	hheat = TH1F("hheat", "hheat",100, 1, 15)
	for elem in arr_heat:
		hheat.Fill(elem)
	hheat.Scale(1./0.14)

	hheat.Draw()
	funcheat.Draw("same")
	print hheat.Integral()
	raw_input()

	#verif fid gamma bckg
	class FidGamma_bckg:
		def __call__( self, x, par ):
			return fFidGamma.Eval(x[0])*par[0]
	funcFidGamma = TF1("h", FidGamma_bckg(), 5, 15, 1)
	funcFidGamma.SetParameter(0,exp_FidGamma/sum_FidGamma)
	print funcFidGamma.Integral(5,15)

	hFidGamma = TH1F("hFidGamma", "hFidGamma",100, 5, 15)
	for elem in arr_FidGamma:
		hFidGamma.Fill(elem)
	hFidGamma.Scale(1./0.1)

	hFidGamma.Draw()
	funcFidGamma.Draw("same")
	print hFidGamma.Integral()

	raw_input()
	
	#verif surf gamma bckg
	class S1Gamma_bckg:
		def __call__( self, x, par ):
			return fS1Gamma.Eval(x[0])*par[0]
	funcS1Gamma = TF1("h", S1Gamma_bckg(), 5, 15, 1)
	funcS1Gamma.SetParameter(0,exp_S1Gamma/sum_S1Gamma)
	print funcS1Gamma.Integral(5,15)

	hS1Gamma = TH1F("hS1Gamma", "hS1Gamma",100, 5, 15)
	for elem in arr_S1Gamma:
		hS1Gamma.Fill(elem)
	hS1Gamma.Scale(1./0.1)

	hS1Gamma.Draw()
	funcS1Gamma.Draw("same")
	print hS1Gamma.Integral()

	raw_input()

	#verif Pb gamma bckg
	class S1Pb_bckg:
		def __call__( self, x, par ):
			return fS1Pb.Eval(x[0])*par[0]
	funcS1Pb = TF1("h", S1Pb_bckg(), 5, 60, 1)
	funcS1Pb.SetParameter(0,exp_S1Pb/sum_S1Pb)
	print funcS1Pb.Integral(5,60)

	hS1Pb = TH1F("hS1Pb", "hS1Pb",100, 5, 60)
	for elem in arr_S1Pb:
		hS1Pb.Fill(elem)
	hS1Pb.Scale(1./0.55)

	hS1Pb.Draw()
	funcS1Pb.Draw("same")
	print hS1Pb.Integral()

	raw_input()

	#verif Beta gamma bckg
	class S1Beta_bckg:
		def __call__( self, x, par ):
			return fS1Beta.Eval(x[0])*par[0]
	funcS1Beta = TF1("h", S1Beta_bckg(), 5, 20, 1)
	funcS1Beta.SetParameter(0,exp_S1Beta/sum_S1Beta)
	print funcS1Beta.Integral(5,20)

	hS1Beta = TH1F("hS1Beta", "hS1Beta",100, 5, 20)
	for elem in arr_S1Beta:
		hS1Beta.Fill(elem)
	hS1Beta.Scale(1./0.15)

	hS1Beta.Draw()
	funcS1Beta.Draw("same")
	print hS1Beta.Integral()

	raw_input()

	#verif Beta gamma bckg
	class S2Beta_bckg:
		def __call__( self, x, par ):
			return fS2Beta.Eval(x[0])*par[0]
	funcS2Beta = TF1("h", S2Beta_bckg(), 5, 20, 1)
	funcS2Beta.SetParameter(0,exp_S2Beta/sum_S2Beta)
	print funcS2Beta.Integral(5,20)

	hS2Beta = TH1F("hS2Beta", "hS2Beta",100, 5, 20)
	for elem in arr_S2Beta:
		hS2Beta.Fill(elem)
	hS2Beta.Scale(1./0.15)

	hS2Beta.Draw()
	funcS2Beta.Draw("same")
	print hS2Beta.Integral()

	raw_input()


def correct_scaling(bolo_name, analysis_type):

	"""Make event selection with usual cuts 
	for true and simulated data. Compare population size.
    
	Detail:
		Use an analysis type such that: no cut on vet, no cut on ion 
		Reproduce the same settings as the standard bckg_scaling no ion cut file.

		Verify that we obtain the same (modulo Poisson) number of events.

		Modify the scaling file accordingly if that is not the case

	Args:
		bolo_name = (str) bolometer name
		analysis_type  = (str) analysis_type

	Returns:
		void

	Raises:
		void
	"""

	BDT_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_better/BDT_" + bolo_name + "/" + analysis_type + "/"

	ttrue,ftrue = PyRPl.open_ROOT_object("../Fond_ERA_merged/" + bolo_name + "_" + analysis_type + "_fond.root", "data")
	tsimu, fsimu = PyRPl.open_ROOT_object(BDT_path +"True_events/ROOT_files/" + bolo_name + "_true_events_tree.root", "t_new1")

	print "all events true", ttrue.GetEntries()
	print "all events simu", tsimu.GetEntries()

	#Load estimators
	d_est = BDT_fh.open_estimator_file(bolo_name)
	
	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name)
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])
	

	l_FidGamma = TEventList("l_FidGamma")
	l_S1Gamma  = TEventList("l_S1Gamma")
	l_S2Gamma  = TEventList("l_S2Gamma")
	
	l_S1Beta   = TEventList("l_S1Beta")
	l_S2Beta   = TEventList("l_S2Beta")
	
	l_S1Pb     = TEventList("l_S1Pb")
	l_S2Pb     = TEventList("l_S2Pb")

	# Triples are not an issue when cutting at 4 sigma on veto !
	ttrue.Draw(">>l_triple"," ((EIA>1.5 && EIB>1.5 && EID>1.5)|| (EIB>1.5 && EIA>1.5 && EIC>1.5) || (EIB>1.5 && EIC>1.5 && EID>1.5) || (EIA>1.5 && EIC>1.5 && EID>1.5) )")
	print "expected triple:" , l_triple.GetN()

	EC_cut = "EC1>1.5 && EC1<15"

	ttrue.Draw(">>l_heatonly" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] )
	ttrue.Draw(">>l_FidGamma" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_FID"]  + ">0.7")
	ttrue.Draw(">>l_S1Gamma" ,EC_cut +" && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" + d_est["Q_S1"] + ">0.65")
	ttrue.Draw(">>l_S2Gamma" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" + d_est["Q_S2"] + ">0.65")
	ttrue.Draw(">>l_S1Beta" ,EC_cut +" && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" +  d_est["Q_S1"] + "<0.65 && " + d_est["Q_S1"] + ">0.2")
	ttrue.Draw(">>l_S2Beta" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_S2"] + "<0.65 && " + d_est["Q_S2"] + ">0.2")
	ttrue.Draw(">>l_S1Pb" ,EC_cut +" && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" +  d_est["Q_S1"] + "<0.15 &&" +  d_est["Q_S1"] + ">0.04")
	ttrue.Draw(">>l_S2Pb" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_S2"] + "<0.15 &&" +  d_est["Q_S2"] + ">0.04")

	list_type = ["heatonly", "FidGamma", "S1Gamma", "S2Gamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb"]
	list_list = [l_heatonly, l_FidGamma, l_S1Gamma, l_S2Gamma, l_S1Beta, l_S2Beta, l_S1Pb, l_S2Pb]
	list_num_evt = [l.GetN() for l in list_list]

	sum_evt = sum(list_num_evt)

	for evt, list_evt in zip(list_type, list_list):
		print evt, list_evt.GetN() #/float(sum_evt)


	lsimu_FidGamma = TEventList("lsimu_FidGamma")
	lsimu_S1Gamma  = TEventList("lsimu_S1Gamma")
	lsimu_S2Gamma  = TEventList("lsimu_S2Gamma")
	
	lsimu_S1Beta   = TEventList("lsimu_S1Beta")
	lsimu_S2Beta   = TEventList("lsimu_S2Beta")
	
	lsimu_S1Pb     = TEventList("lsimu_S1Pb")
	lsimu_S2Pb     = TEventList("lsimu_S2Pb")

	tsimu.Draw(">>lsimu_heatonly" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] )
	tsimu.Draw(">>lsimu_FidGamma" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_FID"]  + ">0.7")
	tsimu.Draw(">>lsimu_S1Gamma" ,EC_cut +" && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" + d_est["Q_S1"] + ">0.65")
	tsimu.Draw(">>lsimu_S2Gamma" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" + d_est["Q_S2"] + ">0.65")
	tsimu.Draw(">>lsimu_S1Beta" ,EC_cut +" && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" +  d_est["Q_S1"] + "<0.65 && " + d_est["Q_S1"] + ">0.2")
	tsimu.Draw(">>lsimu_S2Beta" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_S2"] + "<0.65 && " + d_est["Q_S2"] + ">0.2")
	tsimu.Draw(">>lsimu_S1Pb" ,EC_cut +" && EIA>" + d_std["FWIA"] +" && EIB>" + d_std["FWIB"] +"&& EIC<" + d_std["FWIC"] +"&& EID<" + d_std["FWID"] + " &&" +  d_est["Q_S1"] + "<0.15 &&" +  d_est["Q_S1"] + ">0.04")
	tsimu.Draw(">>lsimu_S2Pb" ,EC_cut +" && EIA<" + d_std["FWIA"] +" && EIB<" + d_std["FWIB"] +"&& EIC>" + d_std["FWIC"] +"&& EID>" + d_std["FWID"] + " &&" +  d_est["Q_S2"] + "<0.15 &&" +  d_est["Q_S2"] + ">0.04")
	
	list_listsimu = [lsimu_heatonly, lsimu_FidGamma, lsimu_S1Gamma, lsimu_S2Gamma, lsimu_S1Beta, lsimu_S2Beta, lsimu_S1Pb, lsimu_S2Pb]

	sum_evt = 0

	for evt, list_evt in zip(list_type, list_listsimu):
		print evt, list_evt.GetN()
		sum_evt+=list_evt.GetN()

def get_scaling(bolo_name, duration, heat_inf):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:
		Compute the expected number of background events
		which pass the standard quality cuts and have EC>1

		Strategy: compute # of events in a given energy band 
		Use normalised bckg model : 
		expected = #observed in band / integral of band

		We get the total number of expected events for the 
		relevant exposure

	Args:
		bolo_name = (str) bolometer name
		duration  = (float) data taking duration in days

	Returns:
		void

	Raises:
		voidt
	"""

	#Generic path
	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")
	
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_PSA_and_fond.root", "data")
	
	
	#Load background files
	arr_heat     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_heatonly.txt", delimiter=",") 
	arr_FidGamma = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_FidGamma.txt", delimiter=",") 
	arr_S1Gamma  = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Gamma.txt", delimiter=",") 
	arr_S2Gamma  = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Gamma.txt", delimiter=",") 
	arr_S1Beta   = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Beta.txt", delimiter=",") 
	arr_S2Beta   = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Beta.txt", delimiter=",") 	
	arr_S1Pb     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Pb.txt", delimiter=",") 
	arr_S2Pb     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Pb.txt", delimiter=",") 
	
	#Load background models, compute integral
	fS1Beta, files_beta = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Beta_spectrum_extrapol.root", "S1Beta_extra")
	fS2Beta = files_beta.Get("S2Beta_extra")
	fS1Pb, files_Pb = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Pb_spectrum_extrapol.root", "S1Pb_extra")
	fS2Pb = files_Pb.Get("S2Pb_extra")

	fFidGamma, file_Gamma = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Gamma_spectrum.root", "FidGamma")
	fS1Gamma = file_Gamma.Get("S1Gamma")
	fS2Gamma = file_Gamma.Get("S2Gamma")

	fheat, file_heat = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_heatonly_spectrum_extrapol.root", "heat_extra")

	#We have used check_scaling above to show that it is reasonable
	#to extrapolate the number of counts of the background based on the rate measured 
	#in the interval xx_inf, xx_sup
	#we now use the same interval xx_inf, xx_sup to measure the expected counts in the 1, 15 interval

	sum_S1Beta, sum_S2Beta   = float(fS1Beta.Integral(ECinf,15)), float(fS1Beta.Integral(ECinf,15))
	sum_S1Pb, sum_S2Pb       = float(fS1Pb.Integral(ECinf,15)), float(fS2Pb.Integral(ECinf,15))
	sum_S1Gamma, sum_S2Gamma = float(fS1Gamma.Integral(ECinf,15)), float(fS2Gamma.Integral(ECinf,15))
	sum_FidGamma, sum_heat   = float(fFidGamma.Integral(ECinf,15)), float(fheat.Integral(ECinf,15))

	Beta_inf, Beta_sup   = 5,6
	Pb_inf, Pb_sup       = 15,40
	Gamma_inf, Gamma_sup = 9,12
	heat_inf, heat_sup   = 0.5,15

	list_S1Beta_band   = [elem for elem in arr_S1Beta if Beta_inf<elem<Beta_sup]
	list_S2Beta_band   = [elem for elem in arr_S2Beta if Beta_inf<elem<Beta_sup]
	
	list_S1Pb_band     = [elem for elem in arr_S1Pb if Pb_inf<elem<Pb_sup]
	list_S2Pb_band     = [elem for elem in arr_S2Pb if Pb_inf<elem<Pb_sup]
	
	list_S1Gamma_band  = [elem for elem in arr_S1Gamma if Gamma_inf<elem<Gamma_sup]
	list_S2Gamma_band  = [elem for elem in arr_S2Gamma if Gamma_inf<elem<Gamma_sup]
	
	list_FidGamma_band = [elem for elem in arr_FidGamma if Gamma_inf<elem<Gamma_sup]
	list_heat_band     = [elem for elem in arr_heat if heat_inf<elem<heat_sup]
		
	exp_S1Beta         = len(list_S1Beta_band)*sum_S1Beta/fS1Beta.Integral(Beta_inf,Beta_sup)
	exp_S2Beta         = len(list_S2Beta_band)*sum_S2Beta/fS2Beta.Integral(Beta_inf,Beta_sup)
	
	exp_S1Pb           = len(list_S1Pb_band)*sum_S1Pb/fS1Pb.Integral(Pb_inf,Pb_sup)
	exp_S2Pb           = len(list_S2Pb_band)*sum_S2Pb/fS2Pb.Integral(Pb_inf,Pb_sup)
	
	exp_S1Gamma        = len(list_S1Gamma_band)*sum_S1Gamma/fS1Gamma.Integral(Gamma_inf,Gamma_sup)
	exp_S2Gamma        = len(list_S2Gamma_band)*sum_S2Gamma/fS2Gamma.Integral(Gamma_inf,Gamma_sup)
	
	exp_FidGamma       = len(list_FidGamma_band)*sum_FidGamma/fFidGamma.Integral(Gamma_inf,Gamma_sup)
	exp_heat           = len(list_heat_band)*sum_heat/fheat.Integral(heat_inf,heat_sup)

	exp_heat = 2*tree.GetEntries(standard_cuts + " && pearsA>0.9 && dtw_valA<100 && 0.5*(EIB+EID)<0 && 0.5*(EC1+EC2)>" + str(ECinf))
	# tree.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts)
	# print exp_heat
	# raw_input()

	#Load duration
	duration = 0
	file_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
	duration_file_name =bolo_name + "_duration_KTH_cut.txt" 
	duration         = 0
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

	with open(file_path_name + bolo_name + "_bckg_scaling_PSA_no_ion_cut_" + str(ECinf) + "_heat_thresh.txt", "w") as fs:
		fs.write("#no ion cut, 0.5*(EC1+EC2)>0, rate in days, " + str(ECinf)  + " to 15 keV in heat\n")
		fs.write("prop_heatonly,prop_FidGamma,prop_S1Gamma,prop_S2Gamma,prop_S1Beta,prop_S2Beta,prop_S1Pb,prop_S2Pb\n")
		fs.write(",".join([str(elem) for elem in list_pct]) + "\n")
		fs.write("rate_heatonly,rate_FidGamma,rate_S1Gamma,rate_S2Gamma,rate_S1Beta,rate_S2Beta,rate_S1Pb,rate_S2Pb\n")
		fs.write(",".join([str(elem) for elem in list_rate]) + "\n")



	# #Verif heat bckg
	# class heat_bckg:
	# 	def __call__( self, x, par ):
	# 		return fheat.Eval(x[0])*par[0]
	# funcheat = TF1("h", heat_bckg(), 1, 15, 1)
	# funcheat.SetParameter(0,exp_heat/sum_heat)
	# print "func", funcheat.Integral(1,15)

	# hheat = TH1F("hheat", "hheat",100, 1, 15)
	# for elem in arr_heat:
	# 	hheat.Fill(elem)
	# hheat.Scale(1./(0.14))

	# hheat.Draw()
	# funcheat.Draw("same")
	# print "hist", hheat.Integral()
	# raw_input()

def get_scaling_with_neutron(bolo_name, duration, heat_inf):

	"""Get the scaling in a [heat_inf, 15] keVee (heat) box
    
	Detail:
		Compute the expected number of background events
		which pass the standard quality cuts and have EC>1

		Strategy: compute # of events in a given energy band 
		Use normalised bckg model : 
		expected = #observed in band / integral of band

		We get the total number of expected events for the 
		relevant exposure

	Args:
		bolo_name = (str) bolometer name
		duration  = (float) data taking duration in days

	Returns:
		void

	Raises:
		voidt
	"""

	#Generic path
	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/" + "Analyse_" + bolo_name
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")
	
	tree, file_tree = PyRPl.open_ROOT_object("../Fond_ERA_merged/"+bolo_name+"_PSA_and_fond.root", "data")
	
	
	#Load background files
	arr_heat     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_heatonly.txt", delimiter=",") 
	arr_FidGamma = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_FidGamma.txt", delimiter=",") 
	arr_S1Gamma  = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Gamma.txt", delimiter=",") 
	arr_S2Gamma  = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Gamma.txt", delimiter=",") 
	arr_S1Beta   = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Beta.txt", delimiter=",") 
	arr_S2Beta   = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Beta.txt", delimiter=",") 	
	arr_S1Pb     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S1Pb.txt", delimiter=",") 
	arr_S2Pb     = np.loadtxt(gen_path+ "/Populations/Pop_for_scaling/" + bolo_name +"_S2Pb.txt", delimiter=",") 
	
	#Load background models, compute integral
	fS1Beta, files_beta = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Beta_spectrum_extrapol.root", "S1Beta_extra")
	fS2Beta = files_beta.Get("S2Beta_extra")
	fS1Pb, files_Pb = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Pb_spectrum_extrapol.root", "S1Pb_extra")
	fS2Pb = files_Pb.Get("S2Pb_extra")

	fFidGamma, file_Gamma = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_Gamma_spectrum.root", "FidGamma")
	fS1Gamma = file_Gamma.Get("S1Gamma")
	fS2Gamma = file_Gamma.Get("S2Gamma")

	fheat, file_heat = PyRPl.open_ROOT_object(gen_path + "/ROOT_files/" + bolo_name + "_heatonly_spectrum_extrapol.root", "heat_extra")

	#We have used check_scaling above to show that it is reasonable
	#to extrapolate the number of counts of the background based on the rate measured 
	#in the interval xx_inf, xx_sup
	#we now use the same interval xx_inf, xx_sup to measure the expected counts in the 1, 15 interval

	sum_S1Beta, sum_S2Beta   = float(fS1Beta.Integral(ECinf,15)), float(fS1Beta.Integral(ECinf,15))
	sum_S1Pb, sum_S2Pb       = float(fS1Pb.Integral(ECinf,15)), float(fS2Pb.Integral(ECinf,15))
	sum_S1Gamma, sum_S2Gamma = float(fS1Gamma.Integral(ECinf,15)), float(fS2Gamma.Integral(ECinf,15))
	sum_FidGamma, sum_heat   = float(fFidGamma.Integral(ECinf,15)), float(fheat.Integral(ECinf,15))

	Beta_inf, Beta_sup   = 5,6
	Pb_inf, Pb_sup       = 15,40
	Gamma_inf, Gamma_sup = 9,12
	heat_inf, heat_sup   = 0.5,15

	list_S1Beta_band   = [elem for elem in arr_S1Beta if Beta_inf<elem<Beta_sup]
	list_S2Beta_band   = [elem for elem in arr_S2Beta if Beta_inf<elem<Beta_sup]
	
	list_S1Pb_band     = [elem for elem in arr_S1Pb if Pb_inf<elem<Pb_sup]
	list_S2Pb_band     = [elem for elem in arr_S2Pb if Pb_inf<elem<Pb_sup]
	
	list_S1Gamma_band  = [elem for elem in arr_S1Gamma if Gamma_inf<elem<Gamma_sup]
	list_S2Gamma_band  = [elem for elem in arr_S2Gamma if Gamma_inf<elem<Gamma_sup]
	
	list_FidGamma_band = [elem for elem in arr_FidGamma if Gamma_inf<elem<Gamma_sup]
	list_heat_band     = [elem for elem in arr_heat if heat_inf<elem<heat_sup]
		
	exp_S1Beta         = len(list_S1Beta_band)*sum_S1Beta/fS1Beta.Integral(Beta_inf,Beta_sup)
	exp_S2Beta         = len(list_S2Beta_band)*sum_S2Beta/fS2Beta.Integral(Beta_inf,Beta_sup)
	
	exp_S1Pb           = len(list_S1Pb_band)*sum_S1Pb/fS1Pb.Integral(Pb_inf,Pb_sup)
	exp_S2Pb           = len(list_S2Pb_band)*sum_S2Pb/fS2Pb.Integral(Pb_inf,Pb_sup)
	
	exp_S1Gamma        = len(list_S1Gamma_band)*sum_S1Gamma/fS1Gamma.Integral(Gamma_inf,Gamma_sup)
	exp_S2Gamma        = len(list_S2Gamma_band)*sum_S2Gamma/fS2Gamma.Integral(Gamma_inf,Gamma_sup)
	
	exp_FidGamma       = len(list_FidGamma_band)*sum_FidGamma/fFidGamma.Integral(Gamma_inf,Gamma_sup)
	exp_heat           = len(list_heat_band)*sum_heat/fheat.Integral(heat_inf,heat_sup)

	exp_heat = 2*tree.GetEntries(standard_cuts + " && pearsA>0.9 && dtw_valA<100 && 0.5*(EIB+EID)<0 && 0.5*(EC1+EC2)>" + str(ECinf))
	# tree.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)", standard_cuts)
	# print exp_heat
	# raw_input()

	#Load duration
	duration = 0
	file_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
	duration_file_name =bolo_name + "_duration_KTH_cut.txt" 
	duration         = 0
	assert(os.path.isfile(file_path_name + duration_file_name) )
	with open(file_path_name + duration_file_name, 'r') as duration_file: 
	    list_duration_lines                                   = duration_file.readlines()[0].rstrip().split(",") 
	    duration                                              = float(list_duration_lines[3])
	    script_utils.print_utility(script_utils.COL("Duration (in days) = "  + str(duration) , "blue"))

	conv_path    = "/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/conv_EE_to_NR.root"
	func_conv, fileconv = PyRPl.open_ROOT_object(conv_path, "conv")
	hneut, fneut = PyRPl.open_ROOT_object("../Analyse_"+ bolo_name + "/ROOT_files/" + bolo_name + "_neutron_simu_hist.root", "hneut")

	lower_bound_ER= func_conv.Eval(ECinf)
	lower_bin, upper_bin = hneut.FindBin(lower_bound_ER), 200 #necessarily integrate up to 20 keVNR =~ 10keVee == bin 200

	import sys
	print   hneut.Integral(lower_bin, upper_bin)
	print   hneut.Integral(lower_bin,hneut.FindBin(7))*(24*200*0.6)/(8150)

	exp_neutron = hneut.Integral(lower_bin, upper_bin)*(duration)/(36.*365.)


	list_exp = [exp_heat, exp_FidGamma, exp_S1Gamma, exp_S2Gamma, exp_S1Beta, exp_S2Beta, exp_S1Pb, exp_S2Pb, exp_neutron]
	list_pct = [float(elem)/float(sum(list_exp)) for elem in list_exp]
	list_rate = [float(elem)/float(duration) for elem in list_exp]
	for index, elem in enumerate(list_exp):
		print elem, list_pct[index], list_rate[index]

	with open(file_path_name + bolo_name + "_bckg_scaling_with_neutron_PSA_no_ion_cut_" + str(ECinf) + "_heat_thresh.txt", "w") as fs:
		fs.write("#no ion cut, 0.5*(EC1+EC2)>0, rate in days, " + str(ECinf)  + " to 15 keV in heat\n")
		fs.write("prop_heatonly,prop_FidGamma,prop_S1Gamma,prop_S2Gamma,prop_S1Beta,prop_S2Beta,prop_S1Pb,prop_S2Pb,prop_neutron\n")
		fs.write(",".join([str(elem) for elem in list_pct]) + "\n")
		fs.write("rate_heatonly,rate_FidGamma,rate_S1Gamma,rate_S2Gamma,rate_S1Beta,rate_S2Beta,rate_S1Pb,rate_S2Pb,rate_neutron\n")
		fs.write(",".join([str(elem) for elem in list_rate]) + "\n")



	# #Verif heat bckg
	# class heat_bckg:
	# 	def __call__( self, x, par ):
	# 		return fheat.Eval(x[0])*par[0]
	# funcheat = TF1("h", heat_bckg(), 1, 15, 1)
	# funcheat.SetParameter(0,exp_heat/sum_heat)
	# print "func", funcheat.Integral(1,15)

	# hheat = TH1F("hheat", "hheat",100, 1, 15)
	# for elem in arr_heat:
	# 	hheat.Fill(elem)
	# hheat.Scale(1./(0.14))

	# hheat.Draw()
	# funcheat.Draw("same")
	# print "hist", hheat.Integral()
	# raw_input()

bolo_name = "FID837"
duration = 1
ECinf = 0.5
# analysis_type = "ana_1.5_min2_50"
analysis_type = "ana_0.5_0_5"

# check_scaling(bolo_name)
get_scaling(bolo_name, duration, ECinf)
# get_scaling_with_neutron(bolo_name, duration, ECinf)
# correct_scaling(bolo_name, analysis_type)