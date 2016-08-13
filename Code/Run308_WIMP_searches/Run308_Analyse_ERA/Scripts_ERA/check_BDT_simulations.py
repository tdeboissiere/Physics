import scipy
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import BDT_file_handler as BDT_fh
import matplotlib.pylab as plt
import statsmodels.api as sm 

def get_BDT_array(bolo_name, analysis_type, mass):

	""" Get the txt file for the 1D distributions of each event type
		Use this intermediate state because of a conflict between root and root_numpy
    
	Detail:

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		mass          = (list) WIMP mass

	Returns:
		void

	Raises:
		void
	"""

	from root_numpy import root2array
	
	BDT_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_better/BDT_" + bolo_name + "/" + analysis_type + "/"
	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	arr_true_events = root2array("../Fond_ERA_merged/" + bolo_name + "_" + analysis_type + "_lowmass_fond.root", "t_merged")
	arr_simu_events = root2array(BDT_path +"True_events/ROOT_files/" + bolo_name + "_true_events_tree.root", "t_new0")

	np.savetxt(pop_path + bolo_name + "_true_events_all.txt", arr_true_events, delimiter = ",")
	np.savetxt(pop_path + bolo_name + "_simu_events_all.txt", arr_simu_events, delimiter = ",")

	simu_heatonly = root2array(BDT_path +"Heatonly/ROOT_files/" + bolo_name + "_heatonly_tree.root", "t_new0")
	simu_FidGamma = root2array(BDT_path +"Gamma/ROOT_files/" + bolo_name + "_FidGamma_tree.root", "t_new0")
	simu_S1Gamma = root2array(BDT_path +"Gamma/ROOT_files/" + bolo_name + "_S1Gamma_tree.root", "t_new0")
	simu_S2Gamma = root2array(BDT_path +"Gamma/ROOT_files/" + bolo_name + "_S2Gamma_tree.root", "t_new0")
	simu_S1Beta = root2array(BDT_path +"Beta_and_Pb/ROOT_files/" + bolo_name + "_S1Beta_tree.root", "t_new0")
	simu_S2Beta = root2array(BDT_path +"Beta_and_Pb/ROOT_files/" + bolo_name + "_S2Beta_tree.root", "t_new0")
	simu_S1Pb = root2array(BDT_path +"Beta_and_Pb/ROOT_files/" + bolo_name + "_S1Pb_tree.root", "t_new0")
	simu_S2Pb = root2array(BDT_path +"Beta_and_Pb/ROOT_files/" + bolo_name + "_S2Pb_tree.root", "t_new0")

	list_evt = ["heatonly", "FidGamma", "S1Gamma", "S2Gamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb"]
	list_simu = [simu_heatonly, simu_FidGamma, simu_S1Gamma, simu_S2Gamma, simu_S1Beta, simu_S2Beta, simu_S1Pb, simu_S2Pb]

	for i in range(8):
		np.savetxt(pop_path + bolo_name + "_BDT_" + list_evt[i] + ".txt", list_simu[i], delimiter = ",")

def check_10keV_distribution(bolo_name):

	""" Compare the 10 keV resolution for heat events
    
	Detail:
		Do KS Test + plot histogram

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		mass          = (list) WIMP mass

	Returns:
		void

	Raises:
		void
	"""

	def KS_signif(len1, len2):
		"""
		95 CL os KS test 
		"""
		return  1.36*sqrt(len1 + len2)/sqrt(len1*len2)


	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}

	#t for true, s for simu

	t_FidGamma = np.loadtxt(pop_path + bolo_name + "_FidGamma_full_info.txt", delimiter=",",  dtype=data_types)
	s_FidGamma = np.loadtxt(pop_path + bolo_name + "_BDT_FidGamma.txt", delimiter=",",  dtype=data_types)

	t_FidGamma = t_FidGamma[np.where(np.logical_and(t_FidGamma["EC1"]>8, t_FidGamma["EC1"]<11) ) ]
	s_FidGamma = s_FidGamma[np.where(np.logical_and(s_FidGamma["EC1"]>8, s_FidGamma["EC1"]<11) ) ]

	#EC1
	t_EC1_FidGamma = t_FidGamma["EC1"]
	s_EC1_FidGamma = s_FidGamma["EC1"]

	print "EC1 ", scipy.stats.ks_2samp(t_EC1_FidGamma, s_EC1_FidGamma),"   ", KS_signif(len(t_EC1_FidGamma) , len(s_EC1_FidGamma))

	hEC1_t_FidGamma = TH1F("hEC1_t_FidGamma", "hEC1_t_FidGamma", 100, 8,11)
	hEC1_s_FidGamma = TH1F("hEC1_s_FidGamma", "hEC1_s_FidGamma", 100, 8,11)

	PyRPl.fill_TH1(hEC1_t_FidGamma, t_EC1_FidGamma)
	PyRPl.fill_TH1(hEC1_s_FidGamma, s_EC1_FidGamma)

	PyRPl.process_TH1(hEC1_t_FidGamma, X_title = "EC1", color = kRed)
	PyRPl.process_TH1(hEC1_s_FidGamma, X_title = "EC1", color = kBlack)

	hEC1_t_FidGamma.Scale(1./hEC1_t_FidGamma.Integral())
	hEC1_s_FidGamma.Scale(1./hEC1_s_FidGamma.Integral())

	hEC1_t_FidGamma.Draw("")
	hEC1_s_FidGamma.Draw("same")
	raw_input()

	#EC2
	t_EC2_FidGamma = t_FidGamma["EC2"]
	s_EC2_FidGamma = s_FidGamma["EC2"]

	print "EC2 ", scipy.stats.ks_2samp(t_EC2_FidGamma, s_EC2_FidGamma),"   ", KS_signif(len(t_EC2_FidGamma) , len(s_EC2_FidGamma))

	hEC2_t_FidGamma = TH1F("hEC2_t_FidGamma", "hEC2_t_FidGamma", 100, 8,11)
	hEC2_s_FidGamma = TH1F("hEC2_s_FidGamma", "hEC2_s_FidGamma", 100, 8,11)

	PyRPl.fill_TH1(hEC2_t_FidGamma, t_EC2_FidGamma)
	PyRPl.fill_TH1(hEC2_s_FidGamma, s_EC2_FidGamma)

	PyRPl.process_TH1(hEC2_t_FidGamma, X_title = "EC2", color = kRed)
	PyRPl.process_TH1(hEC2_s_FidGamma, X_title = "EC2", color = kBlack)

	hEC2_t_FidGamma.Scale(1./hEC2_t_FidGamma.Integral())
	hEC2_s_FidGamma.Scale(1./hEC2_s_FidGamma.Integral())

	hEC2_t_FidGamma.Draw("")
	hEC2_s_FidGamma.Draw("same")

	raw_input()

	#EIB
	t_EIB_FidGamma = t_FidGamma["EIB"]
	s_EIB_FidGamma = s_FidGamma["EIB"]

	print "EIB", scipy.stats.ks_2samp(t_EIB_FidGamma, s_EIB_FidGamma),"   ", KS_signif(len(t_EIB_FidGamma) , len(s_EIB_FidGamma))

	hEIB_t_FidGamma = TH1F("hEIB_t_FidGamma", "hEIB_t_FidGamma", 100, 8,11)
	hEIB_s_FidGamma = TH1F("hEIB_s_FidGamma", "hEIB_s_FidGamma", 100, 8,11)

	PyRPl.fill_TH1(hEIB_t_FidGamma, t_EIB_FidGamma)
	PyRPl.fill_TH1(hEIB_s_FidGamma, s_EIB_FidGamma)

	PyRPl.process_TH1(hEIB_t_FidGamma, X_title = "EIB", color = kRed)
	PyRPl.process_TH1(hEIB_s_FidGamma, X_title = "EIB", color = kBlack)

	hEIB_t_FidGamma.Scale(1./hEIB_t_FidGamma.Integral())
	hEIB_s_FidGamma.Scale(1./hEIB_s_FidGamma.Integral())

	hEIB_t_FidGamma.Draw("")
	hEIB_s_FidGamma.Draw("same")
	raw_input()

	#EID
	t_EID_FidGamma = t_FidGamma["EID"]
	s_EID_FidGamma = s_FidGamma["EID"]

	print "EID", scipy.stats.ks_2samp(t_EID_FidGamma, s_EID_FidGamma),"   ", KS_signif(len(t_EID_FidGamma) , len(s_EID_FidGamma))

	hEID_t_FidGamma = TH1F("hEID_t_FidGamma", "hEID_t_FidGamma", 100, 8,11)
	hEID_s_FidGamma = TH1F("hEID_s_FidGamma", "hEID_s_FidGamma", 100, 8,11)

	PyRPl.fill_TH1(hEID_t_FidGamma, t_EID_FidGamma)
	PyRPl.fill_TH1(hEID_s_FidGamma, s_EID_FidGamma)

	PyRPl.process_TH1(hEID_t_FidGamma, X_title = "EID", color = kRed)
	PyRPl.process_TH1(hEID_s_FidGamma, X_title = "EID", color = kBlack)

	hEID_t_FidGamma.Scale(1./hEID_t_FidGamma.Integral())
	hEID_s_FidGamma.Scale(1./hEID_s_FidGamma.Integral())

	hEID_t_FidGamma.Draw("")
	hEID_s_FidGamma.Draw("same")

	raw_input()

def check_veto_distribution(bolo_name):

	""" Verify the IA and IC distributions by looking at FID + Heat only
    
	Detail:
		Do KS Test + plot histogram

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		mass          = (list) WIMP mass

	Returns:
		void

	Raises:
		void
	"""

	def KS_signif(len1, len2):
		"""
		95 CL os KS test 
		"""
		return  1.36*sqrt(len1 + len2)/sqrt(len1*len2)


	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}

	#t for true, s for simu

	t_heatonly = np.loadtxt(pop_path + bolo_name + "_heatonly_full_info.txt", delimiter=",",  dtype=data_types)
	t_FidGamma = np.loadtxt(pop_path + bolo_name + "_FidGamma_full_info.txt", delimiter=",",  dtype=data_types)

	s_heatonly = np.loadtxt(pop_path + bolo_name + "_BDT_heatonly.txt", delimiter=",",  dtype=data_types)
	s_FidGamma = np.loadtxt(pop_path + bolo_name + "_BDT_FidGamma.txt", delimiter=",",  dtype=data_types)

	#IA
	t_IA_heatonly = t_heatonly["EIA"]
	s_IA_heatonly = s_heatonly["EIA"]

	t_IA_FidGamma = t_FidGamma["EIA"]
	s_IA_FidGamma = s_FidGamma["EIA"]


	print "heat to heat ", scipy.stats.ks_2samp(t_IA_heatonly, s_IA_heatonly),"   ", KS_signif(len(t_IA_heatonly) , len(s_IA_heatonly))
	print "fid to fid ", scipy.stats.ks_2samp(t_IA_FidGamma, s_IA_FidGamma),"   ", KS_signif(len(t_IA_FidGamma) , len(s_IA_FidGamma))

	hIA_t_heatonly = TH1F("hIA_t_heatonly", "hIA_t_heatonly", 50, -2, 2)
	hIA_s_heatonly = TH1F("hIA_s_heatonly", "hIA_s_heatonly", 50, -2, 2)
	hIA_t_FidGamma = TH1F("hIA_t_FidGamma", "hIA_t_FidGamma", 50, -2, 2)
	hIA_s_FidGamma = TH1F("hIA_s_FidGamma", "hIA_s_FidGamma", 50, -2, 2)

	PyRPl.fill_TH1(hIA_t_heatonly, t_IA_heatonly)
	PyRPl.fill_TH1(hIA_s_heatonly, s_IA_heatonly)
	PyRPl.fill_TH1(hIA_t_FidGamma, t_IA_FidGamma)
	PyRPl.fill_TH1(hIA_s_FidGamma, s_IA_FidGamma)

	PyRPl.process_TH1(hIA_t_heatonly, X_title = "IA", color = kRed)
	PyRPl.process_TH1(hIA_s_heatonly, X_title = "IA", color = kBlack)
	PyRPl.process_TH1(hIA_t_FidGamma, X_title = "IA", color = kBlue)
	PyRPl.process_TH1(hIA_s_FidGamma, X_title = "IA", color = kOrange-1)

	hIA_t_heatonly.Scale(1./hIA_t_heatonly.Integral())
	hIA_s_heatonly.Scale(1./hIA_s_heatonly.Integral())
	hIA_t_FidGamma.Scale(1./hIA_t_FidGamma.Integral())
	hIA_s_FidGamma.Scale(1./hIA_s_FidGamma.Integral())

	hIA_t_heatonly.Draw()
	hIA_s_heatonly.Draw("same")
	hIA_t_FidGamma.Draw("same")
	hIA_s_FidGamma.Draw("same")
	raw_input()

	#IC
	t_IC_heatonly = t_heatonly["EIC"]
	s_IC_heatonly = s_heatonly["EIC"]

	t_IC_FidGamma = t_FidGamma["EIC"]
	s_IC_FidGamma = s_FidGamma["EIC"]

	print "heat to heat ", scipy.stats.ks_2samp(t_IC_heatonly, s_IC_heatonly),"   ", KS_signif(len(t_IC_heatonly) , len(s_IC_heatonly))
	print "fid to fid ", scipy.stats.ks_2samp(t_IC_FidGamma, s_IC_FidGamma),"   ", KS_signif(len(t_IC_FidGamma) , len(s_IC_FidGamma))
	print "heat to fid ", scipy.stats.ks_2samp(t_IC_heatonly, s_IC_FidGamma),"   ", KS_signif(len(t_IC_heatonly) , len(s_IC_FidGamma))
	print "fid to heat ", scipy.stats.ks_2samp(t_IC_FidGamma, s_IC_heatonly),"   ", KS_signif(len(t_IC_FidGamma) , len(s_IC_heatonly))

	hIC_t_heatonly = TH1F("hIC_t_heatonly", "hIC_t_heatonly", 50, -2, 2)
	hIC_s_heatonly = TH1F("hIC_s_heatonly", "hIC_s_heatonly", 50, -2, 2)
	hIC_t_FidGamma = TH1F("hIC_t_FidGamma", "hIC_t_FidGamma", 50, -2, 2)
	hIC_s_FidGamma = TH1F("hIC_s_FidGamma", "hIC_s_FidGamma", 50, -2, 2)

	PyRPl.fill_TH1(hIC_t_heatonly, t_IC_heatonly)
	PyRPl.fill_TH1(hIC_s_heatonly, s_IC_heatonly)
	PyRPl.fill_TH1(hIC_t_FidGamma, t_IC_FidGamma)
	PyRPl.fill_TH1(hIC_s_FidGamma, s_IC_FidGamma)

	PyRPl.process_TH1(hIC_t_heatonly, X_title = "IC", color = kRed)
	PyRPl.process_TH1(hIC_s_heatonly, X_title = "IC", color = kBlack)
	PyRPl.process_TH1(hIC_t_FidGamma, X_title = "IC", color = kBlue)
	PyRPl.process_TH1(hIC_s_FidGamma, X_title = "IC", color = kOrange-1)

	hIC_t_heatonly.Scale(1./hIC_t_heatonly.Integral())
	hIC_s_heatonly.Scale(1./hIC_s_heatonly.Integral())
	hIC_t_FidGamma.Scale(1./hIC_t_FidGamma.Integral())
	hIC_s_FidGamma.Scale(1./hIC_s_FidGamma.Integral())

	hIC_t_heatonly.Draw()
	hIC_s_heatonly.Draw("same")
	hIC_t_FidGamma.Draw("same")
	hIC_s_FidGamma.Draw("same")

	raw_input()

def check_collectrode_distribution(bolo_name):

	""" Verify the IB and ID distributions by looking at FID + Heat only
    
	Detail:
		Do KS Test + plot histogram

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		mass          = (list) WIMP mass

	Returns:
		void

	Raises:
		void
	"""

	def KS_signif(len1, len2):
		"""
		95 CL os KS test 
		"""
		return  1.36*sqrt(len1 + len2)/sqrt(len1*len2)


	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}

	#t for true, s for simu

	t_heatonly = np.loadtxt(pop_path + bolo_name + "_heatonly_full_info.txt", delimiter=",",  dtype=data_types)
	t_S2Gamma  = np.loadtxt(pop_path + bolo_name + "_S2Gamma_full_info.txt", delimiter=",",  dtype=data_types)
	t_S2Beta   = np.loadtxt(pop_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)
	t_S1Gamma  = np.loadtxt(pop_path + bolo_name + "_S1Gamma_full_info.txt", delimiter=",",  dtype=data_types)
	t_S1Beta   = np.loadtxt(pop_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)

	s_heatonly = np.loadtxt(pop_path + bolo_name + "_BDT_heatonly.txt", delimiter=",",  dtype=data_types)
	s_S2Gamma = np.loadtxt(pop_path + bolo_name + "_BDT_S2Gamma.txt", delimiter=",",  dtype=data_types)
	s_S2Beta = np.loadtxt(pop_path + bolo_name + "_BDT_S2Beta.txt", delimiter=",",  dtype=data_types)
	s_S1Gamma = np.loadtxt(pop_path + bolo_name + "_BDT_S1Gamma.txt", delimiter=",",  dtype=data_types)
	s_S1Beta = np.loadtxt(pop_path + bolo_name + "_BDT_S1Beta.txt", delimiter=",",  dtype=data_types)
	# t_S2Gamma = t_S2Gamma[np.where(np.logical_and(t_S2Gamma["EC1"]>10,t_S2Gamma["EC1"]< 15))]
	# t_S2Beta = t_S2Beta[np.where(np.logical_and(t_S2Beta["EC1"]>10,t_S2Beta["EC1"]<15)) ]

	#IB
	t_IB_heatonly = t_heatonly["EIB"]
	s_IB_heatonly = s_heatonly["EIB"]

	t_IB_Surf = np.hstack( ( t_S2Gamma["EIB"], t_S2Beta["EIB"] ) )
	s_IB_Surf = np.hstack( ( s_S2Gamma["EIB"], s_S2Beta["EIB"] ) )

	print "heat to heat ", scipy.stats.ks_2samp(t_IB_heatonly, s_IB_heatonly),"   ", KS_signif(len(t_IB_heatonly) , len(s_IB_heatonly))
	print "fid to fid ", scipy.stats.ks_2samp(t_IB_Surf, s_IB_Surf),"   ", KS_signif(len(t_IB_Surf) , len(s_IB_Surf))

	hIB_t_heatonly = TH1F("hIB_t_heatonly", "hIB_t_heatonly", 200, -2, 2)
	hIB_s_heatonly = TH1F("hIB_s_heatonly", "hIB_s_heatonly", 200, -2, 2)
	hIB_t_Surf = TH1F("hIB_t_Surf", "hIB_t_Surf", 200, -2, 2)
	hIB_s_Surf = TH1F("hIB_s_Surf", "hIB_s_Surf", 200, -2, 2)

	PyRPl.fill_TH1(hIB_t_heatonly, t_IB_heatonly)
	PyRPl.fill_TH1(hIB_s_heatonly, s_IB_heatonly)
	PyRPl.fill_TH1(hIB_t_Surf, t_IB_Surf)
	PyRPl.fill_TH1(hIB_s_Surf, s_IB_Surf)

	PyRPl.process_TH1(hIB_t_heatonly, X_title = "IB", color = kRed)
	PyRPl.process_TH1(hIB_s_heatonly, X_title = "IB", color = kBlack)
	PyRPl.process_TH1(hIB_t_Surf, X_title = "IB", color = kBlue)
	PyRPl.process_TH1(hIB_s_Surf, X_title = "IB", color = kOrange-1)

	hIB_t_heatonly.Scale(1./hIB_t_heatonly.Integral())
	hIB_s_heatonly.Scale(1./hIB_s_heatonly.Integral())
	hIB_t_Surf.Scale(1./hIB_t_Surf.Integral())
	hIB_s_Surf.Scale(1./hIB_s_Surf.Integral())

	hIB_t_heatonly.Draw()
	hIB_s_heatonly.Draw("same")
	raw_input()
	hIB_t_Surf.Draw("same")
	hIB_s_Surf.Draw("same")
	raw_input()

	#ID
	t_ID_heatonly = t_heatonly["EID"]
	s_ID_heatonly = s_heatonly["EID"]

	t_ID_Surf = np.hstack( ( t_S1Gamma["EID"], t_S1Beta["EID"] ) )
	s_ID_Surf = np.hstack( ( s_S1Gamma["EID"], s_S1Beta["EID"] ) )

	print "heat to heat ", scipy.stats.ks_2samp(t_ID_heatonly, s_ID_heatonly),"   ", KS_signif(len(t_ID_heatonly) , len(s_ID_heatonly))
	print "fid to fid ", scipy.stats.ks_2samp(t_ID_Surf, s_ID_Surf),"   ", KS_signif(len(t_ID_Surf) , len(s_ID_Surf))

	hID_t_heatonly = TH1F("hID_t_heatonly", "hID_t_heatonly", 200, -2, 2)
	hID_s_heatonly = TH1F("hID_s_heatonly", "hID_s_heatonly", 200, -2, 2)
	hID_t_Surf = TH1F("hID_t_Surf", "hID_t_Surf", 200, -2, 2)
	hID_s_Surf = TH1F("hID_s_Surf", "hID_s_Surf", 200, -2, 2)

	PyRPl.fill_TH1(hID_t_heatonly, t_ID_heatonly)
	PyRPl.fill_TH1(hID_s_heatonly, s_ID_heatonly)
	PyRPl.fill_TH1(hID_t_Surf, t_ID_Surf)
	PyRPl.fill_TH1(hID_s_Surf, s_ID_Surf)

	PyRPl.process_TH1(hID_t_heatonly, X_title = "ID", color = kRed)
	PyRPl.process_TH1(hID_s_heatonly, X_title = "ID", color = kBlack)
	PyRPl.process_TH1(hID_t_Surf, X_title = "ID", color = kBlue)
	PyRPl.process_TH1(hID_s_Surf, X_title = "ID", color = kOrange-1)

	hID_t_heatonly.Scale(1./hID_t_heatonly.Integral())
	hID_s_heatonly.Scale(1./hID_s_heatonly.Integral())
	hID_t_Surf.Scale(1./hID_t_Surf.Integral())
	hID_s_Surf.Scale(1./hID_s_Surf.Integral())

	hID_t_heatonly.Draw()
	hID_s_heatonly.Draw("same")
	raw_input()
	hID_t_Surf.Draw("same")
	hID_s_Surf.Draw("same")
	raw_input()

def check_intrinsic_scattering(bolo_name):

	""" Look for intrinsic scattering differences between simu and data
    
	Detail:

	Args:
		bolo_name     = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}

	#t for true, s for simu

	t_FidGamma = np.loadtxt(pop_path + bolo_name + "_FidGamma_full_info.txt", delimiter=",",  dtype=data_types)
	s_FidGamma = np.loadtxt(pop_path + bolo_name + "_BDT_FidGamma.txt", delimiter=",",  dtype=data_types)
	t_heatonly = np.loadtxt(pop_path + bolo_name + "_heatonly_full_info.txt", delimiter=",",  dtype=data_types)
	s_heatonly = np.loadtxt(pop_path + bolo_name + "_BDT_heatonly.txt", delimiter=",",  dtype=data_types)

	# #FidGamma 1D distribution just for check
	# t_heat_dist = TH1F("t_heat_dist", "t_heat_dist", 100, 6, 12)
	# s_heat_dist = TH1F("s_heat_dist", "s_heat_dist", 100, 6, 12)

	# PyRPl.fill_TH1(t_heat_dist, t_FidGamma["EID"])
	# PyRPl.fill_TH1(s_heat_dist, s_FidGamma["EID"])

	# PyRPl.process_TH1(t_heat_dist, X_title = "EID", color = kRed)
	# PyRPl.process_TH1(s_heat_dist, X_title = "EID", color = kBlack)

	# t_heat_dist.Scale(1./t_heat_dist.Integral())
	# s_heat_dist.Scale(1./s_heat_dist.Integral())

	# s_heat_dist.Draw("")
	# t_heat_dist.Draw("same")
	# raw_input()

	#FidGamma
	l_t = np.where(np.logical_and(t_FidGamma["EC1"]>8 , t_FidGamma["EC2"]>8))
	l_s = np.where(np.logical_and(s_FidGamma["EC1"]>8 , s_FidGamma["EC2"]>8))

	t_FidGamma=t_FidGamma[l_t]
	s_FidGamma=s_FidGamma[l_s]

	true_dist = t_FidGamma["EC1"]-t_FidGamma["EC2"]
	simu_dist = s_FidGamma["EC1"]-s_FidGamma["EC2"]

	t_heat_dist = TH1F("t_heat_dist", "t_heat_dist", 100,-2,2)
	s_heat_dist = TH1F("s_heat_dist", "s_heat_dist", 100,-2,2)

	PyRPl.fill_TH1(t_heat_dist, true_dist)
	PyRPl.fill_TH1(s_heat_dist, simu_dist)

	PyRPl.process_TH1(t_heat_dist, X_title = "distance", color = kRed)
	PyRPl.process_TH1(s_heat_dist, X_title = "distance", color = kBlack)

	t_heat_dist.Scale(1./t_heat_dist.Integral())
	s_heat_dist.Scale(1./s_heat_dist.Integral())

	s_heat_dist.Draw("")
	# s_heat_dist.Fit("gaus")
	t_heat_dist.Draw("same")
	# t_heat_dist.Fit("gaus")
	raw_input()

	#Heatonly
	l_t = np.where(np.logical_and(t_heatonly["EC1"]>2 , t_heatonly["EC2"]>2))
	l_s = np.where(np.logical_and(s_heatonly["EC1"]>2 , s_heatonly["EC2"]>2))

	t_heatonly=t_heatonly[l_t]
	s_heatonly=s_heatonly[l_s]

	# true_dist = 0.414*t_heatonly["EIB"]+(1-0.414)*t_heatonly["EID"]
	# simu_dist = 0.414*s_heatonly["EIB"]+(1-0.414)*s_heatonly["EID"]
	true_dist = t_heatonly["EIB"]-t_heatonly["EID"]
	simu_dist = s_heatonly["EIB"]-s_heatonly["EID"]

	t_heat_dist = TH1F("t_heat_dist", "t_heat_dist", 100,-2,2)
	s_heat_dist = TH1F("s_heat_dist", "s_heat_dist", 100,-2,2)

	PyRPl.fill_TH1(t_heat_dist, true_dist)
	PyRPl.fill_TH1(s_heat_dist, simu_dist)

	PyRPl.process_TH1(t_heat_dist, X_title = "distance", color = kRed)
	PyRPl.process_TH1(s_heat_dist, X_title = "distance", color = kBlack)

	t_heat_dist.Scale(1./t_heat_dist.Integral())
	s_heat_dist.Scale(1./s_heat_dist.Integral())

	s_heat_dist.Draw("")
	s_heat_dist.Fit("gaus")
	t_heat_dist.Draw("same")
	t_heat_dist.Fit("gaus")
	raw_input()

def check_BDT_simulations(bolo_name, analysis_type, mass):

	""" Verify the 1D distributions of each event type
    
	Detail:

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		mass          = (list) WIMP mass

	Returns:
		void

	Raises:
		void
	"""

	plt.ion()

	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
	BDT_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_better/BDT_" + bolo_name + "/" + analysis_type + "/Application/"

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIC", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}
	arr_heatonly = np.loadtxt(pop_path + bolo_name + "_heatonly_full_info.txt", delimiter=",",  dtype=data_types)
	arr_FidGamma = np.loadtxt(pop_path + bolo_name + "_FidGamma_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S1Gamma = np.loadtxt(pop_path + bolo_name + "_S1Gamma_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Gamma = np.loadtxt(pop_path + bolo_name + "_S2Gamma_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S1Beta = np.loadtxt(pop_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(pop_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb   = np.loadtxt(pop_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb   = np.loadtxt(pop_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	arr_all = np.loadtxt(pop_path + bolo_name + "_true_events_all.txt", delimiter=",",  dtype=data_types)
	arr_FidGamma_after_cuts = np.loadtxt(pop_path + bolo_name + "_true_events_FidGamma.txt", delimiter=",",  dtype=data_types)
	arr_heatonly_after_cuts = np.loadtxt(pop_path + bolo_name + "_true_events_heatonly.txt", delimiter=",",  dtype=data_types)

	simu_heatonly = np.loadtxt(pop_path + bolo_name + "_BDT_heatonly.txt", delimiter=",",  dtype=data_types)
	simu_FidGamma = np.loadtxt(pop_path + bolo_name + "_BDT_FidGamma.txt", delimiter=",",  dtype=data_types)
	simu_S1Gamma = np.loadtxt(pop_path + bolo_name + "_BDT_S1Gamma.txt", delimiter=",",  dtype=data_types)
	simu_S2Gamma = np.loadtxt(pop_path + bolo_name + "_BDT_S2Gamma.txt", delimiter=",",  dtype=data_types)
	simu_S1Beta = np.loadtxt(pop_path + bolo_name + "_BDT_S1Beta.txt", delimiter=",",  dtype=data_types)
	simu_S2Beta = np.loadtxt(pop_path + bolo_name + "_BDT_S2Beta.txt", delimiter=",",  dtype=data_types)
	simu_S1Pb   = np.loadtxt(pop_path + bolo_name + "_BDT_S1Pb.txt", delimiter=",",  dtype=data_types)
	simu_S2Pb   = np.loadtxt(pop_path + bolo_name + "_BDT_S2Pb.txt", delimiter=",",  dtype=data_types)

	print simu_FidGamma.shape[0]
	list_FidGamma_index = np.random.randint(0,simu_FidGamma.shape[0], arr_FidGamma.shape[0])
	simu_FidGamma = simu_FidGamma[list_FidGamma_index]
	print simu_FidGamma.shape[0]

	simu_all = np.loadtxt(pop_path + bolo_name + "_simu_events_all.txt", delimiter=",",  dtype=data_types)

	list_channel = ["EC1", "EC2", "EIC", "EIB", "EIC", "EID"]
	list_evt = ["heatonly", "FidGamma", "S1Gamma", "S2Gamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb"]
	list_arr = [arr_heatonly_after_cuts, arr_FidGamma_after_cuts, arr_S1Gamma, arr_S2Gamma, arr_S1Beta, arr_S2Beta, arr_S1Pb, arr_S2Pb]
	list_simu = [simu_heatonly, simu_FidGamma, simu_S1Gamma, simu_S2Gamma, simu_S1Beta, simu_S2Beta, simu_S1Pb, simu_S2Pb]
	list_canvas = [TCanvas("cc"+evt, "cc"+ evt) for evt in list_evt]

	#Check the channels where no signal is expected

	d_hist_arr={}
	d_hist_simu={}
	for evt in list_evt:
		d_hist_arr[evt] = [TH1F(evt + "_" + chan, evt + "_" + chan, 100, -2,2) for chan in list_channel]
		d_hist_simu[evt] = [TH1F(evt + "_simu_" + chan, evt + "_simu_" + chan, 100, -2,2) for chan in list_channel]
	d_list_chan = {}

	d_list_chan["heatonly"] = ["EIC", "EIB", "EIC", "EID"]
	d_list_chan["FidGamma"] = ["EIC", "EIC"]
	for evt in list_evt:
		if "S1" in evt:
			d_list_chan[evt] = ["EIC", "EID"]
		if "S2" in evt:
			d_list_chan[evt] = ["EIC", "EIB"]

	for i in range(len(list_evt)):
		list_canvas[i].Divide(2,2)
		evt_type = list_evt[i]
		arr = list_arr[i]
		arr_simu = list_simu[i]
		for k in range(len(d_list_chan[evt_type])):
			chan = d_list_chan[evt_type][k]
			list_canvas[i].cd(k+1)
			for evt in arr[chan]:
				d_hist_arr[evt_type][k].Fill(evt)
			for evt in arr_simu[chan]:
				d_hist_simu[evt_type][k].Fill(evt)
			d_hist_arr[evt_type][k].Scale(1./d_hist_arr[evt_type][k].Integral())
			d_hist_simu[evt_type][k].Scale(1./d_hist_simu[evt_type][k].Integral())
			PyRPl.process_TH1(d_hist_arr[evt_type][k], X_title = chan, color = kRed)
			PyRPl.process_TH1(d_hist_simu[evt_type][k], X_title = chan, color = kBlack)
			d_hist_arr[evt_type][k].Draw()
			d_hist_simu[evt_type][k].Draw("same")			


	# #Check the channels where signal is expected
	# #FidGamma
	# list_canvas[1].Divide(2,2)
	# evt_type = list_evt[0]
	# arr = list_arr[1]
	# arr_simu = list_simu[1]
	# list_harr = [TH1F(evt_type + "_" + chan, evt_type + "_" + chan, 200, -2,15) for chan in ["EC1", "EC2", "EIB", "EID"]]
	# list_hsimu = [TH1F(evt_type + "_simu_" + chan, evt_type + "_simu_" + chan, 200, -2,15) for chan in ["EC1", "EC2", "EIB", "EID"]]
	# for index, k in enumerate([0,1,3,5]):
	# 	chan = list_channel[k]
	# 	list_canvas[1].cd(index+1)

	# 	for evt in arr[chan]:
	# 		list_harr[index].Fill(evt)
	# 	for evt in arr_simu[chan]:
	# 		list_hsimu[index].Fill(evt)
	# 	list_harr[index].Scale(1./list_harr[index].Integral())
	# 	list_hsimu[index].Scale(1./list_hsimu[index].Integral())
	# 	PyRPl.process_TH1(list_harr[index], X_title = chan, color = kRed)
	# 	PyRPl.process_TH1(list_hsimu[index], X_title = chan, color = kBlack)
	# 	list_harr[index].Draw()
	# 	list_hsimu[index].Draw("same")



	raw_input()

bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
# analysis_type = "ana_1.5_0_50"
analysis_type = "ana_1.5_min2_50"
mass = "5"
get_BDT_array(bolo_name, analysis_type, mass)
# check_BDT_simulations(bolo_name, analysis_type, mass)
# check_veto_distribution(bolo_name)
# check_collectrode_distribution(bolo_name)
# check_10keV_distribution(bolo_name)
# check_intrinsic_scattering(bolo_name)