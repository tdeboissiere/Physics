from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Poisson_file_handler as Poisson_fh
import matplotlib.pylab as plt

def check_Qmodel_surface_validity_Gamma(bolo_name):

	""" Check Q model validity (check there is no straggle)
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		void
	"""

	plt.ion()

	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	#Load the estimator
	d_est             = Poisson_fh.open_estimator_file(bolo_name)

    #Best estimator for heat: coefficients
	coeff_EC1, coeff_EC2 = float(d_est["HEAT"][:5]), 1 - float(d_est["HEAT"][:5])
	coeff_EIB, coeff_EID = float(d_est["FID"][:5]), 1-float(d_est["FID"][:5])

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}
	arr_FidGamma = np.loadtxt(pop_path + bolo_name + "_FidGamma_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S1Gamma = np.loadtxt(pop_path + bolo_name + "_S1Gamma_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Gamma = np.loadtxt(pop_path + bolo_name + "_S2Gamma_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_FidGamma = coeff_EIB*arr_FidGamma["EIB"] + coeff_EID*arr_FidGamma["EID"]
	arr_EC_FidGamma = coeff_EC1*arr_FidGamma["EC1"] + coeff_EC2*arr_FidGamma["EC2"]

	coeff_EIA, coeff_EIB = float(d_est["S1"][:5]), 1-float(d_est["S1"][:5])
	coeff_EIC, coeff_EID = float(d_est["S2"][:5]), 1-float(d_est["S2"][:5])

	arr_EI_S1Gamma, arr_EI_S2Gamma = coeff_EIA*arr_S1Gamma["EIA"] + coeff_EIB*arr_S1Gamma["EIB"], coeff_EIC*arr_S2Gamma["EIC"] + coeff_EID*arr_S2Gamma["EID"]
	arr_EC_S1Gamma, arr_EC_S2Gamma = coeff_EC1*arr_S1Gamma["EC1"] + coeff_EC2*arr_S1Gamma["EC2"], coeff_EC1*arr_S2Gamma["EC1"] + coeff_EC2*arr_S2Gamma["EC2"]

	lS1Gamma, lS2Gamma = np.where(arr_EC_S1Gamma<40), np.where(arr_EC_S2Gamma<40)
	lFidGamma = np.where(arr_EC_FidGamma<40)

	arr_EI_S1Gamma, arr_EC_S1Gamma = arr_EI_S1Gamma[lS1Gamma], arr_EC_S1Gamma[lS1Gamma]
	arr_EI_S2Gamma, arr_EC_S2Gamma = arr_EI_S2Gamma[lS2Gamma], arr_EC_S2Gamma[lS2Gamma]
	arr_EI_FidGamma, arr_EC_FidGamma = arr_EI_FidGamma[lFidGamma], arr_EC_FidGamma[lFidGamma]

	arr_EI_S1Gamma, arr_EC_S1Gamma = np.array(arr_EI_S1Gamma).astype(float), np.array(arr_EC_S1Gamma).astype(float)
	arr_EI_S2Gamma, arr_EC_S2Gamma = np.array(arr_EI_S2Gamma).astype(float), np.array(arr_EC_S2Gamma).astype(float)
	arr_EI_FidGamma, arr_EC_FidGamma = np.array(arr_EI_FidGamma).astype(float), np.array(arr_EC_FidGamma).astype(float)

	arr_Q_S1Gamma = arr_EI_S1Gamma/((1+8./3)*arr_EC_S1Gamma - arr_EI_S1Gamma*5.5/3)
	arr_Q_S2Gamma = arr_EI_S1Gamma/((1+8./3)*arr_EC_S1Gamma - arr_EI_S1Gamma*5.5/3)
	arr_Q_FidGamma = arr_EI_FidGamma/((1+8./3)*arr_EC_FidGamma - arr_EI_FidGamma*8./3)

	gr_S1Gamma, gr_S2Gamma   = TGraph(len(arr_EI_S1Gamma), arr_EC_S1Gamma, arr_EI_S1Gamma),  TGraph(len(arr_EI_S2Gamma), arr_EC_S2Gamma, arr_EI_S2Gamma)
	gr_QS1Gamma, gr_QS2Gamma = TGraph(len(arr_Q_S1Gamma), arr_EC_S1Gamma, arr_Q_S1Gamma),  TGraph(len(arr_Q_S2Gamma), arr_EC_S2Gamma, arr_Q_S2Gamma)
	gr_FidGamma= TGraph(len(arr_EI_FidGamma), arr_EC_FidGamma, arr_EI_FidGamma)


	gr_QFidGamma = TGraph(len(arr_Q_FidGamma), arr_EC_FidGamma, arr_Q_FidGamma)
	PyRPl.process_TGraph(gr_FidGamma, X_title = "Heat", Y_title = "Ion", color=kOrange-3)
	PyRPl.process_TGraph(gr_QFidGamma, X_title = "Heat", Y_title = "Q", color=kOrange-3)

	PyRPl.process_TGraph(gr_S1Gamma, X_title = "Heat", Y_title = "Ion", color=kGreen+2), PyRPl.process_TGraph(gr_S2Gamma, X_title = "Heat", Y_title = "Ion", color=kBlue)
	PyRPl.process_TGraph(gr_QS1Gamma, X_title = "Heat", Y_title = "Q", color=kGreen+2), PyRPl.process_TGraph(gr_QS2Gamma, X_title = "Heat", Y_title = "Q", color=kBlue)


	list_graph = [gr_FidGamma, gr_S1Gamma, gr_S2Gamma]
	list_Qgraph = [gr_QFidGamma, gr_QS1Gamma, gr_QS2Gamma]
	
	Ana_path   = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"

	tt, ff = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_FidGamma_tree_test.root", "t_new")
	ttS1, ffS1 = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_S1Gamma_tree_test.root", "t_new")
	ttS2, ffS2 = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_S2Gamma_tree_test.root", "t_new")

	cc = TCanvas("cc", "cc")
	gr_FidGamma.Draw("A*")
	tt.Draw("EIB:EC1", "", "same*")

	ccS1 = TCanvas("ccS1", "ccS1")
	gr_S1Gamma.Draw("A*")
	ttS1.Draw("EIB:EC1", "", "same*")

	ccS2 = TCanvas("ccS2", "ccS2")
	gr_S2Gamma.Draw("A*")
	ttS2.Draw("EID:EC1", "", "same*")

	raw_input()


bolo_name = "FID837"
check_Qmodel_surface_validity_Gamma(bolo_name)