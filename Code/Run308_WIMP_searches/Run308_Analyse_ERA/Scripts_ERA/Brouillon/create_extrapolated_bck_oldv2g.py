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

import BDT_file_handler as BDT_fh

def fit_surface_relation(bolo_name):

	"""Plot EI vs EH and fit
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		Assertion Error if the initial histograms do not exist
	"""

	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

    #Load the estimator
	d_est = Poisson_fh.open_estimator_file(bolo_name)
    
    #Best estimator for heat: coefficients
	coeff_EC1 = float(d_est["HEAT"][:5])
	coeff_EC2 = 1 - coeff_EC1

    #Best estimator for surf1 ion: coefficients
	coeff_EIA = float(d_est["S1"][:5])
	coeff_EIB = 1 - coeff_EIA

    #Best estimator for surf2 ion: coefficients
	coeff_EIC = float(d_est["S2"][:5])
	coeff_EID = 1 - coeff_EIC

	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Beta = np.loadtxt(gen_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(gen_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Beta = coeff_EIA*arr_S1Beta["EIA"] + coeff_EIB*arr_S1Beta["EIB"]
	arr_EI_S2Beta = coeff_EIC*arr_S2Beta["EIC"] + coeff_EID*arr_S2Beta["EID"]

	arr_EC_S1Beta = coeff_EC1*arr_S1Beta["EC1"] + coeff_EC2*arr_S1Beta["EC2"]
	arr_EC_S2Beta = coeff_EC1*arr_S2Beta["EC1"] + coeff_EC2*arr_S2Beta["EC2"]

	arr_EI_S1Beta = np.array(arr_EI_S1Beta).astype(float)
	arr_EC_S1Beta = np.array(arr_EC_S1Beta).astype(float)

	arr_EI_S2Beta = np.array(arr_EI_S2Beta).astype(float)
	arr_EC_S2Beta = np.array(arr_EC_S2Beta).astype(float)

	gr_S1Beta = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta)
	gr_S2Beta = TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)

	PyRPl.process_TGraph(gr_S1Beta, X_title = "Heat", Y_title = "Ion", color=kRed)
	PyRPl.process_TGraph(gr_S2Beta, X_title = "Heat", Y_title = "Ion", color=kBlue)

	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Pb = np.loadtxt(gen_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt(gen_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Pb = coeff_EIA*arr_S1Pb["EIA"] + coeff_EIB*arr_S1Pb["EIB"]
	arr_EI_S2Pb = coeff_EIC*arr_S2Pb["EIC"] + coeff_EID*arr_S2Pb["EID"]

	arr_EC_S1Pb = coeff_EC1*arr_S1Pb["EC1"] + coeff_EC2*arr_S1Pb["EC2"]
	arr_EC_S2Pb = coeff_EC1*arr_S2Pb["EC1"] + coeff_EC2*arr_S2Pb["EC2"]

	arr_EI_S1Pb = np.array(arr_EI_S1Pb).astype(float)
	arr_EC_S1Pb = np.array(arr_EC_S1Pb).astype(float)

	arr_EI_S2Pb = np.array(arr_EI_S2Pb).astype(float)
	arr_EC_S2Pb = np.array(arr_EC_S2Pb).astype(float)

	gr_S1Pb = TGraph(len(arr_EI_S1Pb), arr_EC_S1Pb, arr_EI_S1Pb)
	gr_S2Pb = TGraph(len(arr_EI_S2Pb), arr_EC_S2Pb, arr_EI_S2Pb)

	PyRPl.process_TGraph(gr_S1Pb, X_title = "Heat", Y_title = "Ion", color=kRed)
	PyRPl.process_TGraph(gr_S2Pb, X_title = "Heat", Y_title = "Ion", color=kBlue)

	list_h= [TH2F("h" + str(i), "h" + str(i), 100, 0,40, 100, 0, 40) for i in range(4)]


	out_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/ROOT_files/"

	fout = TFile(out_path + bolo_name + "_ion_heat_surface_relation.root" , "recreate")

	cc = TCanvas("cc", "cc")
	cc.Divide(2,2)
	cc.cd(1)
	list_h[0].Draw()
	gr_S1Beta.Draw("same*")
	gr_S1Beta.Fit("pol3")
	cc.cd(2)
	list_h[1].Draw()
	gr_S2Beta.Draw("same*")
	gr_S2Beta.Fit("pol3")
	cc.cd(3)
	list_h[2].Draw()
	gr_S1Pb.Draw("same*")
	gr_S1Pb.Fit("pol4")
	cc.cd(4)
	list_h[3].Draw()
	gr_S2Pb.Draw("same*")
	gr_S2Pb.Fit("pol4")
	raw_input()

	fS1Beta = gr_S1Beta.GetFunction("pol3")
	fS1Beta.SetName("S1Beta")
	fS1Beta.Write()

	fS2Beta = gr_S2Beta.GetFunction("pol3")
	fS2Beta.SetName("S2Beta")
	fS2Beta.Write()

	fS1Pb = gr_S1Pb.GetFunction("pol4")
	fS1Pb.SetName("S1Pb")
	fS1Pb.Write()

	fS2Pb = gr_S2Pb.GetFunction("pol4")
	fS2Pb.SetName("S2Pb")
	fS2Pb.Write()

	fout.Close()

def fit_surface_KDE(bolo_name):

	"""KDE estimation of EI vs EH surface
    
	Detail:

	Args:
		bolo_name = (str) bolometer name

	Returns:
		void

	Raises:
		Assertion Error if the initial histograms do not exist
	"""

	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

    #Load the estimator
	d_est = Poisson_fh.open_estimator_file(bolo_name)
    
    #Best estimator for heat: coefficients
	coeff_EC1 = float(d_est["HEAT"][:5])
	coeff_EC2 = 1 - coeff_EC1

    #Best estimator for surf1 ion: coefficients
	coeff_EIA = float(d_est["S1"][:5])
	coeff_EIB = 1 - coeff_EIA

    #Best estimator for surf2 ion: coefficients
	coeff_EIC = float(d_est["S2"][:5])
	coeff_EID = 1 - coeff_EIC

	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}

	######
	# Beta
	######
	arr_S1Beta = np.loadtxt(gen_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(gen_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Beta = coeff_EIA*arr_S1Beta["EIA"] + coeff_EIB*arr_S1Beta["EIB"]
	arr_EI_S2Beta = coeff_EIC*arr_S2Beta["EIC"] + coeff_EID*arr_S2Beta["EID"]

	arr_EC_S1Beta = coeff_EC1*arr_S1Beta["EC1"] + coeff_EC2*arr_S1Beta["EC2"]
	arr_EC_S2Beta = coeff_EC1*arr_S2Beta["EC1"] + coeff_EC2*arr_S2Beta["EC2"]

	arr_EI_S1Beta = np.array(arr_EI_S1Beta).astype(float)
	arr_EC_S1Beta = np.array(arr_EC_S1Beta).astype(float)

	arr_EI_S2Beta = np.array(arr_EI_S2Beta).astype(float)
	arr_EC_S2Beta = np.array(arr_EC_S2Beta).astype(float)

	######
	# Pb
	######
	arr_S1Pb = np.loadtxt(gen_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt(gen_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Pb = coeff_EIA*arr_S1Pb["EIA"] + coeff_EIB*arr_S1Pb["EIB"]
	arr_EI_S2Pb = coeff_EIC*arr_S2Pb["EIC"] + coeff_EID*arr_S2Pb["EID"]

	arr_EC_S1Pb = coeff_EC1*arr_S1Pb["EC1"] + coeff_EC2*arr_S1Pb["EC2"]
	arr_EC_S2Pb = coeff_EC1*arr_S2Pb["EC1"] + coeff_EC2*arr_S2Pb["EC2"]

	arr_EI_S1Pb = np.array(arr_EI_S1Pb).astype(float)
	arr_EC_S1Pb = np.array(arr_EC_S1Pb).astype(float)

	arr_EI_S2Pb = np.array(arr_EI_S2Pb).astype(float)
	arr_EC_S2Pb = np.array(arr_EC_S2Pb).astype(float)


	#boundaries
	nx, ny = 100, 100
	x, y = np.meshgrid(np.linspace(0, 60, nx),np.linspace(0, 60, ny))
	xy = np.vstack([x.ravel(), y.ravel()]).T
	values = np.vstack([arr_EC_S1Pb, arr_EI_S1Pb])

	from sklearn.grid_search import GridSearchCV
	grid = GridSearchCV(KernelDensity(),{'bandwidth': np.linspace(0.1, 5.0, 30)},cv=20) # 20-fold cross-validation
	grid.fit(zip(*values))
	bw= grid.best_params_["bandwidth"]

	# -------------------------------------------------------
	# Perform a kernel density estimate on the data using sklearn.
	kde = KernelDensity(kernel='gaussian', bandwidth=2*bw).fit(zip(*values))
	log_dens = np.exp(kde.score_samples(xy)) 
	
	log_dens = log_dens.reshape(x.shape)

	levels = np.linspace(0, log_dens.max(), 25)
	plt.contourf(x, y, log_dens, levels=levels, cmap=plt.cm.Reds)
	
	plt.show()
	raw_input()


	# # Get KDE value for the point.
	# iso2 = kernel_sk.score_samples([[x1, y1]])
	# print 'iso2 = ', np.exp(iso2[0])

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

	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name)


	#Load the corresponding files
	gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
	data_dir   = gen_path + "Analyse_" + bolo_name +"/ROOT_files/ROOT_scaling/" + bolo_name 
	output_dir = gen_path + "Analyse_" + bolo_name +"/ROOT_files/" + bolo_name 

	print data_dir

	#Check the files exist
	assert(os.path.isfile(data_dir + "_spectrum_FidGamma.root") )
	assert(os.path.isfile(data_dir + "_spectrum_S1Gamma.root") )
	assert(os.path.isfile(data_dir + "_spectrum_S2Gamma.root") )

	hFidGamma, file_FidGamma   = PyRPl.open_ROOT_object(data_dir + "_spectrum_FidGamma.root", "h1")
	hS1Gamma, file_S1Gamma  = PyRPl.open_ROOT_object(data_dir + "_spectrum_S1Gamma.root", "h1")
	hS2Gamma, file_S2Gamma  = PyRPl.open_ROOT_object(data_dir + "_spectrum_S2Gamma.root", "h1")

    ##################
    #
	#   Beta from bolo
	#
	##################

    #Load the estimator
	d_est = Poisson_fh.open_estimator_file(bolo_name)
    
    #Best estimator for heat: coefficients
	coeff_EC1 = float(d_est["HEAT"][:5])
	coeff_EC2 = 1 - coeff_EC1

	surf_bol_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Beta = np.loadtxt(surf_bol_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(surf_bol_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)

	#Only use EC1 in this analysis (it has best resolution)
	arr_EC_S1Beta =arr_S1Beta["EC1"]
	arr_EC_S2Beta =arr_S2Beta["EC1"]

	arr_EC_S1Beta = np.array(arr_EC_S1Beta).astype(float)
	arr_EC_S2Beta = np.array(arr_EC_S2Beta).astype(float)

	#Data driven model: build histo of the bolo data
	hS1Beta_bol = TH1F("hS1Beta_bol", "hS1Beta_bol", 120, 0, 60)
	hS2Beta_bol = TH1F("hS2Beta_bol", "hS2Beta_bol", 120, 0, 60)

	#Fill histos
	for heat in arr_EC_S1Beta:
		hS1Beta_bol.Fill(heat)
	for heat in arr_EC_S2Beta:
		hS2Beta_bol.Fill(heat)

	hS1Beta_bol.Smooth(51)
	sS1Beta_bol = TSpline5(hS1Beta_bol)
	sS1Beta_bol.SetLineColor(kRed)
	hS1Beta_bol.Draw()
	sS1Beta_bol.Draw("same")
	S1Beta_endpt = 5
	slope_S1Beta=sS1Beta_bol.Derivative(S1Beta_endpt)
	offset_S1Beta=sS1Beta_bol.Eval(S1Beta_endpt)-S1Beta_endpt*slope_S1Beta

	class S1Beta_bol:
		def __call__( self, x, par ):
			if x[0]>S1Beta_endpt:
				return sS1Beta_bol.Eval(x[0]) + par[0]
			else:
				return slope_S1Beta*x[0]+offset_S1Beta

	fS1Beta_bol_extra = TF1( "S1Beta_extra", S1Beta_bol(), 0, 60., 1 )
	fS1Beta_bol_extra.SetParameter(0,0)
	file_S1Beta_bol_extra =TFile(output_dir + "_Beta_from_bolo_spectrum_extrapol.root", "recreate")
	fS1Beta_bol_extra.Write()
	file_S1Beta_bol_extra.Close()

	raw_input()

	hS2Beta_bol.Smooth(41)
	sS2Beta_bol = TSpline5(hS2Beta_bol)
	sS2Beta_bol.SetLineColor(kRed)
	hS2Beta_bol.Draw()
	sS2Beta_bol.Draw("same")
	S2Beta_endpt = 5
	slope_S2Beta=sS2Beta_bol.Derivative(S2Beta_endpt)
	offset_S2Beta=sS2Beta_bol.Eval(S2Beta_endpt)-S2Beta_endpt*slope_S2Beta

	class S2Beta_bol:
		def __call__( self, x, par ):
			if x[0]>S2Beta_endpt:
				return sS2Beta_bol.Eval(x[0]) + par[0]
			else:
				return slope_S2Beta*x[0]+offset_S2Beta

	fS2Beta_bol_extra = TF1( "S2Beta_extra", S2Beta_bol(), 0, 60., 1 )
	fS2Beta_bol_extra.SetParameter(0,0)
	file_S2Beta_bol_extra =TFile(output_dir + "_Beta_from_bolo_spectrum_extrapol.root", "update")
	fS2Beta_bol_extra.Write()
	file_S2Beta_bol_extra.Close()

	raw_input()

    ##################
    #
	#   Pb from bolo
	#
	##################

    #Load the estimator
	d_est = Poisson_fh.open_estimator_file(bolo_name)
    
    #Best estimator for heat: coefficients
	coeff_EC1 = float(d_est["HEAT"][:5])
	coeff_EC2 = 1 - coeff_EC1

	surf_bol_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Pb = np.loadtxt(surf_bol_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt(surf_bol_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	#Only use EC1 in this analysis (it has best resolution)
	arr_EC_S1Pb =arr_S1Pb["EC1"]
	arr_EC_S2Pb =arr_S2Pb["EC1"]

	arr_EC_S1Pb = np.array(arr_EC_S1Pb).astype(float)
	arr_EC_S2Pb = np.array(arr_EC_S2Pb).astype(float)

	#Data driven model: build histo of the bolo data
	hS1Pb_bol = TH1F("hS1Pb_bol", "hS1Pb_bol", 30, 5, 40)
	hS2Pb_bol = TH1F("hS2Pb_bol", "hS2Pb_bol", 50, 5, 40)

	#Fill histos
	for heat in arr_EC_S1Pb:
		hS1Pb_bol.Fill(heat)
	for heat in arr_EC_S2Pb:
		hS2Pb_bol.Fill(heat)

	hS1Pb_bol.Draw()
	hS1Pb_bol.Fit("pol3")
	S1Pbfit =  hS1Pb_bol.GetFunction("pol3")
	S1Pb_endpt = 5
	S1Pb_zero_pt =  S1Pbfit.GetX(0,30,50)
	slope_S1Pb=S1Pbfit.Derivative(S1Pb_endpt)
	offset_S1Pb=S1Pbfit.Eval(S1Pb_endpt)-S1Pb_endpt*slope_S1Pb
	raw_input()

	class S1Pb_bol:
		def __call__( self, x, par ):
			if x[0]<=S1Pb_endpt:
				return slope_S1Pb*x[0]+offset_S1Pb + par[0]
			elif S1Pb_zero_pt>x[0]>S1Pb_endpt:
				return S1Pbfit.Eval(x[0])
			elif x[0]>S1Pb_zero_pt:
				return 0

	fS1Pb_bol_extra = TF1( "S1Pb_extra", S1Pb_bol(), 0, 60., 1 )
	fS1Pb_bol_extra.SetParameter(0,0)
	file_S1Pb_bol_extra =TFile(output_dir + "_Pb_from_bolo_spectrum_extrapol.root", "recreate")
	fS1Pb_bol_extra.Write()
	file_S1Pb_bol_extra.Close()


	hS2Pb_bol.Draw()
	hS2Pb_bol.Fit("pol3")
	S2Pbfit =  hS2Pb_bol.GetFunction("pol3")
	S2Pb_endpt = 5
	S2Pb_zero_pt =  S2Pbfit.GetX(0,30,50)
	slope_S2Pb=S2Pbfit.Derivative(S2Pb_endpt)
	offset_S2Pb=S2Pbfit.Eval(S2Pb_endpt)-S2Pb_endpt*slope_S2Pb
	raw_input()

	class S2Pb_bol:
		def __call__( self, x, par ):
			if x[0]<=S2Pb_endpt:
				return slope_S2Pb*x[0]+offset_S2Pb + par[0]
			elif S2Pb_zero_pt>x[0]>S2Pb_endpt:
				return S2Pbfit.Eval(x[0])
			elif x[0]>S2Pb_zero_pt:
				return 0

	fS2Pb_bol_extra = TF1( "S2Pb_extra", S2Pb_bol(), 0, 60., 1 )
	fS2Pb_bol_extra.SetParameter(0,0)
	file_S2Pb_bol_extra =TFile(output_dir + "_Pb_from_bolo_spectrum_extrapol.root", "update")
	fS2Pb_bol_extra.Write()
	file_S2Pb_bol_extra.Close()

	# raw_input()

    ##################
    #
	#   Beta and Pb
	#
	##################

	data_types = {"names": ("RUN", "SN", "EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("i", "i", "f", "f", "f", "f", "f", "f")}
	arr_S1Beta = np.loadtxt(gen_path + "/Text_files/S1Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(gen_path + "/Text_files/S2Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb = np.loadtxt(gen_path + "/Text_files/S1Pb_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt(gen_path + "/Text_files/S2Pb_heatremoved.txt", delimiter=",",  dtype=data_types)

	hS1Beta = TH1F("hS1Beta", "hS1Beta", 300, 0, 60)
	hS2Beta = TH1F("hS2Beta", "hS2Beta", 300, 0, 60)

	hS1Pb = TH1F("hS1Pb", "hS1Pb", 300, 0, 60)
	hS2Pb = TH1F("hS2Pb", "hS2Pb", 300, 0, 60)

	#Fill histos
	for i in range(arr_S1Beta.shape[0]):
		hS1Beta.Fill(0.5*(arr_S1Beta["EC1"][i]+arr_S1Beta["EC2"][i]))
	for i in range(arr_S2Beta.shape[0]):
		hS2Beta.Fill(0.5*(arr_S2Beta["EC1"][i]+arr_S2Beta["EC2"][i]))
	for i in range(arr_S1Pb.shape[0]):
		hS1Pb.Fill(0.5*(arr_S1Pb["EC1"][i]+arr_S1Pb["EC2"][i]))
	for i in range(arr_S2Pb.shape[0]):
		hS2Pb.Fill(0.5*(arr_S2Pb["EC1"][i]+arr_S2Pb["EC2"][i]))


	########
	#S1Beta
	########
	hS1Beta.Smooth(2000)
	sS1Beta = TSpline5(hS1Beta)
	sS1Beta.SetLineColor(kRed)
	hS1Beta.Draw()
	sS1Beta.Draw("same")

	S1Beta_endpt = 5
	slope_S1Beta=sS1Beta.Derivative(S1Beta_endpt)
	offset_S1Beta=sS1Beta.Eval(S1Beta_endpt)-S1Beta_endpt*slope_S1Beta

	class S1Beta:
		def __call__( self, x, par ):
			if x[0]>S1Beta_endpt:
				return sS1Beta.Eval(x[0]) + par[0]
			else:
				return slope_S1Beta*x[0]+offset_S1Beta

	fS1Beta_extra = TF1( "S1Beta_extra", S1Beta(), 0, 60., 1 )
	fS1Beta_extra.SetParameter(0,0)
	fS1Beta_extra.Draw()
	raw_input()

	file_S1Beta_extra =TFile(output_dir + "_Surf_spectrum_extrapol.root", "recreate")
	fS1Beta_extra.Write()
	file_S1Beta_extra.Close()

	########
	#S2Beta
	########
	hS2Beta.Smooth(2000)
	sS2Beta = TSpline5(hS2Beta)
	sS2Beta.SetLineColor(kRed)
	hS2Beta.Draw()
	sS2Beta.Draw("same")

	S2Beta_endpt = 5
	slope_S2Beta=sS2Beta.Derivative(S2Beta_endpt)
	offset_S2Beta=sS2Beta.Eval(S2Beta_endpt)-S2Beta_endpt*slope_S2Beta

	class S2Beta:
		def __call__( self, x, par ):
			if x[0]>S2Beta_endpt:
				return sS2Beta.Eval(x[0]) + par[0]
			else:
				return slope_S2Beta*x[0]+offset_S2Beta

	fS2Beta_extra = TF1( "S2Beta_extra", S2Beta(), 0, 60., 1 )
	fS2Beta_extra.SetParameter(0,0)
	fS2Beta_extra.Draw()
	raw_input()

	file_S2Beta_extra =TFile(output_dir + "_Surf_spectrum_extrapol.root", "update")
	fS2Beta_extra.Write()
	file_S2Beta_extra.Close()

	########
	#S1Pb
	########
	hS1Pb.Smooth(2000)
	sS1Pb = TSpline5(hS1Pb)
	sS1Pb.SetLineColor(kRed)
	hS1Pb.Draw()
	sS1Pb.Draw("same")

	S1Pb_endpt = 15
	slope_S1Pb=sS1Pb.Derivative(S1Pb_endpt)
	offset_S1Pb=sS1Pb.Eval(S1Pb_endpt)-S1Pb_endpt*slope_S1Pb

	class S1Pb:
		def __call__( self, x, par ):
			if x[0]>S1Pb_endpt:
				return sS1Pb.Eval(x[0]) + par[0]
			else:
				return slope_S1Pb*x[0]+offset_S1Pb

	fS1Pb_extra = TF1( "S1Pb_extra", S1Pb(), 0, 60., 1 )
	fS1Pb_extra.SetParameter(0,0)
	fS1Pb_extra.Draw()
	raw_input()

	file_S1Pb_extra =TFile(output_dir + "_Surf_spectrum_extrapol.root", "update")
	fS1Pb_extra.Write()
	file_S1Pb_extra.Close()

	########
	#S2Pb
	########
	hS2Pb.Smooth(2000)
	sS2Pb = TSpline5(hS2Pb)
	sS2Pb.SetLineColor(kRed)
	hS2Pb.Draw()
	sS2Pb.Draw("same")

	S2Pb_endpt = 13
	slope_S2Pb=sS2Pb.Derivative(S2Pb_endpt)
	offset_S2Pb=sS2Pb.Eval(S2Pb_endpt)-S2Pb_endpt*slope_S2Pb

	class S2Pb:
		def __call__( self, x, par ):
			if x[0]>S2Pb_endpt:
				return sS2Pb.Eval(x[0]) + par[0]
			else:
				return slope_S2Pb*x[0]+offset_S2Pb

	fS2Pb_extra = TF1( "S2Pb_extra", S2Pb(), 0, 60., 1 )
	fS2Pb_extra.SetParameter(0,0)
	fS2Pb_extra.Draw()
	raw_input()

	file_S2Pb_extra =TFile(output_dir + "_Surf_spectrum_extrapol.root", "update")
	fS2Pb_extra.Write()
	file_S2Pb_extra.Close()



    ##################
    #
	#      Gamma Fid
	#
	##################

	class Gamma:
		def __call__( self, x, par ):
			flat_part = par[0] 
			peak_10_4 = par[1]*TMath.Gaus(x[0], 10.37, d_std["OWC1"])
			peak_9_66 = 0.097*par[1]*TMath.Gaus(x[0], 9.66, d_std["OWC1"])
			peak_8_98 = par[2]*TMath.Gaus(x[0], 8.98, d_std["OWC1"])
			peak_7_11 = par[3]*TMath.Gaus(x[0], 7.11, d_std["OWC1"])
			peak_6_54 = par[4]*TMath.Gaus(x[0], 6.54, d_std["OWC1"])
			peak_5_99 = par[5]*TMath.Gaus(x[0], 5.99, d_std["OWC1"])
			peak_5_46 = par[6]*TMath.Gaus(x[0], 5.46, d_std["OWC1"])
			peak_4_97 = par[7]*TMath.Gaus(x[0], 4.97, d_std["OWC1"])
			peak_1_3  = 0.1*par[1]*TMath.Gaus(x[0], 1.3, d_std["OWC1"])
			list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
			return flat_part + sum(list_peaks)

	# Call FidGamma function to get derivative
	fFidGamma = TF1( "FidGamma", Gamma(), 0, 40., 8 )
	fFidGamma.SetNpx(500)
	fFidGamma.SetParName(0,"p0")
	fFidGamma.SetParName(1,"A_10.4")
	fFidGamma.SetParName(2,"A_8.98")
	fFidGamma.SetParName(3,"A_7.11")
	fFidGamma.SetParName(4,"A_6.54")
	fFidGamma.SetParName(5,"A_5.99")
	fFidGamma.SetParName(6,"A_5.46")
	fFidGamma.SetParName(7,"A_4.97")

	fFidGamma.SetParameters(1,1,1,1,1,1,1,1)

	for i in range(8):
		fFidGamma.SetParLimits(i, 0, 100)
     
	FidGamma_endpt =0



	hFidGamma.Draw()
	cst = hFidGamma.GetBinWidth(20)*hFidGamma.Integral(hFidGamma.FindBin(11.5),hFidGamma.FindBin(14.5))/(14.5-11.5)
	fFidGamma.FixParameter(0,cst)
	hFidGamma.Fit("FidGamma","LL","",FidGamma_endpt,15)

	raw_input()
	c1.Print("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Figures/" +bolo_name + "_FidGamma_fit.png")

	slope_FidGamma=fFidGamma.Derivative(FidGamma_endpt)
	offset_FidGamma=fFidGamma.Eval(FidGamma_endpt)-FidGamma_endpt*slope_FidGamma

	class FidGamma_bckg_extra:
		def __call__( self, x, par ):
			if x[0]>0.5:
				return fFidGamma.Eval(x[0]) +par[0]
			else: 
				return slope_FidGamma*x[0]+offset_FidGamma

	fFidGamma_extra = TF1( "FidGamma_extra", FidGamma_bckg_extra(), 0, 15., 1 )
	fFidGamma_extra.SetNpx(500)
	fFidGamma_extra.SetParameter(0,0)
	

	file_FidGamma_extra =TFile(output_dir + "_Gamma_spectrum_extrapol.root", "recreate")
	fFidGamma_extra.Write()
	file_FidGamma_extra.Close()


	#Deconv model: multiply ampl by 100 and divide sigma by 100 to conserved correct amplitudes
	class Gamma_deconv:
		def __call__( self, x, par ):
			flat_part = par[0] 
			peak_10_4 = par[1]*TMath.Gaus(x[0], 10.37, 0.05*d_std["OWC1"])
			peak_9_66 = 0.097*par[1]*TMath.Gaus(x[0], 9.66, 0.05*d_std["OWC1"])
			peak_8_98 = par[2]*TMath.Gaus(x[0], 8.98, 0.05*d_std["OWC1"])
			peak_7_11 = par[3]*TMath.Gaus(x[0], 7.11, 0.05*d_std["OWC1"])
			peak_6_54 = par[4]*TMath.Gaus(x[0], 6.54, 0.05*d_std["OWC1"])
			peak_5_99 = par[5]*TMath.Gaus(x[0], 5.99, 0.05*d_std["OWC1"])
			peak_5_46 = par[6]*TMath.Gaus(x[0], 5.46, 0.05*d_std["OWC1"])
			peak_4_97 = par[7]*TMath.Gaus(x[0], 4.97, 0.05*d_std["OWC1"])
			peak_1_3  = 0.1*par[1]*TMath.Gaus(x[0], 1.3, 0.05*d_std["OWC1"])
			list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
			return flat_part + sum(list_peaks)

	fFidGamma_deconv = TF1("FidGamma_deconv", Gamma_deconv(), 0, 15, 8)

	fFidGamma_deconv.SetNpx(5000)
	fFidGamma_deconv.SetParameter(0,cst)
	fFidGamma_deconv.SetParameter(1,20*fFidGamma.GetParameter(1))
	fFidGamma_deconv.SetParameter(2,20*fFidGamma.GetParameter(2))
	fFidGamma_deconv.SetParameter(3,20*fFidGamma.GetParameter(3))
	fFidGamma_deconv.SetParameter(4,20*fFidGamma.GetParameter(4))
	fFidGamma_deconv.SetParameter(5,20*fFidGamma.GetParameter(5))
	fFidGamma_deconv.SetParameter(6,20*fFidGamma.GetParameter(6))
	fFidGamma_deconv.SetParameter(7,20*fFidGamma.GetParameter(7))

	file_FidGamma_deconv =TFile(output_dir + "_Gamma_spectrum_extrapol_deconv.root", "recreate")
	fFidGamma_deconv.Write()
	file_FidGamma_deconv.Close()

    ##################
    #
	#      Gamma S1
	#
	##################


	class SurfGamma:
		def __call__( self, x, par ):
			luke = (1+5.5/3)/(1+8./3)*par[8]
			flat_part = par[0] 
			peak_10_4 = par[1]*TMath.Gaus(x[0], 10.37*luke, d_std["OWC1"])
			peak_9_66 = 0.097*par[1]*TMath.Gaus(x[0], 9.66*luke, d_std["OWC1"])
			peak_8_98 = par[2]*TMath.Gaus(x[0], 8.98*luke, d_std["OWC1"])
			peak_7_11 = par[3]*TMath.Gaus(x[0], 7.11*luke, d_std["OWC1"])
			peak_6_54 = par[4]*TMath.Gaus(x[0], 6.54*luke, d_std["OWC1"])
			peak_5_99 = par[5]*TMath.Gaus(x[0], 5.99*luke, d_std["OWC1"])
			peak_5_46 = par[6]*TMath.Gaus(x[0], 5.46*luke, d_std["OWC1"])
			peak_4_97 = par[7]*TMath.Gaus(x[0], 4.97*luke, d_std["OWC1"])
			peak_1_3  = 0.1*par[1]*TMath.Gaus(x[0], 1.3*luke, d_std["OWC1"])
			list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
			return flat_part + sum(list_peaks)

	fS1Gamma = TF1( "S1Gamma", SurfGamma(), 0, 40., 9 )
	fS1Gamma.SetNpx(500)
	fS1Gamma.SetParName(0,"p0")
	fS1Gamma.SetParName(1,"A_10.4")
	fS1Gamma.SetParName(2,"A_8.98")
	fS1Gamma.SetParName(3,"A_7.11")
	fS1Gamma.SetParName(4,"A_6.54")
	fS1Gamma.SetParName(5,"A_5.99")
	fS1Gamma.SetParName(6,"A_5.46")
	fS1Gamma.SetParName(7,"A_4.97")
	fS1Gamma.SetParName(8,"Luke_shift")

	fS1Gamma.SetParameters(1,1,1,1,1,1,1,1,1)

	for i in range(8):
		fS1Gamma.SetParLimits(i, 0, 100)

	S1Gamma_endpt =0.1

	hS1Gamma.Draw()
	cst = hS1Gamma.GetBinWidth(20)*hS1Gamma.Integral(hS1Gamma.FindBin(11.5),hS1Gamma.FindBin(14.5))/(14.5-11.5)
	fS1Gamma.FixParameter(0,cst)
	hS1Gamma.Draw()
	hS1Gamma.Fit("S1Gamma","LL","",S1Gamma_endpt,15)
	raw_input()
	c1.Print("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Figures/" +bolo_name + "_S1Gamma_fit.png")

	
	slope_S1Gamma=fS1Gamma.Derivative(S1Gamma_endpt)
	offset_S1Gamma=fS1Gamma.Eval(S1Gamma_endpt)-S1Gamma_endpt*slope_S1Gamma

	class S1Gamma_bckg_extra:
		def __call__( self, x, par ):
			if x[0]>0.5:
				return fS1Gamma.Eval(x[0]) +par[0]
			else: 
				return slope_S1Gamma*x[0]+offset_S1Gamma

	fS1Gamma_extra = TF1( "S1Gamma_extra", S1Gamma_bckg_extra(), 0, 15., 1 )
	fS1Gamma_extra.SetNpx(500)
	fS1Gamma_extra.SetParameter(0,0)
	
	file_S1Gamma_extra =TFile(output_dir + "_Gamma_spectrum_extrapol.root", "update")
	fS1Gamma_extra.Write()
	file_S1Gamma_extra.Close()

	#Deconv model: multiply ampl by 100 and divide sigma by 100 to conserved correct amplitudes
	class SurfGamma_deconv:
		def __call__( self, x, par ):
			luke = (1+5.5/3)/(1+8./3)*par[8]
			flat_part = par[0] 
			peak_10_4 = par[1]*TMath.Gaus(x[0], 10.37*luke, 0.05*d_std["OWC1"])
			peak_9_66 = 0.097*par[1]*TMath.Gaus(x[0], 9.66*luke, 0.05*d_std["OWC1"])
			peak_8_98 = par[2]*TMath.Gaus(x[0], 8.98*luke, 0.05*d_std["OWC1"])
			peak_7_11 = par[3]*TMath.Gaus(x[0], 7.11*luke, 0.05*d_std["OWC1"])
			peak_6_54 = par[4]*TMath.Gaus(x[0], 6.54*luke, 0.05*d_std["OWC1"])
			peak_5_99 = par[5]*TMath.Gaus(x[0], 5.99*luke, 0.05*d_std["OWC1"])
			peak_5_46 = par[6]*TMath.Gaus(x[0], 5.46*luke, 0.05*d_std["OWC1"])
			peak_4_97 = par[7]*TMath.Gaus(x[0], 4.97*luke, 0.05*d_std["OWC1"])
			peak_1_3  = 0.1*par[1]*TMath.Gaus(x[0], 1.3*luke, 0.05*d_std["OWC1"])
			list_peaks = [peak_10_4, peak_9_66, peak_8_98, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
			return flat_part + sum(list_peaks)

	fS1Gamma_deconv = TF1("S1Gamma_deconv", SurfGamma_deconv(), 0, 15, 9)

	fS1Gamma_deconv.SetNpx(5000)
	fS1Gamma_deconv.SetParameter(0,cst)
	fS1Gamma_deconv.SetParameter(1,20*fS1Gamma.GetParameter(1))
	fS1Gamma_deconv.SetParameter(2,20*fS1Gamma.GetParameter(2))
	fS1Gamma_deconv.SetParameter(3,20*fS1Gamma.GetParameter(3))
	fS1Gamma_deconv.SetParameter(4,20*fS1Gamma.GetParameter(4))
	fS1Gamma_deconv.SetParameter(5,20*fS1Gamma.GetParameter(5))
	fS1Gamma_deconv.SetParameter(6,20*fS1Gamma.GetParameter(6))
	fS1Gamma_deconv.SetParameter(7,20*fS1Gamma.GetParameter(7))
	fS1Gamma_deconv.SetParameter(8,fS1Gamma.GetParameter(8))

	file_S1Gamma_deconv =TFile(output_dir + "_Gamma_spectrum_extrapol_deconv.root", "update")
	fS1Gamma_deconv.Write()
	file_S1Gamma_deconv.Close()

    ##################
    #
	#      Gamma S2
	#
	##################

	fS2Gamma = TF1( "S2Gamma", SurfGamma(), 0, 40., 9 )
	fS2Gamma.SetNpx(500)
	fS2Gamma.SetParName(0,"p0")
	fS2Gamma.SetParName(1,"A_10.4")
	fS2Gamma.SetParName(2,"A_8.98")
	fS2Gamma.SetParName(3,"A_7.11")
	fS2Gamma.SetParName(4,"A_6.54")
	fS2Gamma.SetParName(5,"A_5.99")
	fS2Gamma.SetParName(6,"A_5.46")
	fS2Gamma.SetParName(7,"A_4.97")
	fS2Gamma.SetParName(8,"luke_shift")

	fS2Gamma.SetParameters(1,1,1,1,1,1,1,1,1)

	for i in range(8):
		fS2Gamma.SetParLimits(i, 0, 100)

	S2Gamma_endpt =0.1

	hS2Gamma.Draw()
	cst = hS2Gamma.GetBinWidth(20)*hS2Gamma.Integral(hS2Gamma.FindBin(11.5),hS2Gamma.FindBin(14.5))/(14.5-11.5)
	fS2Gamma.FixParameter(0,cst)
	hS2Gamma.Draw()
	hS2Gamma.Fit("S2Gamma","LL","",S2Gamma_endpt,15)
	raw_input()
	c1.Print("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Figures/" +bolo_name + "_S2Gamma_fit.png")

	
	slope_S2Gamma=fS2Gamma.Derivative(S2Gamma_endpt)
	offset_S2Gamma=fS2Gamma.Eval(S2Gamma_endpt)-S2Gamma_endpt*slope_S2Gamma

	class S2Gamma_bckg_extra:
		def __call__( self, x, par ):
			if x[0]>0.5:
				return fS2Gamma.Eval(x[0]) +par[0]
			else: 
				return slope_S2Gamma*x[0]+offset_S2Gamma

	fS2Gamma_extra = TF1( "S2Gamma_extra", S2Gamma_bckg_extra(), 0, 15., 1 )
	fS2Gamma_extra.SetNpx(500)
	fS2Gamma_extra.SetParameter(0,0)
	
	file_S2Gamma_extra =TFile(output_dir + "_Gamma_spectrum_extrapol.root", "update")
	fS2Gamma_extra.Write()
	file_S2Gamma_extra.Close()

	#Gamma S2 deconv
	fS2Gamma_deconv = TF1("S2Gamma_deconv", SurfGamma_deconv(), 0, 15, 9)

	fS2Gamma_deconv.SetNpx(5000)
	fS2Gamma_deconv.SetParameter(0,cst)
	fS2Gamma_deconv.SetParameter(1,20*fS2Gamma.GetParameter(1))
	fS2Gamma_deconv.SetParameter(2,20*fS2Gamma.GetParameter(2))
	fS2Gamma_deconv.SetParameter(3,20*fS2Gamma.GetParameter(3))
	fS2Gamma_deconv.SetParameter(4,20*fS2Gamma.GetParameter(4))
	fS2Gamma_deconv.SetParameter(5,20*fS2Gamma.GetParameter(5))
	fS2Gamma_deconv.SetParameter(6,20*fS2Gamma.GetParameter(6))
	fS2Gamma_deconv.SetParameter(7,20*fS2Gamma.GetParameter(7))
	fS2Gamma_deconv.SetParameter(8,fS2Gamma.GetParameter(8))

	file_S2Gamma_deconv =TFile(output_dir + "_Gamma_spectrum_extrapol_deconv.root", "update")
	fS2Gamma_deconv.Write()
	file_S2Gamma_deconv.Close()

    ##################
    #
	#    Heatonly
	#
	##################

	#Naming convention in case of ERA
	hist_name = "h1"

	#Get the TH1F heat data
	heat_path_file = gen_path + "Analyse_" + bolo_name +"/ROOT_files/ROOT_scaling/" + bolo_name + "_spectrum_heatonly.root"

	#Check the file exists then open it
	assert(os.path.isfile(heat_path_file))
	hist_heat, file_heatonly = PyRPl.open_ROOT_object(heat_path_file, hist_name)

	#Define triple exp fit function
	class triple_exp:
		def __call__( self, x, par ):
			return par[0]*TMath.Exp(par[1]*x[0]) + par[2]*TMath.Exp(par[3]*x[0]) + par[4]*TMath.Exp(par[5]*x[0])

	f_triple_exp = TF1("heat_extra", triple_exp(), 1, 15, 6)
	f_triple_exp.SetParameters(1,-1,1,-1, 1, -1)
	f_triple_exp.SetParameters(1,-1, 1, -0.1)
	f_triple_exp.SetParLimits(1,-10,0)
	f_triple_exp.SetParLimits(3,-1,0)
	# f_triple_exp.FixParameter(5,-0.32)
	# f_triple_exp.FixParameter(4,42)


	#Output directory
	output_heat_path = gen_path + "Analyse_" + bolo_name +"/ROOT_files/"
	

	hist_heat.Fit("heat_extra","LL","",1,15)
	gPad.SetLogy()
	# print hist_heat.Integral(hist_heat.FindBin(1.5), hist_heat.FindBin(5))
	# print f_triple_exp.Integral(1.5,5)/hist_heat.GetBinWidth(20)
	raw_input()

	heatonly_endpt = 1.1

	slope_heatonly=f_triple_exp.Derivative(heatonly_endpt)
	offset_heatonly=f_triple_exp.Eval(heatonly_endpt)-heatonly_endpt*slope_heatonly

	# class heatonly_bckg_extra:
	# 	def __call__( self, x, par ):
	# 		if x[0]>1.1:
	# 			return f_triple_exp.Eval(x[0]) +par[0]
	# 		else: 
	# 			return slope_heatonly*x[0]+offset_heatonly

	# fheatonly_extra = TF1( "heatonly_extra", heatonly_bckg_extra(), 0, 15., 1 )
	# fheatonly_extra.SetNpx(500)
	# fheatonly_extra.SetParameter(0,0)

	file_heat = TFile(output_heat_path + bolo_name + "_heatonly_spectrum_extrapol.root", "recreate")
	f_triple_exp.Write()
	file_heat.Close()
	
	##################
	#
	#    Heatonly 2D
	#
	##################
	
	file_tree   = TFile("../Fond_ERA_merged/"+bolo_name+"_lowmass_fond.root")
	tree        = file_tree.Get("t_merged")
	
	#Load standard cuts
	TCut_path_name = script_utils.create_directory("../Cut_files/")  
	TCut_file_name ="TCuts.txt" 
	file_TCut ="" 
	#Add an exception if the file does not exist
	try:
		file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
	except IOError:  
		script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
		sys.exit()
	
	# Load the cut values. 
	list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
	standard_cuts =""
	# Add a boolean flag to check if the bolo has its cuts in the file
	is_bolo_in_file =False
	for line in list_file_TCut_lines:
		if bolo_name == line[0]:
			standard_cuts = line[1]
			is_bolo_in_file = True
	assert(is_bolo_in_file)
	
	standard_cuts = standard_cuts + "&&KTH<1&&KTH>0&&0.5*(EIB+EID)<0"

	hheat = TH2F("heat2D", "heat2D", 200, -2, 15, 200, -2, 15)
	tree.Project("heat2D", "EC2:EC1", standard_cuts)

	hheat.Draw("colz")
	raw_input()

	file_heat2D = TFile(output_heat_path + bolo_name + "_heatonly_2D.root", "recreate")
	hheat.Write()
	file_heat2D.Close()
	

if __name__ == '__main__':

	bolo_name ="FID837"
	# fit_surface_relation(bolo_name)
	# fit_surface_KDE(bolo_name)
	create_extrapolated_bckg(bolo_name)