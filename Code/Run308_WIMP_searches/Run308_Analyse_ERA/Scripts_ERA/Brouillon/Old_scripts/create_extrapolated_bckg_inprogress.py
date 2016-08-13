#! /usr/bin/env python

from ROOT import *
import Poisson_file_handler as Poisson_fh
import PyROOTPlots as PyRPl
import os
import script_utils as script_utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from scipy import interpolate


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

	data_types = {"names": ("EC1_ERA", "EC2_ERA", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Beta = np.loadtxt(gen_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(gen_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Beta = coeff_EIA*arr_S1Beta["EIA"] + coeff_EIB*arr_S1Beta["EIB"]
	arr_EI_S2Beta = coeff_EIC*arr_S2Beta["EIC"] + coeff_EID*arr_S2Beta["EID"]

	arr_EC_S1Beta = coeff_EC1*arr_S1Beta["EC1_ERA"] + coeff_EC2*arr_S1Beta["EC2_ERA"]
	arr_EC_S2Beta = coeff_EC1*arr_S2Beta["EC1_ERA"] + coeff_EC2*arr_S2Beta["EC2_ERA"]

	arr_EI_S1Beta = np.array(arr_EI_S1Beta).astype(float)
	arr_EC_S1Beta = np.array(arr_EC_S1Beta).astype(float)

	arr_EI_S2Beta = np.array(arr_EI_S2Beta).astype(float)
	arr_EC_S2Beta = np.array(arr_EC_S2Beta).astype(float)

	gr_S1Beta = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta)
	gr_S2Beta = TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)

	PyRPl.process_TGraph(gr_S1Beta, X_title = "Heat", Y_title = "Ion", color=kRed)
	PyRPl.process_TGraph(gr_S2Beta, X_title = "Heat", Y_title = "Ion", color=kBlue)

	data_types = {"names": ("EC1_ERA", "EC2_ERA", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Pb = np.loadtxt(gen_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt(gen_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Pb = coeff_EIA*arr_S1Pb["EIA"] + coeff_EIB*arr_S1Pb["EIB"]
	arr_EI_S2Pb = coeff_EIC*arr_S2Pb["EIC"] + coeff_EID*arr_S2Pb["EID"]

	arr_EC_S1Pb = coeff_EC1*arr_S1Pb["EC1_ERA"] + coeff_EC2*arr_S1Pb["EC2_ERA"]
	arr_EC_S2Pb = coeff_EC1*arr_S2Pb["EC1_ERA"] + coeff_EC2*arr_S2Pb["EC2_ERA"]

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
	gr_S1Pb.Fit("pol3")
	cc.cd(4)
	list_h[3].Draw()
	gr_S2Pb.Draw("same*")
	gr_S2Pb.Fit("pol3")
	raw_input()

	fS1Beta = gr_S1Beta.GetFunction("pol3")
	fS1Beta.SetName("S1Beta")
	fS1Beta.Write()

	fS2Beta = gr_S2Beta.GetFunction("pol3")
	fS2Beta.SetName("S2Beta")
	fS2Beta.Write()

	fS1Pb = gr_S1Pb.GetFunction("pol3")
	fS1Pb.SetName("S1Pb")
	fS1Pb.Write()

	fS2Pb = gr_S2Pb.GetFunction("pol3")
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

	data_types = {"names": ("EC1_ERA", "EC2_ERA", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}

	######
	# Beta
	######
	arr_S1Beta = np.loadtxt(gen_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(gen_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Beta = coeff_EIA*arr_S1Beta["EIA"] + coeff_EIB*arr_S1Beta["EIB"]
	arr_EI_S2Beta = coeff_EIC*arr_S2Beta["EIC"] + coeff_EID*arr_S2Beta["EID"]

	arr_EC_S1Beta = coeff_EC1*arr_S1Beta["EC1_ERA"] + coeff_EC2*arr_S1Beta["EC2_ERA"]
	arr_EC_S2Beta = coeff_EC1*arr_S2Beta["EC1_ERA"] + coeff_EC2*arr_S2Beta["EC2_ERA"]

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

	arr_EC_S1Pb = coeff_EC1*arr_S1Pb["EC1_ERA"] + coeff_EC2*arr_S1Pb["EC2_ERA"]
	arr_EC_S2Pb = coeff_EC1*arr_S2Pb["EC1_ERA"] + coeff_EC2*arr_S2Pb["EC2_ERA"]

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
	d_std = Poisson_fh.open_true_event_FWHM_file(bolo_name)


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

	data_types = {"names": ("EC1_ERA", "EC2_ERA", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f", "f", "f")}
	arr_S1Beta = np.loadtxt(surf_bol_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(surf_bol_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EC_S1Beta = coeff_EC1*arr_S1Beta["EC1_ERA"] + coeff_EC2*arr_S1Beta["EC2_ERA"]
	arr_EC_S2Beta = coeff_EC1*arr_S2Beta["EC1_ERA"] + coeff_EC2*arr_S2Beta["EC2_ERA"]

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
			poly_part = par[0] + par[1]*x[0] + par[2]*x[0]*x[0]
			peak_10_4 = par[3]*TMath.Gaus(x[0], par[4], d_std["FWHEAT"])
			peak_8_98 = par[5]*TMath.Gaus(x[0], par[6], d_std["FWHEAT"])
			peak_1_3  = par[7]*TMath.Gaus(x[0], par[8], d_std["FWHEAT"])
			return poly_part + peak_10_4 + peak_8_98 + peak_1_3

	# Call FidGamma function to get derivative
	fFidGamma = TF1( "FidGamma", Gamma(), 0, 40., 9 )
	fFidGamma.SetNpx(500)
	fFidGamma.SetParName(0,"p0")
	fFidGamma.SetParName(1,"p1")
	fFidGamma.SetParName(2,"p2")
	fFidGamma.SetParName(3,"A_10.4")
	fFidGamma.SetParName(4,"mu_10.4")
	fFidGamma.SetParName(5,"A_8.98")
	fFidGamma.SetParName(6,"mu_8.98")
	fFidGamma.SetParName(7,"A_1.3")
	fFidGamma.SetParName(8,"mu_1.3")


	fFidGamma.SetParameter(0,1)
	fFidGamma.FixParameter(1,0)
	fFidGamma.FixParameter(2,0)
	fFidGamma.SetParameter(3,1)
	fFidGamma.SetParameter(5,1)
	fFidGamma.SetParameter(7,1)


	fFidGamma.SetParameter(4, 10.37)
	fFidGamma.SetParLimits(4, 10.2,11)
	
	fFidGamma.SetParameter(6, 9.1)
	fFidGamma.SetParLimits(6, 9,9.5)

	fFidGamma.SetParameter(8, 1.3)
	fFidGamma.SetParLimits(8, 1.1,1.5)
     
	FidGamma_endpt =0
	
	hFidGamma.Draw()
	hFidGamma.Fit("FidGamma","","",FidGamma_endpt,15)
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
	
	#Remove file if it already exists

	file_FidGamma_extra =TFile(output_dir + "_Gamma_spectrum_extrapol.root", "recreate")
	fFidGamma_extra.Write()
	file_FidGamma_extra.Close()

    ##################
    #
	#      Gamma S1
	#
	##################

	# Call S1 function to get derivative
	fS1 = TF1( "S1", Gamma(), 0, 40., 9 )
	fS1.SetNpx(500)
	fS1.SetParName(0,"p0")
	fS1.SetParName(1,"p1")
	fS1.SetParName(2,"p2")
	fS1.SetParName(3,"A_10.4")
	fS1.SetParName(4,"mu_10.4")
	fS1.SetParName(5,"A_8.98")
	fS1.SetParName(6,"mu_8.98")
	fS1.SetParName(7,"A_1.3")
	fS1.SetParName(8,"mu_1.3")

	fS1.SetParameter(0,1)
	fS1.FixParameter(1,0)
	fS1.FixParameter(2,0)
	fS1.SetParameter(3,1)
	fS1.SetParameter(5,1)
	fS1.SetParameter(7,1)

	fS1.SetParLimits(5,0,100)

	fS1.SetParameter(4, 8.2)
	fS1.SetParLimits(4, 8.1,8.5)
	
	fS1.SetParameter(6, 7)
	fS1.SetParLimits(6, 6.5,7.5)

	fS1.SetParameter(8, 1.)
	fS1.SetParLimits(8, 0.8,1.3)

	S1_endpt =0.1
	
	hS1Gamma.Draw()
	hS1Gamma.Fit("S1","","",S1_endpt,15)
	raw_input()
	c1.Print("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Figures/" +bolo_name + "_S1Gamma_fit.png")

	
	slope_S1=fS1.Derivative(S1_endpt)
	offset_S1=fS1.Eval(S1_endpt)-S1_endpt*slope_S1

	class S1_bckg_extra:
		def __call__( self, x, par ):
			if x[0]>0.5:
				return fS1.Eval(x[0]) +par[0]
			else: 
				return slope_S1*x[0]+offset_S1

	fS1_extra = TF1( "S1Gamma_extra", S1_bckg_extra(), 0, 15., 1 )
	fS1_extra.SetNpx(500)
	fS1_extra.SetParameter(0,0)
	
	file_S1_extra =TFile(output_dir + "_Gamma_spectrum_extrapol.root", "update")
	fS1_extra.Write()
	file_S1_extra.Close()


    ##################
    #
	#      Gamma S2
	#
	##################

	# Call S2 function to get derivative
	fS2 = TF1( "S2", Gamma(), 0, 40., 9 )
	fS2.SetNpx(500)
	fS2.SetParName(0,"p0")
	fS2.SetParName(1,"p1")
	fS2.SetParName(2,"p2")
	fS2.SetParName(3,"A_10.4")
	fS2.SetParName(4,"mu_10.4")
	fS2.SetParName(5,"A_8.98")
	fS2.SetParName(6,"mu_8.98")
	fS2.SetParName(7,"A_1.3")
	fS2.SetParName(8,"mu_1.3")

	fS2.SetParameter(0,1)
	fS2.FixParameter(1,0)
	fS2.FixParameter(2,0)
	fS2.SetParameter(3,1)
	fS2.SetParameter(5,1)
	fS2.SetParameter(7,1)

	fS2.SetParLimits(5,0,100)

	fS2.SetParameter(4, 8.2)
	fS2.SetParLimits(4, 8.1,8.5)
	
	fS2.SetParameter(6, 7)
	fS2.SetParLimits(6, 6.5,7.5)

	fS2.SetParameter(8, 1.)
	fS2.SetParLimits(8, 0.8,1.3)

	S2_endpt =0.1
	
	hS2Gamma.Draw()
	hS2Gamma.Fit("S2","","",S2_endpt,15)
	raw_input()
	c1.Print("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Figures/" +bolo_name + "_S2Gamma_fit.png")

	
	slope_S2=fS2.Derivative(S2_endpt)
	offset_S2=fS2.Eval(S2_endpt)-S2_endpt*slope_S2

	class S2_bckg_extra:
		def __call__( self, x, par ):
			if x[0]>0.5:
				return fS2.Eval(x[0]) +par[0]
			else: 
				return slope_S2*x[0]+offset_S2

	fS2_extra = TF1( "S2Gamma_extra", S2_bckg_extra(), 0, 15., 1 )
	fS2_extra.SetNpx(500)
	fS2_extra.SetParameter(0,0)
	
	file_S2_extra =TFile(output_dir + "_Gamma_spectrum_extrapol.root", "update")
	fS2_extra.Write()
	file_S2_extra.Close()


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
	f_triple_exp.SetParLimits(1,-10,0)
	f_triple_exp.SetParLimits(3,-10,0)
	f_triple_exp.SetParLimits(5,-10,0)

	f_triple_exp.FixParameter(4,0)

	#Output directory
	output_heat_path = gen_path + "Analyse_" + bolo_name +"/ROOT_files/"
	

	hist_heat.Fit("heat_extra","","",1,15)
	gPad.SetLogy()
	raw_input()

	heatonly_endpt = 1.1

	slope_heatonly=f_triple_exp.Derivative(heatonly_endpt)
	offset_heatonly=f_triple_exp.Eval(heatonly_endpt)-heatonly_endpt*slope_heatonly

	class heatonly_bckg_extra:
		def __call__( self, x, par ):
			if x[0]>1.1:
				return f_triple_exp.Eval(x[0]) +par[0]
			else: 
				return slope_heatonly*x[0]+offset_heatonly

	fheatonly_extra = TF1( "heatonly_extra", heatonly_bckg_extra(), 0, 15., 1 )
	fheatonly_extra.SetNpx(500)
	fheatonly_extra.SetParameter(0,0)

	file_heat = TFile(output_heat_path + bolo_name + "_heatonly_spectrum_extrapol.root", "recreate")
	f_triple_exp.Write()
	file_heat.Close()


if __name__ == '__main__':

	bolo_name ="FID837"
	# fit_surface_relation(bolo_name)
	fit_surface_KDE(bolo_name)
	# create_extrapolated_bckg(bolo_name)