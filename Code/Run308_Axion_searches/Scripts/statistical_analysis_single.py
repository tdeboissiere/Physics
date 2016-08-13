import utils_ana_jules as ut
import script_utils as script_utils
from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import matplotlib.pyplot as pl
import emcee
import triangle
import scipy.optimize as op
from matplotlib.ticker import MaxNLocator


def log_factorial(k):
	"""Utility to complete log of a factorial with approximate formula
	"""
	if k ==0:
		return 0
	else:
		return 0.5*np.log(2*TMath.Pi()*k) + k*np.log(k) - k + np.log(1+1./(12*k) + 1/(288.*k**2) -139./(51840*k**3)-571./(2488320*k**4) + 163879./(209018880*k**5))

def OOM_approach(bolo_name, exposure):
	"""Order of magnitude approach to estimate the limit
	
    Detail: 
		void

	Args:
		bolo_name = (str) bolo name
		exposure = (float) exposure in kg.d
	
	Returns:
		void

	Raises:
		void
	"""

	h, hfile = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
	fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/" + bolo_name + "_flux.root", bolo_name + "_flux_mass_0")


	class Flux:
		def __call__( self, x, par ):
			return exposure*(par[0]**4)*fflux.Eval(x[0])

	norm_flux = TF1("flux", Flux(), 0, 12, 1)
	norm_flux.SetParameter(0,20E-12)

	h.Draw()
	norm_flux.Draw("same")
	raw_input()


def simple_likelihood(bolo_name, exposure):
	"""A simple 1D frequentist likelihood
	
    Detail: 
		void

	Args:
		bolo_name = (str) bolo name
		exposure = (float) exposure in kg.d
		
	Returns:
		void

	Raises:
		void
	"""

	hdata, file_data = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
	fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", bolo_name + "_FidGamma")
	fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/" + bolo_name + "_flux_perkgdkeV.root", bolo_name + "_flux_mass_0")

	#Histograms parameters
	bin_X, min_X, max_X = hdata.GetNbinsX(), hdata.GetBinCenter(1) - 0.5*hdata.GetBinWidth(0), hdata.GetBinCenter(hdata.GetNbinsX()) + 0.5 * hdata.GetBinWidth(0)

	class Flux:
		def __call__( self, x, par ):
			return exposure*fflux.Eval(x[0]) + par[0]

	#Create bckg model histo in c/keV
	hbckg = TH1F("hbckg", "hbckg", bin_X, min_X, max_X)
	for i in range(1,bin_X+1) :
		hbckg.SetBinContent(i, fmodel.Eval(hbckg.GetBinCenter(i)))

	#Create signal histo in c/keV, gAe set to 1
	hsignal = TH1F("hsignal", "hsignal", bin_X, min_X, max_X)
	norm_flux = TF1("flux", Flux(), min_X, max_X, 1)
	norm_flux.SetParameter(0,0)
	for i in range(1, 301) :
		hsignal.SetBinContent(i, np.power(4E-11,4)*norm_flux.Eval(hsignal.GetBinCenter(i)))

	# hbckg.SetLineColor(kRed)
	# hsignal.SetLineColor(kGreen-3)
	# hbckg.Draw()
	# hdata.Draw("same")
	# hsignal.Draw("same")
	# raw_input()

	#Build likelihood
	class likelihood_hist:
		def __call__( self, x, par ):
			likelihood =0
			for i in range(1,bin_X+1):
				N_expected =hbckg.GetBinContent(i)+np.power(x[0],4)*norm_flux.Eval(hbckg.GetBinCenter(i)) #hsignal.GetBinContent(i)
				N_obs = hdata.GetBinContent(i) + np.power(5E-11,4)*norm_flux.Eval(hbckg.GetBinCenter(i))
				if N_expected>0:
					likelihood +=-N_expected+N_obs*np.log(N_expected)-log_factorial(N_obs)
			return -(likelihood + par[0])

	
	nll = TF1("nll", likelihood_hist(), 49E-12, 51E-12, 1)
	# nll.SetParameter(0,1298.47030762)
	nll.SetNpx(200)
	nll.Draw()
	raw_input()



def frequentist_MC(bolo_name, exposure):
	"""A simple 1D frequentist likelihood
	
    Detail: 
		void

	Args:
		bolo_name = (str) bolo name
		exposure = (float) exposure in kg.d
		
	Returns:
		void

	Raises:
		void
	"""

	hdata, file_data = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
	fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", bolo_name + "_FidGamma")
	fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/" + bolo_name + "_flux.root", bolo_name + "_flux_mass_0")

	#Histograms parameters
	bin_X, min_X, max_X = hdata.GetNbinsX(), hdata.GetBinCenter(1) - 0.5*hdata.GetBinWidth(0), hdata.GetBinCenter(hdata.GetNbinsX()) + 0.5 * hdata.GetBinWidth(0)

	class Flux:
		def __call__( self, x, par ):
			return exposure*fflux.Eval(x[0]) + par[0]
	norm_flux = TF1("flux", Flux(), min_X, max_X, 1)
	norm_flux.SetParameter(0,0)

	#Create bckg model histo in c/keV
	hbckg = TH1F("hbckg", "hbckg", bin_X, min_X, max_X)
	for i in range(1,bin_X+1) :
		hbckg.SetBinContent(i, fmodel.Eval(hbckg.GetBinCenter(i)))

	for gAe in [14, 0.5] :
		list_result = []
		for nsimu in range(100) :
			print gAe, nsimu
			#Simulate fake data 
			hdata_sim = TH1F("hdata_sim"+str(nsimu)+str(gAe), "hdata_sim"+str(nsimu)+str(gAe), bin_X, min_X, max_X)
			for i in range(1, bin_X+1): 
				hdata_sim.SetBinContent(i, hbckg.GetBinContent(i) + np.random.poisson(np.power(gAe*1E-12,4)*norm_flux.Eval(hbckg.GetBinCenter(i))))

			hdata_sim.SetLineColor(kRed)
			# hsignal.SetLineColor(kGreen-3)
			hbckg.Draw()
			hdata_sim.Add(hbckg,-1)
			hdata_sim.Draw("same")
			print hdata_sim.Integral()*hdata_sim.GetBinWidth(10)
			# hsignal.Draw("same")
			raw_input()

			#Build likelihood
			class likelihood_hist:
				def __call__( self, x, par ):
					likelihood =0
					for i in range(1,bin_X+1):
						N_expected =hbckg.GetBinContent(i)+np.power(x[0]*1E-12,4)*norm_flux.Eval(hbckg.GetBinCenter(i)) #hsignal.GetBinContent(i)
						N_obs = hdata_sim.GetBinContent(i)
						if N_expected>0:
							likelihood +=-N_expected+N_obs*np.log(N_expected)-log_factorial(N_obs)
					return -(likelihood + par[0])

			nll = TF1("nll", likelihood_hist(), 1E-2,1E3, 1)
			nll.SetParameter(0,0)
			list_result.append(nll.GetMinimumX())
		pl.hist(list_result, histtype = "step", bins = 30)
		pl.show()
		raw_input()

	# #Create signal histo in c/keV, gAe set to 1
	# hsignal = TH1F("hsignal", "hsignal", bin_X, min_X, max_X)
	# for i in range(1, 301) :
	# 	hsignal.SetBinContent(i, np.power(1E-12,4)*norm_flux.Eval(hsignal.GetBinCenter(i)))

	# hbckg.SetLineColor(kRed)
	# hsignal.SetLineColor(kGreen-3)
	# hbckg.Draw()
	# hdata.Draw("same")
	# hsignal.Draw("same")
	# raw_input()



	

	raw_input()

def MCMC_likelihood(bolo_name, exposure):
	"""2D MCMC likelihood
	
    Detail: 
		void

	Args:
		bolo_name = (str) bolo name
		exposure = (float) exposure in kg.d
		
	Returns:
		void

	Raises:
		void
	"""

	hdata, file_data = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
	fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", bolo_name + "_FidGamma")
	fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/" + bolo_name + "_flux_perkgdkeV.root", bolo_name + "_flux_mass_0")

	#Histograms parameters
	bin_X, min_X, max_X = hdata.GetNbinsX(), hdata.GetBinCenter(1) - 0.5*hdata.GetBinWidth(0), hdata.GetBinCenter(hdata.GetNbinsX()) + 0.5 * hdata.GetBinWidth(0)


	class Flux:
		def __call__( self, x, par ):
			# return exposure*fflux.Eval(x[0]) + par[0]
			return np.power(2.59E-11,4)*fflux.Eval(x[0]) + par[0]

	norm_flux = TF1("flux", Flux(), min_X, max_X, 1)
	h = TH1F("h", "h", 100, 2.5, 7.5)
	h.SetMaximum(1.4)
	h.SetMinimum(0)
	h.Draw()
	norm_flux.SetParameter(0,0)
	norm_flux.Draw("same")
	raw_input()

	#Create bckg, signal and data array in c/keV
	arr_bckg = np.array([fmodel.Eval(hdata.GetBinCenter(i)) for i in range(1, bin_X+1)]).astype(float)
	arr_signal = np.array([norm_flux.Eval(hdata.GetBinCenter(i)) for i in range(1, bin_X+1)]).astype(float)
	arr_data = np.array([hdata.GetBinContent(i) for i in range(1, bin_X+1)]).astype(float)
	arr_log_data = np.array([log_factorial(e) for e in arr_data]).astype(float)

	def lnprior(theta):
		bckg_norm, gAe = theta
		# if 0.9 < bckg_norm < 1.1 and 0.01 < gAe < 100 :
		# 	return 0
		# if 0.9 < bckg_norm < 1.1 and -35 < gAe < -20 :
		# 	return 0
		if  0.01 < gAe < 100:
			return 0.0 - 0.5*(np.log(2*np.pi) + np.power(bckg_norm-1,2))
		return -np.inf

	def lnlike(theta, arr_data, arr_log_data, arr_bckg, arr_signal):
		bckg_norm, gAe = theta
		# exp = bckg_norm*arr_bckg +np.exp(4*gAe)*arr_signal
		exp = bckg_norm*arr_bckg +np.power(gAe*1E-12,4)*arr_signal
		return np.sum(-exp + arr_data*np.log(exp) - arr_log_data)

	def lnprob(theta, arr_data, log_data, arr_bckg, arr_signal):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, arr_data, log_data, arr_bckg, arr_signal)


	nll    = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [1, 1], bounds = ((1E-2,10), (0.01, 100)), args=(arr_data, arr_log_data, arr_bckg, arr_signal))
	m_ml   = result["x"]
	print "Maximum likelihood:", m_ml

	ndim, nwalkers = 2, 100
	pos            = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler        = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(arr_data, arr_log_data, arr_bckg, arr_signal))

	# Clear and run the production chain.
	script_utils.print_utility("  Running MCMC  ")
	sampler.run_mcmc(pos, 500, rstate0 =np.random.get_state())
	script_utils.print_utility("  Done.  ")

	pl.clf()
	fig, axes                                    = pl.subplots(2, 1, sharex=True, figsize=(8, 9))
	axes[0].plot(sampler.chain[:, :, 0].T, color ="k", alpha=0.4)
	axes[0].yaxis.set_major_locator(MaxNLocator(5))
	axes[0].set_ylabel("N bckg")

	axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
	axes[1].yaxis.set_major_locator(MaxNLocator(5))
	axes[1].set_ylabel("gAe")
	axes[1].set_yscale("log")
	axes[1].set_xlabel("step number")
	
	fig.tight_layout()
	fig.savefig("../Figures/line-time.png")

	# Make the triangle plot.
	burnin = 200
	print sampler.chain.shape
	samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
	
	fig = triangle.corner(samples[0:samples.shape[0]:10,:], labels=["N bckg", "gAe"])
	fig.savefig("../Figures/line-triangle.png")
	pl.close("all")
	print np.percentile(samples[0:samples.shape[0]:10,1], 90)
	raw_input()


if __name__ == "__main__":

	list_bolo_name = ["FID837"]
	# OOM_approach("FID837", 73.2)
	simple_likelihood("FID837", 73.2)	
	# MCMC_likelihood("FID837", 73.2)	
	# MCMC_likelihood_1D("FID837", 73.2)
	# frequentist_MC("FID837", 73.2)