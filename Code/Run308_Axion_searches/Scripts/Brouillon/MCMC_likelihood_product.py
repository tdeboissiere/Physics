def MCMC_likelihood_product(list_bolo_name, mass, d_exposure):
	"""2D MCMC likelihood (multiply individual likelihoods)
	
    Detail: 
		void

	Args:
		mass = (float) axion mass
		
	Returns:
		void

	Raises:
		void
	"""

	#Dic to store the data of each relevant bolometer 
	d_data = {}
	for bolo_name in list_bolo_name:

		d_data[bolo_name] = {"arr_bckg": [], "arr_signal" : [], "arr_data": [], "arr_log_data" : []}

		hdata, file_data = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
		fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", bolo_name + "_FidGamma")
		fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/" + bolo_name + "_flux_perkgdkeV.root", bolo_name + "_flux_mass_" + str(mass))	

		#Histograms parameters
		bin_X, min_X, max_X = hdata.GetNbinsX(), hdata.GetBinCenter(1) - 0.5*hdata.GetBinWidth(0), hdata.GetBinCenter(hdata.GetNbinsX()) + 0.5 * hdata.GetBinWidth(0)

		class Flux:
			def __call__( self, x, par ):
				return d_exposure[bolo_name]*fflux.Eval(x[0]) + par[0]

		norm_flux = TF1("flux", Flux(), min_X, max_X, 1)
		norm_flux.SetParameter(0,0)

		#Create bckg, signal and data array in c/keV
		d_data[bolo_name]["arr_bckg"]     = np.array([fmodel.Eval(hdata.GetBinCenter(i)) for i in range(1, bin_X+1)]).astype(float)
		d_data[bolo_name]["arr_signal"]   = np.array([norm_flux.Eval(hdata.GetBinCenter(i)) for i in range(1, bin_X+1)]).astype(float)
		d_data[bolo_name]["arr_data"]     = np.array([hdata.GetBinContent(i) for i in range(1, bin_X+1)]).astype(float)
		d_data[bolo_name]["arr_log_data"] = np.array([log_factorial(e) for e in d_data[bolo_name]["arr_data"]]).astype(float)


	def lnprior(theta):
		gAe = theta
		if  0.01 < gAe < 100:
			return 0.0
		return -np.inf

	def lnlike(theta):
		gAe = theta
		ll_tot = 0
		for bolo_name in list_bolo_name :
			exp_bolo = d_data[bolo_name]["arr_bckg"] +np.power(gAe*1E-12,4)*d_data[bolo_name]["arr_signal"]
			ll_tot += np.sum(-exp_bolo + d_data[bolo_name]["arr_data"]*np.log(exp_bolo+1E-8) - d_data[bolo_name]["arr_log_data"])
		return ll_tot

	def lnprob(theta):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta)

	ndim, nwalkers = 1, 100
	pos            = [10 + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler        = emcee.EnsembleSampler(nwalkers, ndim, lnprob)

	# Clear and run the production chain.
	script_utils.print_utility("  Running MCMC  ")
	sampler.run_mcmc(pos, 500, rstate0 =np.random.get_state())
	script_utils.print_utility("  Done.  ")

	# pl.clf()
	# fig, axes                                    = pl.subplots(2, 1, sharex=True, figsize=(8, 9))
	# axes[0].plot(sampler.chain[:, :, 0].T, color ="k", alpha=0.4)
	# axes[0].yaxis.set_major_locator(MaxNLocator(5))
	# axes[0].set_ylabel("N bckg")

	# axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
	# axes[1].yaxis.set_major_locator(MaxNLocator(5))
	# axes[1].set_ylabel("gAe")
	# axes[1].set_yscale("log")
	# axes[1].set_xlabel("step number")
	
	# fig.tight_layout()
	# fig.savefig("../Figures/line-time.png")

	# Make the triangle plot.
	burnin = 200
	print sampler.chain.shape
	samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
	
	# fig = triangle.corner(samples[0:samples.shape[0]:10,:], labels=["N bckg", "gAe"])
	# fig.savefig("../Figures/line-triangle.png")
	# pl.close("all")
	print np.percentile(samples[0:samples.shape[0]:10], 90)
	raw_input()