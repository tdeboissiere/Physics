	plt.figure()
	arr_EC_S1Beta = np.array([elem for elem in arr_EC_S1Beta if elem<60])

	# from sklearn.grid_search import GridSearchCV
	# grid = GridSearchCV(KernelDensity(),{'bandwidth': np.linspace(0.1, 5.0, 30)},cv=20) # 20-fold cross-validation
	# grid.fit(arr_EC_S1Beta[:, np.newaxis])
	# print grid.best_params_
	# yields 0.75

	X_plot = np.linspace(0, 60, 1000)[:, np.newaxis]
	from astroML.utils import check_random_state
	rng = check_random_state(None)
	# resample the data:
	ind = rng.randint(arr_EC_S1Beta.shape[0], size=(100, int(0.8*arr_EC_S1Beta.shape[0])))
	# print ind.shape, arr_EC_S1Beta.shape[0], ind[0,:]
	log_dens = np.zeros(X_plot.shape[0])
	for i in range(100):
		# kde = KernelDensity(kernel="gaussian", bandwidth=1.9).fit(arr_EC_S1Beta[ind[i]][:, np.newaxis])
		kde = KernelDensity(kernel="gaussian", bandwidth=0.75).fit(arr_EC_S1Beta[ind[i]][:, np.newaxis])
		log_dens+= kde.score_samples(X_plot)
		# plt.plot(X_plot[:, 0], np.exp(log_dens), '-',label ="kernel = '{0}'".format("gaussian"))


	filt = savitzky_golay(np.exp(0.01*log_dens), window_size=401, order=6)
	print (np.diff(filt) / np.diff(X_plot[:,0]))[100], X_plot[:,0][100]
	plt.hist(arr_EC_S1Beta,  bins=100, normed=1, histtype='step')	
	plt.plot(X_plot[:, 0], np.exp(0.01*log_dens), '-',label ="kernel = '{0}'".format("gaussian"))
	plt.plot(X_plot[:, 0], filt, '-')	
	plt.ylim([0,0.12])
	plt.show()



	gr = TGraph(1000, np.array(X_plot[:,0]).astype(float), np.array(filt).astype(float))
	gr.Draw("A*")
	raw_input()

	S1Beta_endpt = 6
	endpt_index = int(1000*S1Beta_endpt/60.)
	slope_S1Beta=(np.diff(filt) / np.diff(X_plot[:,0]))[endpt_index]
	offset_S1Beta=gr.Eval(S1Beta_endpt)-S1Beta_endpt*slope_S1Beta

	class S1Beta_bol:
		def __call__( self, x, par ):
			if x[0]>S1Beta_endpt:
				return gr.Eval(x[0]) + par[0]
			else:
				return slope_S1Beta*x[0]+offset_S1Beta

	fS1Beta_bol_extra = TF1( "S1Beta_extra", S1Beta_bol(), 0, 60., 1 )
	fS1Beta_bol_extra.SetParameter(0,0)
	hS1Beta_bol.Draw()
	fS1Beta_bol_extra.Draw("same")

	file_S1Beta_bol_extra =TFile(output_dir + "_Beta_from_bolo_spectrum_extrapol.root", "recreate")
	fS1Beta_bol_extra.Write()
	file_S1Beta_bol_extra.Close()