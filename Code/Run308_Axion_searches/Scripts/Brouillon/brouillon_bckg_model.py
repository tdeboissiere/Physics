def get_elec_recoil_bckgmodel(list_bolo_name):
	"""Create ROOT file with elec recoil spectrum bckg model TF1
	
    Detail: 
		void

	Args:
		list_bolo_name = (str) list of bolo names
		
	Returns:
		void

	Raises:
		void
	"""

	for bolo_name in list_bolo_name :

		ftrigger, file_trigger = PyRPl.open_ROOT_object("../ROOT_files/Axion/Trigger_eff/" + bolo_name + "_trigger_eff.root", "trigger_eff")
		fcut, file_cut = PyRPl.open_ROOT_object("../ROOT_files/Axion/Cut_eff/" + bolo_name + "_cut_eff.root", "cut_eff")

		#Get data tree
		t, f = PyRPl.open_ROOT_object("../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root", "t_axion")
		#Get cut
		cut  = ut.get_electronic_cut_line(bolo_name)
		#Create hist of FWbest reso and get theoretical 0 keV resolution
		hFWbest = TH1F("hFWbest_" + bolo_name , "hFWbest_" + bolo_name, 100, 0, 15)
		t.Project("hFWbest_" + bolo_name, "FWbest")
		reso_0 = hFWbest.GetMean()/2.3548

		#Get peak amplitude
		peak_amp = np.loadtxt("../Tables/gamma_background_Fid_" + bolo_name + ".txt", skiprows=2).astype(float)
		d_amp = {}
		for i in range(peak_amp.shape[0]):
			d_amp[peak_amp[i][0]]=peak_amp[i][1]

		##############
		#Fit functions
		##############

		class GammaBckg_fit:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = d_amp[0]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = d_amp[10.37]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = d_amp[9.66]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = d_amp[8.98]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				peak_7_71  = d_amp[7.71]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso**2)
				peak_7_11  = d_amp[7.11]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11]
				return flat_part + sum(list_peaks)

		class GammaBckg_full:
			#Use this class to obtain the full background model for statistical analysis
			def __call__( self, x, par ):
				flat_part  = d_amp[0]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = d_amp[10.37]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = d_amp[9.66]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = d_amp[8.98]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				peak_7_71  = d_amp[7.71]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso**2)
				peak_7_11  = d_amp[7.11]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso**2)
				peak_6_54  = d_amp[6.54]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso**2)
				peak_5_99  = d_amp[5.99]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso**2)
				peak_5_46  = d_amp[5.46]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso**2)
				peak_4_97  = d_amp[4.97]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso**2)
				peak_1_3   = d_amp[1.3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.3)**2/reso**2) + d_amp[1.19]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.19)**2/reso**2)
				peak_1_3   += d_amp[1.1]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.1)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return (flat_part + sum(list_peaks))*ftrigger.Eval(x[0])*fcut.Eval(x[0])

		class GammaBckg_deconv:
			#Use this class for simulations
			def __call__( self, x, par ):
				flat_part  = d_amp[0]
				reso_deconv = TMath.Abs(par[1]*x[0]+1E-3)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso_deconv)

				peak_10_37 = d_amp[10.37]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso_deconv**2)
				peak_9_66  = d_amp[9.66]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso_deconv**2)
				peak_8_98  = d_amp[8.98]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso_deconv**2)
				peak_7_71  = d_amp[7.71]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso_deconv**2)
				peak_7_11  = d_amp[7.11]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso_deconv**2)
				peak_6_54  = d_amp[6.54]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso_deconv**2)
				peak_5_99  = d_amp[5.99]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso_deconv**2)
				peak_5_46  = d_amp[5.46]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso_deconv**2)
				peak_4_97  = d_amp[4.97]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso_deconv**2)
				peak_1_3   = d_amp[1.3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.3)**2/reso_deconv**2) + d_amp[1.19]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.19)**2/reso_deconv**2)
				peak_1_3   += d_amp[1.1]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.1)**2/reso_deconv**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				#No trigger efficiency taken into account
				return flat_part + sum(list_peaks)

		h = TH1F("h_" + bolo_name , "h_" + bolo_name, 100, 0, 15)
		t.Project("h_" + bolo_name, "Ebest", cut)
		#Creat fitting function and output file
		fout = TFile("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", "recreate")
		fFidGamma_fit = TF1(bolo_name +  "_FidGamma_fit", GammaBckg_fit(), 0, 15., 2 )
		fFidGamma_fit.SetParName(0,"offset")
		fFidGamma_fit.SetParName(1,"reso")
		fFidGamma_fit.SetParameters(1,0.1)
		#Convert hist to /keV
		h.Scale(1./h.GetBinWidth(5))
		#Draw and fit
		h.Draw()
		h.Fit(bolo_name +  "_FidGamma_fit", "LL", "", 7.5, 11.5)
		raw_input()

		fFidGamma = TF1(bolo_name +  "_FidGamma", GammaBckg_full(), 0, 15., 2 )
		fFidGamma_noeff = TF1(bolo_name +  "_FidGamma_noeff", GammaBckg_noeff(), 0, 15., 2 )
		fFidGamma_deconv = TF1(bolo_name +  "_FidGamma_deconv", GammaBckg_deconv(), 0, 15., 2 )

		fFidGamma.SetNpx(500)
		fFidGamma.SetParameter(0,fFidGamma_fit.GetParameter(0))
		fFidGamma.SetParameter(1,fFidGamma_fit.GetParameter(1))

		fFidGamma_noeff.SetNpx(500)
		fFidGamma_noeff.SetParameter(0,fFidGamma_fit.GetParameter(0))
		fFidGamma_noeff.SetParameter(1,fFidGamma_fit.GetParameter(1))

		fFidGamma_deconv.SetNpx(500)
		fFidGamma_deconv.SetParameter(0,fFidGamma_fit.GetParameter(0))
		fFidGamma_deconv.SetParameter(1,fFidGamma_fit.GetParameter(1))

		print bolo_name +","+ str(reso_0) +","+ str(fFidGamma.GetParameter(1))
		#Save to file
		fFidGamma.Write()
		fFidGamma_noeff.Write()
		fFidGamma_deconv.Write()
		fout.Close()
		del fFidGamma
		del fFidGamma_noeff
		del fFidGamma_deconv
		del fFidGamma_fit
		del fout