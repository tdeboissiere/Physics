import scipy
from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import matplotlib.pylab as plt
import statsmodels.api as sm 

def check_BDT_simulations_slice_KS(bolo_name, analysis_type, mass):

	""" Verify the 2D distributions (slice the EC, EFID plane in bands and compute KS)
    
	Detail:
		Example: take 1<EC<2 and EI in [0,15] and compute KS test on ion distributions

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
	BDT_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_better/BDT_" + bolo_name + "/" + analysis_type + "/"

	ttrue,ftrue = PyRPl.open_ROOT_object("../Fond_ERA_merged/" + bolo_name + "_" + analysis_type + "_lowmass_fond.root", "t_merged")
	tsimu, fsimu = PyRPl.open_ROOT_object(BDT_path +"True_events/ROOT_files/" + bolo_name + "_true_events_tree.root", "t_new0")

	print "true: ", ttrue.GetEntries("EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")
	print "simu: ", tsimu.GetEntries("EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")

	ttrue.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)>>hist(1000,-2,15,1000,-2,15", "")
	tsimu.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)>>hist2(1000,-2,15,1000,-2,15", "")

	# ttrue.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)>>hist(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")
	# tsimu.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)>>hist2(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")
	# ttrue.Draw("0.414*EIB+(1-0.414)*EID:0.574*EC1+(1-0.574)*EC2>>hist(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0 && EIA<2 && EIC<2")
	# tsimu.Draw("0.414*EIB+(1-0.414)*EID:0.574*EC1+(1-0.574)*EC2>>hist2(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")

	# ttrue.Draw("EIB:EID>>hist(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0 && EIA<2 && EIC<2")
	# tsimu.Draw("EIB:EID>>hist2(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")
	# ttrue.Draw("EC1:EC2>>hist(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0 && EIA<2 && EIC<2")
	# tsimu.Draw("EC1:EC2>>hist2(1000,-2,15,1000,-2,15", "EIB>2 && EID>2 && abs(EIB-EID)<1 && abs(EIB-EC1)<4 && EC1>2.0")

	hist.SetMarkerColor(kRed)
	hist.SetMarkerStyle(20)
	hist2.SetMarkerStyle(20)
	hist.Draw()
	hist2.Draw("same")

	raw_input()

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}

	arr_true = np.loadtxt(pop_path + bolo_name + "_true_events_all.txt", delimiter=",",  dtype=data_types)
	arr_simu = np.loadtxt(pop_path + bolo_name + "_simu_events_all.txt", delimiter=",",  dtype=data_types)

	EI_true = 0.5*(arr_true["EIB"]+arr_true["EID"])
	EC_true = 0.5*(arr_true["EC1"]+arr_true["EC2"])

	EI_simu = 0.5*(arr_simu["EIB"]+arr_simu["EID"])
	EC_simu = 0.5*(arr_simu["EC1"]+arr_simu["EC2"])

	h2Darr = TH2F("h2Darr", "h2Darr", 1000, -2, 15, 1000, -2, 15)
	h2Dsimu = TH2F("h2Dsimu", "h2Dsimu", 1000, -2, 15, 1000, -2, 15)

	for i in range(EI_true.shape[0]):
		h2Darr.Fill(EC_true[i], EI_true[i])
	for i in range(EI_simu.shape[0]):
		h2Dsimu.Fill(EC_simu[i],EI_simu[i])

	PyRPl.process_TH2(h2Darr, X_title = "EC", Y_title = "EI", color = kRed)
	PyRPl.process_TH2(h2Dsimu, X_title = "EC", Y_title = "EI", color = kBlack)

	h2Darr.Draw()
	h2Dsimu.Draw("same")

	#Slices on EC
	for EC in range(2,15):
		l_true = np.where(np.logical_and(EC_true>EC-1 , EC_true<EC))
		l_simu = np.where(np.logical_and(EC_simu>EC-1 , EC_simu<EC))

		slice_EI_true = EI_true[l_true]
		slice_EI_simu = EI_simu[l_simu]

		print scipy.stats.ks_2samp(slice_EI_true, slice_EI_simu),"   ", 1.36*sqrt(len(slice_EI_true) + len(slice_EI_simu))/sqrt(len(slice_EI_true) * len(slice_EI_simu))

		true_cdf = sm.distributions.ECDF(slice_EI_true)
		simu_cdf = sm.distributions.ECDF(slice_EI_simu)

		x_true = np.linspace(min(slice_EI_true), max(slice_EI_true))
		x_simu = np.linspace(min(slice_EI_simu), max(slice_EI_simu))
		y_true = true_cdf(x_true)
		y_simu = simu_cdf(x_simu)

		plt.step(x_true, y_true, "r", label = "True IonFid CDF @ EC in [" + str(EC-1) + "," + str(EC) + "]" )
		plt.step(x_simu, y_simu, "k", label = "Simu IonFid CDF @ EC in [" + str(EC-1) + "," + str(EC) + "]")
		plt.legend(loc="upper left", prop={"size":10})

		plt.show()
		raw_input()
		plt.clf()

	#Slices on EI
	for EI in range(1,15):
		l_true = np.where(np.logical_and(EI_true>EI-1 , EI_true<EI))
		l_simu = np.where(np.logical_and(EI_simu>EI-1 , EI_simu<EI))

		slice_EC_true = EC_true[l_true]
		slice_EC_simu = EC_simu[l_simu]

		print scipy.stats.ks_2samp(slice_EC_true, slice_EC_simu),"   ", 1.36*sqrt(len(slice_EC_true) + len(slice_EC_simu))/sqrt(len(slice_EC_true) * len(slice_EC_simu))

		true_cdf = sm.distributions.ECDF(slice_EC_true)
		simu_cdf = sm.distributions.ECDF(slice_EC_simu)

		x_true = np.linspace(min(slice_EC_true), max(slice_EC_true))
		x_simu = np.linspace(min(slice_EC_simu), max(slice_EC_simu))
		y_true = true_cdf(x_true)
		y_simu = simu_cdf(x_simu)

		plt.step(x_true, y_true, "r", label = "True IonFid CDF @ EI in [" + str(EI-1) + "," + str(EI) + "]" )
		plt.step(x_simu, y_simu, "k", label = "Simu IonFid CDF @ EI in [" + str(EI-1) + "," + str(EI) + "]")
		plt.legend(loc="upper left", prop={"size":10})

		plt.show()
		raw_input()
		plt.clf()


bolo_name = "FID837"
analysis_type = "ana_1.5_min2_50"
mass = "5"
check_BDT_simulations_slice_KS(bolo_name, analysis_type, mass)