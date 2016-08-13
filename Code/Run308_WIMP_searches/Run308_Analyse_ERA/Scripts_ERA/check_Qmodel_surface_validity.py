from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Poisson_file_handler as Poisson_fh
import matplotlib.pylab as plt
from sklearn.neighbors import KernelDensity

def fit_surface_relation_FID808():

	"""Plot EI vs EH and fit
	for FID808, the bolo with surface calibration
    
	Detail:

	Args:

	Returns:
		void

	Raises:
		void
	"""


	# FWIA	FWIB	FWIC	FWID	FWC1	FWC2	VFID	VET	
	# 1.09	0.74	1.01	0.89	1.12	1.80	8	5.50

	#Open event files
	data_types = {"names": ("RUN", "SN", "EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("i", "i", "f", "f", "f", "f",  "f", "f")}

	arr_S1Beta = np.loadtxt("../Text_files/S1Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt("../Text_files/S2Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb = np.loadtxt("../Text_files/S1Pb_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt("../Text_files/S2Pb_heatremoved.txt", delimiter=",",  dtype=data_types)

	arr_EC_S1Beta, arr_EI_S1Beta = arr_S1Beta["EC1"], arr_S1Beta["EIB"]
	arr_EC_S2Beta, arr_EI_S2Beta = arr_S2Beta["EC1"], arr_S2Beta["EID"]
	arr_EC_S1Pb, arr_EI_S1Pb = arr_S1Pb["EC1"], arr_S1Pb["EIB"]
	arr_EC_S2Pb, arr_EI_S2Pb = arr_S2Pb["EC1"], arr_S2Pb["EID"]

	lS1Beta, lS2Beta, lS1Pb, lS2Pb = np.where(arr_EC_S1Beta<40), np.where(arr_EC_S2Beta<40), np.where(arr_EC_S1Pb<40), np.where(arr_EC_S2Pb<40)

	arr_EI_S1Beta, arr_EC_S1Beta = arr_EI_S1Beta[lS1Beta], arr_EC_S1Beta[lS1Beta]
	arr_EI_S2Beta, arr_EC_S2Beta = arr_EI_S2Beta[lS2Beta], arr_EC_S2Beta[lS2Beta]
	arr_EI_S1Pb, arr_EC_S1Pb     = arr_EI_S1Pb[lS1Pb], arr_EC_S1Pb[lS1Pb]
	arr_EI_S2Pb, arr_EC_S2Pb     = arr_EI_S2Pb[lS2Pb], arr_EC_S2Pb[lS2Pb]

	arr_EI_S1Beta, arr_EC_S1Beta = np.array(arr_EI_S1Beta).astype(float), np.array(arr_EC_S1Beta).astype(float)
	arr_EI_S2Beta, arr_EC_S2Beta = np.array(arr_EI_S2Beta).astype(float), np.array(arr_EC_S2Beta).astype(float)
	arr_EI_S1Pb, arr_EC_S1Pb = np.array(arr_EI_S1Pb).astype(float), np.array(arr_EC_S1Pb).astype(float)
	arr_EI_S2Pb, arr_EC_S2Pb = np.array(arr_EI_S2Pb).astype(float), np.array(arr_EC_S2Pb).astype(float)

	arr_Q_S1Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
	arr_Q_S2Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
	arr_Q_S1Pb   = arr_EI_S1Pb/((1+8./3)*arr_EC_S1Pb - arr_EI_S1Pb*5.5/3)
	arr_Q_S2Pb   = arr_EI_S2Pb/((1+8./3)*arr_EC_S2Pb - arr_EI_S2Pb*5.5/3)

	gr_S1Beta, gr_S2Beta   = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta),  TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)
	gr_S1Pb, gr_S2Pb       = TGraph(len(arr_EI_S1Pb), arr_EC_S1Pb, arr_EI_S1Pb),  TGraph(len(arr_EI_S2Pb), arr_EC_S2Pb, arr_EI_S2Pb)
	gr_QS1Beta, gr_QS2Beta = TGraph(len(arr_Q_S1Beta), arr_EC_S1Beta, arr_Q_S1Beta),  TGraph(len(arr_Q_S2Beta), arr_EC_S2Beta, arr_Q_S2Beta)
	gr_QS1Pb, gr_QS2Pb     = TGraph(len(arr_Q_S1Pb), arr_EC_S1Pb, arr_Q_S1Pb),  TGraph(len(arr_Q_S2Pb), arr_EC_S2Pb, arr_Q_S2Pb)

	PyRPl.process_TGraph(gr_S1Beta, X_title = "Heat", Y_title = "Ion", color=kOrange-3), PyRPl.process_TGraph(gr_S2Beta, X_title = "Heat", Y_title = "Ion", color=kBlue)
	PyRPl.process_TGraph(gr_S1Pb, X_title = "Heat", Y_title = "Ion", color=kBlack), PyRPl.process_TGraph(gr_S2Pb, X_title = "Heat", Y_title = "Ion", color=kGreen+2)

	PyRPl.process_TGraph(gr_QS1Beta, X_title = "Heat", Y_title = "Q", color=kOrange-3), PyRPl.process_TGraph(gr_QS2Beta, X_title = "Heat", Y_title = "Q", color=kBlue)
	PyRPl.process_TGraph(gr_QS1Pb, X_title = "Heat", Y_title = "Q", color=kBlack), PyRPl.process_TGraph(gr_QS2Pb, X_title = "Heat", Y_title = "Q", color=kGreen+2)

	out_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ROOT_files/"
	fout = TFile(out_path + "FID808_ion_heat_surface_relation.root" , "recreate")

	list_graph = [gr_S1Beta, gr_S2Beta, gr_S1Pb, gr_S2Pb]
	list_Qgraph = [gr_QS1Beta, gr_QS2Beta, gr_QS1Pb, gr_QS2Pb]

	cc = TCanvas("cc", "cc")
	cc.Divide(2,2)
	for i in range(4):
		cc.cd(i+1)
		list_graph[i].Draw("A*")
		list_graph[i].Fit("pol4")

	ccQ = TCanvas("ccQ", "ccQ")
	ccQ.Divide(2,2)
	for i in range(4):
		ccQ.cd(i+1)
		list_Qgraph[i].Draw("A*")

	list_func = [gr.GetFunction("pol4") for gr in list_graph]
	list_names = ["S1Beta", "S2Beta", "S1Pb", "S2Pb"]
	for i in range(4):
		list_func[i].SetName(list_names[i])
		list_func[i].Write()

	fout.Close()
	raw_input()

def check_Qmodel_surface_validity_FID808():

	""" Check Q model validity (look for straggle)
    
	Detail:

	Args:

	Returns:
		void

	Raises:
		void
	"""

	plt.ion()

	# FWIA	FWIB	FWIC	FWID	FWC1	FWC2	VFID	VET	
	# 1.09	0.74	1.01	0.89	1.12	1.80	8	5.50

	#Open event files
	data_types = {"names": ("RUN", "SN", "EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("i", "i", "f", "f", "f", "f",  "f", "f")}

	arr_S1Beta = np.loadtxt("../Text_files/S1Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt("../Text_files/S2Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb = np.loadtxt("../Text_files/S1Pb_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt("../Text_files/S2Pb_heatremoved.txt", delimiter=",",  dtype=data_types)

	arr_EC_S1Beta, arr_EI_S1Beta = arr_S1Beta["EC1"], arr_S1Beta["EIB"]
	arr_EC_S2Beta, arr_EI_S2Beta = arr_S2Beta["EC1"], arr_S2Beta["EID"]
	arr_EC_S1Pb, arr_EI_S1Pb = arr_S1Pb["EC1"], arr_S1Pb["EIB"]
	arr_EC_S2Pb, arr_EI_S2Pb = arr_S2Pb["EC1"], arr_S2Pb["EID"]

	lS1Beta, lS2Beta, lS1Pb, lS2Pb = np.where(arr_EC_S1Beta<40), np.where(arr_EC_S2Beta<40), np.where(arr_EC_S1Pb<40), np.where(arr_EC_S2Pb<40)

	arr_EI_S1Beta, arr_EC_S1Beta = arr_EI_S1Beta[lS1Beta], arr_EC_S1Beta[lS1Beta]
	arr_EI_S2Beta, arr_EC_S2Beta = arr_EI_S2Beta[lS2Beta], arr_EC_S2Beta[lS2Beta]
	arr_EI_S1Pb, arr_EC_S1Pb     = arr_EI_S1Pb[lS1Pb], arr_EC_S1Pb[lS1Pb]
	arr_EI_S2Pb, arr_EC_S2Pb     = arr_EI_S2Pb[lS2Pb], arr_EC_S2Pb[lS2Pb]

	arr_EI_S1Beta, arr_EC_S1Beta = np.array(arr_EI_S1Beta).astype(float), np.array(arr_EC_S1Beta).astype(float)
	arr_EI_S2Beta, arr_EC_S2Beta = np.array(arr_EI_S2Beta).astype(float), np.array(arr_EC_S2Beta).astype(float)
	arr_EI_S1Pb, arr_EC_S1Pb = np.array(arr_EI_S1Pb).astype(float), np.array(arr_EC_S1Pb).astype(float)
	arr_EI_S2Pb, arr_EC_S2Pb = np.array(arr_EI_S2Pb).astype(float), np.array(arr_EC_S2Pb).astype(float)

	arr_Q_S1Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
	arr_Q_S2Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
	arr_Q_S1Pb   = arr_EI_S1Pb/((1+8./3)*arr_EC_S1Pb - arr_EI_S1Pb*5.5/3)
	arr_Q_S2Pb   = arr_EI_S2Pb/((1+8./3)*arr_EC_S2Pb - arr_EI_S2Pb*5.5/3)

	gr_S1Beta, gr_S2Beta   = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta),  TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)
	gr_S1Pb, gr_S2Pb       = TGraph(len(arr_EI_S1Pb), arr_EC_S1Pb, arr_EI_S1Pb),  TGraph(len(arr_EI_S2Pb), arr_EC_S2Pb, arr_EI_S2Pb)
	gr_QS1Beta, gr_QS2Beta = TGraph(len(arr_Q_S1Beta), arr_EC_S1Beta, arr_Q_S1Beta),  TGraph(len(arr_Q_S2Beta), arr_EC_S2Beta, arr_Q_S2Beta)
	gr_QS1Pb, gr_QS2Pb     = TGraph(len(arr_Q_S1Pb), arr_EC_S1Pb, arr_Q_S1Pb),  TGraph(len(arr_Q_S2Pb), arr_EC_S2Pb, arr_Q_S2Pb)

	PyRPl.process_TGraph(gr_S1Beta, X_title = "Heat", Y_title = "Ion", color=kOrange-3), PyRPl.process_TGraph(gr_S2Beta, X_title = "Heat", Y_title = "Ion", color=kBlue)
	PyRPl.process_TGraph(gr_S1Pb, X_title = "Heat", Y_title = "Ion", color=kBlack), PyRPl.process_TGraph(gr_S2Pb, X_title = "Heat", Y_title = "Ion", color=kGreen+2)

	PyRPl.process_TGraph(gr_QS1Beta, X_title = "Heat", Y_title = "Q", color=kOrange-3), PyRPl.process_TGraph(gr_QS2Beta, X_title = "Heat", Y_title = "Q", color=kBlue)
	PyRPl.process_TGraph(gr_QS1Pb, X_title = "Heat", Y_title = "Q", color=kBlack), PyRPl.process_TGraph(gr_QS2Pb, X_title = "Heat", Y_title = "Q", color=kGreen+2)


	#Open fitted Heat/Ion relation file
	out_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ROOT_files/"
	func_convS1Beta, file_conv = PyRPl.open_ROOT_object(out_path + "FID808_ion_heat_surface_relation.root", "S1Beta")
	func_convS2Beta, file_conv = PyRPl.open_ROOT_object(out_path + "FID808_ion_heat_surface_relation.root", "S2Beta")
	func_convS1Pb, file_conv = PyRPl.open_ROOT_object(out_path + "FID808_ion_heat_surface_relation.root", "S1Pb")
	func_convS2Pb, file_conv = PyRPl.open_ROOT_object(out_path + "FID808_ion_heat_surface_relation.root", "S2Pb")

	list_Qgraph = [gr_QS1Beta, gr_QS2Beta, gr_QS1Pb, gr_QS2Pb]
	list_graph = [gr_S1Beta, gr_S2Beta, gr_S1Pb, gr_S2Pb]

	cc = TCanvas("cc", "cc")
	cc.Divide(2,2)
	for i in range(4):
		cc.cd(i+1)
		list_graph[i].Draw("A*")


	ccQ = TCanvas("ccQ", "ccQ")
	ccQ.Divide(2,2)
	for i in range(4):
		ccQ.cd(i+1)
		list_Qgraph[i].Draw("A*")

	# raw_input()
	d_zone_Beta = {}
	for i in range(1,11):
		d_zone_Beta[i] = [i-1, i+1]
	d_zone_Beta[14] = [10,18]
	d_zone_Beta[24] = [18,30]

	d_zone_Pb = {}
	for i in range(4,11):
		d_zone_Pb[i] = [i-1, i+1]
	d_zone_Pb[15] = [10,20]
	d_zone_Pb[20] = [15,25]
	d_zone_Pb[30] = [20,40]

	# d_zone_Beta = {}
	# d_zone_Beta[4] = [2,6]
	# d_zone_Beta[8] = [6,10]
	# d_zone_Beta[12.5] = [10,15]
	# # d_zone_Beta[17.5] = [15,25]
	# d_zone_Beta[40] = [30,50]

	# d_zone_Pb = {}
	# d_zone_Pb[14] = [8,20]
	# d_zone_Pb[30] = [20,40]
	
	list_S1Beta_res, list_S2Beta_res, list_S1Pb_res, list_S2Pb_res = [], [], [], []

	for key in sorted(d_zone_Beta.keys()):
		lS1 = np.where(np.logical_and(d_zone_Beta[key][0]<arr_EC_S1Beta, arr_EC_S1Beta<d_zone_Beta[key][1]))
		lS2 = np.where(np.logical_and(d_zone_Beta[key][0]<arr_EC_S2Beta, arr_EC_S2Beta<d_zone_Beta[key][1]))
		EC_S1, EI_S1 = arr_EC_S1Beta[lS1], arr_EI_S1Beta[lS1]
		EC_S2, EI_S2 = arr_EC_S2Beta[lS2], arr_EI_S2Beta[lS2]

		hS1 = TH1F("hS1", "hS1", 50, -5.0,5.0)
		hS2 = TH1F("hS2", "hS2", 50, -5.0,5.0)
		
		for ec, ei in zip(EC_S1, EI_S1):
			hS1.Fill(ei-func_convS1Beta(ec))

		for ec, ei in zip(EC_S2, EI_S2):
			hS2.Fill(ei-func_convS2Beta(ec))

		ccQ.cd(1)
		hS1.Draw()
		hS1.Fit("gaus")
		list_S1Beta_res.append(hS1.GetFunction("gaus").GetParameter(2))

		ccQ.cd(2)
		hS2.Draw()
		hS2.Fit("gaus")
		# raw_input()
		list_S2Beta_res.append(hS2.GetFunction("gaus").GetParameter(2))

	for key in sorted(d_zone_Pb.keys()):
		lS1 = np.where(np.logical_and(d_zone_Pb[key][0]<arr_EC_S1Pb, arr_EC_S1Pb<d_zone_Pb[key][1]))
		lS2 = np.where(np.logical_and(d_zone_Pb[key][0]<arr_EC_S2Pb, arr_EC_S2Pb<d_zone_Pb[key][1]))
		EC_S1, EI_S1 = arr_EC_S1Pb[lS1], arr_EI_S1Pb[lS1]
		EC_S2, EI_S2 = arr_EC_S2Pb[lS2], arr_EI_S2Pb[lS2]

		hS1 = TH1F("hS1", "hS1", 50, -5.0,5.0)
		hS2 = TH1F("hS2", "hS2", 50, -5.0,5.0)
		
		for ec, ei in zip(EC_S1, EI_S1):
			hS1.Fill(ei-func_convS1Pb(ec))

		for ec, ei in zip(EC_S2, EI_S2):
			hS2.Fill(ei-func_convS2Pb(ec))

		ccQ.cd(1)
		hS1.Draw()
		hS1.Fit("gaus")
		list_S1Pb_res.append(hS1.GetFunction("gaus").GetParameter(2))

		ccQ.cd(2)
		hS2.Draw()
		hS2.Fit("gaus")
		# raw_input()
		list_S2Pb_res.append(hS2.GetFunction("gaus").GetParameter(2))

	# plt.plot(sorted(d_zone_Beta.keys()), list_S1Beta_res, "r")
	# plt.plot(sorted(d_zone_Beta.keys()), list_S2Beta_res, "b")
	# plt.plot(sorted(d_zone_Pb.keys()), list_S1Pb_res, "k")
	# plt.plot(sorted(d_zone_Pb.keys()), list_S2Pb_res, "g")
	# plt.show()
	# raw_input()

	# FWIA	FWIB	FWIC	FWID	FWC1	FWC2	VFID	VET	
	# 1.09	0.74	1.01	0.89	1.12	1.80	8	5.50
	s_IB = 0.74/2.3548
	s_ID = 0.89/2.3548
	s_heat = 1.12/2.3548

	list_S1Beta_strag = [TMath.Sqrt(res**2 -s_IB**2 + (func_convS1Beta.Derivative(0)*s_heat)**2) for res in list_S1Beta_res] 
	list_S2Beta_strag = [TMath.Sqrt(res**2 -s_ID**2 + (func_convS2Beta.Derivative(0)*s_heat)**2) for res in list_S2Beta_res]
	list_S1Pb_strag   = [TMath.Sqrt(res**2 -s_IB**2 + (func_convS1Pb.Derivative(0)*s_heat)**2) for res in list_S1Pb_res]
	list_S2Pb_strag   = [TMath.Sqrt(res**2 -s_ID**2 + (func_convS2Pb.Derivative(0)*s_heat)**2) for res in list_S2Pb_res]

	print list_S1Beta_strag
	print list_S2Beta_strag
	print list_S1Pb_strag
	print list_S2Pb_strag

	plt.plot(sorted(d_zone_Beta.keys()), list_S1Beta_strag, "r")
	plt.plot(sorted(d_zone_Beta.keys()), list_S2Beta_strag, "b")
	plt.plot(sorted(d_zone_Pb.keys()), list_S1Pb_strag, "k")
	plt.plot(sorted(d_zone_Pb.keys()), list_S2Pb_strag, "g")
	plt.show()
	# raw_input()

	x_beta = np.array([0] + sorted(d_zone_Beta.keys())).astype(float)
	x_Pb   = np.array([0] + sorted(d_zone_Pb.keys())).astype(float)

	print  x_beta
	print np.array([list_S1Beta_strag[0]] + list_S1Beta_strag)

	gr_strag_S1Beta = TGraph(len(x_beta), x_beta , np.array([list_S1Beta_strag[0]] + list_S1Beta_strag))  
	gr_strag_S2Beta = TGraph(len(x_beta), x_beta, np.array([list_S2Beta_strag[0]] + list_S2Beta_strag))
	gr_strag_S1Pb   = TGraph(len(x_Pb),x_Pb, np.array([list_S1Pb_strag[0]] + list_S1Pb_strag))  
	gr_strag_S2Pb   = TGraph(len(x_Pb), x_Pb, np.array([list_S2Pb_strag[0]] + list_S2Pb_strag))

	list_strag_graph = [gr_strag_S1Beta, gr_strag_S2Beta, gr_strag_S1Pb, gr_strag_S2Pb]

	for i in range(1):
		ccQ.cd(i+1)
		list_strag_graph[i].Draw("A*")

	raw_input()

	class strag_S1Beta:
	   def __call__( self, x, par ):
	      return gr_strag_S1Beta.Eval(x[0])+par[0]

	class strag_S2Beta:
	   def __call__( self, x, par ):
	      return gr_strag_S2Beta.Eval(x[0])+par[0]

	class strag_S1Pb:
	   def __call__( self, x, par ):
	      return gr_strag_S1Pb.Eval(x[0])+par[0]

	class strag_S2Pb:
	   def __call__( self, x, par ):
	      return gr_strag_S2Pb.Eval(x[0])+par[0]

	fstrag_S1Beta = TF1("S1Beta_strag", strag_S1Beta(), 0, 25, 1)
	fstrag_S2Beta = TF1("S2Beta_strag", strag_S2Beta(), 0, 25, 1)
	fstrag_S1Pb = TF1("S1Pb_strag", strag_S1Pb(), 0, 25, 1)
	fstrag_S2Pb = TF1("S2Pb_strag", strag_S2Pb(), 0, 25, 1)

	list_fstrag = [fstrag_S1Beta, fstrag_S2Beta, fstrag_S1Pb, fstrag_S2Pb]

	fout = TFile(out_path + "FID808_straggle_std.root", "recreate")
	for func in list_fstrag:
		func.Write()
	fout.Close()

def check_Qmodel_surface_validity(bolo_name):

	""" Check Q model validity (look for straggle)
    
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
	d_std_true_events = Poisson_fh.open_true_event_FWHM_file(bolo_name)

    #Best estimator for heat: coefficients
	coeff_EC1, coeff_EC2 = float(d_est["HEAT"][:5]), 1 - float(d_est["HEAT"][:5])
	coeff_EIA, coeff_EIB = float(d_est["S1"][:5]), 1-float(d_est["S1"][:5])
	coeff_EIC, coeff_EID = float(d_est["S2"][:5]), 1-float(d_est["S2"][:5])

	#Open event files
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}
	arr_S1Beta = np.loadtxt(pop_path + bolo_name + "_S1Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt(pop_path + bolo_name + "_S2Beta_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb   = np.loadtxt(pop_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb   = np.loadtxt(pop_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Beta, arr_EI_S2Beta = coeff_EIA*arr_S1Beta["EIA"] + coeff_EIB*arr_S1Beta["EIB"], coeff_EIC*arr_S2Beta["EIC"] + coeff_EID*arr_S2Beta["EID"]
	arr_EC_S1Beta, arr_EC_S2Beta = coeff_EC1*arr_S1Beta["EC1"] + coeff_EC2*arr_S1Beta["EC2"], coeff_EC1*arr_S2Beta["EC1"] + coeff_EC2*arr_S2Beta["EC2"]
	arr_EI_S1Pb, arr_EI_S2Pb     = coeff_EIA*arr_S1Pb["EIA"] + coeff_EIB*arr_S1Pb["EIB"], coeff_EIC*arr_S2Pb["EIC"] + coeff_EID*arr_S2Pb["EID"]
	arr_EC_S1Pb, arr_EC_S2Pb     = coeff_EC1*arr_S1Pb["EC1"] + coeff_EC2*arr_S1Pb["EC2"], coeff_EC1*arr_S2Pb["EC1"] + coeff_EC2*arr_S2Pb["EC2"]

	lS1Beta, lS2Beta, lS1Pb, lS2Pb = np.where(arr_EC_S1Beta<40), np.where(arr_EC_S2Beta<40), np.where(arr_EC_S1Pb<40), np.where(arr_EC_S2Pb<40)

	arr_EI_S1Beta, arr_EC_S1Beta = arr_EI_S1Beta[lS1Beta], arr_EC_S1Beta[lS1Beta]
	arr_EI_S2Beta, arr_EC_S2Beta = arr_EI_S2Beta[lS2Beta], arr_EC_S2Beta[lS2Beta]
	arr_EI_S1Pb, arr_EC_S1Pb     = arr_EI_S1Pb[lS1Pb], arr_EC_S1Pb[lS1Pb]
	arr_EI_S2Pb, arr_EC_S2Pb     = arr_EI_S2Pb[lS2Pb], arr_EC_S2Pb[lS2Pb]

	arr_EI_S1Beta, arr_EC_S1Beta = np.array(arr_EI_S1Beta).astype(float), np.array(arr_EC_S1Beta).astype(float)
	arr_EI_S2Beta, arr_EC_S2Beta = np.array(arr_EI_S2Beta).astype(float), np.array(arr_EC_S2Beta).astype(float)
	arr_EI_S1Pb, arr_EC_S1Pb     = np.array(arr_EI_S1Pb).astype(float), np.array(arr_EC_S1Pb).astype(float)
	arr_EI_S2Pb, arr_EC_S2Pb     = np.array(arr_EI_S2Pb).astype(float), np.array(arr_EC_S2Pb).astype(float)

	arr_Q_S1Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
	arr_Q_S2Beta = arr_EI_S1Beta/((1+8./3)*arr_EC_S1Beta - arr_EI_S1Beta*5.5/3)
	arr_Q_S1Pb   = arr_EI_S1Pb/((1+8./3)*arr_EC_S1Pb - arr_EI_S1Pb*5.5/3)
	arr_Q_S2Pb   = arr_EI_S2Pb/((1+8./3)*arr_EC_S2Pb - arr_EI_S2Pb*5.5/3)

	gr_S1Beta, gr_S2Beta   = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta),  TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)
	gr_S1Pb, gr_S2Pb       = TGraph(len(arr_EI_S1Pb), arr_EC_S1Pb, arr_EI_S1Pb),  TGraph(len(arr_EI_S2Pb), arr_EC_S2Pb, arr_EI_S2Pb)
	gr_QS1Beta, gr_QS2Beta = TGraph(len(arr_Q_S1Beta), arr_EC_S1Beta, arr_Q_S1Beta),  TGraph(len(arr_Q_S2Beta), arr_EC_S2Beta, arr_Q_S2Beta)
	gr_QS1Pb, gr_QS2Pb     = TGraph(len(arr_Q_S1Pb), arr_EC_S1Pb, arr_Q_S1Pb),  TGraph(len(arr_Q_S2Pb), arr_EC_S2Pb, arr_Q_S2Pb)

	PyRPl.process_TGraph(gr_S1Beta, color=kBlack)
	PyRPl.process_TGraph(gr_S2Beta, color=kBlack)
	PyRPl.process_TGraph(gr_S1Pb, color=kBlack)
	PyRPl.process_TGraph(gr_S2Pb, color=kBlack)
	
	hbound = TH2F("hbound", "hbound", 100, 1.5, 15, 100, 0, 15)    
	PyRPl.process_TH2(hbound, X_title = "Combined Heat (keVee)", Y_title = "Combined Surface Ion (keV)")

	PyRPl.process_TGraph(gr_QS1Beta, X_title = "Heat", Y_title = "Q", color=kOrange-3), PyRPl.process_TGraph(gr_QS2Beta, X_title = "Heat", Y_title = "Q", color=kBlack)
	PyRPl.process_TGraph(gr_QS1Pb, X_title = "Heat", Y_title = "Q", color=kRed), PyRPl.process_TGraph(gr_QS2Pb, X_title = "Heat", Y_title = "Q", color=kBlack)


	#Open fitted Heat/Ion relation file
	Ana_path             = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/"
	convert_path = Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_ion_heat_surface_relation.root"
	func_convS1Beta, file_convS1Beta = PyRPl.open_ROOT_object(convert_path, "S1Beta")
	func_convS2Beta, file_convS2Beta = PyRPl.open_ROOT_object(convert_path, "S2Beta")
	func_convS1Pb, file_convS1Pb     = PyRPl.open_ROOT_object(convert_path, "S1Pb")
	func_convS2Pb, file_convS2Pb     = PyRPl.open_ROOT_object(convert_path, "S2Pb")


	list_Qgraph = [gr_QS1Beta, gr_QS2Beta, gr_QS1Pb, gr_QS2Pb]
	list_graph = [gr_S1Beta, gr_S2Beta, gr_S1Pb, gr_S2Pb]
	list_func_conv = [func_convS1Beta, func_convS2Beta, func_convS1Pb, func_convS2Pb]

	# #boundaries
	# nx, ny = 100, 100
	# x, y = np.meshgrid(np.linspace(2, 40, nx),np.linspace(2, 40, ny))
	# xy = np.vstack([x.ravel(), y.ravel()]).T
	# values = np.vstack([arr_EC_S1Beta, arr_EI_S1Beta])

	# # from sklearn.grid_search import GridSearchCV
	# # grid = GridSearchCV(KernelDensity(),{'bandwidth': np.linspace(0.1, 5.0, 30)},cv=10) # 10-fold cross-validation
	# # grid.fit(zip(*values))
	# # bw= grid.best_params_["bandwidth"]
	# # print bw
	# bw = 0.2

	# # -------------------------------------------------------
	# # Perform a kernel density estimate on the data using sklearn.
	# kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(zip(*values))
	# log_dens = np.exp(kde.score_samples(xy)) 
	
	# log_dens = log_dens.reshape(x.shape)

	# # levels = np.linspace(0, log_dens.max(), 25)
	# # plt.contourf(x, y, log_dens, levels=levels, cmap=plt.cm.Reds)
	# # plt.show()
	# # raw_input()

	# len_samples = 50000
	# samples = kde.sample(len_samples)
	# print samples.shape

	# arr_EC_S1Beta = np.array(samples[:,0]).astype(float)
	# arr_EI_S1Beta = np.array(samples[:,1]).astype(float)

	# gr_KDE_S1Beta   = TGraph(len_samples, np.array(samples[:,0]).astype(float), np.array(samples[:,1]).astype(float))
	# PyRPl.process_TGraph(gr_KDE_S1Beta, X_title = "Heat", Y_title = "Ion", color=kBlack)


	tt_S1Beta, ff_S1Beta = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_S1Beta_tree_strag_test.root", "t_new")
	tt_S2Beta, ff_S2Beta = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_S2Beta_tree_strag_test.root", "t_new")
	tt_S1Pb, ff_S1Pb = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_S1Pb_tree_strag_test.root", "t_new")
	tt_S2Pb, ff_S2Pb = PyRPl.open_ROOT_object(Ana_path + "Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_S2Pb_tree_strag_test.root", "t_new")

	list_trees = [tt_S1Beta, tt_S2Beta, tt_S1Pb, tt_S2Pb]
	list_var_tree = ["EIB", "EID", "EIB", "EID"]

	cc = TCanvas("cc", "cc")
	cc.Divide(2,2)
	for i in range(4):
		cc.cd(i+1)
		# gr_KDE_S1Beta.Draw("A*")
		hbound.Draw()
		list_graph[i].Draw("same*")
		list_func_conv[i].Draw("same")
		# list_trees[i].Draw(list_var_tree[i] + ":EC1", "", "same*")

	# cch = TCanvas("cch", "cch")
	# hsimu = TH1F("hsimu", "hsimu", 500, 0,30)
	# hsimu_nores = TH1F("hsimu_nores", "hsimu_nores", 500, 0,30)
	# h = TH1F("h", "h", 100, 0,30)
	# for elem in arr_EC_S1Beta:
	# 	h.Fill(elem)
	# tt.Project("hsimu", "EC2")
	# tt_nores.Project("hsimu_nores", "EC2")
	# h.SetLineColor(kRed)
	# hsimu.Scale(h.Integral()/float(hsimu.Integral()))
	# hsimu_nores.Scale(h.Integral()/float(hsimu_nores.Integral()))
	# # h.Draw()
	# hsimu_nores.SetLineColor(kBlue)
	# hsimu_nores.Draw()
	# hsimu.Draw("same")

	raw_input()
	cc.Print("../Analyse_" + bolo_name + "/Figures/" + bolo_name + "_ion_heat_surf_relation.eps")

	ccQ = TCanvas("ccQ", "ccQ")
	ccQ.Divide(2,2)
	# for i in range(4):
	# 	ccQ.cd(i+1)
	# 	list_Qgraph[i].Draw("A*")


	d_zone_Beta = {}
	d_zone_Beta[4] = [2,6]
	d_zone_Beta[8] = [6,10]
	d_zone_Beta[12.5] = [10,15]
	# d_zone_Beta[17.5] = [15,25]
	d_zone_Beta[40] = [30,50]

	d_zone_Pb = {}
	d_zone_Pb[14] = [5,20]
	d_zone_Pb[30] = [20,40]
	
	list_S1Beta_res, list_S2Beta_res, list_S1Pb_res, list_S2Pb_res = [], [], [], []

	for key in sorted(d_zone_Beta.keys()):
		lS1 = np.where(np.logical_and(d_zone_Beta[key][0]<arr_EC_S1Beta, arr_EC_S1Beta<d_zone_Beta[key][1]))
		lS2 = np.where(np.logical_and(d_zone_Beta[key][0]<arr_EC_S2Beta, arr_EC_S2Beta<d_zone_Beta[key][1]))
		EC_S1, EI_S1 = arr_EC_S1Beta[lS1], arr_EI_S1Beta[lS1]
		EC_S2, EI_S2 = arr_EC_S2Beta[lS2], arr_EI_S2Beta[lS2]

		hS1 = TH1F("hS1", "hS1", 50, -15.0, 15.0)
		hS2 = TH1F("hS2", "hS2", 50, -15.0, 15.0)
		
		for ec, ei in zip(EC_S1, EI_S1):
			hS1.Fill(ei-func_convS1Beta(ec))

		for ec, ei in zip(EC_S2, EI_S2):
			hS2.Fill(ei-func_convS2Beta(ec))

		ccQ.cd(1)
		hS1.Draw()
		hS1.Fit("gaus", "LL")
		list_S1Beta_res.append(hS1.GetFunction("gaus").GetParameter(2))

		ccQ.cd(2)
		hS2.Draw()
		hS2.Fit("gaus", "LL")
		raw_input()
		list_S2Beta_res.append(hS2.GetFunction("gaus").GetParameter(2))

	for key in sorted(d_zone_Pb.keys()):
		lS1 = np.where(np.logical_and(d_zone_Pb[key][0]<arr_EC_S1Pb, arr_EC_S1Pb<d_zone_Pb[key][1]))
		lS2 = np.where(np.logical_and(d_zone_Pb[key][0]<arr_EC_S2Pb, arr_EC_S2Pb<d_zone_Pb[key][1]))
		EC_S1, EI_S1 = arr_EC_S1Pb[lS1], arr_EI_S1Pb[lS1]
		EC_S2, EI_S2 = arr_EC_S2Pb[lS2], arr_EI_S2Pb[lS2]

		hS1 = TH1F("hS1", "hS1", 50, -5,5)
		hS2 = TH1F("hS2", "hS2", 50, -5,5)
		
		for ec, ei in zip(EC_S1, EI_S1):
			hS1.Fill(ei-func_convS1Pb(ec))

		for ec, ei in zip(EC_S2, EI_S2):
			hS2.Fill(ei-func_convS2Pb(ec))

		ccQ.cd(3)
		hS1.Draw()
		hS1.Fit("gaus", "LL")
		list_S1Pb_res.append(hS1.GetFunction("gaus").GetParameter(2))

		ccQ.cd(4)
		hS2.Draw()
		hS2.Fit("gaus", "LL")
		raw_input()
		list_S2Pb_res.append(hS2.GetFunction("gaus").GetParameter(2))

	# raw_input()

	# plt.plot(sorted(d_zone_Beta.keys()), list_S1Beta_res, "r")
	# plt.plot(sorted(d_zone_Beta.keys()), list_S2Beta_res, "b")
	# plt.plot(sorted(d_zone_Pb.keys()), list_S1Pb_res, "k")
	# plt.plot(sorted(d_zone_Pb.keys()), list_S2Pb_res, "g")
	# plt.show()
	# raw_input()

	list_S1Beta_strag = [TMath.Sqrt(res**2 -d_std_true_events["FWS1"]**2 + (func_convS1Beta.Derivative(0)*d_std_true_events["FWHEAT"])**2) for res in list_S1Beta_res] 
	list_S2Beta_strag = [TMath.Sqrt(res**2 -d_std_true_events["FWS2"]**2 + (func_convS2Beta.Derivative(0)*d_std_true_events["FWHEAT"])**2) for res in list_S2Beta_res]
	list_S1Pb_strag   = [TMath.Sqrt(res**2 -d_std_true_events["FWS1"]**2 + (func_convS1Pb.Derivative(0)*d_std_true_events["FWHEAT"])**2) for res in list_S1Pb_res]
	list_S2Pb_strag   = [TMath.Sqrt(res**2 -d_std_true_events["FWS2"]**2 + (func_convS2Pb.Derivative(0)*d_std_true_events["FWHEAT"])**2) for res in list_S2Pb_res]

	plt.plot(sorted(d_zone_Beta.keys()), list_S1Beta_strag, "r")
	plt.plot(sorted(d_zone_Beta.keys()), list_S2Beta_strag, "b")
	plt.plot(sorted(d_zone_Pb.keys()), list_S1Pb_strag, "k")
	plt.plot(sorted(d_zone_Pb.keys()), list_S2Pb_strag, "g")
	plt.ylim([0,5.5])
	plt.show()
	raw_input()

	x_beta = np.array([0] + sorted(d_zone_Beta.keys())).astype(float)
	x_Pb   = np.array([0] + sorted(d_zone_Pb.keys())).astype(float)

	print  x_beta
	print np.array([list_S1Beta_strag[0]] + list_S1Beta_strag)

	gr_strag_S1Beta = TGraph(len(x_beta), x_beta , np.array([list_S1Beta_strag[0]] + list_S1Beta_strag))  
	gr_strag_S2Beta = TGraph(len(x_beta), x_beta, np.array([list_S2Beta_strag[0]] + list_S2Beta_strag))
	#value 0.2 at 0 from FID808 calib
	gr_strag_S1Pb   = TGraph(len(x_Pb),x_Pb, np.array([0.2] + list_S1Pb_strag))  
	gr_strag_S2Pb   = TGraph(len(x_Pb), x_Pb, np.array([0.2] + list_S2Pb_strag))

	list_strag_graph = [gr_strag_S1Beta, gr_strag_S2Beta, gr_strag_S1Pb, gr_strag_S2Pb]

	for i in range(4):
		ccQ.cd(i+1)
		list_strag_graph[i].Draw("A*")
		list_strag_graph[i].Fit("pol3")

	raw_input()

def check_Pb_selection(bolo_name):

	""" Do a heat/FID ion plot to check the selection of Pb events
    
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
	arr_all  = np.loadtxt(pop_path + bolo_name + "_all_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb = np.loadtxt(pop_path + bolo_name + "_S1Pb_full_info.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt(pop_path + bolo_name + "_S2Pb_full_info.txt", delimiter=",",  dtype=data_types)

	arr_EI_S1Pb = coeff_EIB*arr_S1Pb["EIB"]+coeff_EID*arr_S1Pb["EID"]
	arr_EI_S2Pb = coeff_EIB*arr_S2Pb["EIB"]+coeff_EID*arr_S2Pb["EID"]
	arr_EI_all = coeff_EIB*arr_all["EIB"]+coeff_EID*arr_all["EID"]

	arr_EC_S1Pb = coeff_EC1*arr_S1Pb["EC1"]+coeff_EC2*arr_S1Pb["EC2"]
	arr_EC_S2Pb = coeff_EC1*arr_S2Pb["EC1"]+coeff_EC2*arr_S2Pb["EC2"]
	arr_EC_all = coeff_EC1*arr_all["EC1"]+coeff_EC2*arr_all["EC2"]

	hS1Pb = TH2F("hS1Pb", "hS1Pb",400, -2, 40, 400, -2, 40)
	hS2Pb = TH2F("hS2Pb", "hS2Pb",400, -2, 40, 400, -2, 40)
	hall = TH2F("hall", "hall",400, -2, 40, 400, -2, 40)

	PyRPl.fill_TH2(hS1Pb, arr_EC_S1Pb, arr_EI_S1Pb)
	PyRPl.fill_TH2(hS2Pb, arr_EC_S2Pb, arr_EI_S2Pb)
	PyRPl.fill_TH2(hall, arr_EC_all, arr_EI_all)

	PyRPl.process_TH2(hS1Pb, X_title = "EC", Y_title = "EI", color = kRed, marker_style = 20)
	PyRPl.process_TH2(hS2Pb, X_title = "EC", Y_title = "EI", color = kGreen+1, marker_style = 20)
	PyRPl.process_TH2(hall,  X_title = "EC", Y_title = "EI", color = kBlack, marker_style = 20)

	hall.Draw()
	hS1Pb.Draw("same")
	hS2Pb.Draw("same")
	raw_input()

bolo_name = "FID837"
# fit_surface_relation_FID808()
# check_Qmodel_surface_validity_FID808()
check_Qmodel_surface_validity(bolo_name)
# check_Pb_selection(bolo_name)