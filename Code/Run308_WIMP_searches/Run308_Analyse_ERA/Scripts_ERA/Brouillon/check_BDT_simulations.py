from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import BDT_file_handler as BDT_fh
import matplotlib.pylab as plt
from sklearn.neighbors import KernelDensity


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
	arr_true_ion = 0.5*(arr_true_events["EIB"]+ arr_true_events["EID"])
	lheatonly = np.where(arr_true_ion<1)
	lFidGamma = np.where(arr_true_ion>1.5)

	np.savetxt(pop_path + bolo_name + "_true_events_all.txt", arr_true_events, delimiter = ",")
	np.savetxt(pop_path + bolo_name + "_true_events_heatonly.txt", arr_true_events[lheatonly], delimiter = ",")
	np.savetxt(pop_path + bolo_name + "_true_events_FidGamma.txt", arr_true_events[lFidGamma], delimiter = ",")

	arr_simu_events = root2array(BDT_path +"True_events/ROOT_files/" + bolo_name + "_true_events_tree.root", "t_new0")
	arr_simu_ion = 0.5*(arr_simu_events["EIB"]+ arr_simu_events["EID"])
	lheatonly = np.where(arr_simu_ion<1)
	lFidGamma = np.where(arr_simu_ion>1.5)

	np.savetxt(pop_path + bolo_name + "_simu_events_all.txt", arr_simu_events, delimiter = ",")
	np.savetxt(pop_path + bolo_name + "_simu_events_heatonly.txt", arr_simu_events[lheatonly], delimiter = ",")
	np.savetxt(pop_path + bolo_name + "_simu_events_FidGamma.txt", arr_simu_events[lFidGamma], delimiter = ",")


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
	data_types = {"names": ("EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("f", "f", "f", "f",  "f", "f")}
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

	list_channel = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID"]
	list_evt = ["heatonly", "FidGamma", "S1Gamma", "S2Gamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb"]
	list_arr = [arr_heatonly_after_cuts, arr_FidGamma_after_cuts, arr_S1Gamma, arr_S2Gamma, arr_S1Beta, arr_S2Beta, arr_S1Pb, arr_S2Pb]
	list_simu = [simu_heatonly, simu_FidGamma, simu_S1Gamma, simu_S2Gamma, simu_S1Beta, simu_S2Beta, simu_S1Pb, simu_S2Pb]
	list_canvas = [TCanvas("cc"+evt, "cc"+ evt) for evt in list_evt]



	# #Check the channels where no signal is expected

	# d_hist_arr={}
	# d_hist_simu={}
	# for evt in list_evt:
	# 	d_hist_arr[evt] = [TH1F(evt + "_" + chan, evt + "_" + chan, 100, -2,2) for chan in list_channel]
	# 	d_hist_simu[evt] = [TH1F(evt + "_simu_" + chan, evt + "_simu_" + chan, 100, -2,2) for chan in list_channel]
	# d_list_chan = {}

	# d_list_chan["heatonly"] = ["EIA", "EIB", "EIC", "EID"]
	# d_list_chan["FidGamma"] = ["EIA", "EIC"]
	# for evt in list_evt:
	# 	if "S1" in evt:
	# 		d_list_chan[evt] = ["EIC", "EID"]
	# 	if "S2" in evt:
	# 		d_list_chan[evt] = ["EIA", "EIB"]

	# for i in range(len(list_evt)):
	# 	list_canvas[i].Divide(2,2)
	# 	evt_type = list_evt[i]
	# 	arr = list_arr[i]
	# 	arr_simu = list_simu[i]
	# 	for k in range(len(d_list_chan[evt_type])):
	# 		chan = d_list_chan[evt_type][k]
	# 		list_canvas[i].cd(k+1)
	# 		for evt in arr[chan]:
	# 			d_hist_arr[evt_type][k].Fill(evt)
	# 		for evt in arr_simu[chan]:
	# 			d_hist_simu[evt_type][k].Fill(evt)
	# 		d_hist_arr[evt_type][k].Scale(1./d_hist_arr[evt_type][k].Integral())
	# 		d_hist_simu[evt_type][k].Scale(1./d_hist_simu[evt_type][k].Integral())
	# 		PyRPl.process_TH1(d_hist_arr[evt_type][k], X_title = chan, color = kRed)
	# 		PyRPl.process_TH1(d_hist_simu[evt_type][k], X_title = chan, color = kBlack)
	# 		d_hist_arr[evt_type][k].Draw()
	# 		d_hist_simu[evt_type][k].Draw("same")			


	#Check the channels where signal is expected
	#FidGamma
	list_canvas[1].Divide(2,2)
	evt_type = list_evt[0]
	arr = list_arr[1]
	arr_simu = list_simu[1]
	list_harr = [TH1F(evt_type + "_" + chan, evt_type + "_" + chan, 200, -2,15) for chan in ["EC1", "EC2", "EIB", "EID"]]
	list_hsimu = [TH1F(evt_type + "_simu_" + chan, evt_type + "_simu_" + chan, 200, -2,15) for chan in ["EC1", "EC2", "EIB", "EID"]]
	for index, k in enumerate([0,1,3,5]):
		chan = list_channel[k]
		list_canvas[1].cd(index+1)

		for evt in arr[chan]:
			list_harr[index].Fill(evt)
		for evt in arr_simu[chan]:
			list_hsimu[index].Fill(evt)
		list_harr[index].Scale(1./list_harr[index].Integral())
		list_hsimu[index].Scale(1./list_hsimu[index].Integral())
		PyRPl.process_TH1(list_harr[index], X_title = chan, color = kRed)
		PyRPl.process_TH1(list_hsimu[index], X_title = chan, color = kBlack)
		list_harr[index].Draw()
		list_hsimu[index].Draw("same")


	#Plot gamma for true data + simulated data
	# list_canvas[1].Divide(2,2)
	# evt_type = list_evt[1]
	# arr = arr_FidGamma_after_cuts
	# arr_simu = simu_FidGamma_after_cuts
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

	# #Do 2D Plots for FidGamma
	cc = TCanvas("cc", "cc")
	h2Darr = TH2F("h2Darr", "h2Darr", 1000, -2, 15, 1000, -2, 15)
	h2Dsimu = TH2F("h2Dsimu", "h2Dsimu", 1000, -2, 15, 1000, -2, 15)

	arr_all_EC = 0.5*(arr_all["EC1"]+arr_all["EC2"])
	simu_all_EC = 0.5*(simu_all["EC1"]+simu_all["EC2"])

	l_arr = np.where(np.logical_and(arr_all_EC>2 , arr_all_EC<3))
	l_simu = np.where(np.logical_and(simu_all_EC>2 , simu_all_EC<3))

	arr_d = arr_all[l_arr]["EC1"]-arr_all[l_arr]["EC2"]
	simu_d = simu_all[l_simu]["EC1"]-simu_all[l_simu]["EC2"]

	hdarr = TH1F("hdarr", "hdarr", 50 , -1, 0.5)
	hdsimu = TH1F("hdsimu", "hdsimu", 50 , -1, 0.5)

	for elem in arr_d:
		hdarr.Fill(elem)

	for elem in simu_d:
		hdsimu.Fill(elem)


	PyRPl.process_TH1(hdarr, X_title = "EC", Y_title = "EID", color = kRed)
	PyRPl.process_TH1(hdsimu, X_title = "EC", Y_title = "EID", color = kBlack)

	hdarr.Draw()
	hdsimu.Draw("same")
	raw_input()

	for i in range(arr_all["EC2"].shape[0]):
		EI = 0.414*arr_all["EIB"][i]+(1-0.414)*arr_all["EID"][i]
		EC = 0.574*arr_all["EC1"][i]+(1-0.574)*arr_all["EC2"][i]
		h2Darr.Fill(EC, EI)
	for i in range(simu_all["EC2"].shape[0]):
		EI = 0.414*simu_all["EIB"][i]+(1-0.414)*simu_all["EID"][i]
		EC = 0.574*simu_all["EC1"][i]+(1-0.574)*simu_all["EC2"][i]
		h2Dsimu.Fill(EC,EI)


	# for i in range(arr_all["EC2"].shape[0]):
	# 	h2Darr.Fill(arr_all["EC2"][i], arr_all["EC1"][i])
	# for i in range(simu_all["EC2"].shape[0]):
	# 	h2Dsimu.Fill(simu_all["EC2"][i], simu_all["EC1"][i])

	# for i in range(arr_FidGamma["EC2"].shape[0]):
	# 	EI = 0.414*arr_FidGamma["EIB"][i]+(1-0.414)*arr_FidGamma["EID"][i]
	# 	EC = 0.574*arr_FidGamma["EC1"][i]+(1-0.574)*arr_FidGamma["EC2"][i]
	# 	h2Darr.Fill(EC, EI)
	# for i in range(simu_FidGamma_XD["EC2"].shape[0]):
	# 	EI = 0.414*simu_FidGamma_XD["EIB"][i]+(1-0.414)*simu_FidGamma_XD["EID"][i]
	# 	EC = 0.574*simu_FidGamma_XD["EC1"][i]+(1-0.574)*simu_FidGamma_XD["EC2"][i]
	# 	h2Dsimu.Fill(EC,EI)

	# for i in range(arr_FidGamma["EIB"].shape[0]):
	# 	h2Darr.Fill(arr_FidGamma["EIB"][i], arr_FidGamma["EC1"][i])
	# for i in range(simu_FidGamma_XD["EIB"].shape[0]):
	# 	h2Dsimu.Fill(simu_FidGamma_XD["EIB"][i], simu_FidGamma_XD["EC1"][i])

	# for i in range(arr_FidGamma["EC2"].shape[0]):
	# 	h2Darr.Fill(arr_FidGamma["EC2"][i], arr_FidGamma["EID"][i])
	# for i in range(simu_FidGamma["EC2"].shape[0]):
	# 	h2Dsimu.Fill(simu_FidGamma["EC2"][i], simu_FidGamma["EID"][i])

	PyRPl.process_TH2(h2Darr, X_title = "EC", Y_title = "EID", color = kRed, marker_style = 20)
	PyRPl.process_TH2(h2Dsimu, X_title = "EC", Y_title = "EID", color = kBlack, marker_style = 20)
	h2Darr.Draw()
	h2Dsimu.Draw("same")

	raw_input()

bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
analysis_type = "ana_1.5_0_50"
mass = "5"
# get_BDT_array(bolo_name, analysis_type, mass)
check_BDT_simulations(bolo_name, analysis_type, mass)