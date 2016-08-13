from ROOT import *
import numpy as np
import script_utils as script_utils
import os
import scipy.stats
import PyROOTPlots as PyRPl
from math import sqrt as sqrt

def compute_KS_on_surface(bolo_name):

	"""Compute KS test FIDXXX versus FID808
	
	Detail:
		2 tests: ROOT and scipy

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""

	#Get .txt file data for "true" KS test
	pop_path = "../Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"

	arr_S1Beta = np.loadtxt(pop_path + "/" + bolo_name + "_S1Beta.txt", delimiter=",").astype(float)
	arr_S2Beta = np.loadtxt(pop_path + "/" + bolo_name + "_S2Beta.txt", delimiter=",").astype(float)

	arr_S1Pb = np.loadtxt(pop_path + "/" + bolo_name + "_S1Pb.txt", delimiter=",").astype(float)
	arr_S2Pb = np.loadtxt(pop_path + "/" + bolo_name + "_S2Pb.txt", delimiter=",").astype(float)

	#Get reference arrays of FID808

	print pop_path

	arr_S1Beta_ref = np.loadtxt("../Text_files/S1Beta_heatremoved.txt", delimiter=",").astype(float)
	arr_S2Beta_ref = np.loadtxt("../Text_files/S2Beta_heatremoved.txt", delimiter=",").astype(float)

	arr_S1Beta_ref = 0.5*(arr_S1Beta_ref[:,2] + arr_S1Beta_ref[:,3])
	arr_S2Beta_ref = 0.5*(arr_S2Beta_ref[:,2] + arr_S2Beta_ref[:,3])

	arr_S1Pb_ref = np.loadtxt("../Text_files/S1Pb_heatremoved.txt", delimiter=",").astype(float)
	arr_S2Pb_ref = np.loadtxt("../Text_files/S2Pb_heatremoved.txt", delimiter=",").astype(float)

	arr_S1Pb_ref = 0.5*(arr_S1Pb_ref[:,2] + arr_S1Pb_ref[:,3])
	arr_S2Pb_ref = 0.5*(arr_S2Pb_ref[:,2] + arr_S2Pb_ref[:,3])

	# print scipy.stats.ks_2samp(arr_S1Beta_ref, arr_S1Beta_ref)	
	# raw_input()

	np.savetxt("../Text_files/S1Beta_heatremoved_onlyheat.txt", arr_S1Beta_ref)
	np.savetxt("../Text_files/S1Beta_heatremoved_onlyheat.txt", arr_S2Beta_ref)
	np.savetxt("../Text_files/S1Pb_heatremoved_onlyheat.txt", arr_S1Pb_ref)
	np.savetxt("../Text_files/S1Pb_heatremoved_onlyheat.txt", arr_S2Pb_ref)

	binBeta, XBeta, YBeta = 50, 2, 20
	binPb, XPb, YPb = 200, 15, 40




	# #Generate data from the model built with KDE
	# fil = TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name + "_Beta_from_bolo_spectrum_extrapol.root")
	# f = fil.Get("S1Beta_extra")

	# list_ev = [f.GetRandom(2,20) for i in range(50000)]
	# list_ev2 = [f.GetRandom(2,20) for i in range(50000)]
	# arr_S1Beta_simu = np.array(list_ev).astype(float)
	# arr_S1Beta_simu2 = np.array(list_ev2).astype(float)
	# hS1Beta_simu = TH1F("hS1Beta_simu", "hS1Beta_simu", binBeta, XBeta, YBeta)
	# for elem in arr_S1Beta_simu:
	# 	hS1Beta_simu.Fill(elem)
	# hS1Beta_simu.SetLineColor(kBlue)
	# hS1Beta_simu.Scale(hS1Beta.Integral()/float(hS1Beta_simu.Integral()))


	# print hS1Beta.KolmogorovTest(hS1Beta_ref)
	# print hS2Beta.KolmogorovTest(hS2Beta_ref)
	# print hS1Pb.KolmogorovTest(hS1Pb_ref)
	# print hS2Pb.KolmogorovTest(hS2Pb_ref)


	# print scipy.stats.ks_2samp(arr_S1Beta, arr_S1Beta_ref)
	# print scipy.stats.ks_2samp(arr_S2Beta, arr_S2Beta_ref)

	# print scipy.stats.ks_2samp(arr_S1Pb, arr_S1Pb_ref)
	# print scipy.stats.ks_2samp(arr_S2Pb, arr_S2Pb_ref)

	# #Pick events in good "boxes"
	arr_S1Beta = np.array([elem for elem in arr_S1Beta if 4<elem<20])
	arr_S2Beta = np.array([elem for elem in arr_S2Beta if 4<elem<20])
	arr_S1Pb = np.array([elem for elem in arr_S1Pb if 15<elem<40])
	arr_S2Pb = np.array([elem for elem in arr_S2Pb if 15<elem<40])

	arr_S1Beta_ref = np.array([elem for elem in arr_S1Beta_ref if 4<elem<20])
	arr_S2Beta_ref = np.array([elem for elem in arr_S2Beta_ref if 4<elem<20])
	arr_S1Pb_ref = np.array([elem for elem in arr_S1Pb_ref if 15<elem<40])
	arr_S2Pb_ref = np.array([elem for elem in arr_S2Pb_ref if 15<elem<40])

	list_S1Beta_index = np.random.randint(0,arr_S1Beta_ref.shape[0], arr_S1Beta.shape[0])
	list_S2Beta_index = np.random.randint(0,arr_S2Beta_ref.shape[0], arr_S2Beta.shape[0])
	list_S1Pb_index = np.random.randint(0,arr_S1Pb_ref.shape[0], arr_S1Pb.shape[0])
	list_S2Pb_index = np.random.randint(0,arr_S2Pb_ref.shape[0], arr_S2Pb.shape[0])

	arr_S1Beta_ref = arr_S1Beta_ref[list_S1Beta_index]
	arr_S2Beta_ref = arr_S2Beta_ref[list_S2Beta_index]

	arr_S1Pb_ref = arr_S1Pb_ref[list_S1Pb_index]
	arr_S2Pb_ref = arr_S2Pb_ref[list_S2Pb_index]

	#Using root TH1
	#Between 2 and 30 for Beta
	#Between 15 and 40 for Pb
	hS1Beta = TH1F("hS1Beta", "hS1Beta", binBeta, XBeta, YBeta)
	hS2Beta = TH1F("hS2Beta", "hS2Beta", binBeta, XBeta, YBeta)
	
	hS1Pb = TH1F("hS1Pb", "hS1Pb", binPb, XPb, YPb)
	hS2Pb = TH1F("hS2Pb", "hS2Pb", binPb, XPb, YPb)

	for elem in arr_S1Beta:
		hS1Beta.Fill(elem)
	for elem in arr_S2Beta:
		hS2Beta.Fill(elem)

	for elem in arr_S1Pb:
		hS1Pb.Fill(elem)
	for elem in arr_S2Pb:
		hS2Pb.Fill(elem)

	hS1Beta_ref = TH1F("hS1Beta_ref", "hS1Beta_ref", binBeta, XBeta, YBeta)
	hS2Beta_ref = TH1F("hS2Beta_ref", "hS2Beta_ref", binBeta, XBeta, YBeta)
	
	hS1Pb_ref = TH1F("hS1Pb_ref", "hS1Pb_ref", binPb, XPb, YPb)
	hS2Pb_ref = TH1F("hS2Pb_ref", "hS2Pb_ref", binPb, XPb, YPb)

	for elem in arr_S1Beta_ref:
		hS1Beta_ref.Fill(elem)
	for elem in arr_S2Beta_ref:
		hS2Beta_ref.Fill(elem)

	for elem in arr_S1Pb_ref:
		hS1Pb_ref.Fill(elem)
	for elem in arr_S2Pb_ref:
		hS2Pb_ref.Fill(elem)

	list_hist_ref = [hS1Beta_ref, hS2Beta_ref, hS1Pb_ref, hS2Pb_ref]
	list_hist = [hS1Beta, hS2Beta, hS1Pb, hS2Pb]
	for i, elem in enumerate(list_hist_ref):
		PyRPl.process_TH1(elem, color = kRed)
		elem.Scale(float(list_hist[i].Integral())/float(elem.Integral()))


	# arr_S1Beta_simu = np.array([elem for elem in arr_S1Beta_simu if 5<elem<20])
	# arr_S1Beta_simu2 = np.array([elem for elem in arr_S1Beta_simu2 if 5<elem<20])

	# print scipy.stats.ks_2samp(arr_S1Beta_simu, arr_S1Beta_simu2),"   ", 1.36*sqrt(len(arr_S1Beta_simu) + len(arr_S1Beta_simu2))/sqrt(len(arr_S1Beta_simu2) * len(arr_S1Beta_simu))
	# raw_input()

	# print
	# print

	# #1.36 for 95 CL see Wikipedia
	# print scipy.stats.ks_2samp(arr_S1Beta, arr_S1Beta_ref),"   ", 1.63*sqrt(len(arr_S1Beta) + len(arr_S1Beta_ref))/sqrt(len(arr_S1Beta) * len(arr_S1Beta_ref))
	# print scipy.stats.ks_2samp(arr_S2Beta, arr_S2Beta_ref),"   ", 1.63*sqrt(len(arr_S2Beta) + len(arr_S2Beta_ref))/sqrt(len(arr_S2Beta) * len(arr_S2Beta_ref))

	# print scipy.stats.ks_2samp(arr_S1Pb, arr_S1Pb_ref),"   ", 1.63*sqrt(len(arr_S1Pb) + len(arr_S1Pb_ref))/sqrt(len(arr_S1Pb) * len(arr_S1Pb_ref))
	# print scipy.stats.ks_2samp(arr_S2Pb, arr_S2Pb_ref),"   ", 1.63*sqrt(len(arr_S2Pb) + len(arr_S2Pb_ref))/sqrt(len(arr_S2Pb) * len(arr_S2Pb_ref))

	# arr_S1Beta = np.concatenate((arr_S1Beta, arr_S2Beta))
	# arr_S1Beta_ref = np.concatenate((arr_S1Beta_ref, arr_S2Beta_ref))
	# arr_S1Pb = np.concatenate((arr_S1Pb, arr_S2Pb))
	# arr_S1Pb_ref = np.concatenate((arr_S1Pb_ref, arr_S2Pb_ref))

	# raw_input()

	return scipy.stats.ks_2samp(arr_S1Beta, arr_S1Beta_ref)[1], scipy.stats.ks_2samp(arr_S2Beta, arr_S2Beta_ref)[1], scipy.stats.ks_2samp(arr_S1Pb, arr_S1Pb_ref)[1], scipy.stats.ks_2samp(arr_S2Pb, arr_S2Pb_ref)[1]

	# print scipy.stats.ks_2samp(arr_S1Beta, arr_S1Beta_simu),"   ", 1.36*sqrt(len(arr_S1Beta) + len(arr_S1Beta_simu))/sqrt(len(arr_S1Beta) * len(arr_S1Beta_simu))

	# print hS1Beta.KolmogorovTest(hS1Beta_simu)
	# print len(list(arr_S1Beta))
	# print len(list(arr_S1Beta_simu))

	# print len(arr_S1Beta), len(arr_S1Beta_ref), sqrt(len(arr_S1Beta) + len(arr_S1Beta_ref))/sqrt(len(arr_S1Beta) * len(arr_S1Beta_ref))
	# print len(arr_S2Beta), len(arr_S2Beta_ref), sqrt(len(arr_S2Beta) + len(arr_S2Beta_ref))/sqrt(len(arr_S2Beta) * len(arr_S2Beta_ref))
	# print len(arr_S1Pb), len(arr_S1Pb_ref), sqrt(len(arr_S1Pb) + len(arr_S1Pb_ref))/sqrt(len(arr_S1Pb) * len(arr_S1Pb_ref))
	# print len(arr_S2Pb), len(arr_S2Pb_ref), sqrt(len(arr_S2Pb) + len(arr_S2Pb_ref))/sqrt(len(arr_S2Pb) * len(arr_S2Pb_ref))

	# print
	# print

	# # hS1Beta.Draw()
	# # hS1Beta_simu.Draw("same")

	# cc = TCanvas("cc", "cc")
	# cc.Divide(2,2)
	# cc.cd(1)
	# hS1Beta.Draw()
	# hS1Beta_ref.Draw("same")

	# cc.cd(2)
	# hS2Beta.Draw()
	# hS2Beta_ref.Draw("same")

	# cc.cd(3)
	# hS1Pb.Draw()
	# hS1Pb_ref.Draw("same")

	# cc.cd(4)
	# hS2Pb.Draw()
	# hS2Pb_ref.Draw("same")



	# import statsmodels.api as sm # recommended import according to the docs
	# import matplotlib.pylab as plt

	
	# S1Beta_cdf = sm.distributions.ECDF(arr_S1Beta)
	# S1Beta_ref_cdf = sm.distributions.ECDF(arr_S1Beta_ref)

	# x_S1Beta = np.linspace(min(arr_S1Beta), max(arr_S1Beta))
	# x_S1Beta_ref = np.linspace(min(arr_S1Beta_ref), max(arr_S1Beta_ref))
	# y_S1Beta = S1Beta_cdf(x_S1Beta)
	# y_S1Beta_ref = S1Beta_ref_cdf(x_S1Beta_ref)

	# S2Beta_cdf = sm.distributions.ECDF(arr_S2Beta)
	# S2Beta_ref_cdf = sm.distributions.ECDF(arr_S2Beta_ref)

	# x_S2Beta = np.linspace(min(arr_S2Beta), max(arr_S2Beta))
	# x_S2Beta_ref = np.linspace(min(arr_S2Beta_ref), max(arr_S2Beta_ref))
	# y_S2Beta = S2Beta_cdf(x_S2Beta)
	# y_S2Beta_ref = S2Beta_ref_cdf(x_S2Beta_ref)

	# S1Pb_cdf = sm.distributions.ECDF(arr_S1Pb)
	# S1Pb_ref_cdf = sm.distributions.ECDF(arr_S1Pb_ref)

	# x_S1Pb = np.linspace(min(arr_S1Pb), max(arr_S1Pb))
	# x_S1Pb_ref = np.linspace(min(arr_S1Pb_ref), max(arr_S1Pb_ref))
	# y_S1Pb = S1Pb_cdf(x_S1Pb)
	# y_S1Pb_ref = S1Pb_ref_cdf(x_S1Pb_ref)

	# S2Pb_cdf = sm.distributions.ECDF(arr_S2Pb)
	# S2Pb_ref_cdf = sm.distributions.ECDF(arr_S2Pb_ref)

	# x_S2Pb = np.linspace(min(arr_S2Pb), max(arr_S2Pb))
	# x_S2Pb_ref = np.linspace(min(arr_S2Pb_ref), max(arr_S2Pb_ref))
	# y_S2Pb = S2Pb_cdf(x_S2Pb)
	# y_S2Pb_ref = S2Pb_ref_cdf(x_S2Pb_ref)


	# plt.subplot(221)
	# # plt.set_title("Sharing x per column, y per row")
	# plt.step(x_S1Beta, y_S1Beta, "k", label = "FID837 S1Beta")
	# plt.step(x_S1Beta_ref, y_S1Beta_ref, "r", label = "FID808 S1Beta")
	# plt.legend(loc="upper left", prop={"size":10})

	# plt.subplot(222)
	# plt.step(x_S2Beta, y_S2Beta, "k", label = "FID837 S2Beta")
	# plt.step(x_S2Beta_ref, y_S2Beta_ref, "r", label = "FID808 S2Beta")
	# plt.legend(loc="upper left", prop={"size":10})

	# plt.subplot(223)	
	# plt.step(x_S1Pb, y_S1Pb, "k", label = "FID837 S1Pb")
	# plt.step(x_S1Pb_ref, y_S1Pb_ref, "r", label = "FID808 S1Pb")
	# plt.legend(loc="upper left", prop={"size":10})

	# plt.subplot(224)	
	# plt.step(x_S2Pb, y_S2Pb, "k", label = "FID837 S2Pb")
	# plt.step(x_S2Pb_ref, y_S2Pb_ref, "r", label = "FID808 S2Pb")
	# plt.legend(loc="upper left", prop={"size":10})
	# plt.show()

	# raw_input()

bolo_list=["FID837"]
lS1Beta = []
lS2Beta = []
lS1Pb = []
lS2Pb = []
# bolo_name = "FID838"
for i in range(300):
	for bolo_name in bolo_list:
		print i
		a,b,c,d = compute_KS_on_surface(bolo_name)
		lS1Beta.append(a)
		lS2Beta.append(b)
		lS1Pb.append(c)
		lS2Pb.append(d)


# print np.mean(np.array(lS1Beta))
# print np.mean(np.array(lS2Beta))
print np.mean(np.array(lS1Pb))
print np.mean(np.array(lS2Pb))