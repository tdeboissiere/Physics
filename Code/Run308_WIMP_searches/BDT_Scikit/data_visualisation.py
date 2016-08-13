import numpy as np
import scipy as sc
import pandas as pd
import sklearn as sk
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier as GBC
from pandas import read_csv, DataFrame
import os,sys
import pickle
import data_preparation as dp 
import PyROOTPlots as PyRPl
from ROOT import *
import BDT_file_handler as BDT_fh
import script_utils as script_utils
import matplotlib.colors as mcolors
import xgboost as xgb
from sklearn import preprocessing
from sklearn.decomposition import PCA, KernelPCA, RandomizedPCA
from sklearn.lda import LDA

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {"red": [], "green": [], "blue": []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict["red"].append([item, r1, r2])
            cdict["green"].append([item, g1, g2])
            cdict["blue"].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap("CustomMap", cdict)


def plot_results(d_test, d_eval, d_event_dir, classifier_type, WIMP_mass, bolo_name, analysis_type, exposure, bin_X, min_X, max_X):

	"""
	Detail:
		Plot results

	Args:
		d_test (dict)                         = dict with test data
		d_eval (dict)                         = dict with eval data 
		d_event_dir (dict)                    = dict to get the proper directory of each event class 
		classifier_type (str)                 = type of classifier
		WIMP_mass (str)                       = WIMP mass
		bolo_name (str)                       = bolometer name 
		analysis_type (str)                   = type of analysis (which box cut)
		exposure (float)                      = exposure in days 
		bin_X, min_X, max_X (int float float) = TH1F parameters

	Returns:
		void

	Raises:
		void	 
	"""


	#Get scaling dict for data visualisation
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")

	#Load classifier
	pickle_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	clf_file = open(pickle_dir + classifier_type + "_mass_" + str(WIMP_mass) + ".pkl", 'rb')
	clf = pickle.load(clf_file)
	clf_file.close()

	#Get predictions on test sample
	d_pred = {}
	d_hist = {}
	d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
					 "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, "heatonly":kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}

	for event_type in d_test.keys():
		d_pred[event_type] = clf.predict_proba(d_test[event_type].iloc[:,:6].values)
		d_hist[event_type] = TH1F("h" + event_type, "h" + event_type, bin_X, min_X, max_X)
		PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type][:,1])
		PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
		if "WIMP" not in event_type:
			d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))

	#get predictions on data
	hdata = TH1F("hdata", "hdata", bin_X, min_X, max_X)
	PyRPl.fill_TH1(hdata, clf.predict_proba(d_eval["realdata"].iloc[:,:6].values)[:,1])
	d_hist["WIMP_mass_" + WIMP_mass].Scale(hdata.Integral()/d_hist["WIMP_mass_" + WIMP_mass].Integral())

	d_hist["S1Pb"].Add(d_hist["S2Pb"])
	d_hist["S1Beta"].Add(d_hist["S2Beta"])
	d_hist["FidGamma"].Add(d_hist["S1Gamma"])
	d_hist["FidGamma"].Add(d_hist["S2Gamma"])

	list_hist =[d_hist["S1Pb"], d_hist["S1Beta"], d_hist["FidGamma"], d_hist["heatonly"], d_hist["WIMP_mass_" + WIMP_mass]]
	hs=THStack("hs", "hs")
	for hist in list_hist:
		hs.Add(hist)

	cc = TCanvas("cc", "cc")
	h1=TH1F("h1","h1", bin_X, min_X, max_X)
	PyRPl.process_TH1(h1, X_title="BDT ouput", min_Y = 1E-1, max_Y = 20000)
	
	gPad.SetLogy()
	h1.Draw()
	hs.Draw("same")
	hdata.Draw("sameE1")

	leg = TLegend(0.14,0.50,0.33,0.87)
	leg.AddEntry(d_hist["S1Pb"].GetName(),"Lead" ,"f")
	leg.AddEntry(d_hist["S1Beta"].GetName(),"Beta", "f")
	leg.AddEntry(d_hist["FidGamma"].GetName(),"Gamma", "f")
	leg.AddEntry(d_hist["heatonly"].GetName(),"Heat-only", "f")
	leg.AddEntry(d_hist["WIMP_mass_" + WIMP_mass].GetName(),"WIMP " + WIMP_mass+ " GeV","f")
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)
	leg.Draw("same")
	
	print d_hist["WIMP_mass_" + WIMP_mass].Integral(d_hist["WIMP_mass_" + WIMP_mass].FindBin(0.9), 1)/d_hist["WIMP_mass_" + WIMP_mass].Integral()

	raw_input()

	for key in d_hist.keys(): del d_hist[key]
	del d_hist
	del h1


def plot_control(d_test, d_eval, d_event_dir, classifier_type, WIMP_mass, bolo_name, analysis_type, exposure, bin_X, min_X, max_X):

	"""
	Detail:
		Plot results

	Args:
		d_test (dict)                         = dict with test data
		d_eval (dict)                         = dict with eval data 
		d_event_dir (dict)                    = dict to get the proper directory of each event class 
		classifier_type (str)                 = type of classifier
		WIMP_mass (str)                       = WIMP mass
		bolo_name (str)                       = bolometer name 
		analysis_type (str)                   = type of analysis (which box cut)
		exposure (float)                      = exposure in days 
		bin_X, min_X, max_X (int float float) = TH1F parameters

	Returns:
		void

	Raises:
		void	 
	"""

	#Load classifier
	pickle_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	clf_file = open(pickle_dir + classifier_type + "_mass_" + str(WIMP_mass) + ".pkl", 'rb')
	clf = pickle.load(clf_file)
	clf_file.close()

	#Get predictions
	d_pred = {}
	for event_type in d_test.keys():
		d_pred[event_type] = clf.predict_proba(d_test[event_type].iloc[:,:6].values)
	d_pred["realdata"] = clf.predict_proba(d_eval["realdata"].iloc[:,:6].values)
	
	#Get color map
	c = mcolors.ColorConverter().to_rgb
	rvb = make_colormap([c("red"), c("red"), 0.33, c("red"), c("green"), 0.86, c("green")])
	
	#Build list of test data of interest
	list_event_type = ["FidGamma", "WIMP_mass_" + WIMP_mass]

	#Output dir 
	fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/")

	#Compute new columns for data frames
	d_eval["realdata"]["EC"] = 0.5*(d_eval["realdata"]["EC1"] + d_eval["realdata"]["EC2"])
	d_eval["realdata"]["EFID"] = 0.5*(d_eval["realdata"]["EIB"] + d_eval["realdata"]["EID"])
	d_eval["realdata"]["MVA"] = d_pred["realdata"][:,1]
	for event_type in list_event_type:
		d_test[event_type]["EC"] = 0.5*(d_test[event_type]["EC1"] + d_test[event_type]["EC2"])
		d_test[event_type]["EFID"] = 0.5*(d_test[event_type]["EIB"] + d_test[event_type]["EID"])
		d_test[event_type]["MVA"] = d_pred[event_type][:,1]


	#Plots for real data
	plt.figure()
	plt.scatter(d_test["WIMP_mass_" + WIMP_mass]["EC"].values, d_test["WIMP_mass_" + WIMP_mass]["EFID"].values, color ="0.6", s=1 )
	plt.scatter(d_eval["realdata"]["EC"], d_eval["realdata"]["EFID"], s =20, c=d_eval["realdata"]["MVA"].values, cmap=rvb, vmin = -0.5, vmax = 1)

	# d_test["WIMP_mass_" + WIMP_mass].plot(kind='scatter', x='EC', y='EFID', color="0.6", s = 1)
	# ax=d_eval["realdata"].plot(kind='scatter', x='EC', y='EFID', c=d_eval["realdata"]["MVA"], cmap=rvb, vmin = 0, vmax = 1) 
	
	cbar = plt.colorbar()
	cbar.set_label("MVA output", labelpad = 15, fontsize= 18)
	plt.xlabel("Heat (keV)", fontsize = 20)
	plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
	plt.ylim([0,5])
	plt.xlim([0.5,5])
	plt.grid(True)
	plt.savefig(fig_dir + bolo_name + "_real_data_WIMP_mass_" + WIMP_mass  + "_" + classifier_type + ".png")
	plt.close("all")


	for event_type in list_event_type:
		plt.figure()

		plt.scatter(d_test["WIMP_mass_" + WIMP_mass]["EC"].values, d_test["WIMP_mass_" + WIMP_mass]["EFID"].values, color ="0.6", s=1 )
		plt.scatter(d_test[event_type]["EC"], d_test[event_type]["EFID"], s =20, c=d_test[event_type]["MVA"].values, cmap=rvb, vmin = -0.5, vmax = 1)

		# d_test["WIMP_mass_" + WIMP_mass].plot(kind='scatter', x='EC', y='EFID', color="0.6", s = 1)
		# d_test[event_type].plot(kind='scatter', x='EC', y='EFID', c=d_test[event_type]["MVA"], cmap=rvb, vmin = 0, vmax = 1)

		cbar = plt.colorbar()
		cbar.set_label("BDT output", labelpad = 15, fontsize= 18)
		plt.xlabel("Heat (keV)", fontsize = 20)
		plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
		plt.ylim([0,5])
		plt.xlim([0.5,5])
		plt.grid(True)
		plt.savefig(fig_dir + bolo_name + "_" + event_type + "_WIMP_mass_" + WIMP_mass  + "_" + classifier_type + ".png")
		plt.close("all")


def plot_results_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):
	"""
	Detail:
		Plot results

	Args:
		d_test (dict)                         = dict with test data
		d_eval (dict)                         = dict with eval data 
		d_event_dir (dict)                    = dict to get the proper directory of each event class 
		WIMP_mass (str)                       = WIMP mass
		bolo_name (str)                       = bolometer name 
		analysis_type (str)                   = type of analysis (which box cut)
		MVA_tag (str)                         = indicates which scaling file to use
		exposure (float)                      = exposure in days 
		bin_X, min_X, max_X (int float float) = TH1F parameters

	Returns:
		void

	Raises:
		void	 
	"""

	try:
		kwargs["weight_dir"]
	except KeyError:
		sys.exit()

	#Get scaling dict to set the weights
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

	# #Load PCA
	# pickle_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	# pca_file = open(pickle_dir + "pca_classifier_mass_" + str(WIMP_mass) + ".pkl", 'rb')
	# pca = pickle.load(pca_file)
	# pca_file.close()

	key_heat = ""
	for key in d_test.keys():
		if "heat" in key:
			key_heat = key

	# Get classifier
	model_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
	modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
	if kwargs.has_key("classifier_name"):
		modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
	bst = xgb.Booster({'nthread':16}, model_file = modelfile)

	#Get predictions on test sample
	d_pred = {}
	d_hist = {}
	d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
					 "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, key_heat: kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}

	for event_type in d_test.keys():
		d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values) )
		d_hist[event_type] = TH1F("h" + event_type+ WIMP_mass, "h" + event_type+ WIMP_mass, bin_X, min_X, max_X)
		PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
		PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
		if "WIMP" not in event_type:
			d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
			print "Event type:", event_type, "\tExpected #:", float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure

	#get predictions on data
	hdata = TH1F("hdata" + WIMP_mass, "hdata" + WIMP_mass, bin_X, min_X, max_X)
	PyRPl.fill_TH1(hdata, bst.predict( xgb.DMatrix(d_eval["realdata"].iloc[:,:].values) ) )
	d_hist["WIMP_mass_" + WIMP_mass].Scale(hdata.Integral()/d_hist["WIMP_mass_" + WIMP_mass].Integral())

	d_hist["S1Pb"].Add(d_hist["S2Pb"])
	d_hist["S1Beta"].Add(d_hist["S2Beta"])
	d_hist["FidGamma"].Add(d_hist["S1Gamma"])
	d_hist["FidGamma"].Add(d_hist["S2Gamma"])

	list_hist =[d_hist["S1Pb"], d_hist["S1Beta"], d_hist["FidGamma"], d_hist[key_heat], d_hist["WIMP_mass_" + WIMP_mass]]
	hs=THStack("hs" + WIMP_mass, "hs" + WIMP_mass)
	for hist in list_hist:
		hs.Add(hist)

	list_hist =[d_hist["S1Pb"], d_hist["S1Beta"], d_hist["FidGamma"], d_hist[key_heat], d_hist["WIMP_mass_" + WIMP_mass]]
	hsum_bckg=TH1F("hsum_bckg","hsum_bckg", bin_X, min_X, max_X)
	for i in range(1, bin_X+1): 
		hsum_bckg.SetBinContent(i, sum([h.GetBinContent(i) for h in list_hist[:-1]]))

	# print "Chi2: ", hdata.Chi2Test(hsum_bckg, "P")
	del hsum_bckg

	cc = TCanvas("cc", "cc")
	h1=TH1F("h1" + WIMP_mass,"h1" + WIMP_mass, bin_X, min_X, max_X)
	PyRPl.process_TH1(h1, X_title="BDT ouput", min_Y = 1E-1, max_Y = 20000)
	
	gPad.SetLogy()
	h1.Draw()
	hs.Draw("same")
	hdata.Draw("sameE1")

	leg = TLegend(0.14,0.50,0.33,0.87)
	leg.AddEntry(d_hist["S1Pb"].GetName(),"Lead" ,"f")
	leg.AddEntry(d_hist["S1Beta"].GetName(),"Beta", "f")
	leg.AddEntry(d_hist["FidGamma"].GetName(),"Gamma", "f")
	leg.AddEntry(d_hist[key_heat].GetName(),"Heat-only", "f")
	leg.AddEntry(d_hist["WIMP_mass_" + WIMP_mass].GetName(),"WIMP " + WIMP_mass+ " GeV","f")
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)
	leg.Draw("same")
	
	raw_input()
	fig_dir = script_utils.create_directory("./Figures/"+ bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
	cc.Print(fig_dir + bolo_name + "_BDT_mass_" + str(WIMP_mass) + ".eps")

def plot_control_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, list_variables, **kwargs):
	"""
	Detail:
		Plot results control

	Args:
		d_test (dict)                         = dict with test data
		d_eval (dict)                         = dict with eval data 
		d_event_dir (dict)                    = dict to get the proper directory of each event class 
		WIMP_mass (str)                       = WIMP mass
		bolo_name (str)                       = bolometer name 
		analysis_type (str)                   = type of analysis (which box cut)
		MVA_tag (str)                         = indicates which scaling file to use
		exposure (float)                      = exposure in days 
		bin_X, min_X, max_X (int float float) = TH1F parameters

	Returns:
		void

	Raises:
		void	 
	"""


	try:
		kwargs["weight_dir"]
	except KeyError:
		sys.exit()

	#Get scaling dict to set the weights
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

	key_heat = ""
	for key in d_test.keys():
		if "heat" in key:
			key_heat = key

	##################
	# Temporary to check effect of EC versus 0.5(EC1+EC2)
	# Effect small
	####################
	# d_eval = {}
	# data_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Eval_data/" + bolo_name + "/" + analysis_type + "/")
	# d_eval["realdata"] = pd.read_csv(data_dir + bolo_name + "_" + analysis_type +  "_fond.csv", usecols = ["EC1","EC2","EIA","EIB","EIC","EID","EC","EFID", "HR"])
	# temp_eval_EC = d_eval["realdata"]["EC"]
	# d_eval["realdata"] = d_eval["realdata"][list_variables]

	# Get classifier
	model_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
	modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
	if kwargs.has_key("classifier_name"):
		modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
	bst = xgb.Booster({'nthread':16}, model_file = modelfile)

	#Get predictions
	d_pred = {}
	for event_type in d_test.keys():
		d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values))
	d_pred["realdata"] = bst.predict( xgb.DMatrix(d_eval["realdata"].iloc[:,:].values) ) 
	
	#Get color map
	c = mcolors.ColorConverter().to_rgb
	rvb = make_colormap([c("red"), c("red"), 0.33, c("red"), c("green"), 0.86, c("green")])
	
	#Build list of test data of interest
	list_event_type = ["FidGamma", "heatonly",  "WIMP_mass_" + WIMP_mass]
	list_event_type = ["WIMP_mass_" + WIMP_mass]
	# list_event_type = ["heatonly"]

	#Output dir 
	fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/" + kwargs["weight_dir"] + "/")

	#Compute new columns for data frames
	d_eval["realdata"]["EC"] = 0.5*(d_eval["realdata"]["EC1"] + d_eval["realdata"]["EC2"]) #  temp_eval_EC 
	if "EIB" in list(d_eval["realdata"].columns.values) and "EID" in list(d_eval["realdata"].columns.values):
		d_eval["realdata"]["EFID"] = 0.5*(d_eval["realdata"]["EIB"] + d_eval["realdata"]["EID"])
	d_eval["realdata"]["MVA"] = d_pred["realdata"]

	# plt.hist(d_eval["realdata"]["MVA"], bins=100)
	# plt.show()
	# raw_input()

	for event_type in list_event_type:
		d_test[event_type]["EC"] = 0.5*(d_test[event_type]["EC1"] + d_test[event_type]["EC2"])
		if "EIB" in list(d_test[event_type].columns.values) and "EID" in list(d_test[event_type].columns.values):
			d_test[event_type]["EFID"] = 0.5*(d_test[event_type]["EIB"] + d_test[event_type]["EID"])
		d_test[event_type]["MVA"] = d_pred[event_type]

	# l = np.where(d_eval["realdata"]["MVA"]>0)
	# d_eval["realdata"] = d_eval["realdata"].iloc[l]
	# l = np.where(np.logical_and(d_eval["realdata"]["EIA"]<1, d_eval["realdata"]["EIC"]<1))
	# d_eval["realdata"] = d_eval["realdata"].iloc[l]



	# l = np.where(d_test["WIMP_mass_" + WIMP_mass]["MVA"]>3)
	# d_test["WIMP_mass_" + WIMP_mass] = d_test["WIMP_mass_" + WIMP_mass].iloc[l]


	#Plots for real data
	plt.figure()
	plt.scatter(d_test["WIMP_mass_" + WIMP_mass]["EC"].values, d_test["WIMP_mass_" + WIMP_mass]["EFID"].values, color ="0.6", s=1 )
	plt.scatter(d_eval["realdata"]["EC"], d_eval["realdata"]["EFID"], s =20, c=d_eval["realdata"]["MVA"].values, cmap=rvb, vmin = -10, vmax = 10)

	cbar = plt.colorbar()
	cbar.set_label("BDT output", labelpad = 15, fontsize= 18)
	plt.xlabel("Heat (keV)", fontsize = 20)
	plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
	plt.ylim([0,5])
	plt.xlim([0.5,5])
	plt.grid(True)
	plt.savefig(fig_dir + bolo_name + "_real_data_WIMP_mass_" + WIMP_mass +  ".png")
	plt.close("all")


	# for event_type in list_event_type:
	# 	plt.figure()

	# 	# plt.scatter(d_test["WIMP_mass_" + WIMP_mass]["EC"].values, d_test["WIMP_mass_" + WIMP_mass]["EFID"].values, color ="0.6", s=1 )
	# 	plt.scatter(d_test[event_type]["EC"], d_test[event_type]["EFID"], s =20, c=d_test[event_type]["MVA"].values, cmap=rvb, vmin = -13, vmax = 6)

	# 	cbar = plt.colorbar()
	# 	cbar.set_label("BDT output", labelpad = 15, fontsize= 18)
	# 	plt.xlabel("Heat (keV)", fontsize = 20)
	# 	plt.ylabel("Fiducial Ionisation (keV)", fontsize = 20)
	# 	plt.ylim([0,5])
	# 	plt.xlim([0.5,5])
	# 	plt.grid(True)
	# 	plt.savefig(fig_dir + bolo_name + "_" + event_type + "_WIMP_mass_" + WIMP_mass  + ".png")
	# 	plt.close("all")



def plot_PCA_stuff(d_test, d_eval, d_event_dir, classifier_type, WIMP_mass, bolo_name, analysis_type, exposure, bin_X, min_X, max_X, pca_index):

	"""
	Detail:
		Plot PCA results

	Args:
		d_test (dict)                         = dict with test data
		d_eval (dict)                         = dict with eval data 
		d_event_dir (dict)                    = dict to get the proper directory of each event class 
		classifier_type (str)                 = type of classifier
		WIMP_mass (str)                       = WIMP mass
		bolo_name (str)                       = bolometer name 
		analysis_type (str)                   = type of analysis (which box cut)
		exposure (float)                      = exposure in days 
		bin_X, min_X, max_X (int float float) = TH1F parameters

	Returns:
		void

	Raises:
		void	 
	"""

	#Get scaling dict for data visualisation
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")

	#Load PCA
	pickle_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	pca_file = open(pickle_dir + "pca_classifier_mass_" + str(WIMP_mass) + ".pkl", 'rb')
	pca = pickle.load(pca_file)
	pca_file.close()

	# Get classifier
	model_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
	bst = xgb.Booster({'nthread':16}, model_file = modelfile)

	#Get predictions on test sample
	d_hist = {}
	d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
					 "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, "heatonly":kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}


	print pca.transform(d_test["FidGamma"].iloc[:,:-2].values).shape
	print pca.transform(d_eval["realdata"].iloc[:,:].values).shape

	for event_type in ["S1Gamma", "S2Gamma", "FidGamma", "S1Beta", "S2Beta", "S1Pb", "S2Pb", "heatonly"]:
		d_hist[event_type] = TH1F("h" + event_type, "h" + event_type, bin_X, min_X, max_X)
		PyRPl.fill_TH1(d_hist[event_type], bst.predict( xgb.DMatrix(pca.transform(d_test[event_type].iloc[:,:-2].values) )))
		PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
		d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))

	#get predictions on data
	hdata = TH1F("hdata", "hdata", bin_X, min_X, max_X)
	PyRPl.fill_TH1(hdata, bst.predict( xgb.DMatrix(pca.transform(d_eval["realdata"].iloc[:,:].values))))

	d_hist["S1Pb"].Add(d_hist["S2Pb"])
	d_hist["S1Beta"].Add(d_hist["S2Beta"])
	d_hist["FidGamma"].Add(d_hist["S1Gamma"])
	d_hist["FidGamma"].Add(d_hist["S2Gamma"])

	list_hist =[d_hist["S1Pb"], d_hist["S1Beta"], d_hist["FidGamma"], d_hist["heatonly"]]
	hs=THStack("hs", "hs")
	for hist in list_hist:
		hs.Add(hist)

	cc = TCanvas("cc", "cc")
	h1=TH1F("h1","h1", bin_X, min_X, max_X)
	PyRPl.process_TH1(h1, X_title="PCA var", min_Y = 1E-1, max_Y = 20000)
	
	gPad.SetLogy()
	h1.Draw()
	hs.Draw("same")
	hdata.Draw("sameE1")

	leg = TLegend(0.14,0.50,0.33,0.87)
	leg.AddEntry(d_hist["S1Pb"].GetName(),"Lead" ,"f")
	leg.AddEntry(d_hist["S1Beta"].GetName(),"Beta", "f")
	leg.AddEntry(d_hist["FidGamma"].GetName(),"Gamma", "f")
	leg.AddEntry(d_hist["heatonly"].GetName(),"Heat-only", "f")
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)
	leg.Draw("same")
	

	raw_input()