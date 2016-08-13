import numpy as np
import scipy as sc
import pandas as pd
import sklearn as sk
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier as GBC
from pandas import read_csv, DataFrame
import os,sys
sys.path.append("../")
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



def plot_BDT_at_limit_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):
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

	#Load number of WIMP events 
	file_lim = np.loadtxt("./Text_files/" + bolo_name + "_" + analysis_type + "_NWIMP_at_limit_MCMC.txt", delimiter = ",")
	d_lim_norm = {}
	for i in range(file_lim.shape[0]):
		d_lim_norm[str(int(file_lim[i][0]))] = np.mean(file_lim[i][1:])

	print d_lim_norm

	key_heat = ""
	for key in d_test.keys():
		if "heat" in key:
			key_heat = key

	# Get classifier
	model_dir = script_utils.create_directory("../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
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
	d_hist["WIMP_mass_" + WIMP_mass].Scale(d_lim_norm[WIMP_mass]/d_hist["WIMP_mass_" + WIMP_mass].Integral())
	# d_hist["WIMP_mass_" + WIMP_mass].Scale(500./d_hist["WIMP_mass_" + WIMP_mass].Integral())

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
	cc.Print(fig_dir + bolo_name + "_BDT_at_limit_mass_" + str(WIMP_mass) + ".eps")


def launch_classification():

	"""
	Detail:
		Main script to launch classification

	Args:
		bolo_name (str)        = bolometer name
		analysis_type (str)    = type of analysis (which cuts)
		bool_train (bool)      = boolean to pick test/training sample
		list_event_type (list) = list of event class

	Returns:
		void

	Raises:
		void	 
	"""

	bolo_name = "FID837"
	data_dir = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
	analysis_type = "ana_0.5_0_5"
	exposure = 65.
	bin_X, min_X, max_X = 200, -20, 20
	MVA_tag = ""
	list_mass = ["4", "5", "6", "7","10", "25"]
	list_mass = ["4"]

	d_list_variables = {}
	for mass in list_mass:
		d_list_variables[mass] = ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]

	#Loop over masses
	for WIMP_mass in list_mass:
		d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
						 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

		script_utils.print_utility("Processing WIMP mass: " + str(WIMP_mass) + " GeV")

		#Load data
		d_train = dp.get_data_array(bolo_name, 0, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[WIMP_mass], datasplit=0.5)
		d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[WIMP_mass], datasplit=1)
		d_eval  = dp.get_eval_array(bolo_name, analysis_type, "", d_list_variables[WIMP_mass])

		# setup parameters for xgboost
		param = {}
		param["objective"] = "binary:logitraw"
		param["eta"] = 0.1
		param["gamma"] = 0
		param["max_depth"] = 5
		if int(WIMP_mass) >8:
			param["max_depth"] = 7
		param["eval_metric"] = "auc"
		param["silent"] = 1
		param["subsample"] = 0.9
		num_rounds = 150

		#Plot results 
		plot_BDT_at_limit_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, weight_dir = "With_weights")


launch_classification()