from ROOT import *
import numpy as np
import script_utils as script_utils
import os, sys
sys.path.append("../")
import scipy.stats
import PyROOTPlots as PyRPl
from math import sqrt as sqrt
import pandas as pd
import data_preparation as dp 
import xgboost as xgb
import matplotlib.pylab as plt
import BDT_file_handler as BDT_fh
import matplotlib.gridspec as gridspec


def get_KS_training_test(bolo_name, analysis_type, exposure, MVA_tag, **kwargs):
	"""
	Detail:
		Plot results

	Args:
		bolo_name (str)     = bolometer name 
		analysis_type (str) = type of analysis (which box cut)
		MVA_tag (str)       = indicates which scaling file to use
		exposure (float)    = exposure in days

	Returns:
		void

	Raises:
		void	 
	"""

	try:
		kwargs["weight_dir"]
	except KeyError:
		sys.exit()

	list_mass = ["5", "6", "7","10", "25"]
	list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test"]
	#Get scaling dict for data visualisation
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")
	
	d_event_dir = {}
	#Loop over masses
	for WIMP_mass in list_mass:
		d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
						 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

		#Load data
		d_train = dp.get_data_array(bolo_name, 0, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit=1)
		d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit=1)

		# Get classifier
		model_dir = script_utils.create_directory("../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
		modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
		if kwargs.has_key("classifier_name"):
			modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
		bst = xgb.Booster({'nthread':16}, model_file = modelfile)

		#Get predictions on test sample
		d_pred_test = {}
		d_pred_train = {}

		for event_type in d_test.keys():
			d_pred_test[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-3].values) )
			d_pred_train[event_type] = bst.predict( xgb.DMatrix(d_train[event_type].iloc[:,:-3].values) )


		arr_train =  np.concatenate([d_pred_train[event_type] for event_type in d_test.keys()])
		arr_test =  np.concatenate([d_pred_test[event_type] for event_type in d_test.keys()])
		print scipy.stats.ks_2samp(arr_train, arr_test),"   ", 1.36*sqrt(len(arr_train) + len(arr_test))/sqrt(len(arr_train) * len(arr_test))


	# bin_X, min_X, max_X = 200, -20, 20

	# #Get predictions on test sample
	# d_hist_train = {}
	# d_hist_test = {}
	# d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
	# 				 "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10,"heatonly": kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}

	# for event_type in d_test.keys():
	# 	d_hist_train[event_type] = TH1F("htrain" + event_type, "htrain" + event_type, bin_X, min_X, max_X)
	# 	d_hist_test[event_type] = TH1F("htest" + event_type, "htest" + event_type, bin_X, min_X, max_X)
	# 	PyRPl.fill_TH1(d_hist_train[event_type], d_pred_train[event_type])
	# 	PyRPl.fill_TH1(d_hist_test[event_type], d_pred_test[event_type])
	# 	PyRPl.process_TH1(d_hist_train[event_type], use_fill_bool = True, color = d_color[event_type] )
	# 	PyRPl.process_TH1(d_hist_test[event_type], use_fill_bool = True, color = d_color[event_type] )
	# 	if "WIMP" not in event_type:
	# 		d_hist_train[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist_train[event_type].Integral()))
	# 		d_hist_test[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist_test[event_type].Integral()))


	# d_hist_train["S1Pb"].Add(d_hist_train["S2Pb"])
	# d_hist_train["S1Beta"].Add(d_hist_train["S2Beta"])
	# d_hist_train["FidGamma"].Add(d_hist_train["S1Gamma"])
	# d_hist_train["FidGamma"].Add(d_hist_train["S2Gamma"])

	# d_hist_test["S1Pb"].Add(d_hist_test["S2Pb"])
	# d_hist_test["S1Beta"].Add(d_hist_test["S2Beta"])
	# d_hist_test["FidGamma"].Add(d_hist_test["S1Gamma"])
	# d_hist_test["FidGamma"].Add(d_hist_test["S2Gamma"])

	# list_hist_train =[d_hist_train["S1Pb"], d_hist_train["S1Beta"], d_hist_train["FidGamma"], d_hist_train["heatonly"], d_hist_train["WIMP_mass_" + WIMP_mass]]
	# list_hist_test =[d_hist_test["S1Pb"], d_hist_test["S1Beta"], d_hist_test["FidGamma"], d_hist_test["heatonly"], d_hist_test["WIMP_mass_" + WIMP_mass]]
	# hstrain=THStack("hstrain", "hstrain")
	# hstest=THStack("hstest", "hstest")
	# for hist in list_hist_train:
	# 	hstrain.Add(hist)

	# for hist in list_hist_test:
	# 	hstest.Add(hist)

	# cc = TCanvas("cc", "cc")
	# cc.Divide(2,1)
	# h1=TH1F("h1","h1", bin_X, min_X, max_X)
	# PyRPl.process_TH1(h1, X_title="BDT ouput", min_Y = 1E-1, max_Y = 20000)
	
	# cc.cd(1)
	# gPad.SetLogy()
	# h1.Draw()
	# hstrain.Draw("same")

	# cc.cd(2)
	# gPad.SetLogy()
	# h1.Draw()
	# hstest.Draw("same")


	# raw_input()

	
	# #Also plot data for visual check
	# list_ev = ["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly", "WIMP_mass_" + WIMP_mass]
	# d_color = {"S1Pb":"peru", "S2Pb":"saddlebrown", "S1Beta":"olive", "S2Beta":"forestgreen", "S1Gamma":"royalblue", "S2Gamma":"dodgerblue", 
 #                "FidGamma":"deepskyblue", "heatonly":"red", "WIMP_mass_" + WIMP_mass:"dimgray", "neutron":"magenta"}
	# d_index_train, d_index_test, d_weight_test, d_weight_train = {}, {}, {}, {}
	# df_train = pd.concat( [d_train[event_type] for event_type in d_train.keys() ], ignore_index = True )
	# df_test = pd.concat( [d_test[event_type] for event_type in d_test.keys() ], ignore_index = True )

	# for event_type in d_train.keys():
	# 	d_index_train[event_type] = df_train.ix[df_train["EventID"]==event_type]
	# 	d_index_test[event_type] = df_test.ix[df_test["EventID"]==event_type]
	# 	if "WIMP" not in event_type:
	# 		d_weight_train[event_type] = float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure*np.ones(len(d_index_train[event_type]))/len(d_index_train[event_type])
	# 		d_weight_test[event_type] = float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure*np.ones(len(d_index_test[event_type]))/len(d_index_test[event_type])
	# 	else:
	# 		d_weight_train[event_type] = 8000*np.ones(len(d_index_train[event_type]))/len(d_index_train[event_type])
	# 		d_weight_test[event_type] = 8000*np.ones(len(d_index_test[event_type]))/len(d_index_test[event_type])


	# list_weights_train = [d_weight_train[event_type] for event_type in list_ev]
	# list_weights_test = [d_weight_test[event_type] for event_type in list_ev]
	# list_color = [d_color[event_type] for event_type in list_ev]
	# gs = gridspec.GridSpec(2,1)

	# plt.subplot(gs[0]).hist([d_index_train[event_type][val] for event_type in list_ev], weights =list_weights_train, color=list_color , stacked=True, histtype="stepfilled", bins=200)
	# plt.subplot(gs[1]).hist([d_index_test[event_type][val] for event_type in list_ev], weights =list_weights_test, color=list_color , stacked=True, histtype="stepfilled", bins=200)
	# for i in range(2):
	# 	plt.subplot(gs[index]).set_yscale("log")
	# 	plt.subplot(gs[index]).set_xlabel(val)
	# plt.tight_layout()
	# plt.show()



	raw_input()


bolo_name = "FID837"
analysis_type = "ana_1.5_0_5"
exposure = 65.
get_KS_training_test(bolo_name, analysis_type, exposure, "", weight_dir = "With_weights")