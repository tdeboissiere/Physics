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
import matplotlib.gridspec as gridspec

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


def plot_scatter():
	path_train = "/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Training_data/FID837/ana_0.5_0_5/FID837_heatonly.csv"
	path_true = "/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Eval_data/FID837/ana_0.5_0_5/FID837_ana_0.5_0_5_fond.csv"


	df_train = pd.read_csv(path_train, usecols = ["EC1","EC2","EIA","EIB","EIC","EID", "HR"])
	df_true = pd.read_csv(path_true, usecols = ["EC1","EC2","EIA","EIB","EIC","EID", "HR"])

	# ax = df_train[:40000].plot(x="EC1", y="EC2", kind="scatter", color = "r", s=1)
	# df_true[:40000].plot(x="EC1", y="EC2", kind="scatter", color = "DarkBlue", s=1, ax=ax)
	# plt.show() 
	# raw_input()

	list_weights = [(1./len(df_train))*np.ones(len(df_train)), (1./len(df_true))*np.ones(len(df_true))]
	list_color = ["r", "b"]
	gs = gridspec.GridSpec(2,1)
	
	plt.subplot(gs[0]).hist([df_train["EC1"], df_true["EC1"]], weights=list_weights, color=list_color , histtype="step", bins=400)
	plt.subplot(gs[0]).set_yscale("log")
	plt.subplot(gs[0]).set_xlabel("EC1")
	plt.subplot(gs[1]).hist([df_train["EC2"], df_true["EC2"]], weights=list_weights, color=list_color , histtype="step", bins=400)
	plt.subplot(gs[1]).set_yscale("log")
	plt.subplot(gs[1]).set_xlabel("EC2")
	plt.tight_layout()
	plt.show()
	raw_input()



def plot_control_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, **kwargs):
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
	rvb = make_colormap([c("red"), c("red"), 0.33, c("red"), c("green"), 0.95, c("green")])
	
	#Build list of test data of interest
	list_event_type = ["FidGamma", "heatonly",  "WIMP_mass_" + WIMP_mass]
	list_event_type = ["heatonly"]

	#Output dir 
	fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/" + kwargs["weight_dir"] + "/")

	#Compute new columns for data frames
	d_eval["realdata"]["MVA"] = d_pred["realdata"]

	for event_type in list_event_type:
		d_test[event_type]["MVA"] = d_pred[event_type]

	# plt.scatter(d_eval["realdata"]["EC2"],d_eval["realdata"]["MVA"], c = "b")
	# plt.scatter(d_test["heatonly"]["EC2"][:35000],d_test["heatonly"]["MVA"][:35000], c="r", alpha="0.4")
	# plt.show()
	# raw_input()

	#Plots for real data
	plt.figure()
	plt.scatter(d_test["WIMP_mass_" + WIMP_mass]["EC1"].values, d_test["WIMP_mass_" + WIMP_mass]["EC2"].values, color ="0.6", s=1 )
	plt.scatter(d_eval["realdata"]["EC1"], d_eval["realdata"]["EC2"], c=d_eval["realdata"]["MVA"].values, cmap=rvb, vmin = -13, vmax = 0, s=10)
	cbar = plt.colorbar()
	cbar.set_label("MVA output", labelpad = 15, fontsize= 18)
	plt.xlabel("EC1 (keV)", fontsize = 20)
	plt.ylabel("EC2 (keV)", fontsize = 20)
	plt.ylim([0,1.5])
	plt.xlim([0,1.5])
	plt.grid(True)
	plt.savefig(fig_dir + bolo_name + "_real_data_WIMP_mass_" + WIMP_mass +  ".png")
	plt.close("all")


	for event_type in list_event_type:
		plt.figure()

		plt.scatter(d_test["WIMP_mass_" + WIMP_mass]["EC1"].values, d_test["WIMP_mass_" + WIMP_mass]["EC2"].values, color ="0.6", s=1 )
		plt.scatter(d_test[event_type]["EC1"][:35000], d_test[event_type]["EC2"][:35000], c=d_test[event_type]["MVA"].values[:35000], cmap=rvb, vmin = -13, vmax = 0, s=10)

		# d_test["WIMP_mass_" + WIMP_mass].plot(kind='scatter', x='EC', y='EFID', color="0.6", s = 1)
		# d_test[event_type].plot(kind='scatter', x='EC', y='EFID', c=d_test[event_type]["MVA"], cmap=rvb, vmin = 0, vmax = 1)

		cbar = plt.colorbar()
		cbar.set_label("BDT output", labelpad = 15, fontsize= 18)
		plt.xlabel("EC1 (keV)", fontsize = 20)
		plt.ylabel("EC2 (keV)", fontsize = 20)
		plt.ylim([0,1.5])
		plt.xlim([0,1.5])
		plt.grid(True)
		plt.savefig(fig_dir + bolo_name + "_" + event_type + "_WIMP_mass_" + WIMP_mass  + ".png")
		plt.close("all")

def launch_control_plots():

	"""
	Detail:
		launch control plots

	Args:

	Returns:
		void

	Raises:
		void	 
	"""

	bolo_name = "FID837"
	analysis_type = "ana_0.5_0_5"
	exposure = 65.
	bin_X, min_X, max_X = 500, -20, 20
	MVA_tag = ""
	list_mass = ["3", "4", "5", "6", "7","10", "25"]
	list_mass = ["3"]

	list_variables = ["EC1", "EC2"]

	#Loop over masses
	for WIMP_mass in list_mass:
		d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
						 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

		script_utils.print_utility("Processing WIMP mass: " + str(WIMP_mass) + " GeV")

		#Load data
		d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit=1)
		d_eval  = dp.get_eval_array(bolo_name, analysis_type, "", list_variables)

		plot_control_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, weight_dir = "With_weights")

launch_control_plots() 
# plot_scatter()