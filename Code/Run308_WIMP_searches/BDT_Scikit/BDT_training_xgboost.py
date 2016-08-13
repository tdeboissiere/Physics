#!/usr/bin/python
# this is the example script to use xgboost to train
import numpy as np
import scipy as sc
import pandas as pd
import sklearn as sk
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier as GBC
from pandas import read_csv, DataFrame
import os, sys
import pickle
import data_preparation as dp 
import data_visualisation as dv 
import PyROOTPlots as PyRPl
from ROOT import *
import BDT_file_handler as BDT_fh
import script_utils as script_utils
from pysmac.optimize import fmin as smac_fmin
from sklearn.cross_validation import cross_val_score
from sklearn.grid_search import GridSearchCV
from pysmac.optimize import fmin as smac_fmin
from sklearn import grid_search
from hyperopt import fmin as hyperopt_fmin
from hyperopt import tpe, hp, STATUS_OK, space_eval
import xgboost as xgb
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.decomposition import PCA, KernelPCA, RandomizedPCA
from sklearn.lda import LDA
import matplotlib.cm as cm
from scipy.stats import gaussian_kde
import Analysis_utilities as Ana_ut

def get_pca(d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type):

	"""
	Detail:
		Get PCA for test data set

	Args:
		d_test (dict)      = dict with testing data
		d_event_dir (dict)  = dict to get the proper directory of each event class 
		WIMP_mass (str)     = WIMP mass
		bolo_name (str)     = bolometer name 
		analysis_type (str) = type of analysis (which box cut)

	Returns:
		void

	Raises:
		void	 
	"""

	#Clean classification dir
	output_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	try :
		os.remove(output_dir + "pca_classifier_mass_" + str(WIMP_mass) + ".pkl")
	except OSError:
		pass

	x_test, weight_test, y_test = dp.prepare_data_for_scikit(d_test)
	pca = PCA()
	x_test = pca.fit_transform(x_test)
	#put it to file to load it later 
	pickle_out = open(output_dir + "pca_classifier_mass_" + str(WIMP_mass) + ".pkl", "wb")
	pickle.dump(pca, pickle_out)
	pickle_out.close()

def xgboost_classification_unweighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_rounds, **kwargs):

	"""
	Detail:
		Classify data with xgboost no weights

	Args:
		d_train (dict)      = dict with training data
		d_event_dir (dict)  = dict to get the proper directory of each event class 
		WIMP_mass (str)     = WIMP mass
		bolo_name (str)     = bolometer name 
		analysis_type (str) = type of analysis (which box cut)
		param (dict)        = dict of xgboost parameters
		num_rounds (int)    = number of boosting rounds

	Returns:
		void

	Raises:
		void	 
	"""

	#Clean classification dir
	output_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/No_weights/")

	x_train, weight_train, y_train = dp.prepare_data_for_scikit(d_train)
	x_test, weight_test, y_test = dp.prepare_data_for_scikit(d_test)

	# construct xgboost.DMatrix from numpy array
	xgmat = xgb.DMatrix( x_train, label=y_train)

	plst = list(param.items())+[("eval_metric", "auc")]
	watchlist = [ (xgmat,"train-unweighted") ]
	watchlist = []
	bst = xgb.train( plst, xgmat, num_rounds, watchlist )

	# save out model
	if kwargs.has_key("classifier_name"):
		try :
			os.remove(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model")
		except OSError:
			pass
		bst.save_model(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model")
	else :
		try :
			os.remove(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model")
		except OSError:
			pass
		bst.save_model(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model")

	print bst.get_fscore()
	print ("finish training")
	
	# Construct matrix for test set
	xgmat_test = xgb.DMatrix(x_test, weight = weight_test)
	y_test_pred = bst.predict(xgmat_test)
	
	AUC= roc_auc_score(y_test, y_test_pred, sample_weight = weight_test)
	print "AUC =", AUC

def xgboost_classification_weighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_rounds, **kwargs):

	"""
	Detail:
		Classify data with xgboost with weights

	Args:
		d_train (dict)      = dict with training data
		d_event_dir (dict)  = dict to get the proper directory of each event class 
		WIMP_mass (str)     = WIMP mass
		bolo_name (str)     = bolometer name 
		analysis_type (str) = type of analysis (which box cut)
		param (dict)        = dict of xgboost parameters
		num_rounds (int)    = number of boosting rounds

	Returns:
		void

	Raises:
		void	 
	"""

	#Clean classification dir
	output_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/With_weights/")

	x_train, weight_train, y_train = dp.prepare_data_for_scikit(d_train)
	x_test, weight_test, y_test = dp.prepare_data_for_scikit(d_test)

	# construct xgboost.DMatrix from numpy array
	xgmat = xgb.DMatrix( x_train, label=y_train, weight = weight_train)

	sum_wpos = sum(weight_train[y_train == 1])
	sum_wneg = sum(weight_train[y_train == 0])
	param['scale_pos_weight'] = sum_wneg / sum_wpos
	plst = list(param.items())+[("eval_metric", "auc")]
	watchlist = [ (xgmat,"train-weighted") ]		
	watchlist = [  ]		
	bst = xgb.train( plst, xgmat, num_rounds, watchlist )

	# save out model
	if kwargs.has_key("classifier_name"):
		try :
			os.remove(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model")
		except OSError:
			pass
		bst.save_model(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model")
	else :
		try :
			os.remove(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model")
		except OSError:
			pass
		bst.save_model(output_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model")

	print bst.get_fscore()
	print ("finish training")
	
	# Construct matrix for test set
	xgmat_test = xgb.DMatrix(x_test, weight = weight_test)
	y_test_pred = bst.predict(xgmat_test)
	
	AUC= roc_auc_score(y_test, y_test_pred, sample_weight = weight_test)
	print "WIMP_mass", WIMP_mass, "AUC =", AUC

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
	analysis_type = "ana_1.5_0_5"
	exposure = 65.
	bin_X, min_X, max_X = 100, -20, 20
	MVA_tag = ""
	# list_mass = ["3", "4", "5", "6", "7", "10", "25"]
	list_mass = ["4", "5", "6", "7","10", "25"]
	list_mass = ["5", "6", "7","10", "25"]
	# list_mass = ["4"]

	# #Prepare data
	# d_data_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
	# 					 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "realdata": ""}
	# for WIMP_mass in list_mass:
	# 	d_data_dir["WIMP_mass_" + WIMP_mass] = "WIMP"

	# for event_type in d_data_dir.keys():
	# 	#O for training, 1 for test
	# 	dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 0)
	# 	dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 1)

	# sys.exit()

	# list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test", "HR"]
	list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID"]
	# list_variables = ["EIA", "EIB","EIC", "EID"]
	# list_variables = ["EC1", "EC2"]
	
	d_list_variables = {}
	d_list_variables["3"] =  ["EC1", "EC2", "EIA", "EIC", "EFID"]
	d_list_variables["4"] =  ["EC1", "EC2", "EIA", "EIC", "EFID"]
	d_list_variables["5"] = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
	d_list_variables["6"] = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
	d_list_variables["7"] =  ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
	d_list_variables["10"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]
	d_list_variables["25"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test", "prod"]
	for mass in list_mass:
		d_list_variables[mass] =  ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test", "HR"]


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

		#Launch classification
		# xgboost_classification_weighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_rounds)

		#Plot results 
		dv.plot_results_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, weight_dir = "With_weights")
		# for pca_index in range(7):
		# 	dv.plot_PCA_stuff(d_test, d_eval, d_event_dir, classifier_type, WIMP_mass, bolo_name, analysis_type, exposure, bin_X, min_X, max_X, pca_index)
		# dv.plot_control_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, list_variables, weight_dir = "With_weights")

def launch_classification_heat_rate_cut():

	"""
	Detail:
		Main script to launch classification using heat with heat rate cut

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
	data_dir_heat_rate = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr_heatrate_cut/"
	analysis_type = "ana_1.5_0_5"
	bin_X, min_X, max_X = 200, -40, 40
	# list_mass = ["5", "6", "7", "10", "25"]
	list_mass = ["5"]

	#Prepare data
	d_data_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
						 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma"}
	for WIMP_mass in list_mass:
		d_data_dir["WIMP_mass_" + WIMP_mass] = "WIMP"

	for event_type in d_data_dir.keys():
		#O for training, 1 for test
		dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 0)
		dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 1)

	for heat_fraction in ["0.1","0.2","0.3","0.4","0.5","0.8","1"]:
		d_data_dir["realdata"] = "_heat_fraction_" + heat_fraction
		d_data_dir["heatonly_heat_fraction_" + heat_fraction] = "Heatonly"
		dp.root_to_csv(bolo_name, data_dir_heat_rate, analysis_type, "realdata", d_data_dir, 0)
		dp.root_to_csv(bolo_name, data_dir_heat_rate, analysis_type, "heatonly_heat_fraction_" + heat_fraction, d_data_dir, 0)
		dp.root_to_csv(bolo_name, data_dir_heat_rate, analysis_type, "heatonly_heat_fraction_" + heat_fraction, d_data_dir, 1)


	list_variables = ["EIA", "EIB", "EID", "EIC","EC1", "EC2",  "test"]


	# for heat_fraction in ["0.2","0.3","0.4","0.5","0.8","1"]: 
	for heat_fraction in ["0.3"]: 

		MVA_tag = "heat_fraction_" + heat_fraction
		exposure = Ana_ut.get_exposure_in_days(bolo_name, "duration_KTH_heat_cut_" + heat_fraction)
		print exposure

		#Loop over masses
		for WIMP_mass in list_mass:
			d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
							 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly_heat_fraction_" + heat_fraction:"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

			#Load data
			d_train = dp.get_data_array(bolo_name, 0, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit = .2)
			d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit = 1)
			d_eval  = dp.get_eval_array(bolo_name, analysis_type, "_heat_fraction_" + heat_fraction, list_variables)

			# setup parameters for xgboost
			param = {}
			param["objective"] = "binary:logitraw"
			param["eta"] = 0.1
			param["gamma"] = 0
			param["max_depth"] = 3
			param["eval_metric"] = "auc"
			param["silent"] = 1
			param["subsample"] = 0.9
			num_rounds = 80

			#Launch classification
			# xgboost_classification_weighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_rounds, classifier_name=MVA_tag)

			#Plot results 
			dv.plot_results_xgboost(d_test, d_eval, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, classifier_name=MVA_tag, weight_dir = "With_weights")
			# dv.plot_control_xgboost(d_test, d_eval, d_event_dir, classifier_type, WIMP_mass, bolo_name, analysis_type, exposure, bin_X, min_X, max_X)


launch_classification()
# launch_classification_heat_rate_cut()
