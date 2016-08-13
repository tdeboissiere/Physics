import numpy as np
import scipy as sc
import pandas as pd
import sklearn as sk
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier as GBC
from pandas import read_csv, DataFrame
import os
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
from sklearn import grid_search, linear_model
from hyperopt import fmin as hyperopt_fmin
from hyperopt import tpe, hp, STATUS_OK, space_eval
from sklearn.metrics import roc_curve, auc, roc_auc_score


def gbc_classification(d_train, d_event_dir, WIMP_mass, bolo_name, analysis_type):

	"""
	Detail:
		Classify data with GBC

	Args:
		d_train (dict)      = dict with training data
		d_event_dir (dict)  = dict to get the proper directory of each event class 
		WIMP_mass (str)     = WIMP mass
		bolo_name (str)     = bolometer name 
		analysis_type (str) = type of analysis (which box cut)

	Returns:
		void

	Raises:
		void	 
	"""

	x_train,y_train = dp.prepare_data_for_scikit(d_train)

	# ################
	# ##Optimisation 
	# ###############
	# def objective_function(x_int):
	# 	objective_function.n_iterations += 1
	# 	n_estimators, max_depth         = x_int
	# 	n_estimators                    = int(n_estimators)
	# 	max_depth                       = int(max_depth)
	# 	clf                             = GBC(max_depth=max_depth, n_estimators=n_estimators)
	# 	clf.fit(x_train,y_train)
	# 	scores                          = cross_val_score(clf, x_train, y_train, cv=2, scoring='accuracy')
	# 	print objective_function.n_iterations, \
	# 		": n_estimators = ", n_estimators, \
	# 		"\tmax_depth    = ", max_depth, \
	# 		"accuracy       = ", np.mean(scores)
	# 	return 1 - np.mean(scores)


	# objective_function.n_iterations = 0
	# best = hyperopt_fmin(objective_function,
	#     space=(hp.qloguniform('n_estimators', np.log(10), np.log(1000), 10), 
	#            hp.qloguniform('max_depth', np.log(2), np.log(100), 1)),
	#     algo=tpe.suggest,
	#     max_evals=2)

	# print best

	# objective_function.n_iterations = 0
	# xmin, fval = smac_fmin(objective_function, x0_int =(1,1), xmin_int=(1, 1), xmax_int=(10, 2), max_evaluations =2)
	# print "hello"
	# raw_input()

	#5 GeV : 250,6
	#25 GeV : 250,6

	clf = GBC(n_estimators=250, max_depth=6,min_samples_leaf=200,max_features=6,verbose=1)
	clf.fit(x_train,y_train) 

	#put it to file to load it later 
	output_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	pickle_out = open(output_dir + "gbc_classifier_mass_" + str(WIMP_mass) + ".pkl", "wb")
	pickle.dump(clf, pickle_out)
	pickle_out.close()


def logit_classification(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, **kwargs):

	"""
	Detail:
		Classify data with GBC

	Args:
		d_train (dict)      = dict with training data
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
	output_dir = script_utils.create_directory("./Classifier_files/" + bolo_name + "/" + analysis_type + "/With_weights/")

	x_train, weight_train, y_train = dp.prepare_data_for_scikit(d_train)
	x_test, weight_test, y_test = dp.prepare_data_for_scikit(d_test)

	logreg = linear_model.LogisticRegression(C=0.1, max_iter=1000)
	logreg.fit(x_train, y_train)

	# Construct matrix for test set
	y_test_pred = logreg.predict(x_test)
	
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
	analysis_type = "ana_0.5_0_5"
	exposure = 65.
	bin_X, min_X, max_X = 500, -20, 20
	MVA_tag = ""
	# list_mass = ["3", "4", "5", "6", "7", "10", "25"]
	list_mass = ["3", "4", "5", "6", "7","10", "25"]
	# list_mass = ["5", "6", "7","10", "25"]
	list_mass = ["25"]

	#Prepare data
	d_data_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
						 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "realdata": ""}
	for WIMP_mass in list_mass:
		d_data_dir["WIMP_mass_" + WIMP_mass] = "WIMP"

	for event_type in d_data_dir.keys():
		#O for training, 1 for test
		dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 0)
		dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 1)

	# sys.exit()

	# list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test", "HR"]
	list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID"]
	# list_variables = ["EIA", "EIB","EIC", "EID"]
	# list_variables = ["EC1", "EC2"]
	
	d_list_variables = {}
	d_list_variables["3"] =  ["EC1","EC2", "EIA", "EIC", "EFID"]
	d_list_variables["4"] =  ["EC1","EC2", "EIA", "EIC", "EFID"]
	d_list_variables["5"] = ["EC1","EC2", "EIA", "EIC", "EFID"]
	d_list_variables["6"] = ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
	d_list_variables["7"] =  ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
	d_list_variables["10"] =  ["EC1","EC2", "EIA", "EIC", "EFID", "test"]
	d_list_variables["25"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]	
	d_list_variables["25"] =  ["EC","EFID"]	
	for mass in list_mass:
		d_list_variables[mass] = ["EC1", "EC2", "EIA", "EIB","EIC", "EID"]


	#Loop over masses
	for WIMP_mass in list_mass:
		d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
						 "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

		script_utils.print_utility("Processing WIMP mass: " + str(WIMP_mass) + " GeV")

		#Load data
		d_train = dp.get_data_array(bolo_name, 0, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[WIMP_mass], datasplit=0.5)
		d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[WIMP_mass], datasplit=1)
		d_eval  = dp.get_eval_array(bolo_name, analysis_type, "", d_list_variables[WIMP_mass])

		#Launch classification
		logit_classification(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type)




launch_classification()