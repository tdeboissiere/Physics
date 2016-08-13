import numpy as np
import scipy as sp
import sys, os
sys.path.append("../")
import xgboost as xgb
import sklearn.cross_validation as cv
import itertools
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Normalizer, Binarizer
import data_preparation as dp 
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pylab as plt
import BDT_file_handler as BDT_fh
import BDT_CV_xgboost_for_visual as xgCV
import plotting_CV_functions as plfunc
import BDT_sensitivity as BDT_sensi
import Analysis_utilities as Ana_ut
import script_utils as script_utils


def launch_BDT_sensitivity_plots():

	"""
	Detail:
		Get BDT sensitivity plots

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
	analysis_type = "ana_1.5_0_5"
	bin_X, min_X, max_X = 2000, -40, 40
	exposure = 66
	list_mass = ["5", "6", "7", "10", "25"]
	list_mass = ["5"]
	d_bounds_1_5 = {"3":[0, 10], "4":[0, 10], "5":[0, 11], "6":[0, 10], "7":[0, 10], "10":[0, 10], "25":[0, 10]}

	list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID", "test", "HR"]
	#Loop over masses
	for WIMP_mass in list_mass:
		d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
							"S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "WIMP_mass_" + WIMP_mass: "WIMP"}

		# try :
		# 	eff_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
		# 	os.remove(eff_dir + bolo_name + "_WIMP_mass_" + WIMP_mass + "_BDT_cut_eff.root")
		# except OSError:
		# 	pass

		# try :
		# 	eff_dir = script_utils.create_directory("./ROOT_files/" + bolo_name + "/" + analysis_type + "/")
		# 	os.remove(eff_dir + bolo_name + "_WIMP_mass_" + WIMP_mass + "_BDT_sensitivity.root")
		# except OSError:
		# 	pass

		# BDT_sensi.get_plot_BDT_sensitivity_curve(WIMP_mass, bolo_name, analysis_type, bin_X, min_X, max_X, ["0.3","0.4","0.5","0.8","1"])
		# 
		#Single plot for PhD
		MVA_tag = ""
		BDT_sensi.get_plot_BDT_sensitivity_curve_PhD(bolo_name, WIMP_mass, analysis_type, MVA_tag, bin_X, min_X, max_X, exposure, list_variables, d_bounds_1_5, weight_dir = "With_weights")

		sys.exit()

		for heat_fraction in ["0.3","0.4","0.5","0.8","1"]: 
		# for heat_fraction in ["0.3", "1"]: 

			for key in d_event_dir.keys():
				if "heat" in key:
					d_event_dir.pop(key, None)
			d_event_dir["heatonly_heat_fraction_" + heat_fraction] = "Heatonly"

			MVA_tag = "heat_fraction_" + heat_fraction
			exposure = Ana_ut.get_exposure_in_days(bolo_name, "duration_KTH_heat_cut_" + heat_fraction)
			print "Exposure ", exposure


			#Load data
			d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit = 1)

			#Get sensitivity curves
			BDT_sensi.get_BDT_sensitivity_curve(d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, classifier_name=MVA_tag, weight_dir = "With_weights")

			#Get efficiency curves
			# BDT_sensi.get_BDT_efficiency_curves(d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, MVA_tag, exposure, bin_X, min_X, max_X, classifier_name=MVA_tag)

launch_BDT_sensitivity_plots()