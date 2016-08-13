from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import script_utils as script_utils
import math,sys
from array import array
import BDT_file_handler as BDT_fh
import xgboost as xgb
sys.path.append("../../")
import data_preparation as dp

def get_events_passing_cuts(bolo_name, WIMP_mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure):

	"""Simulate data and find the number of events that pass
	the BDT cut

	
	Detail:
		void

	Args:
		bolo_name           = (str) bolometer name
		mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
		analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
		bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
		exposure            = (float) exposure of the simulated data

	Returns:
		void

	Raises:
		void
	"""       

	d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
							"S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly"}
	WIMP_mass = str(WIMP_mass)
	d_event_dir["WIMP_mass_" + WIMP_mass] = "WIMP"

	d_test = dp.get_data_array(bolo_name, 1, analysis_type, d_event_dir.keys(), exposure, "fulldata")

	#Load cut value on BDT output
	cut_val = 0     
	with open ("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + ".txt", "r") as fcut:
		stuff = [elem.rstrip().split(",") for elem in fcut.readlines()]
		for elem in stuff:
			mass_val = elem[0]
			if int(WIMP_mass) ==int(mass_val):
				cut_val = float(elem[1])
	
	#Get scaling dict for data visualisation
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")

	# Get classifier
	model_dir = script_utils.create_directory("../../Classifier_files/" + bolo_name + "/" + analysis_type + "/")
	modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
	bst       = xgb.Booster({'nthread':16}, model_file = modelfile)

	#Get predictions on test sample
	d_pred = {}
	d_hist = {}
	d_color = {"S1Pb":kOrange-8, "S2Pb":kOrange-9, "S1Beta":kGreen+2, "S2Beta":kGreen-3,
                     "S1Gamma":kBlue-7, "S2Gamma":kBlue, "FidGamma":kAzure+10, "heatonly":kRed, "WIMP_mass_" + WIMP_mass:kGray, "neutron":kMagenta}

	for event_type in d_test.keys():
		d_pred[event_type] = bst.predict( xgb.DMatrix(d_test[event_type].iloc[:,:-2].values) )
		d_hist[event_type] = TH1F("h" + event_type, "h" + event_type, bin_X, min_X, max_X)
		PyRPl.fill_TH1(d_hist[event_type], d_pred[event_type])
		PyRPl.process_TH1(d_hist[event_type], use_fill_bool = True, color = d_color[event_type] )
		if "WIMP" not in event_type:
			d_hist[event_type].Scale(float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/float(d_hist[event_type].Integral()))
	d_hist["WIMP_mass_" + WIMP_mass].Scale(1./d_hist["WIMP_mass_" + WIMP_mass].Integral())

	list_hist_bckg =[d_hist["S1Pb"], d_hist["S2Pb"], d_hist["S1Beta"], d_hist["S2Beta"], d_hist["S1Gamma"], d_hist["S2Gamma"], d_hist["FidGamma"], d_hist["heatonly"]]

	hsum =TH1F("hsum","hsum", bin_X, min_X, max_X)
	for i in range(1,bin_X+1):
		sumcontent = sum([h.GetBinContent(i) for h in list_hist_bckg])
		hsum.SetBinContent(i, sumcontent)

	#Run Poisson simulations
	list_event_pass_cut=[]
	for nsimu in range(100):
		hdatasimu = TH1F("hdatasimu","hdatasimu", bin_X, min_X, max_X)
		for i in range(1,bin_X+1):
			hdatasimu.SetBinContent(i, np.random.poisson(hsum.GetBinContent(i)))
		bin_cut = hdatasimu.FindBin(cut_val)
		num_entry_cut = int(hdatasimu.Integral(bin_cut, max_X))
		list_event_pass_cut.append(str(num_entry_cut))
		del hdatasimu

	return list_event_pass_cut


def get_true_events_passing_cuts(bolo_name, WIMP_mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure, list_variables, **kwargs):

	"""Simulate data and find the number of events that pass
	the BDT cut

	
	Detail:
		void

	Args:
		bolo_name           = (str) bolometer name
		mass                = (int) WIMP mass
        d_cut               = (dict) analysis cut dict
		analysis_type       = (str) name of analysis (name indicates which ion cut, which resolution...)
		bin_X, min_X, max_X = (int, float, float) = settings for BDT histogram
		exposure            = (float) exposure of the simulated data
        list_variables (list)      = list of variables to retain for BDT


	Returns:
		void

	Raises:
		void
	"""       

	# Get classifier
	model_dir = script_utils.create_directory("../../Classifier_files/" + bolo_name + "/" + analysis_type + "/"+ kwargs["weight_dir"] + "/")
	modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) +".model"
	if kwargs.has_key("classifier_name"):
		modelfile = model_dir + "xgboost_classifier_mass_" + str(WIMP_mass) + "_" + kwargs["classifier_name"] + ".model"
	bst = xgb.Booster({'nthread':16}, model_file = modelfile)

	d_eval  = dp.get_eval_array(bolo_name, analysis_type,"", list_variables)

	#Load cut value on BDT output
	cut_val = 0     
	with open ("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis_type + ".txt", "r") as fcut:
		stuff = [elem.rstrip().split(",") for elem in fcut.readlines()]
		for elem in stuff:
			mass_val = elem[0]
			if int(WIMP_mass) ==int(mass_val):
				cut_val = float(elem[1])
	
	#Get predictions on data
	# hdata = TH1F("hdata", "hdata", bin_X, min_X, max_X)
	# PyRPl.fill_TH1(hdata, bst.predict( xgb.DMatrix(d_eval["realdata"].iloc[:,:].values) ) )

	vec_pred = np.array(bst.predict( xgb.DMatrix(d_eval["realdata"].iloc[:,:].values)))
	l = np.where(vec_pred > cut_val)
	return [str(len(vec_pred[l]))]

	# #Run Poisson simulations
	# list_event_pass_cut=[]
	# bin_cut = hdata.FindBin(cut_val)
	# num_entry_cut = int(hdata.Integral(bin_cut, max_X))
	# list_event_pass_cut.append(str(num_entry_cut))

	# return list_event_pass_cut


bolo_name = "FID837"
list_mass = ["3", "4", "5", "6", "7", "10", "25"]
# list_mass = ["5", "6", "7", "10", "25"]
# list_mass = ["5"]
analysis_type = "ana_0.5_0_5"
bin_X, min_X, max_X = 200, -20, 20
exposure = 66
d_cut       = {"ECinf": 0.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}
list_variables = ["EC1", "EC2", "EIA", "EIB","EIC", "EID",  "test", "HR"]
d_list_variables = {}
d_list_variables["3"] =  ["EC1", "EC2", "EIA", "EIC", "EFID"]
d_list_variables["4"] =  ["EC1", "EC2", "EIA", "EIC", "EFID"]
d_list_variables["5"] = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
d_list_variables["6"] = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
d_list_variables["7"] =  ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
d_list_variables["10"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"]
d_list_variables["25"] =  ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test", "prod"]
# for mass in list_mass:
# 		d_list_variables[mass] = ["EC1", "EC2", "EIA", "EIB","EIC", "EID"]    

# with open("./Text_files/events_passing_cuts_" + analysis_type + "_" + str(exposure) + ".txt", "w") as fev:
# with open("./Text_files/events_passing_cuts_" + analysis_type +  ".txt", "w") as fev:
# 	for WIMP_mass in list_mass:
# 		print WIMP_mass
# 		list_ev = get_events_passing_cuts(bolo_name, WIMP_mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure)
# 		fev.write(WIMP_mass + "," + ",".join(list_ev) + "\n")

with open("./Text_files/true_events_passing_cuts_" + analysis_type + ".txt", "w") as fev:
	for WIMP_mass in list_mass:
		print WIMP_mass
		list_ev = get_true_events_passing_cuts(bolo_name, WIMP_mass, d_cut, analysis_type, bin_X, min_X, max_X, exposure, d_list_variables[WIMP_mass], weight_dir = "With_weights")
		fev.write(WIMP_mass + "," + ",".join(list_ev) + "\n")