from ROOT import *
# from root_numpy import root2array
import script_utils as script_utils
import numpy as np
import pandas as pd
import BDT_file_handler as BDT_fh
import pickle

def root_to_csv_KDE():

	"""
	Detail:
		Transform root TTree to .csv files to store the data

	Args:
		bolo_name (str)     = bolometer name
		data_dir (str)      = data directory (ROOT TTree files)
		analysis_type (str) = type of analysis (which cuts)
		event_type (str)    = event class
		d_event_dir (dict)  = dict to get the proper directory of each event class 
		bool_train (bool)   = boolean to pick test/training sample

	Returns:
		void

	Raises:
		void	 
	"""

	bolo_name = "FID837"
	data_dir = "./ROOT_files/ana_0.5_0_5_for_KDE/"
	arr = root2array(data_dir + bolo_name + "_WIMP_mass_10_tree.root", "t_new0")
	np.savetxt(data_dir + bolo_name + "_WIMP_mass_10.csv", arr, delimiter = ",", fmt = "%.5f", header = "EC1,EC2,EIA,EIB,EIC,EID,HR", comments="")

# root_to_csv_KDE()

def root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_event_dir, bool_train):

	"""
	Detail:
		Transform root TTree to .csv files to store the data

	Args:
		bolo_name (str)     = bolometer name
		data_dir (str)      = data directory (ROOT TTree files)
		analysis_type (str) = type of analysis (which cuts)
		event_type (str)    = event class
		d_event_dir (dict)  = dict to get the proper directory of each event class 
		bool_train (bool)   = boolean to pick test/training sample

	Returns:
		void

	Raises:
		void	 
	"""

	if event_type == "realdata":
		data_dir = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/"
		arr = root2array(data_dir + bolo_name + "_" + analysis_type + d_event_dir[event_type] + "_fond.root", "data")
		# arrbis = []
		# for i in range(arr.shape[0]):
			# fid = 0.5*(arr["EIB"][i]+arr["EID"][i])
			# if fid<0.7:
				# arrbis.append([arr["EC1"][i],arr["EC2"][i],arr["EIA"][i],arr["EIB"][i],arr["EIC"][i],arr["EID"][i],arr["HR"][i]])
		out_dir = script_utils.create_directory("./Eval_data/" + bolo_name + "/" + analysis_type + "/")
		np.savetxt(out_dir + bolo_name + "_" + analysis_type  + d_event_dir[event_type] + "_fond.csv", arr, delimiter = ",", fmt = "%.5f", header = "EC1,EC2,EIA,EIB,EIC,EID,EC,EFID,HR,RUN,SN", comments = "")		
		# np.savetxt(out_dir + bolo_name + "_heatonly_from_fond.csv", arrbis, delimiter = ",", fmt = "%.5f", header = "EC1,EC2,EIA,EIB,EIC,EID,HR", comments = "")		

	else:
		data_dir += "BDT_" + bolo_name + "/" + analysis_type + "/" + d_event_dir[event_type] + "/ROOT_files/"
		arr = root2array(data_dir + bolo_name + "_" + event_type + "_tree.root", "t_new" + str(bool_train))

		if bool_train == 0:
			out_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Training_data/" + bolo_name + "/" + analysis_type + "/")
			np.savetxt(out_dir + bolo_name + "_" + event_type + ".csv", arr, delimiter = ",", fmt = "%.5f", header = "EC1,EC2,EIA,EIB,EIC,EID,HR", comments="")
		else:
			out_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Test_data/" + bolo_name + "/" + analysis_type + "/")
			np.savetxt(out_dir + bolo_name + "_" + event_type + ".csv", arr, delimiter = ",", fmt = "%.5f", header = "EC1,EC2,EIA,EIB,EIC,EID,HR", comments="")


def get_data_array(bolo_name, bool_train, analysis_type, MVA_tag, list_event_type, exposure, list_variables, **kwargs):

	"""
	Detail:
		Load .csv data as pandas tables store to dict 

	Args:
		bolo_name (str)            = bolometer name
		bool_train (bool)          = boolean to pick test/training sample
		analysis_type (str)        = type of analysis (which cuts)
		MVA_tag (str)              = tag to select the MVA scaling files
		list_event_type (list)     = list of event class
		exposure (float)           = exposure duration in days
		datasplit (str, kwarg)     = fraction of data to be randomly selected
		d_num_events (dict, kwarg) = number of training event dict per class
		list_variables (list)      = list of variables to retain for BDT
	Returns:
		d_array (dict) = a dict with the training or test data for each class

	Raises:
		d_array	 
	"""

	#Get scaling dict to set the weights
	d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, MVA_tag)

	d_array = {}
	data_dir = ""
	if not bool_train:
		data_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Training_data/" + bolo_name + "/" + analysis_type + "/")
	else:
		data_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Test_data/" + bolo_name + "/" + analysis_type + "/")

	list_bckg_type =  [elem for elem in list_event_type if "WIMP" not in elem]
	list_signal_type  = [elem for elem in list_event_type if "WIMP" in elem]
	
	# kde_dir = script_utils.create_directory("./KDE_test/ana_0.5_0_5_for_KDE/" )
	# kde_file = open(kde_dir + bolo_name + "_KDE_heatonly.pkl", "rb")
	# kde = pickle.load(kde_file)
	# kde_file.close()

	for event_type in list_bckg_type :
		d_array[event_type] = pd.read_csv(data_dir + bolo_name + "_" + event_type + ".csv", usecols = ["EC1","EC2","EIA","EIB","EIC","EID", "HR"])
		col_EC    = 0.5*(d_array[event_type]["EC1"] + d_array[event_type]["EC2"])
		col_EI    = 0.5*(d_array[event_type]["EIA"] + d_array[event_type]["EIB"] + d_array[event_type]["EIC"] + d_array[event_type]["EID"])
		col_EIFID = 0.5*(d_array[event_type]["EIB"] + d_array[event_type]["EID"])
		col_EIS1  = 0.5*(d_array[event_type]["EIA"] + d_array[event_type]["EIB"])
		col_EIS2  = 0.5*(d_array[event_type]["EIC"] + d_array[event_type]["EID"])
		col_EI_for_ER = ( 1.5*d_array[event_type]["EIA"] + 4*d_array[event_type]["EIB"] + 1.5*d_array[event_type]["EIC"] + 4*d_array[event_type]["EID"])

		# d_array[event_type]["ECdiff"] = np.log(2+d_array[event_type]["EC1"] - d_array[event_type]["EC2"])
		d_array[event_type]["ECprod"] = (d_array[event_type]["EC1"]*d_array[event_type]["EC2"])
		d_array[event_type]["EFIDprod"] = np.log(5+d_array[event_type]["EIB"]*d_array[event_type]["EID"])
		d_array[event_type]["Evetprod"] = (d_array[event_type]["EIA"]*d_array[event_type]["EIC"])
		# d_array[event_type]["EIprod"] = (d_array[event_type]["EFIDprod"]*d_array[event_type]["Evetprod"])
		d_array[event_type]["ECdiff"] = (d_array[event_type]["EC1"] - d_array[event_type]["EC2"]) /(d_array[event_type]["EC1"] + d_array[event_type]["EC2"])
		d_array[event_type]["EC"] = 0.5*(d_array[event_type]["EC1"] + d_array[event_type]["EC2"] )
		d_array[event_type]["EFID"] = 0.5*(d_array[event_type]["EIB"] + d_array[event_type]["EID"])
		# d_array[event_type]["EFIDdiff"] = np.log(1+np.abs((d_array[event_type]["EIB"] - d_array[event_type]["EID"])/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])))
		d_array[event_type]["EFIDdiff"] = (d_array[event_type]["EIB"] - d_array[event_type]["EID"]) #/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])
		d_array[event_type]["sum_ion"] = (d_array[event_type]["EIB"] + d_array[event_type]["EID"] + d_array[event_type]["EIA"] + d_array[event_type]["EIC"]) #/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])

		# d_array[event_type]["EFIDdiff"] = (d_array[event_type]["EIB"] - d_array[event_type]["EID"])/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])

		temp_df = pd.concat([col_EIS1, col_EIS2,col_EIFID], axis=1, keys = ["EIS1", "EIS2", "EFID"])
		d_array[event_type]["max_ion"] = temp_df[["EIS1", "EIS2", "EFID"]].max(axis=1)
		d_array[event_type]["max_vet"] = temp_df[["EIS1", "EIS2"]].max(axis=1)
		d_array[event_type]["max_vet"] = d_array[event_type][["EIA", "EIC"]].max(axis=1)

		# d_array[event_type]["test"] = np.log(1+ np.abs(col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18)))
		d_array[event_type]["test"] = col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18)
		d_array[event_type]["testEC"] = (d_array[event_type]["EC1"] + d_array[event_type]["EC2"])*(col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18))
		# d_array[event_type]["test"] = np.log(10+(d_array[event_type]["EIB"] - 0.16*np.power(abs((1+8./3)*col_EC - 0.333*col_EI_for_ER), 1.18))*(d_array[event_type]["EID"] - 0.16*np.power(abs((1+8./3)*col_EC - 0.333*col_EI_for_ER), 1.18)))
		d_array[event_type]["ER"] = (1+8./3)*col_EC - 0.333*col_EI_for_ER
		d_array[event_type]["prod"] = d_array[event_type]["EC1"]*d_array[event_type]["EC2"]*d_array[event_type]["EIB"]*d_array[event_type]["EID"]/(1+d_array[event_type]["EIA"]*d_array[event_type]["EIC"])
		d_array[event_type]["to1_3"] =  np.power(d_array[event_type]["EC1"] -1.1,2) +  np.power(d_array[event_type]["EC2"] -1.1,2)
		d_array[event_type]["to1_3"]+=  np.power(d_array[event_type]["EIA"],2) +  np.power(d_array[event_type]["EIB"] -1.1,2)
		d_array[event_type]["to1_3"]+=  np.power(d_array[event_type]["EIC"],2) +  np.power(d_array[event_type]["EID"] -1.1,2)
		col_Q = np.abs(col_EI/d_array[event_type]["ER"])
		d_array[event_type]["Q"] = np.amin( np.concatenate( (np.reshape(col_Q, (col_Q.shape[0], 1)), 2*np.ones((col_Q.shape[0],1))), axis=1 ), axis=1)
		#####################
		#####################
		#####################
		# print event_type, len(d_array[event_type])


		d_array[event_type] = d_array[event_type][list_variables]

		if kwargs.has_key("datasplit"):
			split = kwargs["datasplit"]
			len_data = d_array[event_type].shape[0]
			d_array[event_type] = d_array[event_type][:int(split*len_data)]
		if kwargs.has_key("d_num_events"):
			d_array[event_type] = d_array[event_type][:kwargs["d_num_events"][event_type]]
		d_array[event_type]["EventID"] = event_type
		d_array[event_type]["weight"] = float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure/(len(d_array[event_type]))*np.ones((len(d_array[event_type])))
		d_array[event_type]["tag"] = np.zeros((len(d_array[event_type])))

		# array = d_array[event_type][["EC1", "EC2", "EIA", "EIB", "EIC", "EID"]].values
		# array = kde.score_samples(array)
		# d_array[event_type]["prob"] = array
		# d_array[event_type] = d_array[event_type][list_variables + ["prob", "EventID", "weight", "tag"]]


	#Get sum of bckg weights
	sum_bckg_weights = sum( [ float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure for event_type in list_bckg_type ] )

	for event_type in list_signal_type :
		d_array[event_type] = pd.read_csv(data_dir + bolo_name + "_" + event_type + ".csv", usecols = ["EC1","EC2","EIA","EIB","EIC","EID", "HR"])
		col_EC    = 0.5*(d_array[event_type]["EC1"] + d_array[event_type]["EC2"])
		col_EI    = 0.5*(d_array[event_type]["EIA"] + d_array[event_type]["EIB"] + d_array[event_type]["EIC"] + d_array[event_type]["EID"])
		col_EIFID = 0.5*(d_array[event_type]["EIB"] + d_array[event_type]["EID"])
		col_EIS1  = 0.5*(d_array[event_type]["EIA"] + d_array[event_type]["EIB"])
		col_EIS2  = 0.5*(d_array[event_type]["EIC"] + d_array[event_type]["EID"])
		col_EI_for_ER = ( 1.5*d_array[event_type]["EIA"] + 4*d_array[event_type]["EIB"] + 1.5*d_array[event_type]["EIC"] + 4*d_array[event_type]["EID"])

		# d_array[event_type]["ECdiff"] = np.log(2+d_array[event_type]["EC1"] - d_array[event_type]["EC2"])
		d_array[event_type]["ECprod"] = (d_array[event_type]["EC1"]*d_array[event_type]["EC2"])
		d_array[event_type]["EFIDprod"] = np.log(5+d_array[event_type]["EIB"]*d_array[event_type]["EID"])
		d_array[event_type]["Evetprod"] = (d_array[event_type]["EIA"]*d_array[event_type]["EIC"])
		# d_array[event_type]["EIprod"] = (d_array[event_type]["EFIDprod"]*d_array[event_type]["Evetprod"])
		d_array[event_type]["ECdiff"] = (d_array[event_type]["EC1"] - d_array[event_type]["EC2"]) /(d_array[event_type]["EC1"] + d_array[event_type]["EC2"])
		d_array[event_type]["EC"] = 0.5*(d_array[event_type]["EC1"] + d_array[event_type]["EC2"] )
		d_array[event_type]["EFID"] = 0.5*(d_array[event_type]["EIB"] + d_array[event_type]["EID"])
		# d_array[event_type]["EFIDdiff"] = np.log(1+np.abs((d_array[event_type]["EIB"] - d_array[event_type]["EID"])/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])))
		d_array[event_type]["EFIDdiff"] = (d_array[event_type]["EIB"] - d_array[event_type]["EID"]) #/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])
		d_array[event_type]["sum_ion"] = (d_array[event_type]["EIB"] + d_array[event_type]["EID"] + d_array[event_type]["EIA"] + d_array[event_type]["EIC"]) #/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])

		# d_array[event_type]["EFIDdiff"] = (d_array[event_type]["EIB"] - d_array[event_type]["EID"])/(d_array[event_type]["EIB"] + d_array[event_type]["EID"])

		temp_df = pd.concat([col_EIS1, col_EIS2,col_EIFID], axis=1, keys = ["EIS1", "EIS2", "EFID"])
		d_array[event_type]["max_ion"] = temp_df[["EIS1", "EIS2", "EFID"]].max(axis=1)
		d_array[event_type]["max_vet"] = temp_df[["EIS1", "EIS2"]].max(axis=1)
		d_array[event_type]["max_vet"] = d_array[event_type][["EIA", "EIC"]].max(axis=1)

		# d_array[event_type]["test"] = np.log(1+ np.abs(col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18)))
		d_array[event_type]["test"] = col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18)
		d_array[event_type]["testEC"] = (d_array[event_type]["EC1"] + d_array[event_type]["EC2"])*(col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18))
		# d_array[event_type]["test"] = np.log(10+(d_array[event_type]["EIB"] - 0.16*np.power(abs((1+8./3)*col_EC - 0.333*col_EI_for_ER), 1.18))*(d_array[event_type]["EID"] - 0.16*np.power(abs((1+8./3)*col_EC - 0.333*col_EI_for_ER), 1.18)))
		d_array[event_type]["ER"] = (1+8./3)*col_EC - 0.333*col_EI_for_ER
		d_array[event_type]["prod"] = d_array[event_type]["EC1"]*d_array[event_type]["EC2"]*d_array[event_type]["EIB"]*d_array[event_type]["EID"]/(1+d_array[event_type]["EIA"]*d_array[event_type]["EIC"])
		d_array[event_type]["to1_3"] =  np.power(d_array[event_type]["EC1"] -1.1,2) +  np.power(d_array[event_type]["EC2"] -1.1,2)
		d_array[event_type]["to1_3"]+=  np.power(d_array[event_type]["EIA"],2) +  np.power(d_array[event_type]["EIB"] -1.1,2)
		d_array[event_type]["to1_3"]+=  np.power(d_array[event_type]["EIC"],2) +  np.power(d_array[event_type]["EID"] -1.1,2)
		col_Q = np.abs(col_EI/d_array[event_type]["ER"])
		d_array[event_type]["Q"] = np.amin( np.concatenate( (np.reshape(col_Q, (col_Q.shape[0], 1)), 2*np.ones((col_Q.shape[0],1))), axis=1 ), axis=1)

		d_array[event_type] = d_array[event_type][list_variables]

		if kwargs.has_key("datasplit"):
			split = kwargs["datasplit"]
			len_data = d_array[event_type].shape[0]
			d_array[event_type] = d_array[event_type][:int(split*len_data)]

		if kwargs.has_key("d_num_events"):
			d_array[event_type] = d_array[event_type][:kwargs["d_num_events"][event_type]]
		d_array[event_type]["EventID"] = event_type
		d_array[event_type]["weight"] = sum_bckg_weights/(len(d_array[event_type]))*np.ones((len(d_array[event_type])))
		# d_array[event_type]["weight"] = 1./(len(d_array[event_type]))*np.ones((len(d_array[event_type])))
		d_array[event_type]["tag"] = np.ones((len(d_array[event_type])))

		# array = d_array[event_type][["EC1", "EC2", "EIA", "EIB", "EIC", "EID"]].values
		# array = kde.score_samples(array)
		# d_array[event_type]["prob"] = array
		# d_array[event_type] = d_array[event_type][list_variables + ["prob", "EventID", "weight", "tag"]]

	return d_array

def get_eval_array(bolo_name, analysis_type, realdata_tag, list_variables):

	"""
	Detail:
		Load .csv data as pandas tables store to dict 

	Args:
		bolo_name (str)     = bolometer name
		analysis_type (str) = type of analysis (which cuts)
		realdata_tag (str)  = tag to identify heat cut files
		list_variables (list)      = list of variables to retain for BDT


	Returns:
		d_array (dict) = a dict with the eval data

	Raises:
		d_array	 
	"""

	d_array = {}
	data_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/BDT_Scikit/Eval_data/" + bolo_name + "/" + analysis_type + "/")
	d_array["realdata"] = pd.read_csv(data_dir + bolo_name + "_" + analysis_type +  realdata_tag + "_fond.csv", usecols = ["EC1","EC2","EIA","EIB","EIC","EID","EC","EFID", "HR"])
	col_EC = 0.5*(d_array["realdata"]["EC1"] + d_array["realdata"]["EC2"])
	col_EI = 0.5*( d_array["realdata"]["EIA"] + d_array["realdata"]["EIB"] + d_array["realdata"]["EIC"] + d_array["realdata"]["EID"])
	col_EI_for_ER = ( 1.5*d_array["realdata"]["EIA"] + 4*d_array["realdata"]["EIB"] + 1.5*d_array["realdata"]["EIC"] + 4*d_array["realdata"]["EID"])

	col_EIFID = 0.5*(d_array["realdata"]["EIB"] + d_array["realdata"]["EID"])
	col_EIS1  = 0.5*(d_array["realdata"]["EIA"] + d_array["realdata"]["EIB"])
	col_EIS2  = 0.5*(d_array["realdata"]["EIC"] + d_array["realdata"]["EID"])

	temp_df = pd.concat([col_EIS1, col_EIS2,col_EIFID], axis=1, keys = ["EIS1", "EIS2", "EFID"])
	d_array["realdata"]["max_ion"] = temp_df[["EIS1", "EIS2", "EFID"]].max(axis=1)
	d_array["realdata"]["EFIDdiff"] = (d_array["realdata"]["EIB"] - d_array["realdata"]["EID"])

	d_array["realdata"]["to1_3"] =  np.power(d_array["realdata"]["EC1"] -1.1,2) +  np.power(d_array["realdata"]["EC2"] -1.1,2)
	d_array["realdata"]["to1_3"]+=  np.power(d_array["realdata"]["EIA"],2) +  np.power(d_array["realdata"]["EIB"] -1.1,2)
	d_array["realdata"]["to1_3"]+=  np.power(d_array["realdata"]["EIC"],2) +  np.power(d_array["realdata"]["EID"] -1.1,2)
	d_array["realdata"]["ER"] = (1+8./3)*col_EC - 0.333*col_EI_for_ER
	d_array["realdata"]["test"] = col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18)
	d_array["realdata"]["testEC"] = (d_array["realdata"]["EC1"] + d_array["realdata"]["EC2"])*(col_EIFID -0.16*np.power(abs((1+8./3)*col_EC-0.333*col_EI_for_ER), 1.18))
	d_array["realdata"]["prod"] = d_array["realdata"]["EC1"]*d_array["realdata"]["EC2"]*d_array["realdata"]["EIB"]*d_array["realdata"]["EID"]/(1+d_array["realdata"]["EIA"]*d_array["realdata"]["EIC"])

	d_array["realdata"] = d_array["realdata"][list_variables]


	return d_array


def prepare_data_for_scikit(dict_df):

	"""
	Detail:
		Convert the pandas data frame in dict_df 
		into numpy arrays
		Shuffle the arrays

	Args:
		dict_df (dict)                         = dict with training data as pandas dataframe
		Optional kwarg: d_num_events (dict) = dict to indicate how many training events per class

	Returns:
		x_train, y_train (the feature / class vectors)

	Raises:
		void	 
	"""

	#Merge data frames of different event types. 
	df = pd.concat( [dict_df[event_type] for event_type in dict_df.keys() ], ignore_index = True )

	#Shuffle data frames
	df = df.ix[np.random.permutation(df.index)]

	#Prepare training and test samples
	x      = df.iloc[:,:-3].values
	weight = df.iloc[:,-2].values
	y      = df.iloc[:,-1].values

	return x, weight, y


