from ROOT import *
import numpy as np
import script_utils as script_utils
import os

def arrange_FID_event_file(file_name, bolo_name):
	"""Arrange FIDXXX Pb and Beta files
	
	Detail:
		Get variables of interest: Q, Ion, Heat

	Args:
		file_name (str) = name of event file
		bolo_name (str) = bolometer_name
		
	Returns:
		void

	Raises:
		void
	"""

	#Load estimators
	estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
	estimator_file_name =bolo_name + "_estimators.txt" 
	d_estimator         ={}
	assert(os.path.isfile(estimator_path_name + estimator_file_name) )
	with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
		list_estimator_lines = [elem.rstrip().split(",") for elem in estimator_file.readlines()]
		for line in list_estimator_lines:
			d_estimator[line[0]] = line[1]



	data_types = {"names": 
					("EC1_ERA", "EC2_ERA", "EIA", "EIB", "EIC", "EID"), 
				"formats": 
					("f", "f", "f", "f", "f", "f")}

	pop_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Populations/Pop_for_scaling/"
	arr = np.loadtxt(pop_path + bolo_name + file_name + "_full_info.txt", delimiter=",",  dtype=data_types)

	if "S1" in file_name:
		with open("./Populations/" +bolo_name + "/" + bolo_name + file_name + "_arranged.txt", "w" ) as f:
			for i in range(arr.shape[0]):
				EC1_ERA, EC2_ERA, EIA, EIB = arr["EC1_ERA"][i], arr["EC2_ERA"][i], arr["EIA"][i], arr["EIB"][i]
				est_heat = float(d_estimator["HEAT"][:5])*EC1_ERA+(1-float(d_estimator["HEAT"][:5]))*EC2_ERA
				est_ion  = float(d_estimator["S1"][:5])*EIB+(1-float(d_estimator["S1"][:5]))*EIA
				est_ER   = est_heat*(1+8/3.)-est_ion*0.3333*5.5
				est_Q    = est_ion/est_ER
				f.write(str(est_heat) + "," + str(est_ion) + "," + str(est_ER) + "," + str(est_Q) + "\n")

	if "S2" in file_name:
		with open("./Populations/" +bolo_name + "/" + bolo_name + file_name + "_arranged.txt", "w" ) as f:
			for i in range(arr.shape[0]):
				EC1_ERA, EC2_ERA, EIC, EID = arr["EC1_ERA"][i], arr["EC2_ERA"][i], arr["EIC"][i], arr["EID"][i]
				est_heat = float(d_estimator["HEAT"][:5])*EC1_ERA+(1-float(d_estimator["HEAT"][:5]))*EC2_ERA
				est_ion  = float(d_estimator["S2"][:5])*EID+(1-float(d_estimator["S1"][:5]))*EIC
				est_ER   = est_heat*(1+8/3.)-est_ion*0.3333*5.5
				est_Q    = est_ion/est_ER
				f.write(str(est_heat) + "," + str(est_ion) + "," + str(est_ER) + "," + str(est_Q) + "\n")

# file_name = "S1Beta_heatremoved" 
file_name_list = ["_S1Beta", "_S2Beta", "_S1Pb", "_S2Pb"]
bolo_name_list = ["FID837"]
for bolo_name in bolo_name_list:
	for file_name in file_name_list:
		arrange_FID_event_file(file_name, bolo_name)
