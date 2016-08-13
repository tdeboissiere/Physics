from ROOT import *
import numpy as np

def arrange_FID808_event_file(file_name):
	"""Arrange FID808 Pb and Beta files
	
	Detail:
		Get variables of interest: Q, Ion, Heat

	Args:
		file_name (str) = name of event file
		
	Returns:
		void

	Raises:
		void
	"""

	data_types = {"names": 
					("RUN", "SN", "EC1", "EC2", "EIA", "EIB", "EIC", "EID"), 
				"formats": 
					("i", "i", "f", "f", "f", "f", "f", "f")}

	pop_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Text_files/"
	arr = np.loadtxt(pop_path + file_name + ".txt", delimiter=",",  dtype=data_types)

	if "S1" in file_name:
		with open("./Populations/FID808_Run305/" + file_name + "_arranged.txt", "w" ) as f:
			for i in range(arr.shape[0]):
				EC1, EC2, EIA, EIB = arr["EC1"][i], arr["EC2"][i], arr["EIA"][i], arr["EIB"][i]
				est_heat = 0.72*EC1+(1-0.72)*EC2
				est_ion  = 0.68*EIB+(1-0.68)*EIA
				est_ER   = (0.72*EC1+(1-0.72)*EC2)*(1+8/3.)-(0.68*EIB+(1-0.68)*EIA)*0.3333*5.5
				est_Q    = est_ion/est_ER
				f.write(str(est_heat) + "," + str(est_ion) + "," + str(est_ER) + "," + str(est_Q) + "\n")

	if "S2" in file_name:
		with open("./Populations/FID808_Run305/" + file_name + "_arranged.txt", "w" ) as f:
			for i in range(arr.shape[0]):
				EC1, EC2, EIC, EID = arr["EC1"][i], arr["EC2"][i], arr["EIC"][i], arr["EID"][i]
				est_heat = 0.72*EC1+(1-0.72)*EC2
				est_ion  = 0.56*EID+(1-0.56)*EIC
				est_ER   = (0.72*EC1+(1-0.72)*EC2)*(1+8/3.)-(0.56*EID+(1-0.56)*EIC)*0.3333*5.5
				est_Q    = est_ion/est_ER
				f.write(str(est_heat) + "," + str(est_ion) + "," + str(est_ER) + "," + str(est_Q) + "\n")

# file_name = "S1Beta_heatremoved" 
file_name_list = ["S1Beta_heatremoved", "S2Beta_heatremoved", "S1Pb_heatremoved", "S2Pb_heatremoved"]
for file_name in file_name_list:
	arrange_FID808_event_file(file_name)
