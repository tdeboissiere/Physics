import glob
import conversion_ANA_SAMBA as conv
import numpy as np
import script_utils as script_utils

def get_list_runs(bolo_name):
	"""Quick description
	
	Detail:
		get_list_runs

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""

	text_path = script_utils.create_directory("./Text_files/" + bolo_name + "/")

	list_runs_samba_prelim = glob.glob("../Amp_files/FID837/*Chal1*")
	list_runs_ana = sorted([conv.samba2ana(elem[33:41]) for elem in list_runs_samba_prelim])
	list_runs_samba = [conv.ana2samba(int(elem)) for elem in list_runs_ana]

	np.savetxt(text_path + bolo_name + "_runs_ana.txt", np.array(list_runs_ana), fmt="%i")	
	np.savetxt(text_path + bolo_name + "_runs_samba.txt", np.array(list_runs_samba), fmt="%s")


get_list_runs("FID837")