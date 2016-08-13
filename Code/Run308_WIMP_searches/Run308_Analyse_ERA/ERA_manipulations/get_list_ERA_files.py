import glob
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

	text_path = script_utils.create_directory("./" + bolo_name + "/Text_files/")
	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*.root*"))
	np.savetxt(text_path + bolo_name + "_ERA_files.txt", np.array(list_ERA_files), fmt="%s")	

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*og*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_og.txt", np.array(list_ERA_files), fmt="%s")

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*oh*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_oh.txt", np.array(list_ERA_files), fmt="%s")

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*oi*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_oi.txt", np.array(list_ERA_files), fmt="%s")

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*oj*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_oj.txt", np.array(list_ERA_files), fmt="%s")

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*ok*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_ok.txt", np.array(list_ERA_files), fmt="%s")

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*ol*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_ol.txt", np.array(list_ERA_files), fmt="%s")

	list_ERA_files = sorted(glob.glob("/sps/edelweis/EdwRootAna/Run308/FID837/Traces/*pa*"))
	np.savetxt(text_path + bolo_name + "_ERA_files_pa.txt", np.array(list_ERA_files), fmt="%s")


get_list_runs("FID837")