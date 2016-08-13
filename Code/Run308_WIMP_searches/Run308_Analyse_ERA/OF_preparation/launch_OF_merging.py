#!/usr/bin/env python

import os
import script_utils as script_utils
import download_data as download_data
import tree_manipulation as tmanip
import get_calibration_ERA_heat as calib
import BaselineEvol as baseline

def get_data(bolo_name, lyon_dir, lyon_ANA_dir):

	"""Download and merged data
	
	Detail:
		Merge the Amp files from ERA

	Args:
		bolo_name    = (str) bolo name 
		lyon_dir     = (str) directory from which ERA files are downloaded
		lyon_ANA_dir = (str) directory from which ANA files are downloaded

		
	Returns:
		void

	Raises:
		void
	"""

	download_data.download_ANA_data(bolo_name, lyon_ANA_dir)
	tmanip.merge_ANA_files(bolo_name)
	# download_data.download_txt_files(bolo_name, lyon_dir)	
	# download_data.download_amp_files(bolo_name, lyon_dir)


def launch_OF_merging(bolo_name, interactif):

	"""Download and merged data
	
	Detail:
		Merge the Amp files from ERA

	Args:
		bolo_name = (str) bolo name 
		interactif = (0/1) 1 to see baseline evol plots
	
	Returns:
		void

	Raises:
		void
	"""

	# #Merge downloaded ERA files 
	# tmanip.merge_amp_files(bolo_name)
	
	# #reliable pop for calibration
	# calib.get_calib_pop(bolo_name)
	# # 1st estimation of calib coeff
	# calib.get_calibration_coefficient(bolo_name)

	# #Baseline evol
	# baseline.get_baseline_evol(bolo_name, "Chal1", interactif)
	# baseline.get_baseline_evol(bolo_name, "Chal2", interactif)

	# # First round of calibration
	# calib_file = "./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff.txt"
	# baseline_file1 = "./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_Chal1.txt"
	# baseline_file2 = "./Text_files/"  + bolo_name + "/" + bolo_name + "_baseline_evol_Chal2.txt"
	# tmanip.combine_merged_files_keep_trigged_calibrate(bolo_name, calib_file)
	# tmanip.add_FWHM(bolo_name, baseline_file1, baseline_file2)

	# # Include ERA to ANA data
	# tmanip.get_all_event_file(bolo_name)
	# tmanip.add_ERA_heat_to_ANA(bolo_name)
	
	# #Refine calibration
	# calib.refine_calibration_coefficient(bolo_name)

	# # Update baselines
	# baseline.refine_baseline(bolo_name)
	# baseline.refine_baseline(bolo_name)

	# # Second round of calibration
	# calib_file = "./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff_refined.txt"
	# baseline_file1 = "./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_Chal1_refined.txt"
	# baseline_file2 = "./Text_files/"  + bolo_name + "/" + bolo_name + "_baseline_evol_Chal2_refined.txt"
	# tmanip.combine_merged_files_keep_trigged_calibrate(bolo_name, calib_file)
	# tmanip.add_FWHM(bolo_name, baseline_file1, baseline_file2)

	# # Include ERA to ANA data
	# tmanip.get_all_event_file(bolo_name)
	# tmanip.add_ERA_heat_to_ANA(bolo_name)

	# #Add PSA information
	# tmanip.add_wavelet_PSA_correct_ion(bolo_name)

bolo_name = "FID837"
lyon_dir  = "/sps/edelweis/EdwRootAna/Run308/" + bolo_name + "/"
lyon_ANA_dir  = "/sps/edelweis/anarun308/flat/anai-merged/"
interactif = 0
get_data(bolo_name, lyon_dir, lyon_ANA_dir)
# launch_OF_merging( bolo_name, interactif)

