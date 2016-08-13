#!/usr/bin/env python

import os
import script_utils as script_utils
import download_data as download_data
import OF_utilities as OF_ut
import get_run_by_run_simu_event_tree as get_simu_events
import tree_manipulation as tmanip

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

	# download_data.download_ANA_data(bolo_name, lyon_ANA_dir)
	# download_data.download_txt_files(bolo_name, lyon_dir)	
	download_data.download_spectra_files(bolo_name, lyon_dir)


def launch_OF_simulations(bolo_name):

	"""OF simulations
	
	Detail:
		Launch pulse simulations
		apply OF and wavelet processing

	Args:
		bolo_name = (str) bolo name 
	
	Returns:
		void

	Raises:
		void
	"""

	#Load calibration coefficients
	calib_file = "../OF_preparation/Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff.txt"
	d_channel_calib = {}
	assert(os.path.isfile(calib_file))
	with open(calib_file, "r") as fc:
		lines = fc.readlines()
		calibcoeff1, calibcoeff2 = lines[1].rstrip().split(",")[1], lines[1].rstrip().split(",")[2] 
		d_channel_calib["Chal1"] = float(calibcoeff1)
		d_channel_calib["Chal2"] = float(calibcoeff2)

	#Get run start and end index (for notrigged files)
	# OF_ut.get_run_start_and_end_index(bolo_name,"all_notrigged")

	# #Launch pulse/event simulations
	# get_simu_events.create_event_simu_tree(d_channel_calib)

	# #Merge simulated ERA trees
	# tmanip.merge_simu_amp_files(bolo_name)
	
	# # Combine Chal1, Chal2, basic and calibrate 
	# tmanip.combine_merged_files_calibrate(bolo_name, d_channel_calib)

	# Add PSA information 
	tmanip.add_wavelet_PSA_simu(bolo_name)

bolo_name = "FID837"
lyon_dir  = "/sps/edelweis/EdwRootAna/Run308/" + bolo_name + "/"
lyon_ANA_dir  = "/sps/edelweis/EdwRootAna/Run308/" + bolo_name + "/"
# get_data(bolo_name, lyon_dir, lyon_ANA_dir)
launch_OF_simulations( bolo_name)

