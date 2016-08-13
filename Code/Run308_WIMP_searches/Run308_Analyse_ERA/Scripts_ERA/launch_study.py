#!/usr/bin/env python

import script_utils as script_utils
import paramiko, getpass, os, sys
from os import listdir
from os.path import isfile, join
import get_polar as get_polar
import get_run_start_end_time as get_run_start_end_time
import get_threshold_hist as get_threshold_hist
import get_run_duration as get_run_duration
import get_standard_cuts as get_standard_cuts
import get_cut_file as get_cut_file
import get_FWHM as get_FWHM
import get_estimators as get_estimators
import delete_specific_files as delete_files


def download_ANA_data(bolo_name, lyon_ANA_dir):

	"""Quick description
	
	Detail:
		Download ANA files

	Args:
		bolo_name = (str) bolo name 
		lyon_ANA_dir  = (str) directory from which ANA files are downloaded

		
	Returns:
		void

	Raises:
		void
	"""

	data_dir = "../Fond/"

	ssh = paramiko.SSHClient() 
	ssh.load_host_keys(os.path.expanduser(os.path.join('~', '.ssh', 'known_hosts')))
	pwd = getpass.getpass('Password please:   ')  # command line prompt without echo
	ssh.connect('ccage.in2p3.fr', username='tmain', password=pwd)
	stdin, stdout, stderr = ssh.exec_command("cd " + lyon_ANA_dir + " && ls *"+ bolo_name[3:]+ "*")

	#list containing the file names
	liste=stdout.read().splitlines()	

	#print  stderr.read.splitlines()
	sftp = ssh.open_sftp()
	for bolo_file in liste:
	    if ('-fond-' in bolo_file and 'root' in bolo_file and 'empty' not in bolo_file and 'bias' not in bolo_file): 
			script_utils.print_utility(script_utils.COL("Adding file " + bolo_file, "blue"))
			sftp.get(lyon_ANA_dir+bolo_file,data_dir+bolo_file)
	sftp.close()
	ssh.close()

def merge_files(bolo_name):

	"""Quick description
	
	Detail:
		Merge ANA files

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""
	
	#Define non merged file directory
	data_dir="../Fond/"

    #List of files in data_dir
	list_files_all = [ f for f in listdir(data_dir) if isfile(join(data_dir,f)) ]        

	list_files_bolo =sorted([file_bolo for file_bolo in list_files_all if bolo_name.lower() in file_bolo])
	chain=TChain("data","eionbis")
	for partition in list_files_bolo:
		print partition
		chain.AddFile(data_dir+partition)

	fusion_file_path = '../Fond_ERA_merged/'
	fusion_file_name = fusion_file_path+bolo_name.upper()+"_fond.root"
	chain.Merge(fusion_file_name)

def launch_study(bolo_name, data_dir, tree_name = "data"):

	"""Launch the analysis scripts
	
    Detail:
		Launch the general analysis for the given detector
		Get polar, baselines, estimators, duration and general cuts

	Args:
		bolo_name = (str) bolometer type
		data_dir  = (str) the ROOT data tree directory
		tree_name = (str) the ROOT data tree name

	Returns:
		void

	Raises:
		void
	"""

	# download_ANA_data(bolo_name, lyon_ANA_dir)
	
	# source_dir ="/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name
	# list_pattern_to_delete = [".eps", ".png", ".root", ".txt"]
	
	# # delete_files.delete_files(source_dir,  list_pattern_to_delete)
	
	# #Create directory for cut files if it does not exist
	# script_utils.create_directory("../Cut_files/")
	
	# get_polar.get_polar(bolo_name, data_dir, tree_name )
	# get_standard_cuts.get_standard_cuts(bolo_name, data_dir, tree_name)
	# get_cut_file.get_cut_file(bolo_name, data_dir, tree_name)
	# get_run_start_end_time.get_run_start_end_time(bolo_name, data_dir, tree_name)
	# get_threshold_hist.get_threshold_hist(bolo_name, data_dir, tree_name)
	# get_threshold_hist.get_threshold_hist_for_plot(bolo_name, data_dir, tree_name)
	# get_run_duration.get_run_duration_KTH_cut(bolo_name, data_dir, tree_name)
	# for heat_fraction in [0.1,0.2,0.3,0.4,0.5,0.8,1]:
	# 	get_run_duration.get_run_duration_KTH_heat_cut(bolo_name, data_dir, tree_name, heat_fraction)
	# get_FWHM.get_FWHM(bolo_name, data_dir, tree_name)
	# get_FWHM.get_FWHM_hist(bolo_name, data_dir, tree_name)
	get_FWHM.get_FWHM_time_hist_for_BDT(bolo_name, data_dir, tree_name)
	# get_FWHM.get_FWHM_time_hist_for_heatonly_for_BDT(bolo_name, data_dir, tree_name)
	# get_estimators.get_estimators(bolo_name)
	return

bolo_name = "FID837"
data_dir  = "../Fond_ERA_merged/"
lyon_ANA_dir  = "/sps/edelweis/anarun308/flat/anai-merged/"
launch_study( bolo_name, data_dir)
