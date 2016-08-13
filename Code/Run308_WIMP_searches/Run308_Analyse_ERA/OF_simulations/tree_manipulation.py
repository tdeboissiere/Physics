from os import listdir
from os.path import isfile, join
from ROOT import *
import script_utils as script_utils
import conversion_ANA_SAMBA as conv
import PyROOTPlots as PyRPl
from array import array
import sys
import OF_utilities as OF_ut
from root_numpy import tree2rec
import numpy as np


def merge_simu_amp_files(bolo_name):

	"""Quick description
	
	Detail:
		Merge simulated ERA Amp files

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""
	
	#Define non merged file directory
	data_dir="../Amp_files/" + bolo_name + "/Amp_simu/"

	#List of files in data_dir
	list_files_all = [ f for f in listdir(data_dir) if isfile(join(data_dir,f)) ]        
	
	list_files_Chal1 =sorted([file_Chal1 for file_Chal1 in list_files_all if "Chal1" in file_Chal1])
	list_files_Chal2 =sorted([file_Chal2 for file_Chal2 in list_files_all if "Chal2" in file_Chal2])
	list_files_basic =sorted([file_basic for file_basic in list_files_all if "basic" in file_basic])

	chain_Chal1=TChain("wienerntp_FID837_Chal1","eionbis")
	for partition in list_files_Chal1:
		chain_Chal1.AddFile(data_dir+partition)

	chain_Chal2=TChain("wienerntp_FID837_Chal2","eionbis")
	for partition in list_files_Chal2:
		chain_Chal2.AddFile(data_dir+partition)

	chain_basic=TChain("basicntp_FID837","eionbis")
	for partition in list_files_basic:
		chain_basic.AddFile(data_dir+partition)

	merge_file_path = "../Amp_files_merged/" + bolo_name + "/Amp_simu/"
    #Create directory for data if does not exist
	script_utils.create_directory(merge_file_path)

	merge_file_name_Chal1 = merge_file_path+bolo_name+"_wiener_Chal1.root"
	merge_file_name_Chal2 = merge_file_path+bolo_name+"_wiener_Chal2.root"
	merge_file_name_basic = merge_file_path+bolo_name+"_wiener_basic.root"


	chain_Chal1.Merge(merge_file_name_Chal1)
	chain_Chal2.Merge(merge_file_name_Chal2)
	chain_basic.Merge(merge_file_name_basic)


def combine_merged_files_calibrate(bolo_name, d_channel_calib):

	"""Quick description
	
	Detail:
		Combine Chal1, Chal2, basic files
		Also calibrate

	Args:
		bolo_name       = (str) bolo name 
		d_channel_calib = (str) calib coeff dict
		
	Returns:
		void

	Raises:
		AssertionError()
	"""

	#Get calibration coefficients
	calibcoeff1, calibcoeff2 = float(d_channel_calib["Chal1"]), float(d_channel_calib["Chal2"])
	
	#Get merged files
	path_amp = "../Amp_files_merged/" + bolo_name + "/Amp_simu/" + bolo_name + "_wiener_"
	tChal1, fChal1 = PyRPl.open_ROOT_object(path_amp + "Chal1.root", "wienerntp_" + bolo_name + "_Chal1")
	tChal2, fChal2 = PyRPl.open_ROOT_object(path_amp + "Chal2.root", "wienerntp_" + bolo_name + "_Chal2")
	tbasic, fbasic = PyRPl.open_ROOT_object(path_amp + "basic.root", "basicntp_" + bolo_name )
	
	#new tree
	output_dir = script_utils.create_directory("../Fond_ERA_merged/Fond_simu/")
	outputfile = output_dir + bolo_name + "_merged_simu_tree.root"
	filou=TFile(outputfile, "recreate")
	tree = TTree("t_ERA0", "t_ERA0")

	ORI_HEAT  =  array("f", [ 0 ] )
	EC1_ERA   =  array("f", [ 0 ] )
	EC2_ERA   =  array("f", [ 0 ] )
	CHIC1_ERA =  array("f", [ 0 ] )
	CHIC2_ERA =  array("f", [ 0 ] )
	
	RUN       =  array("i", [ 0 ] )
	SN        =  array("i", [ 0 ] )
	DateSec   =  array("L", [ 0 ] )
	
	tree.Branch("ORI_HEAT", ORI_HEAT, "ORI_HEAT/F")
	tree.Branch("EC1_ERA", EC1_ERA, "EC1_ERA/F")
	tree.Branch("EC2_ERA", EC2_ERA, "EC2_ERA/F")
	tree.Branch("CHIC1_ERA", CHIC1_ERA, "CHIC1_ERA/F")
	tree.Branch("CHIC2_ERA", CHIC2_ERA, "CHIC2_ERA/F")
	tree.Branch("RUN", RUN, "RUN/I")
	tree.Branch("SN", SN, "SN/I")
	tree.Branch("DateSec", DateSec, "DateSec/l")

	nEntries =tChal1.GetEntries()

	for i in range(nEntries):
    
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()
		
		tChal1.GetEntry(i)
		tChal2.GetEntry(i)
		tbasic.GetEntry(i)

		ORI_HEAT[0] = tChal1.OriHeat
		EC1_ERA[0], EC2_ERA[0] = (tChal1.WienerAmpl)/calibcoeff1, (tChal2.WienerAmpl)/calibcoeff2
		CHIC1_ERA[0], CHIC2_ERA[0] = tChal1.WienerChi2, tChal2.WienerChi2
		RUN[0], SN[0] = int(conv.samba2ana(tbasic.Run[:-1])), int(tbasic.SambaNum)
		DateSec[0] = tbasic.DateSec
		tree.Fill()


	tree.Write()
	filou.Close()
	print



def add_wavelet_PSA_simu(bolo_name):

	"""Add wavelet PSA to simu event tree
	
	Detail:
		void

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""

	#Get merged ERA files
	path_ERA = "../Fond_ERA_merged/Fond_simu/"
	t_ERA, f_ERA = PyRPl.open_ROOT_object(path_ERA + bolo_name + "_merged_simu_tree.root", "t_ERA0")

	#Load event dictionary
	d_PSA_events = OF_ut.load_PSA_event_dict(bolo_name, "simu_event")	

	#Create new output tree
	output_file_name = path_ERA + bolo_name + "_PSA_merged_simu_tree.root"
	filou=TFile(output_file_name, "recreate")
	tree = TTree("t_merged", "t_merged")

	ORI_HEAT  =  array("f", [ 0 ] )
	EC1_ERA   =  array("f", [ 0 ] )
	EC2_ERA   =  array("f", [ 0 ] )
	CHIC1_ERA =  array("f", [ 0 ] )
	CHIC2_ERA =  array("f", [ 0 ] )
	
	RUN       =  array("i", [ 0 ] )
	SN        =  array("i", [ 0 ] )
	DateSec   =  array("L", [ 0 ] )
	
	pearsA    =  array("f", [ 0 ] )
	dtw_valA  =  array("f", [ 0 ] )

	tree.Branch("ORI_HEAT", ORI_HEAT, "ORI_HEAT/F")
	tree.Branch("EC1_ERA", EC1_ERA, "EC1_ERA/F")
	tree.Branch("EC2_ERA", EC2_ERA, "EC2_ERA/F")
	tree.Branch("CHIC1_ERA", CHIC1_ERA, "CHIC1_ERA/F")
	tree.Branch("CHIC2_ERA", CHIC2_ERA, "CHIC2_ERA/F")
	tree.Branch("RUN", RUN, "RUN/I")
	tree.Branch("SN", SN, "SN/I")
	tree.Branch("DateSec", DateSec, "DateSec/l")

	tree.Branch("pearsA", pearsA, "pearsA/F")
	tree.Branch("dtw_valA", dtw_valA, "dtw_valA/F")

	nEntries = t_ERA.GetEntries()

	for i in range(nEntries):
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		t_ERA.GetEntry(i)
		RUN[0], SN[0] = int(t_ERA.RUN), int(t_ERA.SN)

		try:
			pearsA[0]    = d_PSA_events[str(RUN[0]) + "_" + str(SN[0])][0]
			dtw_valA[0]  = d_PSA_events[str(RUN[0]) + "_" + str(SN[0])][1]
			
			ORI_HEAT[0]  =  t_ERA.ORI_HEAT
			EC1_ERA[0]   =  t_ERA.EC1_ERA
			EC2_ERA[0]   =  t_ERA.EC2_ERA
			CHIC1_ERA[0] =  t_ERA.CHIC1_ERA
			CHIC2_ERA[0] =  t_ERA.CHIC2_ERA			
			DateSec[0]   =  t_ERA.DateSec			

			tree.Fill()

		except KeyError:
			pass

	print
	tree.Write()
	filou.Close()