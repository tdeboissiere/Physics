from os import listdir
from os.path import isfile, join
from ROOT import *
import script_utils as script_utils
import new_conversion_ANA_SAMBA as conv
import PyROOTPlots as PyRPl
from array import array
import sys
import OF_utilities as OF_ut
from root_numpy import tree2rec
import numpy as np


def merge_ANA_files(bolo_name):

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

	fusion_file_path = '../Fond_merged/'
	fusion_file_name = fusion_file_path+bolo_name.upper()+"_fond.root"
	chain.Merge(fusion_file_name)



def merge_amp_files(bolo_name):

	"""Quick description
	
	Detail:
		Merge ERA Amp files

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""
	
	#Define non merged file directory
	data_dir="../Amp_files/" + bolo_name + "/"

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

	merge_file_path = "../Amp_files_merged/" + bolo_name + "/"
    #Create directory for data if does not exist
	script_utils.create_directory(merge_file_path)

	merge_file_name_Chal1 = merge_file_path+bolo_name+"_wiener_Chal1.root"
	merge_file_name_Chal2 = merge_file_path+bolo_name+"_wiener_Chal2.root"
	merge_file_name_basic = merge_file_path+bolo_name+"_wiener_basic.root"


	chain_Chal1.Merge(merge_file_name_Chal1)
	chain_Chal2.Merge(merge_file_name_Chal2)
	chain_basic.Merge(merge_file_name_basic)

def combine_merged_files_keep_trigged_calibrate(bolo_name, calib_file):

	"""Quick description
	
	Detail:
		Combine Chal1, Chal2, basic files
		Keep trigged events only

	Args:
		bolo_name = (str) bolo name
		calib_file = (str) calib file name 
		
	Returns:
		void

	Raises:
		AssertionError()
	"""

	#Get calibration coefficients
	calibcoeff1, calibcoeff2 = 0,0
	with open(calib_file, "r") as fc:
		lines = fc.readlines()
		calibcoeff1, calibcoeff2 = float(lines[1].rstrip().split(",")[1]), float(lines[1].rstrip().split(",")[2])
	
	#Get merged files
	path_amp = "../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_"
	tChal1, fChal1 = PyRPl.open_ROOT_object(path_amp + "Chal1.root", "wienerntp_" + bolo_name + "_Chal1")
	tChal2, fChal2 = PyRPl.open_ROOT_object(path_amp + "Chal2.root", "wienerntp_" + bolo_name + "_Chal2")
	tbasic, fbasic = PyRPl.open_ROOT_object(path_amp + "basic.root", "basicntp_" + bolo_name )
	
	#new tree
	outputfile = path_amp + "merged_trigged_tree.root"
	filou=TFile(outputfile, "recreate")
	tree = TTree("t_ERA0", "t_ERA0")
	
	EC1_ERA   =  array("f", [ 0 ] )
	EC2_ERA   =  array("f", [ 0 ] )
	CHIC1_ERA =  array("f", [ 0 ] )
	CHIC2_ERA =  array("f", [ 0 ] )
	
	Run       =  array("i", [ 0 ] )
	SambaNum  =  array("i", [ 0 ] )
	DateSec   =  array("L", [ 0 ] )
	
	tree.Branch("EC1_ERA", EC1_ERA, "EC1_ERA/F")
	tree.Branch("EC2_ERA", EC2_ERA, "EC2_ERA/F")
	tree.Branch("CHIC1_ERA", CHIC1_ERA, "CHIC1_ERA/F")
	tree.Branch("CHIC2_ERA", CHIC2_ERA, "CHIC2_ERA/F")
	tree.Branch("Run", Run, "Run/I")
	tree.Branch("SambaNum", SambaNum, "SambaNum/I")
	tree.Branch("DateSec", DateSec, "DateSec/l")

	nEntries =tChal1.GetEntries()

	for i in range(nEntries):
    
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()
		
		tChal1.GetEntry(i)
		tChal2.GetEntry(i)
		tbasic.GetEntry(i)

		try:
			assert(tbasic.IsBoloTrigger)
			EC1_ERA[0], EC2_ERA[0] = (tChal1.WienerAmpl)/calibcoeff1, (tChal2.WienerAmpl)/calibcoeff2
			CHIC1_ERA[0], CHIC2_ERA[0] = tChal1.WienerChi2, tChal2.WienerChi2
			Run[0], SambaNum[0] = int(conv.samba2ana(tbasic.Run[:-1])), int(tbasic.SambaNum)
			DateSec[0] = tbasic.DateSec
			tree.Fill()

		except AssertionError:
			pass


	tree.Write()
	filou.Close()
	print




def add_FWHM(bolo_name, baseline_file1, baseline_file2):

	"""Quick description
	
	Detail:
		Add the FWHM (heat) to the merged files 
		(already calibrated, only trigged events)

	Args:
		bolo_name = (str) bolo name 
		baseline_fileXX = (str) baseline file name 
		
	Returns:
		void

	Raises:
		void
	"""


	#Get required dictionaries
	d_baseline1 = OF_ut.load_baseline_dict(bolo_name, baseline_file1, "Chal1")
	d_baseline2 = OF_ut.load_baseline_dict(bolo_name, baseline_file2, "Chal2")
	d_periods   = OF_ut.load_period_dict(bolo_name)

	#Get merged files
	path_amp = "../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_"
	t_ERA_old, fERA_old = PyRPl.open_ROOT_object(path_amp + "merged_trigged_tree.root", "t_ERA0")
	
	#new tree
	outputfile = path_amp + "merged_trigged_fwhm_tree.root"
	filou=TFile(outputfile, "recreate")
	tree = TTree("t_ERA1", "t_ERA1")
	
	EC1_ERA   =  array("f", [ 0 ] )
	EC2_ERA   =  array("f", [ 0 ] )
	CHIC1_ERA =  array("f", [ 0 ] )
	CHIC2_ERA =  array("f", [ 0 ] )
	FWC1_ERA  =  array("f", [ 0 ] )
	FWC2_ERA  =  array("f", [ 0 ] )

	Run       =  array("i", [ 0 ] )
	SambaNum  =  array("i", [ 0 ] )
	DateSec   =  array("L", [ 0 ] )
	
	tree.Branch("EC1_ERA", EC1_ERA, "EC1_ERA/F")
	tree.Branch("EC2_ERA", EC2_ERA, "EC2_ERA/F")
	tree.Branch("CHIC1_ERA", CHIC1_ERA, "CHIC1_ERA/F")
	tree.Branch("CHIC2_ERA", CHIC2_ERA, "CHIC2_ERA/F")
	tree.Branch("FWC1_ERA", FWC1_ERA, "FWC1_ERA/F")
	tree.Branch("FWC2_ERA", FWC2_ERA, "FWC2_ERA/F")
	tree.Branch("Run", Run, "Run/I")
	tree.Branch("SambaNum", SambaNum, "SambaNum/I")
	tree.Branch("DateSec", DateSec, "DateSec/l")

	nEntries =t_ERA_old.GetEntries()

	for i in range(nEntries):
    
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()
		
		t_ERA_old.GetEntry(i)

		EC1_ERA[0], EC2_ERA[0]     = t_ERA_old.EC1_ERA, t_ERA_old.EC2_ERA
		CHIC1_ERA[0], CHIC2_ERA[0] = t_ERA_old.CHIC1_ERA, t_ERA_old.CHIC2_ERA
		Run[0], SambaNum[0]        = t_ERA_old.Run, t_ERA_old.SambaNum
		DateSec[0]                 = t_ERA_old.DateSec
		run_samba                  = conv.ana2samba(Run[0])
		FWC1_ERA[0], FWC2_ERA[0]   = OF_ut.find_event_baseline(bolo_name, run_samba, DateSec[0], d_periods, d_baseline1, d_baseline2)
		tree.Fill()

	tree.Write()
	filou.Close()




def get_all_event_file(bolo_name):

	"""Get all event file for given bolo
	
	Detail:
		Store the index, Run, SambaNum of the event

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""

	#Get merged files
	path_amp = "../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_"
	path_fond = "../Fond_merged/" + bolo_name 
	
	tERA, fERA = PyRPl.open_ROOT_object(path_amp + "merged_trigged_fwhm_tree.root" , "t_ERA1")
	tANA, fANA = PyRPl.open_ROOT_object(path_fond + "_fond.root", "data")

	# file_ERA = path_amp + "merged_trigged_fwhm_tree.root"
	# arr_ERA  = root2array(file_ERA, "t_ERA1")
	arr_ERA  = tree2rec(tERA, branches = ["Run", "SambaNum"])

	# file_ANA = path_fond + "_fond.root"
	# arr_ANA  = root2array(file_ANA, "data")
	arr_ANA  = tree2rec(tANA, branches = ["RUN", "SN"])

	np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name + "_all_event_ANA.txt", arr_ANA, fmt = "%i", delimiter = ",")
	np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name + "_all_event_ERA.txt", arr_ERA, fmt = "%i", delimiter = ",")

def add_ERA_heat_to_ANA(bolo_name):

	"""Add ERA heat to ANA NTP
	
	Detail:
		void

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""

	#Get merged files
	path_amp = "../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_"
	tERA, fERA = PyRPl.open_ROOT_object(path_amp + "merged_trigged_fwhm_tree.root", "t_ERA1")

	path_fond = "../Fond_merged/" + bolo_name 
	tANA, fANA = PyRPl.open_ROOT_object(path_fond + "_fond.root", "data")

	#Load event dictionaries
	d_proc_ERA = OF_ut.load_processing_dict(bolo_name, "ERA")	

	# print d_proc_ERA.keys()[:10], d_proc_ERA.values()[:10]
	# raw_input()

	#new tree
	path_output = script_utils.create_directory("../Fond_ERA_merged/")
	output_file_name = path_output + bolo_name + "_fond.root"
	filou=TFile(output_file_name, "recreate")
	tree = TTree("t_merged", "t_merged")
	
	EC1_ERA   =  array("f", [ 0 ] )
	EC2_ERA   =  array("f", [ 0 ] )
	CHIC1_ERA =  array("f", [ 0 ] )
	CHIC2_ERA =  array("f", [ 0 ] )
	FWC1_ERA  =  array("f", [ 0 ] )
	FWC2_ERA  =  array("f", [ 0 ] )
	
	EC1       =  array("f", [ 0 ] )
	EC2       =  array("f", [ 0 ] )
	EIA       =  array("f", [ 0 ] )
	EIB       =  array("f", [ 0 ] )
	EIC       =  array("f", [ 0 ] )
	EID       =  array("f", [ 0 ] )
	CHIC1     =  array("f", [ 0 ] )
	CHIC2     =  array("f", [ 0 ] )
	CHIA      =  array("f", [ 0 ] )
	CHIB      =  array("f", [ 0 ] )
	CHIC      =  array("f", [ 0 ] )
	CHID      =  array("f", [ 0 ] )
	FWC1      =  array("f", [ 0 ] )
	FWC2      =  array("f", [ 0 ] )
	FWIA      =  array("f", [ 0 ] )
	FWIB      =  array("f", [ 0 ] )
	FWIC      =  array("f", [ 0 ] )
	FWID      =  array("f", [ 0 ] )
	
	RCIA      =  array("f", [ 0 ] )
	RCIB      =  array("f", [ 0 ] )
	RCIC      =  array("f", [ 0 ] )
	RCID      =  array("f", [ 0 ] )
	OTR       =  array("f", [ 0 ] )

	TAFT      =  array("f", [ 0 ] )
	TBEF      =  array("f", [ 0 ] )
	SDEL      =  array("f", [ 0 ] )
	KTH       =  array("f", [ 0 ] )
	VOLT      =  array("f", [ 0 ] )
	VVET      =  array("f", [ 0 ] )
	
	RUN       =  array("i", [ 0 ] )
	SN        =  array("i", [ 0 ] )
	DateSec   =  array("L", [ 0 ] )
	
	#ERA variables
	tree.Branch("EC1_ERA", EC1_ERA, "EC1_ERA/F")
	tree.Branch("EC2_ERA", EC2_ERA, "EC2_ERA/F")
	tree.Branch("CHIC1_ERA", CHIC1_ERA, "CHIC1_ERA/F")
	tree.Branch("CHIC2_ERA", CHIC2_ERA, "CHIC2_ERA/F")
	tree.Branch("FWC1_ERA", FWC1_ERA, "FWC1_ERA/F")
	tree.Branch("FWC2_ERA", FWC2_ERA, "FWC2_ERA/F")
	
	#ANA variables
	tree.Branch("EC1", EC1, "EC1/F")
	tree.Branch("EC2", EC2, "EC2/F")
	tree.Branch("EIA", EIA, "EIA/F")
	tree.Branch("EIB", EIB, "EIB/F")
	tree.Branch("EIC", EIC, "EIC/F")
	tree.Branch("EID", EID, "EID/F")
	tree.Branch("CHIC1", CHIC1, "CHIC1/F")
	tree.Branch("CHIC2", CHIC2, "CHIC2/F")
	tree.Branch("CHIA", CHIA, "CHIA/F")
	tree.Branch("CHIB", CHIB, "CHIB/F")
	tree.Branch("CHIC", CHIC, "CHIC/F")
	tree.Branch("CHID", CHID, "CHID/F")
	tree.Branch("FWC1", FWC1, "FWC1/F")
	tree.Branch("FWC2", FWC2, "FWC2/F")
	tree.Branch("FWIA", FWIA, "FWIA/F")
	tree.Branch("FWIB", FWIB, "FWIB/F")
	tree.Branch("FWIC", FWIC, "FWIC/F")
	tree.Branch("FWID", FWID, "FWID/F")

	tree.Branch("RCIA", RCIA, "RCIA/F")
	tree.Branch("RCIB", RCIB, "RCIB/F")
	tree.Branch("RCIC", RCIC, "RCIC/F")
	tree.Branch("RCID", RCID, "RCID/F")
	tree.Branch("OTR", OTR, "OTR/F")

	tree.Branch("TBEF", TBEF, "TBEF/F")
	tree.Branch("TAFT", TAFT, "TAFT/F")
	tree.Branch("SDEL", SDEL, "SDEL/F")
	tree.Branch("KTH", KTH, "KTH/F")
	tree.Branch("VOLT", VOLT, "VOLT/F")
	tree.Branch("VVET", VVET, "VVET/F")
	
	tree.Branch("RUN", RUN, "RUN/I")
	tree.Branch("SN", SN, "SN/I")
	tree.Branch("DateSec", DateSec, "DateSec/l")	

	nEntries = tANA.GetEntries()

	for i in range(nEntries):
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		tANA.GetEntry(i)
		RUN[0], SN[0] = int(tANA.RUN), int(tANA.SN)

		try:
			ERA_index = d_proc_ERA[str(RUN[0]) + "_" + str(SN[0])]

			tERA.GetEntry(ERA_index)
			
			EC1_ERA[0]   =  tERA.EC1_ERA
			EC2_ERA[0]   =  tERA.EC2_ERA
			CHIC1_ERA[0] =  tERA.CHIC1_ERA
			CHIC2_ERA[0] =  tERA.CHIC2_ERA
			FWC1_ERA[0]  =  tERA.FWC1_ERA
			FWC2_ERA[0]  =  tERA.FWC2_ERA
			
			EC1[0]       =  tANA.EC1
			EC2[0]       =  tANA.EC2
			EIA[0]       =  tANA.EIA
			EIB[0]       =  tANA.EIB
			EIC[0]       =  tANA.EIC
			EID[0]       =  tANA.EID
			CHIC1[0]     =  tANA.CHIC1
			CHIC2[0]     =  tANA.CHIC2
			CHIA[0]      =  tANA.CHIA
			CHIB[0]      =  tANA.CHIB
			CHIC[0]      =  tANA.CHIC
			CHID[0]      =  tANA.CHID
			FWC1[0]      =  tANA.FWC1
			FWC2[0]      =  tANA.FWC2
			FWIA[0]      =  tANA.FWIA
			FWIB[0]      =  tANA.FWIB
			FWIC[0]      =  tANA.FWIC
			FWID[0]      =  tANA.FWID
			
			RCIA[0]      =  tANA.RCIA
			RCIB[0]      =  tANA.RCIB
			RCIC[0]      =  tANA.RCIC
			RCID[0]      =  tANA.RCID
			OTR[0]       = tANA.OTR
			
			TBEF[0]      =  tANA.TBEF
			TAFT[0]      =  tANA.TAFT
			SDEL[0]      =  tANA.SDEL
			KTH[0]       =  tANA.KTH
			VOLT[0]      =  tANA.VOLT
			VVET[0]      =  tANA.VVET
			
			DateSec[0]   =  tERA.DateSec			

			tree.Fill()

		except KeyError:
			pass

	print
	tree.Write()
	filou.Close()


def add_wavelet_PSA_correct_ion(bolo_name):

	"""Add wavelet PSA to merged ERA ANA tree and correct ionisation
	
	Detail:
		void

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""

	#Get ERA/ANA merged _ files
	path_ERA_ANA = "../Fond_ERA_merged/"
	t_ERA_ANA, f_ERA_ANA = PyRPl.open_ROOT_object(path_ERA_ANA + bolo_name + "_fond.root", "t_merged")

	# #Get Wavelet PSA tree
	# path_wave = "../Wavelet/ROOT_files/" + bolo_name + "/" 
	# t_wave, f_wave = PyRPl.open_ROOT_object(path_wave +  bolo_name + "_all_event_wavelet_PSA_tree.root", "t_wavelet")

	#Get the correction to ionisation
	corr_file = "./Text_files/" + bolo_name + "/" + bolo_name + "_ion_calibcorrection.txt"
	corr_IA, corr_IB, corr_IC, corr_ID=1,1,1,1
	assert(isfile(corr_file))
	with open(corr_file, "r") as fc:
		lines = fc.readlines()
		corr_IA, corr_IB = float(lines[1].rstrip().split(",")[1]), float(lines[1].rstrip().split(",")[2])
		corr_IC, corr_ID = float(lines[1].rstrip().split(",")[3]), float(lines[1].rstrip().split(",")[4])

	#Load event dictionary
	d_PSA_events = OF_ut.load_PSA_event_dict(bolo_name, "all_event")	

	#Create new output tree
	path_output = script_utils.create_directory("../Fond_ERA_merged/")
	output_file_name = path_output + bolo_name + "_lowmass_fond.root"
	filou=TFile(output_file_name, "recreate")
	tree = TTree("t_merged", "t_merged")

	EC1_ERA   =  array("f", [ 0 ] )
	EC2_ERA   =  array("f", [ 0 ] )
	CHIC1_ERA =  array("f", [ 0 ] )
	CHIC2_ERA =  array("f", [ 0 ] )
	FWC1_ERA  =  array("f", [ 0 ] )
	FWC2_ERA  =  array("f", [ 0 ] )
	
	EC1       =  array("f", [ 0 ] )
	EC2       =  array("f", [ 0 ] )
	EIA       =  array("f", [ 0 ] )
	EIB       =  array("f", [ 0 ] )
	EIC       =  array("f", [ 0 ] )
	EID       =  array("f", [ 0 ] )
	CHIC1     =  array("f", [ 0 ] )
	CHIC2     =  array("f", [ 0 ] )
	CHIA      =  array("f", [ 0 ] )
	CHIB      =  array("f", [ 0 ] )
	CHIC      =  array("f", [ 0 ] )
	CHID      =  array("f", [ 0 ] )
	FWC1      =  array("f", [ 0 ] )
	FWC2      =  array("f", [ 0 ] )
	FWIA      =  array("f", [ 0 ] )
	FWIB      =  array("f", [ 0 ] )
	FWIC      =  array("f", [ 0 ] )
	FWID      =  array("f", [ 0 ] )
	
	RCIA      =  array("f", [ 0 ] )
	RCIB      =  array("f", [ 0 ] )
	RCIC      =  array("f", [ 0 ] )
	RCID      =  array("f", [ 0 ] )
	OTR       =  array("f", [ 0 ] )
	
	TAFT      =  array("f", [ 0 ] )
	TBEF      =  array("f", [ 0 ] )
	SDEL      =  array("f", [ 0 ] )
	KTH       =  array("f", [ 0 ] )
	VOLT      =  array("f", [ 0 ] )
	VVET      =  array("f", [ 0 ] )
	
	RUN       =  array("i", [ 0 ] )
	SN        =  array("i", [ 0 ] )
	DateSec   =  array("L", [ 0 ] )
	
	pearsA    = array("f", [ 0 ])
	dtw_valA = array("f", [ 0 ])

	#ERA variables
	tree.Branch("EC1_ERA", EC1_ERA, "EC1_ERA/F")
	tree.Branch("EC2_ERA", EC2_ERA, "EC2_ERA/F")
	tree.Branch("CHIC1_ERA", CHIC1_ERA, "CHIC1_ERA/F")
	tree.Branch("CHIC2_ERA", CHIC2_ERA, "CHIC2_ERA/F")
	tree.Branch("FWC1_ERA", FWC1_ERA, "FWC1_ERA/F")
	tree.Branch("FWC2_ERA", FWC2_ERA, "FWC2_ERA/F")
	
	#ANA variables
	tree.Branch("EC1", EC1, "EC1/F")
	tree.Branch("EC2", EC2, "EC2/F")
	tree.Branch("EIA", EIA, "EIA/F")
	tree.Branch("EIB", EIB, "EIB/F")
	tree.Branch("EIC", EIC, "EIC/F")
	tree.Branch("EID", EID, "EID/F")
	tree.Branch("CHIC1", CHIC1, "CHIC1/F")
	tree.Branch("CHIC2", CHIC2, "CHIC2/F")
	tree.Branch("CHIA", CHIA, "CHIA/F")
	tree.Branch("CHIB", CHIB, "CHIB/F")
	tree.Branch("CHIC", CHIC, "CHIC/F")
	tree.Branch("CHID", CHID, "CHID/F")
	tree.Branch("FWC1", FWC1, "FWC1/F")
	tree.Branch("FWC2", FWC2, "FWC2/F")
	tree.Branch("FWIA", FWIA, "FWIA/F")
	tree.Branch("FWIB", FWIB, "FWIB/F")
	tree.Branch("FWIC", FWIC, "FWIC/F")
	tree.Branch("FWID", FWID, "FWID/F")

	tree.Branch("RCIA", RCIA, "RCIA/F")
	tree.Branch("RCIB", RCIB, "RCIB/F")
	tree.Branch("RCIC", RCIC, "RCIC/F")
	tree.Branch("RCID", RCID, "RCID/F")
	tree.Branch("OTR", OTR, "OTR/F")

	tree.Branch("TBEF", TBEF, "TBEF/F")
	tree.Branch("TAFT", TAFT, "TAFT/F")
	tree.Branch("SDEL", SDEL, "SDEL/F")
	tree.Branch("KTH", KTH, "KTH/F")
	tree.Branch("VOLT", VOLT, "VOLT/F")
	tree.Branch("VVET", VVET, "VVET/F")
	
	tree.Branch("RUN", RUN, "RUN/I")
	tree.Branch("SN", SN, "SN/I")
	tree.Branch("DateSec", DateSec, "DateSec/l")	

	tree.Branch("pearsA", pearsA, "pearsA/F")
	tree.Branch("dtw_valA", dtw_valA, "dtw_valA/F")

	nEntries = t_ERA_ANA.GetEntries()

	for i in range(nEntries):
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		t_ERA_ANA.GetEntry(i)
		RUN[0], SN[0] = int(t_ERA_ANA.RUN), int(t_ERA_ANA.SN)

		try:
			pearsA[0] = d_PSA_events[str(RUN[0]) + "_" + str(SN[0])][0]
			dtw_valA[0] = d_PSA_events[str(RUN[0]) + "_" + str(SN[0])][1]
			
			EC1_ERA[0]   =  t_ERA_ANA.EC1_ERA
			EC2_ERA[0]   =  t_ERA_ANA.EC2_ERA
			CHIC1_ERA[0] =  t_ERA_ANA.CHIC1_ERA
			CHIC2_ERA[0] =  t_ERA_ANA.CHIC2_ERA
			FWC1_ERA[0]  =  t_ERA_ANA.FWC1_ERA
			FWC2_ERA[0]  =  t_ERA_ANA.FWC2_ERA
			
			EC1[0]       =  t_ERA_ANA.EC1
			EC2[0]       =  t_ERA_ANA.EC2
			EIA[0]       =  (t_ERA_ANA.EIA)/corr_IA
			EIB[0]       =  (t_ERA_ANA.EIB)/corr_IB
			EIC[0]       =  (t_ERA_ANA.EIC)/corr_IC
			EID[0]       =  (t_ERA_ANA.EID)/corr_ID
			CHIC1[0]     =  t_ERA_ANA.CHIC1
			CHIC2[0]     =  t_ERA_ANA.CHIC2
			CHIA[0]      =  t_ERA_ANA.CHIA
			CHIB[0]      =  t_ERA_ANA.CHIB
			CHIC[0]      =  t_ERA_ANA.CHIC
			CHID[0]      =  t_ERA_ANA.CHID
			FWC1[0]      =  t_ERA_ANA.FWC1
			FWC2[0]      =  t_ERA_ANA.FWC2
			FWIA[0]      =  (t_ERA_ANA.FWIA)/corr_IA
			FWIB[0]      =  (t_ERA_ANA.FWIB)/corr_IB
			FWIC[0]      =  (t_ERA_ANA.FWIC)/corr_IC
			FWID[0]      =  (t_ERA_ANA.FWID)/corr_ID
			
			RCIA[0]      =  t_ERA_ANA.RCIA
			RCIB[0]      =  t_ERA_ANA.RCIB
			RCIC[0]      =  t_ERA_ANA.RCIC
			RCID[0]      =  t_ERA_ANA.RCID
			OTR[0]       = t_ERA_ANA.OTR
			
			TBEF[0]      =  t_ERA_ANA.TBEF
			TAFT[0]      =  t_ERA_ANA.TAFT
			SDEL[0]      =  t_ERA_ANA.SDEL
			KTH[0]       =  t_ERA_ANA.KTH
			VOLT[0]      =  t_ERA_ANA.VOLT
			VVET[0]      =  t_ERA_ANA.VVET
			
			DateSec[0]   =  t_ERA_ANA.DateSec			

			tree.Fill()

		except KeyError:
			pass

	print
	tree.Write()
	filou.Close()