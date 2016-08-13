
import glob
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

def merge_wavelet_PSA_tree(bolo_name):

	""" Merge wavelet PSA tree
	
	Detail:
		Merge wavelet PSA tree

	Args:
		bolo_name = (str) bolo name 
		
	Returns:
		void

	Raises:
		void
	"""
	
	#Define non merged file directory
	data_dir="./ROOT_files/" + bolo_name + "/"

    #List of files in data_dir
	list_files = glob.glob(data_dir + "*event_wavelet_PSA_tree*")
	list_files = sorted(list_files)
	list_files = [elem for elem in list_files if ("all" not in elem and "simu" not in elem)]

	chain=TChain("twavelet","eionbis")
	for partition in list_files:
		print partition
		chain.AddFile(partition)

	fusion_file_name = data_dir + bolo_name + "_all_event_wavelet_PSA_tree.root"
	chain.Merge(fusion_file_name)


def add_wavelet_PSA(bolo_name, analysis_type):

	"""Add wavelet PSA to merged ERA ANA tree
	
	Detail:
		void

	Args:
		bolo_name = (str) bolo name 
		analysis_type = (str) type of analysis
		
	Returns:
		void

	Raises:
		void
	"""

	#Get ERA/ANA merged _ files
	path_ERA_ANA = "../Fond_ERA_merged/"
	# t_noPSA, f_noPSA = PyRPl.open_ROOT_object(path_ERA_ANA + bolo_name + "_" + analysis_type+ "_fond.root", "data")
	t_noPSA, f_noPSA = PyRPl.open_ROOT_object(path_ERA_ANA + bolo_name + "_fond.root", "data")

	# #Get Wavelet PSA tree
	# path_wave = "../Wavelet/ROOT_files/" + bolo_name + "/" 
	# t_wave, f_wave = PyRPl.open_ROOT_object(path_wave +  bolo_name + "_all_event_wavelet_PSA_tree.root", "t_wavelet")

	#Load event dictionary
	d_PSA_events = OF_ut.load_PSA_event_dict(bolo_name, "all_event")	

	#Create new output tree
	path_output = script_utils.create_directory("./ROOT_files/")
	# output_file_name = path_output + bolo_name + "_" + analysis_type + "_PSA_and_fond.root"
	output_file_name = path_output + bolo_name +  "_PSA_and_fond.root"
	filou=TFile(output_file_name, "recreate")
	tree = TTree("data", "data")

	EC1      =  array("f", [ 0 ] )
	EC2      =  array("f", [ 0 ] )
	EIA      =  array("f", [ 0 ] )
	EIB      =  array("f", [ 0 ] )
	EIC      =  array("f", [ 0 ] )
	EID      =  array("f", [ 0 ] )

	XOC1     =  array("f", [ 0 ] )
	XOC2     =  array("f", [ 0 ] )
	CHIA     =  array("f", [ 0 ] )
	CHIB     =  array("f", [ 0 ] )
	CHIC     =  array("f", [ 0 ] )
	CHID     =  array("f", [ 0 ] )
	OWC1     =  array("f", [ 0 ] )
	OWC2     =  array("f", [ 0 ] )
	FWIA     =  array("f", [ 0 ] )
	FWIB     =  array("f", [ 0 ] )
	FWIC     =  array("f", [ 0 ] )
	FWID     =  array("f", [ 0 ] )
	
	RCIA     =  array("f", [ 0 ] )
	RCIB     =  array("f", [ 0 ] )
	RCIC     =  array("f", [ 0 ] )
	RCID     =  array("f", [ 0 ] )
	
	TAFT     =  array("f", [ 0 ] )
	TBEF     =  array("f", [ 0 ] )
	SDEL     =  array("f", [ 0 ] )
	KTH      =  array("f", [ 0 ] )
	VOLT     =  array("f", [ 0 ] )
	VVET     =  array("f", [ 0 ] )
	
	RUN      =  array("i", [ 0 ] )
	SN       =  array("i", [ 0 ] )
	APAT     =  array("f", [ 0 ] )
	MULT     =  array("f", [ 0 ] )
	UT1      =  array("f", [ 0 ] )
	UT2      =  array("f", [ 0 ] )
	
	pearsA   = array("f", [ 0 ])
	dtw_valA = array("f", [ 0 ])

	#ANA variables
	tree.Branch("EC1", EC1, "EC1/F")
	tree.Branch("EC2", EC2, "EC2/F")
	tree.Branch("EIA", EIA, "EIA/F")
	tree.Branch("EIB", EIB, "EIB/F")
	tree.Branch("EIC", EIC, "EIC/F")
	tree.Branch("EID", EID, "EID/F")

	tree.Branch("EIA", EIA, "EIA/F")
	tree.Branch("EIB", EIB, "EIB/F")
	tree.Branch("EIC", EIC, "EIC/F")
	tree.Branch("EID", EID, "EID/F")
	tree.Branch("XOC1", XOC1, "XOC1/F")
	tree.Branch("XOC2", XOC2, "XOC2/F")
	tree.Branch("CHIA", CHIA, "CHIA/F")
	tree.Branch("CHIB", CHIB, "CHIB/F")
	tree.Branch("CHIC", CHIC, "CHIC/F")
	tree.Branch("CHID", CHID, "CHID/F")
	tree.Branch("OWC1", OWC1, "OWC1/F")
	tree.Branch("OWC2", OWC2, "OWC2/F")
	tree.Branch("FWIA", FWIA, "FWIA/F")
	tree.Branch("FWIB", FWIB, "FWIB/F")
	tree.Branch("FWIC", FWIC, "FWIC/F")
	tree.Branch("FWID", FWID, "FWID/F")

	tree.Branch("RCIA", RCIA, "RCIA/F")
	tree.Branch("RCIB", RCIB, "RCIB/F")
	tree.Branch("RCIC", RCIC, "RCIC/F")
	tree.Branch("RCID", RCID, "RCID/F")

	tree.Branch("TBEF", TBEF, "TBEF/F")
	tree.Branch("TAFT", TAFT, "TAFT/F")
	tree.Branch("SDEL", SDEL, "SDEL/F")
	tree.Branch("KTH", KTH, "KTH/F")
	tree.Branch("VOLT", VOLT, "VOLT/F")
	tree.Branch("VVET", VVET, "VVET/F")
	
	tree.Branch("RUN", RUN, "RUN/I")
	tree.Branch("SN", SN, "SN/I")
	tree.Branch("UT1", UT1, "UT1/F")	
	tree.Branch("UT2", UT2, "UT2/F")	
	tree.Branch("APAT", APAT, "APAT/F")	
	tree.Branch("MULT", MULT, "MULT/F")	

	tree.Branch("pearsA", pearsA, "pearsA/F")
	tree.Branch("dtw_valA", dtw_valA, "dtw_valA/F")

	nEntries = t_noPSA.GetEntries()
	nmatch = 0

	for i in range(nEntries):
		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		t_noPSA.GetEntry(i)
		RUN[0], SN[0] = int(t_noPSA.RUN), int(t_noPSA.SN)

		try:
			pearsA[0] = d_PSA_events[str(RUN[0]) + "_" + str(SN[0])][0]
			dtw_valA[0] = d_PSA_events[str(RUN[0]) + "_" + str(SN[0])][1]

			nmatch+=1

			EC1[0]       =  t_noPSA.EC1
			EC2[0]       =  t_noPSA.EC2
			EIA[0]       =  t_noPSA.EIA
			EIB[0]       =  t_noPSA.EIB
			EIC[0]       =  t_noPSA.EIC
			EID[0]       =  t_noPSA.EID
						
			XOC1[0]     =  t_noPSA.XOC1
			XOC2[0]     =  t_noPSA.XOC2
			CHIA[0]      =  t_noPSA.CHIA
			CHIB[0]      =  t_noPSA.CHIB
			CHIC[0]      =  t_noPSA.CHIC
			CHID[0]      =  t_noPSA.CHID
			OWC1[0]      =  t_noPSA.OWC1
			OWC2[0]      =  t_noPSA.OWC2
			FWIA[0]      =  t_noPSA.FWIA
			FWIB[0]      =  t_noPSA.FWIB
			FWIC[0]      =  t_noPSA.FWIC
			FWID[0]      =  t_noPSA.FWID
			
			RCIA[0]      =  t_noPSA.RCIA
			RCIB[0]      =  t_noPSA.RCIB
			RCIC[0]      =  t_noPSA.RCIC
			RCID[0]      =  t_noPSA.RCID
			
			TBEF[0]      =  t_noPSA.TBEF
			TAFT[0]      =  t_noPSA.TAFT
			SDEL[0]      =  t_noPSA.SDEL
			KTH[0]       =  t_noPSA.KTH
			VOLT[0]      =  t_noPSA.VOLT
			VVET[0]      =  t_noPSA.VVET

			UT1[0]      =  t_noPSA.UT1
			UT2[0]      =  t_noPSA.UT2	
			APAT[0]		=t_noPSA.APAT
			MULT[0]		=t_noPSA.MULT

			tree.Fill()

		except KeyError:
			pass

	print nmatch 
	print nEntries

	print
	tree.Write()
	filou.Close()

bolo_name = "FID837"
analysis_type = "ana_0.5_0_5"
add_wavelet_PSA(bolo_name, analysis_type)
# merge_wavelet_PSA_tree(bolo_name)