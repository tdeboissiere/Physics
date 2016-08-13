from ROOT import *
import BDT_file_handler as BDT_fh
import script_utils as script_utils
from array import array
import Analysis_utilities as Ana_ut


def get_ion_corr_tree(bolo_name, d_cut, data_dir, tree_name = "data"):

	""" Skim true event tree so that it passes the standard cut + heat < 0.5 keV cut
	
	Detail:
		The resulting tree will be used to compute the ionisation correlation matrix

	Args:
		bolo_name     = (str) bolometer name
		d_cut         = (dict) indicates the cuts (inf/sup) on heat and ion 
		data_dir      = (str) data directory

	Returns:
		void

	Raises:
		AssertionError 
	"""
	
	file_tree     = TFile(data_dir+bolo_name+"_fond.root")
	tree          = file_tree.Get(tree_name)
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	#Load estimators
	d_est = BDT_fh.open_estimator_file(bolo_name, "")
	
	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name, "")
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	heat_cut = str(d_cut["ECsup"]) + ">" + d_est["HEAT"] 
	veto_cut = "EIA<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIA"] + "&&" + "EIC<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIC"]
	
	all_cuts = "&&".join([standard_cuts, heat_cut, veto_cut])
	print all_cuts.split("&&")

	#Event list
	l_event    = TEventList("l_event")
	tree.Draw(">>l_event",all_cuts )
	pop_len = l_event.GetN()


 	file_skim = TFile(data_dir+bolo_name + "_for_ion_corr_fond.root", "recreate")
	skim_tree          = TTree("data", "data")

	EIA =  array("f", [ 0 ] )
	EIB =  array("f", [ 0 ] )
	EIC =  array("f", [ 0 ] )
	EID =  array("f", [ 0 ] )

	skim_tree.Branch("EIA", EIA, "EIA/F")
	skim_tree.Branch("EIB", EIB, "EIB/F")
	skim_tree.Branch("EIC", EIC, "EIC/F")
	skim_tree.Branch("EID", EID, "EID/F")

	for k in range(pop_len):
		counter = l_event.GetEntry(k)
		tree.GetEntry(counter)

		EIA[0] = tree.EIA
		EIB[0] = tree.EIB
		EIC[0] = tree.EIC
		EID[0] = tree.EID

		skim_tree.Fill()
	
	skim_tree.Write()
	file_skim.Close()

def skim_true_event_tree_with_heatrate(bolo_name, d_cut, analysis_type, data_dir, tree_name = "data"):

	""" Skim true event tree so that it passes the standard cut + analysis cuts
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
		d_cut         = (dict) indicates the cuts (inf/sup) on heat and ion 
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		data_dir      = (str) data directory

	Returns:
		void

	Raises:
		AssertionError 
	"""
	
	file_tree     = TFile(data_dir+bolo_name+"_fond.root")
	tree          = file_tree.Get(tree_name)
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	#Load estimators
	d_est = BDT_fh.open_estimator_file(bolo_name, "")
	
	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name, "")
	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
		d_std[key] = str(2.7*d_std[key])
	
	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	heat_cut = str(d_cut["ECinf"]) + "<EC" + "&&" + str(d_cut["ECsup"]) + ">EC" 
	# ion_cut  = str(d_cut["EIinf"]) + "<" + d_est["FID"] + "&&" + str(d_cut["EIsup"]) + ">" + d_est["FID"]
	ion_cut  = str(d_cut["EIinf"]) + "<0.5*(EIB+EID) &&" + str(d_cut["EIsup"]) + ">0.5*(EIB+EID)"
	veto_cut = "EIA<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIA"] + "&&" + "EIC<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIC"]
	
	all_cuts = "&&".join([standard_cuts, heat_cut, ion_cut, veto_cut])
	print all_cuts.split("&&")

	#Event list
	l_event       = TEventList("l_event")
	tree.Draw(">>l_event",all_cuts )
	
	#Get the histogram of the heatonly rate over time
	file_heatrate = TFile("../Analyse_" + bolo_name + "/ROOT_files/" + bolo_name +  "_heatonly_rate_over_time_" + analysis_type + ".root", "read") 
	heatrate      = file_heatrate.Get("hheat")

	file_skim = TFile(data_dir+bolo_name + "_" + analysis_type  +"_fond.root", "recreate")
	skim_tree = TTree("data", "data")

	EC1  =  array("f", [ 0 ] )
	EC2  =  array("f", [ 0 ] )
	EIA  =  array("f", [ 0 ] )
	EIB  =  array("f", [ 0 ] )
	EIC  =  array("f", [ 0 ] )
	EID  =  array("f", [ 0 ] )
	EC   = array("f", [ 0 ] )
	EFID = array("f", [ 0 ] )
	HR   =  array("f", [ 0 ] )
	RUN  = array("i", [0])
	SN   = array("i", [0])

	skim_tree.Branch("EC1", EC1, "EC1/F")
	skim_tree.Branch("EC2", EC2, "EC2/F")
	skim_tree.Branch("EIA", EIA, "EIA/F")
	skim_tree.Branch("EIB", EIB, "EIB/F")
	skim_tree.Branch("EIC", EIC, "EIC/F")
	skim_tree.Branch("EID", EID, "EID/F")
	skim_tree.Branch("EC", EC, "EC/F")
	skim_tree.Branch("EFID", EFID, "EFID/F")
	skim_tree.Branch("HR", HR, "HR/F")
	skim_tree.Branch("RUN", RUN, "RUN/I")
	skim_tree.Branch("SN", SN, "SN/I")

	for k in range(l_event.GetN()):
		counter = l_event.GetEntry(k)
		tree.GetEntry(counter)

		time = tree.JOUR

		EC1[0]  = tree.EC1
		EC2[0]  = tree.EC2
		EIA[0]  = tree.EIA
		EIB[0]  = tree.EIB
		EIC[0]  = tree.EIC
		EID[0]  = tree.EID
		EC[0]   = tree.EC
		EFID[0] = tree.EFID
		HR[0]   = heatrate.GetBinContent(heatrate.FindBin(time))
		RUN[0]  = int(tree.RUN)
		SN[0]   = int(tree.SN)

		skim_tree.Fill()
	
	skim_tree.Write()
	file_skim.Close()



bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
# analysis_type = "ana_1.5_min2_50"
analysis_type = "ana_0.5_0_5"
#1 sigma on heat only band=  EFid>0.2
d_cut       = {"ECinf": 0.5, "ECsup": 15, "EIinf": 0, "EIsup": 15, "sigma_vet": 5}
# d_cut       = {"ECinf": 1.5, "ECsup": 15, "EIinf": -2, "EIsup": 15, "sigma_vet": 50}
data_dir = "../Fond_ERA_merged/"
skim_true_event_tree_with_heatrate(bolo_name, d_cut, analysis_type, data_dir)

# #######################################
# ## Get correlation tree for ionisation
# ######################################
# bolo_name = "FID837"
# data_dir = "../Fond_ERA_merged/"
# d_cut       = {"ECinf": -5, "ECsup": 0.5, "EIinf": 0*0.2, "EIsup": 15, "sigma_vet": 5}
# get_ion_corr_tree(bolo_name, d_cut, data_dir, tree_name = "data")







# def skim_true_event_tree(bolo_name, d_cut, analysis_type, data_dir, tree_name = "data"):

# 	""" Skim true event tree so that it passes the standard cut + analysis cuts
	
# 	Detail:
# 		Detailed description

# 	Args:
# 		bolo_name     = (str) bolometer name
# 		d_cut         = (dict) indicates the cuts (inf/sup) on heat and ion 
# 		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
# 		data_dir      = (str) data directory

# 	Returns:
# 		void

# 	Raises:
# 		AssertionError 
# 	"""
	
# 	file_tree     = TFile(data_dir+bolo_name+"_fond.root")
# 	tree          = file_tree.Get(tree_name)
	
# 	#Load standard cuts
# 	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

# 	#Load estimators
# 	d_est = BDT_fh.open_estimator_file(bolo_name, "")
	
# 	#Load FWHM
# 	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name, "")
# 	for key in ["OWC1", "OWC2", "FWIA", "FWIB", "FWIC", "FWID"]:
# 		d_std[key] = str(2.7*d_std[key])
	
# 	#Load standard cuts
# 	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

# 	heat_cut = str(d_cut["ECinf"]) + "<" + d_est["HEAT"] + "&&" + str(d_cut["ECsup"]) + ">" + d_est["HEAT"] 
# 	ion_cut  = str(d_cut["EIinf"]) + "<" + d_est["FID"] + "&&" + str(d_cut["EIsup"]) + ">" + d_est["FID"]
# 	veto_cut = "EIA<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIA"] + "&&" + "EIC<(1./2.7)*" + str(d_cut["sigma_vet"]) + "*" + d_std["FWIC"]
	
# 	all_cuts = "&&".join([standard_cuts, heat_cut, ion_cut, veto_cut])
# 	print all_cuts.split("&&")

# 	#Event list
# 	l_event    = TEventList("l_event")
# 	tree.Draw(">>l_event",all_cuts )
# 	pop_len       = l_event.GetN()

#  	file_skim = TFile(data_dir+bolo_name + "_" + analysis_type  +"_fond.root", "recreate")
# 	skim_tree          = TTree("data", "data")

# 	EC1 =  array("f", [ 0 ] )
# 	EC2 =  array("f", [ 0 ] )
# 	EIA =  array("f", [ 0 ] )
# 	EIB =  array("f", [ 0 ] )
# 	EIC =  array("f", [ 0 ] )
# 	EID =  array("f", [ 0 ] )
# 	ENR =  array("f", [ 0 ] )
# 	RUN = array("i", [0])
# 	SN = array("i", [0])

# 	skim_tree.Branch("EC1", EC1, "EC1/F")
# 	skim_tree.Branch("EC2", EC2, "EC2/F")
# 	skim_tree.Branch("EIA", EIA, "EIA/F")
# 	skim_tree.Branch("EIB", EIB, "EIB/F")
# 	skim_tree.Branch("EIC", EIC, "EIC/F")
# 	skim_tree.Branch("EID", EID, "EID/F")
# 	skim_tree.Branch("ENR", ENR, "ENR/F")
# 	skim_tree.Branch("RUN", RUN, "RUN/I")
# 	skim_tree.Branch("SN", SN, "SN/I")

# 	for k in range(pop_len):
# 		counter = l_event.GetEntry(k)
# 		tree.GetEntry(counter)

# 		EC1[0] = tree.EC1
# 		EC2[0] = tree.EC2
# 		EIA[0] = tree.EIA
# 		EIB[0] = tree.EIB
# 		EIC[0] = tree.EIC
# 		EID[0] = tree.EID
# 		ENR[0] = 0
# 		RUN[0] = int(tree.RUN)
# 		SN[0] = int(tree.SN)

# 		skim_tree.Fill()
	
# 	skim_tree.Write()
# 	file_skim.Close()