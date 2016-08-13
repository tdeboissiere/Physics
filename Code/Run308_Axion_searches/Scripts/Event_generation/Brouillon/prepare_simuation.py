import numpy as np
from ROOT import*
import PyROOTPlots as PyRPl
import matplotlib.pylab as plt
import script_utils as script_utils

def get_correlation_matrix(bolo_name, tree_name = "data"):
	""" Get ionisation correlation matrix
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
	"""

	filename = "../../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root"
	tree,f = PyRPl.open_ROOT_object(filename, "t_axion")

	arr_EIA = []
	arr_EIB = []
	arr_EIC = []
	arr_EID = []

	# Loop over entries 
	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		if tree.EC < 0.5 :
			arr_EIA.append(tree.EIA)
			arr_EIB.append(tree.EIB)
			arr_EIC.append(tree.EIC)
			arr_EID.append(tree.EID)

	arr_EIA = np.array(arr_EIA).astype(float)
	arr_EIB = np.array(arr_EIB).astype(float)
	arr_EIC = np.array(arr_EIC).astype(float)
	arr_EID = np.array(arr_EID).astype(float)

	arr_stacked = np.vstack((arr_EIA, arr_EIB))
	arr_stacked = np.vstack((arr_stacked, arr_EIC))
	arr_stacked = np.vstack((arr_stacked, arr_EID))

	L = np.linalg.cholesky(np.cov(arr_stacked))
	np.savetxt("./Text_files/" + bolo_name + "_Lmatrix_coeff.txt", np.ndarray.flatten(L))

def get_histogram_for_simu(bolo_name):
	"""
	Get FWHM histograms

	Detail:

	Arguments:
	bolo_name (str) the bolometer name

	Outputs:

	void
	"""
	
	#Load the data
	filename = "../../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root"
	tree,file_tree   = PyRPl.open_ROOT_object(filename, "t_axion")
	
	l_events    = TEventList("l_events")
	tree.Draw(">>l_events" )
	pop_len     = l_events.GetN()
	list_events = []
	for k in range(pop_len):
		counter = l_events.GetEntry(k)
		tree.GetEntry(counter)
		list_events.append(tree.JOUR)

	#Get min and max time, + safety margin
	tmin, tmax = min(list_events)-5, max(list_events)+5
		
	# 2 D threshold histogram
	binning = 2000
	hist_thresh = TH2F("hist_thresh","hist_thresh",binning,tmin, tmax,1000,0,2)
	hist_FWIA   = TH2F("hist_FWIA","hist_FWIA",binning,tmin, tmax,1000,0,2)
	hist_FWIC   = TH2F("hist_FWIC","hist_FWIC",binning,tmin, tmax,1000,0,2)
	hist_FWFID  = TH2F("hist_FWFID","hist_FWFID",binning,tmin, tmax,1000,0,2)
	hist_FWC    = TH2F("hist_FWC","hist_FWC",binning,tmin, tmax,1000,0,2)
	
	tree.Project( "hist_thresh", "KTH:JOUR")
	tree.Project( "hist_FWIA",  "FWIA:JOUR")
	tree.Project( "hist_FWIC",  "FWIC:JOUR")
	tree.Project( "hist_FWFID",  "FWF:JOUR")
	tree.Project( "hist_FWC",  "FWC:JOUR")
	
	hprojX_KTH    = hist_thresh.ProjectionX()
	
	hprojY_KTH    = hist_thresh.ProjectionY()
	hprojY_FWIA = hist_FWIA.ProjectionY()
	hprojY_FWIC = hist_FWIC.ProjectionY()
	hprojY_FWFID = hist_FWFID.ProjectionY()
	hprojY_FWC = hist_FWC.ProjectionY()
	
	#Fill time and threshold lists used to build the 1D threshold histhograms
	list_time, list_fwia, list_fwic, list_fwfid,  list_fwc, list_jour = [], [], [], [], [],[]
	for ix in range(1,binning) : 
		temp_fwia, temp_fwic,temp_fwfid, temp_fwc =0,0,0,0
		for iy in range(1, hprojY_KTH.GetNbinsX()): 
			fwia = hprojY_FWIA.GetBinCenter(iy)
			fwic = hprojY_FWIC.GetBinCenter(iy)
			fwfid = hprojY_FWFID.GetBinCenter(iy)
			fwc = hprojY_FWC.GetBinCenter(iy)

			if ( hist_FWIA.GetBinContent(ix,iy) !=0) :
				temp_fwia =fwia

			if ( hist_FWIC.GetBinContent(ix,iy) !=0) :
				temp_fwic =fwic

			if ( hist_FWFID.GetBinContent(ix,iy) !=0) :
				temp_fwfid =fwfid

			if ( hist_FWC.GetBinContent(ix,iy) !=0) :
				temp_fwc =fwc

		list_time.append(hprojX_KTH.GetBinCenter(ix))
		if (temp_fwia !=0 and temp_fwic!=0 and temp_fwfid!=0 and temp_fwc!=0):
			list_fwia.append(temp_fwia)
			list_fwic.append(temp_fwic)
			list_fwfid.append(temp_fwfid)
			list_fwc.append(temp_fwc)
			list_jour.append(1)
		else:
			list_fwia.append(0)
			list_fwic.append(0)
			list_fwfid.append(0)
			list_fwc.append(0)
			list_jour.append(0)         

	# ROOT says there is an issue with increasing value for the bins represented by time
	# I disagree.
	# Solution: use a graph to build the 1D threshold histhogram

	gr_fwia = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwia))
	hfwia   = TH1F("hfwia", "hfwia",np.size(list_time),min(list_time), max(list_time))
	for k in range(1, np.size(list_time)+1):
		hfwia.SetBinContent(k,gr_fwia.Eval(list_time[k-1]))

	gr_fwic = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwic))
	hfwic   = TH1F("hfwic", "hfwic",np.size(list_time),min(list_time), max(list_time))
	for k in range(1, np.size(list_time)+1):
		hfwic.SetBinContent(k,gr_fwic.Eval(list_time[k-1]))

	gr_fwfid = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwfid))
	hfwfid   = TH1F("hfwfid", "hfwfid",np.size(list_time),min(list_time), max(list_time))
	for k in range(1, np.size(list_time)+1):
		hfwfid.SetBinContent(k,gr_fwfid.Eval(list_time[k-1]))

	gr_fwc = TGraph(np.size(list_time), np.array(list_time), np.array(list_fwc))
	hfwc   = TH1F("hfwc", "hfwc",np.size(list_time),min(list_time), max(list_time))
	for k in range(1, np.size(list_time)+1):
		hfwc.SetBinContent(k,gr_fwc.Eval(list_time[k-1]))

	hist_jour = TH1F("hist_jour", "hist_jour", np.size(list_time),min(list_time), max(list_time))
	for k in range(1, np.size(list_time)+1):
		hist_jour.SetBinContent(k,list_jour[k-1])
	
	root_path = "./ROOT_files/" + bolo_name + "_fwhm_time_hist.root"
	fout = TFile(root_path, "recreate")
	hist_jour.Write()
	hfwia.Write()
	hfwic.Write()
	hfwfid.Write()
	hfwc.Write()
	fout.Close()

list_bolo_name = ["FID824", "FID825", "FID827", "FID837", "FID838", "FID839", "FID841", "FID842"]
for bolo_name in list_bolo_name :
	script_utils.print_utility("Processing bolo " + bolo_name)
	# get_correlation_matrix(bolo_name)
	get_histogram_for_simu(bolo_name)

