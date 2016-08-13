from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl

def superimpose_data_to_source(file_name, bolo_name):

	"""Superimpose source and data 
	beta/Pb histograms (Qplots)
	
	Detail:
		Go figure

    Args:     
		file_name (str)        = name of event file
		bolo_name (str)        = bolometer_name
		
	Returns:
		void

	Raises:
		void
	"""
       
	data_types = {"names": 
					("EC", "EI", "ER", "Q"), 
				"formats":   
					("f", "f", "f", "f")}


	arr_src  = np.loadtxt("./Populations/FID808_Run305/" + file_name + "_heatremoved_arranged.txt", delimiter=",",  dtype=data_types)
	arr_data =  np.loadtxt("./Populations/" +bolo_name + "/" + bolo_name + "_" + file_name + "_arranged.txt", delimiter=",",  dtype=data_types)

	heat_src  = arr_src["EC"]	
	heat_data = arr_data["EC"]
	Q_src  = arr_src["Q"]	
	Q_data = arr_data["Q"]


	h_src  = TH2F("h_src", "h_src", 200, 0, 50, 100, 0, 1.2)	
	h_data = TH2F("h_data", "h_data", 200, 0, 50, 100, 0, 1.2)

	for index, elem in enumerate(heat_src):
		h_src.Fill(elem, Q_src[index])

	for index, elem in enumerate(heat_data):
		h_data.Fill(elem,Q_data[index])

	PyRPl.process_TH2(h_src, X_title = "Heat (keV)", Y_title = "Q", color = kRed)	
	PyRPl.process_TH2(h_data, X_title = "Heat (keV)", Y_title = "Q", color = kBlack)

	cc = TCanvas("cc", "cc")

	h_src.SetMarkerSize(4)
	h_src.Draw()
	h_data.Draw("same")
	
	cc.Print("./Figures/" +bolo_name + "/" + bolo_name + "_" + file_name + "_Qplot.png")

file_name_list = ["S1Beta", "S2Beta", "S1Pb", "S2Pb"]
bolo_name_list = ["FID837"]
for bolo_name in bolo_name_list:
	for file_name in file_name_list:
		superimpose_data_to_source(file_name, bolo_name)