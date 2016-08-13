from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl

def superimpose_data_to_source(file_name, bolo_name):

	"""Superimpose source and data 
	beta/Pb histograms
	
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

	h_src  = TH1F("h_src", "h_src", 100, 0, 50)	
	h_data = TH1F("h_data", "h_data", 100, 0, 50)

	for elem in heat_src:
		h_src.Fill(elem)

	for elem in heat_data:
		h_data.Fill(elem)

	PyRPl.process_TH1(h_src, X_title = "Heat (keV)", Y_title = "Arbitrary counts", color = kRed)	
	PyRPl.process_TH1(h_data, X_title = "Heat (keV)", Y_title = "Arbitrary counts", color = kBlack)

	cc = TCanvas("cc", "cc")

	h_data.SetMaximum(1.5*h_data.GetMaximum())
	h_src.Scale(float(h_data.Integral())/float(h_src.Integral()))

	h_data.Draw()
	h_src.Draw("same")
	cc.Print("./Figures/" +bolo_name + "/" + bolo_name + "_" + file_name + ".png")

file_name_list = ["S1Beta", "S2Beta", "S1Pb", "S2Pb"]
bolo_name_list = ["FID837"]
for bolo_name in bolo_name_list:
	for file_name in file_name_list:
		superimpose_data_to_source(file_name, bolo_name)