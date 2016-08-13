from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl

def superimpose_Q_values(file_name, bolo_name):

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

	Q_src  = arr_src["Q"]	
	Q_data = arr_data["Q"]

	h_src_highE  = TH1F("h_src_highE", "h_src_highE", 100, 0, 1)	
	h_data_highE = TH1F("h_data_highE", "h_data_highE", 100, 0, 1)

	for index, elem in enumerate(Q_src):
		if (10<heat_src[index]<30):
			h_src_highE.Fill(elem)

	for index, elem in enumerate(Q_data):
		if (10<heat_data[index]<30):
			h_data_highE.Fill(elem)

	h_src_lowE  = TH1F("h_src_lowE", "h_src_lowE", 100, 0, 1)	
	h_data_lowE = TH1F("h_data_lowE", "h_data_lowE", 100, 0,1)

	for index, elem in enumerate(Q_src):
		if (heat_src[index]<8):
			h_src_lowE.Fill(elem)

	for index, elem in enumerate(Q_data):
		if (heat_data[index]<8):
			h_data_lowE.Fill(elem)

	PyRPl.process_TH1(h_src_highE, X_title = "Heat (keV)", Y_title = "Arbitrary counts", color = kRed)	
	PyRPl.process_TH1(h_data_highE, X_title = "Heat (keV)", Y_title = "Arbitrary counts", color = kBlack)

	PyRPl.process_TH1(h_src_lowE, X_title = "Heat (keV)", Y_title = "Arbitrary counts", color = kRed)	
	PyRPl.process_TH1(h_data_lowE, X_title = "Heat (keV)", Y_title = "Arbitrary counts", color = kBlack)

	cc = TCanvas("cc", "cc")

	h_data_highE.SetMaximum(1.5*h_data_highE.GetMaximum())
	h_src_highE.Scale(float(h_data_highE.Integral())/float(h_src_highE.Integral()))

	h_data_highE.Draw()
	h_src_highE.Draw("same")
	cc.Print("./Figures/" +bolo_name + "/" + bolo_name + "_" +  file_name + "_Q_distri_highE.png")

	cc = TCanvas("cc", "cc")

	h_data_lowE.SetMaximum(1.5*h_data_lowE.GetMaximum())
	h_src_lowE.Scale(float(h_data_lowE.Integral())/float(h_src_lowE.Integral()))

	h_data_lowE.Draw()
	h_src_lowE.Draw("same")
	cc.Print("./Figures/" +bolo_name + "/" + bolo_name + "_" + file_name + "_Q_distri_lowE.png")

file_name_list = ["S1Beta", "S2Beta", "S1Pb", "S2Pb"]
bolo_name_list = ["FID827", "FID828", "FID837", "FID838"]
for bolo_name in bolo_name_list:
	for file_name in file_name_list:
		superimpose_Q_values(file_name, bolo_name)