from ROOT import *
import PyROOTPlots as PyRPl
import sys
import script_utils as script_utils
import Analysis_utilities as Ana_ut


def wavelet_plots(bolo_name):

	"""Create various plots to illustrate wavelet
	
	Detail:
		void
              
	Args:
		bolo_name        = (type) bolo name
		
	Returns:
		void
          
	Raises:
		void
	"""     
          

	#heat cut
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")
	cheat = standard_cuts + "&&0.5*(EIB+EID)<0"

	#Get data
	ERA_path     = "../Fond_ERA_merged/"
	tree, fERA   = PyRPl.open_ROOT_object(ERA_path + bolo_name + "_PSA_and_fond.root", "data")

	# #Get efficiency function
	# feff, file_eff = PyRPl.open_ROOT_object("./ROOT_files/" + bolo_name + "/PSA_cut_eff.root", "PSA_eff")

	#Draw heatonly spectrum before and after PSA
	h1D     = TH1F("h1D", "h1D", 1000, 0.5, 5)
	h1D_PSA = TH1F("h1D_PSA", "h1D_PSA", 1000, 0.5, 5)	

	tree.Project("h1D","0.5*(EC1+EC2)", cheat)
	tree.Project("h1D_PSA","0.5*(EC1+EC2)", cheat + "&&pearsA>0.91 && dtw_valA<100")

	PyRPl.process_TH1(h1D, X_title="Heat (keV)", Y_title = "Counts (all periods)", color=kBlue)
	PyRPl.process_TH1(h1D_PSA, X_title="Heat (keV)", Y_title = "Counts (all periods)", color=kRed)

	h1D.SetMaximum(1.2*max(h1D.GetMaximum(), h1D_PSA.GetMaximum()))

	cc = TCanvas("cc", "cc")
	# # gPad.SetLogy()
	# hANA.Draw()
	# hANA_PSA.Draw("same")
	# raw_input()
	# cc.Print("./Figures/" + bolo_name + "/heat_spectrum_ANA_PSA.png")

	h1D.Draw()
	h1D_PSA.Draw("same")
	raw_input()
	cc.Print("./Figures/" + bolo_name + "/heat_spectrum_ERA_PSA.png")



bolo_name = "FID837"
wavelet_plots(bolo_name)
# wavelet_cut_eff_NR(bolo_name, 8.)