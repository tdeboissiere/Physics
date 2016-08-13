from ROOT import *
import PyROOTPlots as PyRPl
import sys
import script_utils as script_utils
import Analysis_utilities as Ana_ut


def optimize_wavelet_cut(bolo_name):

	"""optimize wavelet cut
	
	Detail:
		try to eliminate as much heat only as possible
		while keeping a good enough efficiency
              
	Args:
		bolo_name = (type) bolo name
		
	Returns:
		void
          
	Raises:
		void
	"""     
	
	

	simu_path  = "../Fond_ERA_merged/Fond_simu/"
	tsimu, fsimu =  PyRPl.open_ROOT_object(simu_path + bolo_name + "_PSA_merged_simu_tree.root", "t_merged")

	ERA_path     = "../Fond_ERA_merged/"
	tERA, fERA   = PyRPl.open_ROOT_object(ERA_path + bolo_name + "_PSA_and_fond.root", "t_merged")
		
	#1D plot 
	h2D = TH2F("h2D", "h2D",1000, 0, 5, 1000, 0, 1.2)
	h2Dbis = TH2F("h2Dbis", "h2Dbis",1000, 0, 5, 1000, 0, 500)
	tsimu.Project("h2D", "pearsA:EC1_ERA")
	tsimu.Project("h2Dbis", "dtw_valA:EC1_ERA")
	h2D.Draw()
	raw_input()
	h2Dbis.Draw()
	raw_input()

	#Efficiency histogram
	hsimunoPSA= TH1F("hsimunoPSA", "hsimunoPSA", 100,0,5)
	hsimuPSA= TH1F("hsimuPSA", "hsimuPSA", 100,0,5)

	tsimu.Project("hsimunoPSA","0.5*(EC1_ERA+EC2_ERA)")
	tsimu.Project("hsimuPSA","0.5*(EC1_ERA+EC2_ERA)", "pearsA>0.91 && dtw_valA<100")

	# hsimunoPSA.Draw()
	# hsimuPSA.Draw("same")
	# raw_input()

	h= TH1F("h", "h", 100,0,2)
	PyRPl.process_TH1(h, X_title = "Heat (keVee)")
	h.Draw()

	teff1 = TEfficiency(hsimuPSA, hsimunoPSA)
	teff1.SetName("teff1")
	teff1.SetLineColor(kRed)
	teff1.Draw("same")

	c1.Print("./Figures/" + bolo_name + "/cut_eff.png")

	leg = TLegend(0.1934673,0.5813953,0.5037688,0.7043189,"","brNDC")
	leg.AddEntry(teff1.GetName(), "PSA cut", "leg")
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)
	leg.Draw("same")
	
	raw_input()

	hsimuPSA.Divide(hsimunoPSA)
	hsimuPSA.Smooth(50)
	hsimuPSA.Draw()
	s = TSpline3(hsimuPSA)
	s.SetLineColor(kRed)
	s.Draw("same")
	raw_input()


	class Spline_fit:
		def __call__( self, x, par ):
			if x[0]>4:
				return s.Eval(4)
			elif x[0]<0.05:
				return s.Eval(0.05)
			else:
				return s.Eval(x[0]) + par[0]

	ffit_extra = TF1("PSA_eff", Spline_fit(), 0, 32,1)
	ffit_extra.SetParameter(0,0)

	f = TFile("./ROOT_files/" + bolo_name + "/PSA_cut_eff.root", "recreate")
	ffit_extra.SetNpx(1000)
	ffit_extra.Write()
	f.Close()

def wavelet_cut_eff_NR(bolo_name, fiducial_voltage):

	"""Convert the efficiency to keVNR
	
	Detail:
		void
              
	Args:
		bolo_name        = (type) bolo name
		fiducial_voltage = (float) Vfid
		
	Returns:
		void
          
	Raises:
		void
	"""     
          

	#Get efficiency function
	feff, file_eff = PyRPl.open_ROOT_object("./ROOT_files/" + bolo_name + "/PSA_cut_eff.root", "PSA_eff")

	class NR_to_EE:
	   def __call__( self, x, par ):
			ENR=x[0]
			Vf=par[0]
			Q=0.16*pow(ENR,0.18)
			return ENR*(1+Q*float(Vf)/3.)/(1+float(Vf)/3.)

	fNR_to_EE = TF1("r", NR_to_EE(), 0,40,1)
	fNR_to_EE.SetParameter(0,fiducial_voltage)
	fNR_to_EE.SetNpx(200)

	class PSA_eff_NR:
		def __call__( self, x, par ):

			Eee = fNR_to_EE(x[0])
			eff = feff.Eval(Eee) + par[0]
			if x[0]>1:
				return eff
			else:
				return 0

	feff_NR = TF1("PSA_eff_NR" , PSA_eff_NR(), 0,31,1 )
	feff_NR.SetParameter(0,0)
	feff_NR.SetNpx(500)
	fout = TFile("./ROOT_files/" + bolo_name + "/PSA_cut_eff_NR.root", "recreate")
	feff_NR.Write()
	fout.Close()


# optimize_wavelet_cut("FID837")
wavelet_cut_eff_NR("FID837", 8.)