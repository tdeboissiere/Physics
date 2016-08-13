import utils_ana_jules as ut
import script_utils as script_utils
from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import os

def get_axio_electric():
	"""Get axio electric cross section
	
	Detail:
		void

	Args:
		void

	Returns:
		void

	Raises:
		void
	"""
	
	cross_PE = np.loadtxt("../Text_files/Ge_cross.txt",delimiter=",").astype(float)
	#Define a list of axion masses + constants
	list_mass = [0, 0.01, 0.05, 0.1, 0.15, 0.2]
	alpha = 1/137.0		
	pi    = 3.14159265359	
	me    = 511.
	#Out file
	fout = TFile("../ROOT_files/Axion/Axio_electric/axio_electric_cross_section_cm2perkg.root", "recreate")
	for mass in list_mass :
		beta, betafac = np.zeros(cross_PE.shape[0]), np.zeros(cross_PE.shape[0])
		for i in range(cross_PE.shape[0]):
			if cross_PE[i,0]<= mass :
				beta[i] = 0
				betafac[i] = 0
			else :
			   beta[i]=pow(1-mass**2/(cross_PE[i,0]*cross_PE[i,0]),0.5)
			   betafac[i]=(1-pow(beta[i],2./3.)/3.0)/beta[i]				

		# cross_PE is in b/atom => divide by 120.3 to get cm2/g 
		# multiply by 1000 to get cm2 per kg
		cross_ax = (1000./120.3)*cross_PE[:,1]*cross_PE[:,0]**2*betafac*3/(16*pi*alpha*me**2)
		#graph for interpolation
		gr_cross_ax = TGraph(int(cross_ax.shape[0]), cross_PE[:,0].astype(float), cross_ax.astype(float))

		class cross_axio:
			def __call__( self, x, par ):

				return par[0] + gr_cross_ax.Eval(x[0])

		fcross_axio = TF1("cross_axio_mA_" + str(mass), cross_axio(), 0,30, 1)
		fcross_axio.SetParameter(0,0)
		fcross_axio.SetNpx(5000)
		fcross_axio.Write()
		del fcross_axio

	fout.Close()



def get_CBRD_func():
	"""Interpolate J. Redondo data points to get signal flux
	
	Detail:
		void

	Args:
		void

	Returns:
		void

	Raises:
		void
	"""
	
	flux_arr = np.loadtxt("../Text_files/CBRD_points.txt", skiprows=9).astype(float)
	gAe=1.
	flux_arr[:,1] = 1E19*flux_arr[:,1]*(gAe/(0.511*1E-10))**2
	gr_CBRD = TGraph(int(flux_arr.shape[0]), flux_arr[:,0].astype(float), flux_arr[:,1].astype(float))

	class CBRD:
		def __call__( self, x, par ):

			return par[0] + gr_CBRD.Eval(x[0])

	fCBRD = TF1("CBRD", CBRD(), 0.1,12, 1)
	fCBRD.SetParameter(0,0)
	fCBRD.SetNpx(5000)

	fout = TFile("../ROOT_files/Axion/CBRD/CBRD_flux_per_cm2daykeV_gAeto1.root", "recreate")
	fCBRD.Write()
	fout.Close()


def get_cross_CBRD_prod():
	"""Multiply flux by cross section
	
	Detail:
		void

	Args:
		void

	Returns:
		void
          
	Raises:
		void
	"""
	
	file_CBRD = TFile("../ROOT_files/Axion/CBRD/CBRD_flux_per_cm2daykeV_gAeto1.root", "read")
	func_CBRD = file_CBRD.Get("CBRD")

	
	try :
		os.remove("../ROOT_files/Axion/CBRD/CBRD_cross_section.root")
	except OSError:
		pass

	list_mass = [0, 0.01, 0.05, 0.1, 0.15, 0.2]
	for mass in list_mass :

		file_cross_ax = TFile("../ROOT_files/Axion/Axio_electric/axio_electric_cross_section_cm2perkg.root", "read")
		fcross_ax     = file_cross_ax.Get("cross_axio_mA_" + str(mass))

		class CBRD_cross:
			def __call__( self, x, par ):
				return par[0] + func_CBRD.Eval(x[0])*fcross_ax.Eval(x[0])

		fCBRD_cross = TF1("CBRD_cross_mass_" + str(mass), CBRD_cross(), 0.1,12, 1)
		fCBRD_cross.SetParameter(0,0)
		fCBRD_cross.SetNpx(5000)

		fout = TFile("../ROOT_files/Axion/CBRD/CBRD_cross_section.root", "update")
		fCBRD_cross.Write()
		fout.Close()

def get_trigger_efficiency_function(list_bolo_name):
	"""Get trigger efficiency function
	
	Detail:
		void

    Args:
        list_bolo_name  = (str) list of bolo names

	Returns:
		void

	Raises:
		void
	"""

	for bolo_name in list_bolo_name :
		eff_arr = np.loadtxt("../Tables/trigger_efficiency_" + bolo_name + ".txt", skiprows=2).astype(float)

		gr_trigger = TGraph(int(eff_arr.shape[0]), eff_arr[:,0].astype(float), eff_arr[:,1].astype(float))
		class trigger:
			def __call__( self, x, par ):

				if x[0]<5 :
					return par[0] + gr_trigger.Eval(x[0])
				else : 
					return 1 + par[0]

		ftrigger = TF1("trigger_eff", trigger(), 0, 30, 1)
		ftrigger.SetParameter(0,0)
		ftrigger.SetNpx(5000)

		fout = TFile("../ROOT_files/Axion/Trigger_eff/" + bolo_name + "_trigger_eff.root", "recreate")
		ftrigger.Write()
		fout.Close()


def get_cut_efficiency_bckg(list_bolo_name):
	"""Get trigger efficiency function
	
	Detail:
		void

    Args:
        list_bolo_name  = (str) list of bolo names

	Returns:
		void

	Raises:
		void
	"""

	for bolo_name in list_bolo_name :

		t, f = PyRPl.open_ROOT_object("./Event_generation/ROOT_files/" + bolo_name + "_FidGamma_tree.root", "t_new")
		cut = ut.get_electronic_cut_line(bolo_name)

		hcut = TH1F("hcut" + bolo_name, "hcut" + bolo_name, 200, 0, 10)
		hnocut = TH1F("hnocut" + bolo_name, "hnocut" + bolo_name, 200, 0, 10)
		heff = TH1F("heff" + bolo_name, "heff" + bolo_name, 200, 0, 10)

		t.Project("hcut" + bolo_name, "Ebest", cut)
		t.Project("hnocut" + bolo_name, "Ebest")

		for i in range(1, 201):
			if hnocut.GetBinContent(i) > 0 :
				heff.SetBinContent(i, hcut.GetBinContent(i)/float(hnocut.GetBinContent(i)))
			else :
				heff.SetBinContent(i, 0)

		heff.Smooth(50)
		# heff.Draw()
		s = TSpline3(heff)
		s.SetLineColor(kRed)
		# s.Draw("same")

		class cut_eff:
			def __call__( self, x, par ):

				if x[0]<2:
					return s.Eval(x[0]) + par[0]
				else :
					return s.Eval(2)

		fcuteff = TF1("cut_eff", cut_eff(), 0, 30, 1)
		fcuteff.SetParameter(0,0)
		fcuteff.SetNpx(5000)

		fout = TFile("../ROOT_files/Axion/Cut_eff/" + bolo_name + "_cut_eff.root", "recreate")
		fcuteff.Write()
		fout.Close()

def get_cut_efficiency_signal(list_bolo_name):
	"""Get trigger efficiency function
	
	Detail:
		void

    Args:
        list_bolo_name  = (str) list of bolo names

	Returns:
		void

	Raises:
		void
	"""

	for bolo_name in list_bolo_name :
		try :
			os.remove("../ROOT_files/Axion/Cut_eff_signal/" + bolo_name + "_cut_eff_signal.root")
		except OSError:
			pass
		script_utils.print_utility("Computing for bolo " + bolo_name)
		for mass in [0, 0.01, 0.05, 0.1, 0.15, 0.2]:
			t, f = PyRPl.open_ROOT_object("./Event_generation/ROOT_files/" + bolo_name + "_signal_tree.root", "t_new" + str(mass))
			cut = ut.get_electronic_cut_line(bolo_name)

			hcut = TH1F("hcut" + bolo_name, "hcut" + bolo_name, 200, 0, 10)
			hnocut = TH1F("hnocut" + bolo_name, "hnocut" + bolo_name, 200, 0, 10)
			heff = TH1F("heff" + bolo_name, "heff" + bolo_name, 200, 0, 10)

			t.Project("hcut" + bolo_name, "Ebest", cut)
			t.Project("hnocut" + bolo_name, "Ebest")

			for i in range(1, 201):
				if hnocut.GetBinContent(i) > 0 :
					heff.SetBinContent(i, hcut.GetBinContent(i)/float(hnocut.GetBinContent(i)))
				else :
					heff.SetBinContent(i, 0)

			hcut.SetMaximum(1.1*hnocut.GetMaximum())
			hcut.SetLineColor(kRed)
			hcut.Draw()
			hnocut.Draw("same")
			raw_input()

			heff.Smooth(50)
			# heff.Draw()
			s = TSpline3(heff)
			s.SetLineColor(kRed)
			s.Draw("same")
			raw_input()

			class cut_eff:
				def __call__( self, x, par ):

					if x[0]<2:
						return s.Eval(x[0]) + par[0]
					else :
						return s.Eval(2)

			fcuteff = TF1("cut_eff_signal_mass_" + str(mass), cut_eff(), 0, 30, 1)
			fcuteff.SetParameter(0,0)
			fcuteff.SetNpx(5000)

			fout = TFile("../ROOT_files/Axion/Cut_eff_signal/" + bolo_name + "_cut_eff_signal.root", "update")
			fcuteff.Write()
			fout.Close()

def get_skimmed_axion_tree(list_bolo_name):
    """Skim ROOT TTree to get relevant data for axion analysis
    
    Detail:
        void

    Args:
        list_bolo_name  = (str) list of bolo names
        
    Returns:
        void

    Raises:
        void
    """

    for bolo_name in list_bolo_name :
		script_utils.print_utility("Processing bolo " + bolo_name)
		ut.get_skimmed_axion(bolo_name)


def get_elec_recoil_spectrum(list_bolo_name):
    """Create ROOT file with elec recoil spectrum TH1F
    
    Detail:
        void

    Args:
        list_bolo_name  = (str) list of bolo names
        
    Returns:
        void

    Raises:
        void
    """

    for bolo_name in list_bolo_name :
    	#Get data tree
		t, f = PyRPl.open_ROOT_object("../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root", "t_axion")
		#Get cut
		cut = ut.get_electronic_cut_line(bolo_name)
		#Create hist and hist file
		fout = TFile("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "recreate")
		h = TH1F("h_" + bolo_name , "h_" + bolo_name, 200, 0, 30)
		#Save to file
		t.Project("h_" + bolo_name, "Ebest", cut)
		h.Scale(1./h.GetBinWidth(5))
		h.Write()
		fout.Close()
		del h 
		del fout


def refine_Eric_bckgmodel(list_bolo_name):
	"""Create ROOT file with elec recoil spectrum bckg model TF1
	
    Detail: 
		void

	Args:
		list_bolo_name = (str) list of bolo names
		
	Returns:
		void

	Raises:
		void
	"""

	for bolo_name in list_bolo_name :

		ftrigger, file_trigger = PyRPl.open_ROOT_object("../ROOT_files/Axion/Trigger_eff/" + bolo_name + "_trigger_eff.root", "trigger_eff")
		fcut, file_cut = PyRPl.open_ROOT_object("../ROOT_files/Axion/Cut_eff/" + bolo_name + "_cut_eff.root", "cut_eff")

		#Get data tree
		t, f = PyRPl.open_ROOT_object("../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root", "t_axion")
		#Get cut
		cut  = ut.get_electronic_cut_line(bolo_name)
		#Create hist of FWbest reso and get theoretical 0 keV resolution
		hFWbest = TH1F("hFWbest_" + bolo_name , "hFWbest_" + bolo_name, 100, 0, 15)
		t.Project("hFWbest_" + bolo_name, "FWbest")
		reso_0 = hFWbest.GetMean()/2.3548

		##############
		#Fit functions
		##############
		
		class Flat_fit:
			#Use this class to fit the flat part
			def __call__( self, x, par ):
				return par[0]

		class Reso_fit:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[4]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = par[2]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = 0.097*par[2]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98]
				return flat_part + sum(list_peaks)

		class full_fit:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[2]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = 0.097*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				peak_7_71  = par[5]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso**2)
				peak_7_11  = par[6]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso**2)
				peak_6_54  = par[7]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso**2)
				peak_5_99  = par[8]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso**2)
				peak_5_46  = par[9]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso**2)
				peak_4_97  = par[10]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso**2)
				peak_1_3   = 0.1175*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.297)**2/reso**2)
				peak_1_3   += 0.119*par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.0961)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return (sum(list_peaks) + flat_part)*ftrigger.Eval(x[0])*fcut.Eval(x[0])

		class full_fit_noeff:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[2]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = 0.097*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				peak_7_71  = par[5]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso**2)
				peak_7_11  = par[6]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso**2)
				peak_6_54  = par[7]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso**2)
				peak_5_99  = par[8]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso**2)
				peak_5_46  = par[9]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso**2)
				peak_4_97  = par[10]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso**2)
				peak_1_3   = 0.1175*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.297)**2/reso**2)
				peak_1_3   += 0.119*par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.0961)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return (sum(list_peaks) + flat_part)

		class full_fit_deconv:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[2]
				reso_deconv = TMath.Abs(par[1]*x[0])+1E-3
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso_deconv)

				peak_10_37 = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso_deconv**2)
				peak_9_66  = 0.097*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso_deconv**2)
				peak_8_98  = par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso_deconv**2)
				peak_7_71  = par[5]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso_deconv**2)
				peak_7_11  = par[6]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso_deconv**2)
				peak_6_54  = par[7]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso_deconv**2)
				peak_5_99  = par[8]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso_deconv**2)
				peak_5_46  = par[9]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso_deconv**2)
				peak_4_97  = par[10]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso_deconv**2)
				peak_1_3   = 0.1175*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.297)**2/reso_deconv**2)
				peak_1_3   += 0.119*par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.0961)**2/reso_deconv**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return sum(list_peaks) + flat_part


		h = TH1F("h_" + bolo_name , "h_" + bolo_name, 200, 0, 15)
		t.Project("h_" + bolo_name, "Ebest", cut)
		h.Scale(1./h.GetBinWidth(5))

		fflat = TF1("flat" + bolo_name, Flat_fit(), 2, 4, 1)
		h.Fit("flat" + bolo_name, "LL", "", 2, 4)

		fReso = TF1("Reso_fit" + bolo_name, Reso_fit(), 8,11, 5)
		fReso.SetParameters(1,1,1,1,1)
		fReso.FixParameter(4,fflat.GetParameter(0))
		h.Fit("Reso_fit" + bolo_name, "", "", 8,11)

		ffull = TF1(bolo_name+"_FidGamma", full_fit(), 0, 15, 11)
		ffull.SetNpx(300)
		ffull.SetParameters(1,1,1,1,1,1,1,1,1,1,1)
		for i in range(3, 11) :
			ffull.SetParLimits(i, 0, 10000)
		ffull.FixParameter(2, fflat.GetParameter(0))
		ffull.FixParameter(1, fReso.GetParameter(1))

		h.Fit(bolo_name+"_FidGamma", "LL", "", 2, 11)
		ffull.Draw("same")

		ffull_deconv = TF1(bolo_name + "_FidGamma_deconv", full_fit_deconv(), 0, 15, 11)
		ffull_deconv.SetNpx(300)
		for i in range(11):
			ffull_deconv.SetParameter(i,ffull.GetParameter(i))


		ffull_noeff = TF1(bolo_name + "_FidGamma_noeff", full_fit_noeff(), 0, 15, 11)
		ffull_noeff.SetNpx(300)
		for i in range(11):
			ffull_noeff.SetParameter(i,ffull.GetParameter(i))

		fout = TFile("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", "recreate")
		ffull.Write()
		ffull_noeff.Write()
		ffull_deconv.Write()
		fout.Close()

def refine_Eric_bckgmodel_with_tritium(list_bolo_name, d_option):
	"""Create ROOT file with elec recoil spectrum bckg model TF1
	
    Detail: 
		void

	Args:
		list_bolo_name = (str) list of bolo names
		d_option = (dict) dict of fitting option per bolo
		
	Returns:
		void

	Raises:
		void
	"""

	for bolo_name in list_bolo_name :

		ftrigger, file_trigger = PyRPl.open_ROOT_object("../ROOT_files/Axion/Trigger_eff/" + bolo_name + "_trigger_eff.root", "trigger_eff")
		fcut, file_cut = PyRPl.open_ROOT_object("../ROOT_files/Axion/Cut_eff/" + bolo_name + "_cut_eff.root", "cut_eff")

		#Get data tree
		t, f = PyRPl.open_ROOT_object("../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root", "t_axion")
		#Get cut
		cut  = ut.get_electronic_cut_line(bolo_name)
		#Create hist of FWbest reso and get theoretical 0 keV resolution
		hFWbest = TH1F("hFWbest_" + bolo_name , "hFWbest_" + bolo_name, 100,-2, 2)
		t.Project("hFWbest_" + bolo_name, "FWbest")
		reso_0 = hFWbest.GetMean()/2.3548

		##############
		#Fit functions
		##############

		class Beta_tritium_noeff:
			def __call__( self, x, par ):
				m_e = 9.10938215E-31 #in kg
				c = 299792458 # in m / s
				Z = 2 # for tritium
				alpha = 1/137. # fine structure constant
				E_0 = 18.6 * 1E3 * 1.6E-19 # tritium endpoint in Joule
				T = x[0]*1E3*1.6E-19 #kinetic energy in Joule
				E = T + m_e * np.power(c,2)
				p = np.sqrt( np.power(T,2) + 2*T*m_e*np.power(c,2) )
				eta = 2 * np.pi * alpha * Z * E / (p*c+1E-10)
				F = 2 * np.pi * eta / (1 - np.exp(-2*np.pi * eta))
				if T <E_0 :
					return par[0]*(1/4.37353532216e-56)*F*p*E*np.power(E_0-T,2) #this constant to have unit integral when par[0] = 1
				else :
					return 0

		class Beta_tritium:
			def __call__( self, x, par ):
				m_e = 9.10938215E-31 #in kg
				c = 299792458 # in m / s
				Z = 2 # for tritium
				alpha = 1/137. # fine structure constant
				E_0 = 18.6 * 1E3 * 1.6E-19 # tritium endpoint in Joule
				T = x[0]*1E3*1.6E-19 #kinetic energy in Joule
				E = T + m_e * np.power(c,2)
				p = np.sqrt( np.power(T,2) + 2*T*m_e*np.power(c,2) )
				eta = 2 * np.pi * alpha * Z * E / (p*c+1E-10)
				F = 2 * np.pi * eta / (1 - np.exp(-2*np.pi * eta))
				if T <E_0 :
					return par[0]*(1/4.37353532216e-56)*F*p*E*np.power(E_0-T,2)*ftrigger.Eval(x[0])*fcut.Eval(x[0]) #this constant to have unit integral when par[0] = 1
				else :
					return 0

		ftritium_noeff = TF1("tritium_noeff", Beta_tritium_noeff(), 0, 20, 1)
		ftritium = TF1("tritium", Beta_tritium(), 0, 20, 1)
		ftritium_noeff.SetNpx(200)
		ftritium.SetNpx(200)
		ftritium_noeff.SetParameter(0,1)		
		ftritium.SetParameter(0,1)

		class Flat_fit:
			#Use this class to fit the flat part
			def __call__( self, x, par ):
				return par[0]

		class Reso_fit:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[4]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = par[2]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = 0.097*par[2]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98]
				return flat_part + sum(list_peaks)

		class full_fit:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):

				tritium_part = par[11]*ftritium_noeff.Eval(x[0])
				flat_part  = par[2]
				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = 0.097*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				peak_7_71  = par[5]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso**2)
				peak_7_11  = par[6]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso**2)
				peak_6_54  = par[7]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso**2)
				peak_5_99  = par[8]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso**2)
				peak_5_46  = par[9]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso**2)
				peak_4_97  = par[10]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso**2)
				peak_1_3   = 0.1175*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.297)**2/reso**2)
				peak_1_3   += 0.119*par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.0961)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return (sum(list_peaks) + flat_part + tritium_part)*ftrigger.Eval(x[0])*fcut.Eval(x[0])

		class full_fit_noeff:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[2]
				tritium_part = par[11]*ftritium_noeff.Eval(x[0])

				reso = np.sqrt(reso_0**2 + (par[1]*x[0])**2)
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso)

				peak_10_37 = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso**2)
				peak_9_66  = 0.097*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso**2)
				peak_8_98  = par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso**2)
				peak_7_71  = par[5]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso**2)
				peak_7_11  = par[6]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso**2)
				peak_6_54  = par[7]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso**2)
				peak_5_99  = par[8]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso**2)
				peak_5_46  = par[9]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso**2)
				peak_4_97  = par[10]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso**2)
				peak_1_3   = 0.1175*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.297)**2/reso**2)
				peak_1_3   += 0.119*par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.0961)**2/reso**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return (sum(list_peaks) + flat_part + tritium_part)

		class full_fit_deconv:
			#Use this class to fit the offset and the resolution broadening
			def __call__( self, x, par ):
				flat_part  = par[2]
				tritium_part = par[11]*ftritium_noeff.Eval(x[0])

				reso_deconv = TMath.Abs(par[1]*x[0])+1E-3
				gauss_fact = 1/(np.sqrt(2*TMath.Pi())*reso_deconv)

				peak_10_37 = par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*10.37)**2/reso_deconv**2)
				peak_9_66  = 0.097*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*9.66)**2/reso_deconv**2)
				peak_8_98  = par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*8.98)**2/reso_deconv**2)
				peak_7_71  = par[5]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.71)**2/reso_deconv**2)
				peak_7_11  = par[6]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*7.11)**2/reso_deconv**2)
				peak_6_54  = par[7]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*6.54)**2/reso_deconv**2)
				peak_5_99  = par[8]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.99)**2/reso_deconv**2)
				peak_5_46  = par[9]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*5.46)**2/reso_deconv**2)
				peak_4_97  = par[10]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*4.97)**2/reso_deconv**2)
				peak_1_3   = 0.1175*par[3]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.297)**2/reso_deconv**2)
				peak_1_3   += 0.119*par[4]*gauss_fact*np.exp(-0.5*(x[0]- par[0]*1.0961)**2/reso_deconv**2)
				list_peaks = [peak_10_37, peak_9_66, peak_8_98, peak_7_71, peak_7_11, peak_6_54, peak_5_99, peak_5_46, peak_4_97, peak_1_3]
				return sum(list_peaks) + flat_part + tritium_part


		h = TH1F("h_" + bolo_name , "h_" + bolo_name, 200, 0, 30)
		t.Project("h_" + bolo_name, "Ebest", cut)
		h.Scale(1./h.GetBinWidth(5))

		fflat = TF1("flat" + bolo_name, Flat_fit(), 19, 30, 1)
		h.Fit("flat" + bolo_name, "LL", "", 19, 30)
		raw_input()

		fReso = TF1("Reso_fit" + bolo_name, Reso_fit(), 10,11, 5)
		fReso.SetParameters(1,1,1,1,1)
		fReso.FixParameter(4,fflat.GetParameter(0))
		h.Fit("Reso_fit" + bolo_name, d_option[bolo_name], "", 8,11)
		raw_input()

		ffull = TF1(bolo_name+"_FidGamma", full_fit(), 0, 30, 12)
		ffull.SetNpx(300)
		# ffull.SetParameters(1,1,1,1,1,1,1,1,1,1,1,1)
		for i in range(3, 12) :
			ffull.SetParLimits(i, 0, 10000)
		ffull.FixParameter(2, fflat.GetParameter(0))
		ffull.FixParameter(1, fReso.GetParameter(1))
		ffull.FixParameter(0, fReso.GetParameter(0))

		h.Fit(bolo_name+"_FidGamma", "LL", "", 2, 11)
		raw_input()
		ffull.Draw("same")

		ffull_deconv = TF1(bolo_name + "_FidGamma_deconv", full_fit_deconv(), 0, 30, 12)
		ffull_deconv.SetNpx(300)
		for i in range(12):
			ffull_deconv.SetParameter(i,ffull.GetParameter(i))


		ffull_noeff = TF1(bolo_name + "_FidGamma_noeff", full_fit_noeff(), 0, 30, 12)
		ffull_noeff.SetNpx(300)
		for i in range(12):
			ffull_noeff.SetParameter(i,ffull.GetParameter(i))

		fout = TFile("../ROOT_files/Axion/Model/" + bolo_name + "_model_with_tritium_perkeV.root", "recreate")
		ffull.Write()
		ffull_noeff.Write()
		ffull_deconv.Write()
		#Also save tritium function 
		ftritium_noeff.SetParameter(0, ffull.GetParameter(11))
		ftritium.SetParameter(0, ffull.GetParameter(11))
		ftritium_noeff.Write()
		ftritium.Write()
		fout.Close()

def stack_background(list_bolo_name) :
	"""Create ROOT file with stacked elec recoil spectrum TH1F
	
	Detail:
		void

	Args:
		list_bolo_name = (str) list of bolo names

	Returns:
		void

    Raises: 
		void
	"""
	
	hbckg = TH1F("hstacked", "hstacked", 200, 0, 30)
	
	for bolo_name in list_bolo_name :
    	#Get data tree
		h, f = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
		for i in range(1, 301):
			hbckg.SetBinContent(i, hbckg.GetBinContent(i) + h.GetBinContent(i))
		del h 
		del f

	fout = TFile("../ROOT_files/Axion/Spec/Stacked_spec_perkeV.root", "recreate")
	hbckg.Write()
	fout.Close()

def stack_bckg_model(list_bolo_name, bool_tritium):
	"""Create ROOT file with stacked background model TF1
	
	Detail:
		void

	Args:
		list_bolo_name = (str) list of bolo names
		bool_tritium = (bool) True if use tritium model

	Returns:
		void

    Raises: 
		void
	"""

	list_file_func = []
	list_func_tritium = []
	fstacked_tritium = ""
	if bool_tritium :
		list_file_func = [PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/%s_model_with_tritium_perkeV.root" % bolo_name, "%s_FidGamma" % bolo_name) for bolo_name in list_bolo_name]
		list_func_tritium = [ff[1].Get("tritium") for ff in list_file_func]
		class Stacked_tritium:
			def __call__( self, x, par ):
				stack = sum([tritium.Eval(x[0]) for tritium in list_func_tritium]) + par[0]
				return stack
		fstacked_tritium = TF1("tritium_stacked", Stacked_tritium(), 0, 30, 1)
		fstacked_tritium.SetParameter(0,0)
		fstacked_tritium.SetNpx(500)
	else :
		list_file_func = [PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", bolo_name + "_FidGamma") for bolo_name in list_bolo_name]

	list_func = [ff[0] for ff in list_file_func]
	class Stacked_GammaBckg:
		def __call__( self, x, par ):
			stack = sum([model.Eval(x[0]) for model in list_func]) + par[0]
			return stack

	fstacked = TF1("model_stacked", Stacked_GammaBckg(), 0, 30, 1)
	fstacked.SetParameter(0,0)
	fstacked.SetNpx(500)
	if bool_tritium :
		fout = TFile("../ROOT_files/Axion/Model/Stacked_model_with_tritium_perkeV.root", "recreate")
		fstacked.Write()
		fstacked_tritium.Write()
		fout.Close()	
	else :
		fout = TFile("../ROOT_files/Axion/Model/Stacked_model_perkeV.root", "recreate")
		fstacked.Write()
		fout.Close()	

def stack_signal(list_bolo_name, d_exposure, list_mass):
	"""Create ROOT file with stacked expected signal TF1
	
	Detail:
		void

	Args:
		list_bolo_name = (str) list of bolo names
		
	Returns:
		void

    Raises: 
		void
	"""

	try :
		os.remove("../ROOT_files/Axion/CBRD_convolved/Stacked_flux_perkeV.root")
	except OSError :
		pass

	for mass in list_mass :
		list_file_func = [PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/" + bolo_name + "_flux_perkgdkeV.root", bolo_name + "_flux_mass_" + str(mass)) for bolo_name in list_bolo_name]
		list_func = [ff[0] for ff in list_file_func]
		class Stacked_signal:
			def __call__( self, x, par ):
				stack = sum([d_exposure[list_bolo_name[i]]*list_func[i].Eval(x[0]) for i in range(len(list_func))]) + par[0]
				return stack


		fstacked = TF1("flux_stacked_mass_" + str(mass), Stacked_signal(), 0, 15, 1)
		fstacked.SetParameter(0,0)
		fstacked.SetNpx(500)
		fout = TFile("../ROOT_files/Axion/CBRD_convolved/Stacked_flux_perkeV.root", "update")
		fstacked.Write()
		fout.Close()	

if __name__ == "__main__":

	d_exposure = {}
	with open("../Text_files/exposure.txt", "r") as f:
		stuff = f.readlines()
		for line in stuff :
			line = line.rstrip().split(",")
			d_exposure[line[0]] = float(line[1])

	list_mass = [0, 0.01, 0.05, 0.1, 0.15, 0.2]

	list_bolo_name = ["FID824", "FID825", "FID827", "FID837", "FID838", "FID839", "FID841", "FID842"]
	#Use different fitting options to ensure proper fit
	d_option = {}
	for bolo_name in list_bolo_name :
		d_option[bolo_name] = "LL"
	d_option["FID837"] = ""
	d_option["FID824"] = ""
	d_option["FID838"] = ""
	d_option["FID839"] = ""
	# refine_Eric_bckgmodel(list_bolo_name)
	# refine_Eric_bckgmodel_with_tritium(list_bolo_name, d_option)
	# get_axio_electric()
	# get_CBRD_func()
	# get_cross_CBRD_prod()
	# get_trigger_efficiency_function(list_bolo_name)
	# get_cut_efficiency_bckg(list_bolo_name)
	# get_cut_efficiency_signal(list_bolo_name)
	# get_skimmed_axion_tree(list_bolo_name)
	# get_elec_recoil_spectrum(list_bolo_name)	
	# stack_background(list_bolo_name)	
	stack_bckg_model(list_bolo_name, True)	
	# stack_signal(list_bolo_name, d_exposure, list_mass)
