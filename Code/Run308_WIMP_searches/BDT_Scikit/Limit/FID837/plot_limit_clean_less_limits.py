from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Poisson_90CL as Poisson90


def get_percentile(arr, WIMP_mass,analysis, exposure, detector_mass = 0.6):

	"""See detail
	
	Detail:
		From an array input of 90 CL: 
		extract the percentiles to get the expected 
		95 CL and 68 CL on the 90 CL result

	Args:
		arr           = (np array) array of 90 CL values on the number of WIMPs
		WIMP_mass     = (str) WIMP mass in GeV
		analysis      = (str) type of analysis cut
		exposure      = (float) exposure in days
		detector_mass = (float) fiducial mass
		
	Returns:
		d (dict) = the percentile dictionnary

	Raises:
		void
	"""

	#Convert the obs counts in arr to a 90CL Poisson limit
	for i in range(len(arr)):
		arr[i] = Poisson90.compute_90CL_limit(arr[i])

	#Fill a conversion dictionnary to convert Ncounts to limit on cross section
	d ={}
	count_arr = np.loadtxt("./Text_files/WIMP_counts_for_1kgday_" + analysis + ".txt", delimiter = ",")

	d_conv={}
	for mass, Nevents in count_arr:
		d_conv[str(int(mass))] = float(Nevents)

	d["5"]=np.percentile(arr, 2.5)*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass)
	d["95"]=np.percentile(arr, 97.5)*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass)
	d["32"]=np.percentile(arr, 16)*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass)
	d["68"]=np.percentile(arr, 84)*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass)

	return d

def get_limit(bolo_name, arr, WIMP_mass,analysis, exposure, detector_mass = 0.6):

	"""See detail
	
	Detail:
		From a single observation, get the 90 CL limit

	Args:
		bolo_name = (str) bolometer name
		arr           = (np array) array of 90 CL values on the number of WIMPs
		WIMP_mass     = (str) WIMP mass in GeV
		analysis      = (str) type of analysis cut
		exposure      = (float) exposure in days
		detector_mass = (float) fiducial mass
		
	Returns:
		d (dict) = the percentile dictionnary

	Raises:
		void
	"""

	# #Original number of simulated events
	# d_N_ori         = {"5":15000000, "6":900000, "7":200000, "10":40000, "25":10000}

	# # Load WIMP trees to compute the WIMP events passing the cuts
	# WIMP_path = "../../BDT_" + bolo_name + "/" + analysis + "/WIMP/ROOT_files/"
	# tWIMP5, fWIMP5 = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_5_tree.root", "t_new0")
	# tWIMP6, fWIMP6 = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_6_tree.root", "t_new0")
	# tWIMP7, fWIMP7 = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_7_tree.root", "t_new0")
	# tWIMP10, fWIMP10 = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_10_tree.root", "t_new0")
	# tWIMP25, fWIMP25 = PyRPl.open_ROOT_object(WIMP_path + bolo_name + "_WIMP_mass_25_tree.root", "t_new0")

	# d_N_after       = {}
	# d_N_after["5"]  = tWIMP5.GetEntries()
	# d_N_after["6"]  = tWIMP6.GetEntries()
	# d_N_after["7"]  = tWIMP7.GetEntries()
	# d_N_after["10"] = tWIMP10.GetEntries()
	# d_N_after["25"] = tWIMP25.GetEntries()


	#Convert the obs counts in arr to a 90CL Poisson limit
	for i in range(len(arr)):
		arr[i] = Poisson90.compute_90CL_limit(arr[i])

	#Get the array with the efficiency
	count_eff = np.loadtxt("./Text_files/" + bolo_name + "_BDT_cut_and_eff_" + analysis + ".txt", delimiter = ",")
	d_eff={}
	for mass,truc, eff in count_eff:
		d_eff[str(int(mass))] = float(eff)


	count_arr = np.loadtxt("./Text_files/WIMP_counts_for_1kgday_" + analysis + "_2D.txt", delimiter = ",")
	#Fill a conversion dictionnary to convert Ncounts to limit on cross section
	# #Normalise with 1D PDF and with generated events
	# d_conv={}
	# for mass, Nevents in count_arr:
	# 	d_conv[str(int(mass))] = float(Nevents)

	# return arr[0]*1E-5*d_N_ori[str(int(WIMP_mass))]/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass * d_N_after[str(int(WIMP_mass))])

	# Normalise with 2D PDF
	d_conv={}
	for mass, Nevents in count_arr:
		d_conv[str(int(mass))] = float(Nevents)

	# print arr[0], WIMP_mass, d_conv[str(int(WIMP_mass))]

	print WIMP_mass, arr[0], d_conv[str(int(WIMP_mass))] , exposure , detector_mass, arr[0]*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass)

	return arr[0]*1E-5/(d_conv[str(int(WIMP_mass))] * exposure * detector_mass* d_eff[str(int(WIMP_mass))] )

def get_expected_band_limit(list_mass, analysis, exposure, detector_mass = 0.6):

	"""Get expected sensitivity confidence bands
	
	Detail:
		Return filled area graph for the 68 and 95 CI on the 
		expected sensitivity given a choice of analysis

	Args:
		list_mass = (list) list of WIMP masses
		analysis      = (str) type of analysis cut
		exposure = (float) exposure in days
		detector_mass = (float) fiducial mass

	Returns:
		gr_68, gr_95 (TGraph) = the 2 filled area graphs

	Raises:
		void
	"""

	d_mass_percentile={}
	list_arr_obs = np.loadtxt("./Text_files/events_passing_cuts_" + analysis + ".txt", delimiter = ",")
	for i in range(list_arr_obs.shape[0]):
		WIMP_mass = list_arr_obs[i][0]
		arr_counts = list_arr_obs[i][1:]
		d_mass_percentile[int(WIMP_mass)] = get_percentile(arr_counts, WIMP_mass, analysis, exposure, detector_mass)

	# #Get an array of mass for contour ( in the form: 6,7,8,8,7,6 for instance)
	list_mass_gr =np.array(list_mass + list_mass[::-1]).astype(float)

	list_95 = [d_mass_percentile[mass]["5"] for mass in list_mass] + [d_mass_percentile[mass]["95"] for mass in list_mass[::-1]]
	list_95 = np.array(list_95).astype(float)
	gr_95 = TGraph(len(list_mass_gr), list_mass_gr, list_95)

	list_68 = [d_mass_percentile[mass]["32"] for mass in list_mass] + [d_mass_percentile[mass]["68"] for mass in list_mass[::-1]]
	list_68 = np.array(list_68).astype(float)
	gr_68 = TGraph(len(list_mass_gr), list_mass_gr, list_68)

	# print list_68

	return gr_68, gr_95

def get_true_event_limit(bolo_name, list_mass, analysis, exposure, detector_mass = 0.6):

	"""Get expected sensitivity confidence bands
	
	Detail:
		Return filled area graph for the 68 and 95 CI on the 
		expected sensitivity given a choice of analysis

	Args:
		bolo_name = (str) bolometer name
		list_mass = (list) list of WIMP masses
		analysis      = (str) type of analysis cut
		exposure = (float) exposure in days
		detector_mass = (float) fiducial mass

	Returns:
		gr_68, gr_95 (TGraph) = the 2 filled area graphs

	Raises:
		void
	"""

	d_mass_limit={}
	list_arr_obs = np.loadtxt("./Text_files/true_events_passing_cuts_" + analysis  +".txt", delimiter = ",")
	for i in range(list_arr_obs.shape[0]):
		WIMP_mass = list_arr_obs[i][0]
		arr_counts = list_arr_obs[i][1:]
		d_mass_limit[int(WIMP_mass)] = get_limit(bolo_name, arr_counts, WIMP_mass, analysis, exposure, detector_mass)

	list_lim = [d_mass_limit[int(mass)] for mass in list_mass]

	arr_lim = np.array([[list_mass[i], list_lim[i]] for i in range(len(list_mass))])
	np.savetxt("./Text_files/" + bolo_name + "_" + analysis + "_limit.txt", arr_lim, delimiter = " ")

	gr_lim = TGraph(len(list_mass), np.array(list_mass).astype(float), np.array(list_lim).astype(float))

	# # print list_mass, list_lim
	for i in range(len(list_mass)):
		print list_mass[i], list_lim[i]

	return gr_lim

def get_contour_graph(file_name,  color):

	"""Open the file to get the graph of the contour limit
	
	Detail:
		Return the limit graph of given file_name 
		Also modify color

	Args:
		file_name (str)    = the name of the limit file 
		line_width (int)   = the line width 
		color (ROOT color) = the color for the line
		
	Returns:
		void

	Raises:
		void
	"""


	if "cogent_contour" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E36*np.array(list(arr[:,1])))
		gr.SetFillColor(color)
	elif "cresst" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-6*np.array(list(arr[:,1])))
		gr.SetFillColor(color)
	elif "dama" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-7*np.array(list(arr[:,1])))
		gr.SetFillColor(color)

	elif "cdms" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		for i in range(n_points):
			arr[i,0] = np.power(10, arr[i,0])
			arr[i,1] = 1E36*np.power(10, arr[i,1])
			# print i, arr[i,0], arr[i,1]
		gr = TGraph(n_points, np.array(list(arr[:,0])), np.array(list(arr[:,1])))
		gr.SetFillColor(color)

	elif "cns" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-14*np.array(list(arr[:,1])))
		gr.SetFillColor(color)

	else:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), np.array(list(arr[:,1])))
		gr.SetFillColor(color)

	return gr

def get_limit_graph(file_name, line_width, color):

	"""Open the file to get the graph of the limit 
	for the given file name
	
	Detail:
		Return the limit graph of given file_name 
		Also modify line_width and color

	Args:
		file_name (str)    = the name of the limit file 
		line_width (int)   = the line width 
		color (ROOT color) = the color for the line
		
	Returns:
		void

	Raises:
		void
	"""

	if "cdms" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-10*np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	elif "coupp" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E36*np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	elif "cdmlite" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-3*np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	elif "pico" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-7*np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	elif "xenon10s2" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E36*np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	elif "xenon100" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-10*np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	elif "cns" in file_name:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), 1E-14*np.array(list(arr[:,1])))
		print arr
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	else:
		arr = np.loadtxt(file_name, delimiter = ",") 
		n_points                  =int(arr[:,0].shape[0])
		gr = TGraph(n_points, np.array(list(arr[:,0])), np.array(list(arr[:,1])))
		gr.SetLineColor(color)
		gr.SetLineWidth(line_width)

	return gr

def plot_limit(bolo_name, list_mass, analysis_type, exposure, detector_mass = 0.6):

	"""Plot all the desired limits
	
	Detail:
		Modify the code to add more limits if needed

	Args:
		bolo_name = (str) bolometer name
		list_mass = (list) list of WIMP masses
		analysis_type      = (str) type of analysis cut
		exposure      = (float) exposure in days
		detector_mass = (float) fiducial mass
		
	Returns:
		void
	Raises:
		void
	"""


	# gr_68, gr_95            = get_expected_band_limit(list_mass, analysis_type, exposure, detector_mass)

	gr_lim = get_true_event_limit(bolo_name, list_mass, analysis_type, exposure, detector_mass)
	# gr_lim10 = get_true_event_limit(bolo_name, list_mass, analysis_type, 10*exposure, detector_mass)
	PyRPl.process_TGraph(gr_lim, color = kBlue)

	gr_edw_moriond = get_limit_graph("./Text_files/edw3_ana_1.5_0_5_moriond.txt", 2, kRed)
	gr_edw_poisson = get_limit_graph("./Text_files/edw3_ana_1.5_0_5_poisson.txt", 2, kBlack)
	gr_edw_quentin = get_limit_graph("./Text_files/edw3_quentin_0.5keVthresh.txt", 2, kAzure+9)
	gr_edw_low = get_limit_graph("./Text_files/Published_limits/edw_lowmass_2012.txt", 2, kRed)
	gr_edw_moriond.SetLineStyle(7)
	gr_cns = get_limit_graph("./Text_files/Published_limits/cns.txt", 3, kBlack)
	gr_cns_contour = get_contour_graph("./Text_files/Published_limits/cns_contour.txt",  kGray)
	gr_cdms = get_limit_graph("./Text_files/Published_limits/cdms_limit.txt", 3, kBlue-3)
	gr_cdmlite = get_limit_graph("./Text_files/Published_limits/cdmlite2013.txt", 3, kBlue-3)
	gr_cdmlite.SetLineStyle(7)
	gr_cdms_contour = get_contour_graph("./Text_files/Published_limits/cdms_contour_90CL.txt", kBlue-10)
	gr_dama_contour = get_contour_graph("./Text_files/Published_limits/dama2009.txt", kRed-10)
	gr_cogent_contour = get_contour_graph("./Text_files/Published_limits/cogent_contour_2013.txt", kOrange+1)
	gr_cressta = get_contour_graph("./Text_files/Published_limits/cresst11a.txt", kGreen+2)
	gr_cresstb = get_contour_graph("./Text_files/Published_limits/cresst11b.txt", kGreen+2)
	gr_lux = get_limit_graph("./Text_files/Published_limits/lux2013.txt", 3, kBlack)
	gr_xenon100 = get_limit_graph("./Text_files/Published_limits/xenon100.txt", 2, kRed-5)
	gr_simple = get_limit_graph("./Text_files/Published_limits/simple_2012.txt", 2, kOrange)
	gr_coupp = get_limit_graph("./Text_files/Published_limits/coupp_2012.txt", 2, kGray)
	gr_cresst = get_limit_graph("./Text_files/Published_limits/cresst2014.txt", 2, kGreen+2)
	gr_pico = get_limit_graph("./Text_files/Published_limits/pico2015.txt", 2, kCyan-7)

	h = TH1F("h", "", 100, 1,25)
	PyRPl.process_TH1(h, X_title = "Mass (GeV)", Y_title = "#sigma (pb)", X_title_size = .06, Y_title_size = .06, X_title_offset = .98, Y_title_offset = .95)

	h.SetMinimum(5E-8)
	h.SetMaximum(1E-1)
	h.GetXaxis().SetTitleOffset(1.17)

	# gr_95.SetFillColor(kBlue-9)
	# gr_68.SetFillColor(kBlue-4)

	# gr_95strict.SetFillColor(kOrange-4)
	# gr_68strict.SetFillColor(kOrange-9)

	gr_edw_low.SetName("gr_edw_low")
	gr_edw_poisson.SetName("gr_edw_poisson")
	gr_edw_moriond.SetName("gr_edw_moriond")
	gr_edw_quentin.SetName("gr_edw_quentin")
	gr_cns.SetName("gr_cns")
	gr_cns_contour.SetName("gr_cns_contour")
	gr_cdms.SetName("gr_cdms")
	gr_cdmlite.SetName("gr_cdmlite")
	gr_cdms_contour.SetName("gr_cdms_contour")
	gr_dama_contour.SetName("gr_dama_contour")
	gr_cogent_contour.SetName("gr_cogent_contour")
	gr_cressta.SetName("gr_cressta")
	gr_cresstb.SetName("gr_cresstb")
	gr_cresst.SetName("gr_cresst")
	gr_pico.SetName("gr_pico")
	gr_lux.SetName("gr_lux")
	gr_xenon100.SetName("gr_xenon100")
	gr_simple.SetName("gr_simple")
	gr_coupp.SetName("gr_coupp")
	# gr_68.SetName("gr_68")
	# gr_95.SetName("gr_95")
	gr_lim.SetName("gr_lim")

	gr_lim.SetLineColor(kRed)
	# gr_lim10.SetLineColor(kRed)
	gr_lim.SetLineWidth(5)

	cc = TCanvas("cc", "cc")
	gPad.SetLogy()
	gPad.SetLogx()
	h.GetXaxis().SetTickLength(-0.05)
	h.SetMaximum(1E-1)
	h.SetMinimum(1E-10)
	h.Draw()

	# gr_68.Draw("FC")
	# raw_input()

	# gr_95.Draw("sameF")
	# gr_68.Draw("sameF")
	gr_cdms_contour.Draw("sameF")
	gr_dama_contour.Draw("sameF")
	gr_cogent_contour.Draw("sameF")
	gr_cressta.Draw("sameF")
	gr_cresstb.Draw("sameF")
	gr_cns_contour.Draw("sameF")

	# gr_95strict.Draw("sameF")
	# gr_68strict.Draw("sameF")

	gr_cdms.Draw("sameC")
	gr_cdmlite.Draw("sameC")
	gr_cresst.Draw("sameC")
	# gr_pico.Draw("sameC")
	gr_lux.Draw("sameC")
	gr_xenon100.Draw("sameC")
	# gr_simple.Draw("sameC")
	# gr_coupp.Draw("sameC")
	# gr_edw_poisson.Draw("sameC")
	gr_edw_moriond.SetLineStyle(1)
	gr_edw_moriond.SetLineWidth(4)
	# gr_edw_moriond.Draw("sameL")
	# gr_edw_quentin.Draw("sameC")
	gr_edw_low.SetLineStyle(7)
	gr_edw_low.Draw("sameC")
	gr_lim.Draw("sameC")
	gr_cns.Draw("sameC")
	# gr_lim10.Draw("sameL")

	leg =TLegend(0.564,0.584,0.83,0.857)
	leg.AddEntry("gr_dama_contour", "DAMA" , "f")
	leg.AddEntry("gr_cdms", "SCDMS" , "l")
	leg.AddEntry("gr_cdmlite", "CDMSLite" , "l")
	leg.AddEntry("gr_cdms_contour", "CDMS Si" , "f")
	leg.AddEntry("gr_cogent_contour", "CoGeNT" , "f")
	leg.AddEntry("gr_lux", "LUX" , "l")
	leg.AddEntry("gr_xenon100", "XENON 100 " , "l")
	leg.AddEntry("gr_simple", "SIMPLE" , "l")
	leg.AddEntry("gr_coupp", "COUPP" , "l")
	leg.AddEntry("gr_cresst", "CRESST" , "l")
	leg.AddEntry("gr_pico", "PICO 2L" , "l")
	leg.AddEntry("gr_edw_low", "EDW II" , "l")
	leg.AddEntry("gr_edw_quentin", "EDW III Quentin" , "l")
	leg.AddEntry("gr_edw_poisson", "EDW III Poisson" , "l")

	# leg.AddEntry("gr_95", "EDW III FID837 expected @ 95% CL" , "f")
	leg.AddEntry("gr_lim", "EDW III FID837" , "l")
	leg.SetFillColor(kWhite)
	leg.SetLineColor(kWhite)
	# leg.Draw()


	t_cdms = TText(0.3,1E-10,"truc")
	t_cdms.SetTextSize(0.03)
	t_cdms.SetTextColor(kBlue-3)
	t_cdms.DrawText(8.9,3.1E-7,"CDMS")

	t_dama = TText(0.3,1E-10,"truc")
	t_dama.SetTextSize(0.03)
	t_dama.SetTextColor(kRed-3)
	t_dama.DrawText(10.9,2E-4,"DAMA")

	t_cresst = TText(0.3,1E-10,"truc")
	t_cresst.SetTextSize(0.03)
	t_cresst.SetTextColor(kGreen+2)
	t_cresst.DrawText(19.9,1E-5,"CRESST")


	t_CoGeNT = TText(0.3,1E-10,"truc")
	t_CoGeNT.SetTextSize(0.03)
	t_CoGeNT.SetTextColor(kOrange+7)
	t_CoGeNT.DrawText(6.4,6.5E-5,"CoGeNT")

	t_lux = TText(0.3,1E-10,"truc")
	t_lux.SetTextSize(0.03)
	t_lux.SetTextColor(kBlack)
	t_lux.DrawText(5.8,3.3E-7,"LUX")

	t_edw2 = TText(0.3,1E-10,"truc")
	t_edw2.SetTextSize(0.03)
	t_edw2.SetTextColor(kRed)
	t_edw2.DrawText(7,7E-4,"EDW-II")

	t_edw3 = TText(0.3,1E-10,"truc")
	t_edw3.SetTextSize(0.03)
	t_edw3.SetTextColor(kRed)
	t_edw3.DrawText(9.9,2.16E-6,"EDW-III")

	raw_input()

	cc.Print("./Figures/FID837_" + analysis_type + "_forBeamer_limit.eps")

analysis_type = "ana_1.5_0_5"
exposure = 66
list_mass = [3, 4, 5, 6, 7, 10, 25]
list_mass = [5, 6, 7, 10, 25]
# list_mass = [5]
bolo_name = "FID837"

plot_limit(bolo_name, list_mass, analysis_type, exposure)




