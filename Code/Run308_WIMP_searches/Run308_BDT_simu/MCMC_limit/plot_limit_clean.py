from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl

def get_percentile(arr, WIMP_mass, exposure):

	"""See details
	
	Detail:
		We have simulated the experiment XX times and extracted
		a 90 CL values on the Number of WIMPs 
		This array is given as argument to the function. 
		We compute the appropriate percentile values to get the expected 68 
		and 95 CL on the sensitivity of the experiment 
		The values of the array are converted appropriately to have a limit 
		on the cross section, not the Number of WIMPs

	Args:
		arr       = (np array) array of 90 CL values on the number of WIMPs
		WIMP_mass = (str) WIMP mass in GeV
		exposure  = (float) exposure in kgdays
		
	Returns:
		d (dict) = the percentile dictionnary

	Raises:
		void
	"""

	#Fill a conversion dictionnary to convert Ncounts to limit on cross section
	d ={}
	count_arr = np.loadtxt("./Text_files/WIMP_counts_for_1kgday.txt", delimiter = ",")

	d_conv={}
	for mass, Nevents in count_arr:
		d_conv[str(int(mass))] = float(Nevents)

	d["5"]=np.percentile(arr, 2.5)*1E-5/(d_conv[str(WIMP_mass)] * exposure)
	d["95"]=np.percentile(arr, 97.5)*1E-5/(d_conv[str(WIMP_mass)] * exposure)
	d["32"]=np.percentile(arr, 16)*1E-5/(d_conv[str(WIMP_mass)] * exposure)
	d["68"]=np.percentile(arr, 84)*1E-5/(d_conv[str(WIMP_mass)] * exposure)

	return d

def get_expected_band_limit(analysis_type, exposure, ERA_name):

	"""Get expected sensitivity confidence bands
	
	Detail:
		Return filled area graph for the 68 and 95 CI on the 
		expected sensitivity given a choice of analysis

	Args:
		analysis_type = (str) name of analysis used
		exposure      = (float) exposure in kgdays
		ERA_name      = (str)  "" or _ERA
		
	Returns:
		gr_68, gr_95 (TGraph) = the 2 filled area graphs

	Raises:
		void
	"""

	list_mass =[6,7,10,25]
	arr_6GeV  = np.loadtxt("./Text_files/6GeV" + "_" + analysis_type + "_simu_90CL_NWIMPs" + ERA_name + ".txt")
	arr_7GeV  = np.loadtxt("./Text_files/7GeV" + "_" + analysis_type + "_simu_90CL_NWIMPs" + ERA_name + ".txt")
	arr_10GeV = np.loadtxt("./Text_files/10GeV" + "_" + analysis_type + "_simu_90CL_NWIMPs" + ERA_name + ".txt")
	arr_25GeV = np.loadtxt("./Text_files/25GeV" + "_" + analysis_type + "_simu_90CL_NWIMPs" + ERA_name + ".txt")
	
	d_6GeV = get_percentile(arr_6GeV, 6,exposure)
	d_7GeV = get_percentile(arr_7GeV, 7,exposure)
	d_10GeV = get_percentile(arr_10GeV, 10,exposure)
	d_25GeV = get_percentile(arr_25GeV, 25,exposure)

	list_mass =np.array([6,7,10,25, 25, 10, 7,6]).astype(float)

	list_95   = np.array([d_6GeV["5"], d_7GeV["5"], d_10GeV["5"], d_25GeV["5"], d_25GeV["95"], d_10GeV["95"], d_7GeV["95"], d_6GeV["95"]]).astype(float)
	gr_95     = TGraph(8, list_mass, list_95)
	list_68   = np.array([d_6GeV["32"], d_7GeV["32"], d_10GeV["32"], d_25GeV["32"], d_25GeV["68"], d_10GeV["68"], d_7GeV["68"], d_6GeV["68"]]).astype(float)
	gr_68     = TGraph(8, list_mass, list_68)

	return gr_68, gr_95

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

	arr = np.loadtxt(file_name, delimiter = ",") 
	n_points                  =int(arr[:,0].shape[0])
	gr = TGraph(n_points, np.array(list(arr[:,0])), np.array(list(arr[:,1])))
	gr.SetLineColor(color)
	gr.SetLineWidth(line_width)

	return gr

def plot_limit(analysis_type, exposure):

	"""Plot all the desired limits
	
	Detail:
		Modify the code to add more limits if needed

	Args:
		analysis_type = (str) name of analysis used
		exposure      = (float) exposure in kgdays
		ERA_name      = (str)  "" or _ERA
		
	Returns:
		void
	Raises:
		void
	"""


	gr_68_ERA, gr_95_ERA = get_expected_band_limit(analysis_type, exposure, "_ERA")
	gr_68, gr_95 = get_expected_band_limit(analysis_type, exposure, "")

	gr_edw_low = get_limit_graph("./Text_files/edw_lowmass_2012.txt", 5, kRed)
	gr_cdms = get_limit_graph("./Text_files/cdms_limit.txt", 5, kBlue)

	h = TH1F("h", "", 100, 5,25)
	PyRPl.process_TH1(h, X_title = "Mass (GeV)", Y_title = "#sigma (pb)", X_title_size = .06, Y_title_size = .06, X_title_offset = .98, Y_title_offset = .95)

	h.SetMinimum(5E-8)
	h.SetMaximum(1E-1)

	gr_95_ERA.SetFillColor(kRed-4)
	gr_68_ERA.SetFillColor(kRed-9)

	gr_95.SetFillColor(kBlue-4)
	gr_68.SetFillColor(kBlue-9)

	gr_edw_low.SetName("gr_edw_low")
	gr_cdms.SetName("gr_cdms")
	gr_68_ERA.SetName("gr_68_ERA")

	cc = TCanvas("cc", "cc")
	gPad.SetLogy()
	gPad.SetLogx()
	h.Draw()


	# gr_95.Draw("sameF")
	gr_68.Draw("sameF")

	# gr_95_ERA.Draw("sameF")
	gr_68_ERA.Draw("sameF")

	# gr_edw_low.Draw("sameC")
	# gr_cdms.Draw("sameC")

	leg =TLegend(0.564,0.584,0.83,0.857)
	# leg.AddEntry(gr_edw_low.GetName(),"EDW II low mass ","e")
	# leg.AddEntry(gr_cdms.GetName(),"SCDMS LT (577 kg.d) ","e")
	# leg.AddEntry(gr_68_ERA.GetName(),"","f")

	leg.SetFillColor(kWhite)
	leg.SetLineColor(kWhite)
	leg.Draw()

	raw_input()

	cc.Print("./Figures/expected_limit.png")

analysis_type = "standard_resolution_no_ion_cut"
exposure = 1

plot_limit(analysis_type, exposure)




