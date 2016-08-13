from ROOT import *
import script_utils as script_utils
import sys,os, math
import numpy as np

class exponential:
   def __call__( self, x, par ):
		return par[1]*math.exp(float(x[0])*float(par[0]))

def get_heat_spec_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts):


	file_tree   = TFile(data_dir+bolo_name+"_fond.root")
	tree        = file_tree.Get(tree_name)
	
	cut_line='&&'.join(list_cuts)

	if tree.GetEntries(cut_line) != 0:

		h= TH1F(h_name + "_" + bolo_name, h_name + "_" + bolo_name, bin_X, min_X, max_X)
		tree.Project(h_name + "_" + bolo_name,channel_X, cut_line)
		
		h.SetStats(0)
		
		h.GetXaxis().SetTitle(channel_X_title)
		h.GetXaxis().CenterTitle(kTRUE)
		h.GetXaxis().SetTitleSize(0.06)
		h.GetXaxis().SetTitleOffset(0.8)
		
		h.GetYaxis().SetTitle(channel_Y_title)
		h.GetYaxis().CenterTitle(kTRUE)
		h.GetYaxis().SetTitleSize(0.06)
		h.GetYaxis().SetTitleOffset(0.8)
		# h.GetXaxis().SetTimeDisplay(1)
		
		
		#Do the plots
		
		fexp = TF1("fexp", exponential(), 0., 10., 2)
		fexp.SetParameter(0,-3.)
		fexp.SetParLimits(0,-20.,0)
		
		fexp.SetParameter(1,1.)
		fexp.SetParLimits(1, 0 ,10000000)
		


		# h.Draw()
		# h.Fit("fexp", "")
		# raw_input()
		
		return fexp.GetParameter(0)

	return 0


def launch_scaled_plot(bolo_name, data_dir, tree_name = "data"):

    
	# Load start and end times and duration
	path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
	file_name= bolo_name + "_polar_start_and_end_time.txt"  
	file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")
	
	#Load estimators
	estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
	estimator_file_name =bolo_name + "_estimators.txt" 
	d_estimator         ={}
	assert(os.path.isfile(estimator_path_name + estimator_file_name) )
	with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
		list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
		for line in list_estimator_lines:
			d_estimator[line[0]] = line[1]

	#Load standard cuts
	TCut_path_name = script_utils.create_directory('../Cut_files/')  
	TCut_file_name ="TCuts.txt" 
	file_TCut      ="" 
	#Add an exception if the file does not exist
	try:
		file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
	except IOError:
		script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
		sys.exit()

	# Load the cut values. 
	list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
	standard_cuts        =""
	# Add a boolean flag to check if the bolo has its cuts in the file
	is_bolo_in_file      =False
	for line in list_file_TCut_lines:
		if bolo_name == line[0]:
			standard_cuts = line[1]
			is_bolo_in_file = True
	assert(is_bolo_in_file)

	file_tree   = TFile(data_dir+bolo_name+"_fond.root")
	tree        = file_tree.Get(tree_name)
	
	list_cuts =[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_estimator["HEAT"] + ">2"]+ ["abs(EC1-EC2)<2"] 
	cut_line='&&'.join(list_cuts)

	h1= TH1F("h1", "h1", 40, 2, 10)
	h2= TH1F("h2", "h2", 40, 2, 10)
	h3= TH1F("h3", "h3", 40, 2, 10)
	tree.Project("h1",d_estimator["HEAT"], cut_line + "&& 1E6*UT1+UT2>1408.8E6 && 1E6*UT1+UT2<1409.8E6")
	tree.Project("h2",d_estimator["HEAT"], cut_line + "&&1E6*UT1+UT2>1410.3E6 && 1E6*UT1+UT2<1410.8E6")
	tree.Project("h3",d_estimator["HEAT"], cut_line + "&& 1E6*UT1+UT2>1411E6 && 1E6*UT1+UT2<1412E6")

	list_hist=[h1, h2, h3]

	cc = TCanvas("cc", "cc")
	gPad.SetLogy()
	h = TH1F("heat_rate", "heat_rate", 10, 2,10)
	h.SetStats(0)
	h.SetMaximum(1.3*min([elem.GetMaximum() for elem in list_hist]))
	
	h.GetXaxis().SetTitle("Energy (keV)")
	h.GetXaxis().CenterTitle(kTRUE)
	h.GetXaxis().SetTitleSize(0.06)
	h.GetXaxis().SetTitleOffset(0.8)
	
	h.GetYaxis().SetTitle("Rate (hists rescaled to same Integral)")
	h.GetYaxis().CenterTitle(kTRUE)
	h.GetYaxis().SetTitleSize(0.06)
	h.GetYaxis().SetTitleOffset(0.8)
	# h.GetXaxis().SetTimeDisplay(1)
	
	h.Draw()
	scale = min(list_hist[0].Integral(), list_hist[1].Integral(), list_hist[2].Integral())
	for histo in list_hist: histo.Scale(float(scale)/float(histo.Integral()))

	list_hist[0].SetLineColor(kRed)
	list_hist[1].SetLineColor(kBlue)
	list_hist[2].SetLineColor(kBlack)

	for elem in list_hist: elem.Draw("sameE1")

	leg = TLegend(0.43, 0.72, 0.88, 0.96)
	leg.AddEntry("h1","Calm period " ,"f")
	leg.AddEntry("h2","After 4th Sept period, many events" ,"f")
	leg.AddEntry("h3","Longer time after 4th Sept" ,"f")
	leg.SetFillColor(kWhite)
	leg.SetBorderSize(0)
	leg.Draw()

	raw_input()
	cc.Print('../Analyse_' + bolo_name + '/Figures/' + bolo_name + "_heatonly_spectrum_error.eps")



def launch_fit(bolo_name, data_dir, tree_name = "data"):

    
	# Load start and end times and duration
	path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
	file_name= bolo_name + "_polar_start_and_end_time.txt"  
	file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")
	
	#Load estimators
	estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
	estimator_file_name =bolo_name + "_estimators.txt" 
	d_estimator         ={}
	assert(os.path.isfile(estimator_path_name + estimator_file_name) )
	with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
		list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
		for line in list_estimator_lines:
			d_estimator[line[0]] = line[1]

	#Load standard cuts
	TCut_path_name = script_utils.create_directory('../Cut_files/')  
	TCut_file_name ="TCuts.txt" 
	file_TCut      ="" 
	#Add an exception if the file does not exist
	try:
		file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
	except IOError:
		script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
		sys.exit()

	# Load the cut values. 
	list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
	standard_cuts        =""
	# Add a boolean flag to check if the bolo has its cuts in the file
	is_bolo_in_file      =False
	for line in list_file_TCut_lines:
		if bolo_name == line[0]:
			standard_cuts = line[1]
			is_bolo_in_file = True
	assert(is_bolo_in_file)

    #Get the plots for given bolo
	fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "heat_spec", "heat_spec", "", "Rate","EC2", 100,  2,10
	list_time=[[1408.8E6,1409.8E6],[1410.3E6, 1410.8E6],[1411E6, 1412E6]]
	list_time=[[1409.2E6,1409.75E6],[1410.3E6, 1410.8E6],[1411.2E6, 1411.6E6]]

	with open('../Analyse_' + bolo_name + "/Text_files/" + bolo_name + "_heat_slope_distribution.txt", "w") as heat_slope_file:
		for i in range(len(list_time)):
			tinf      = str(list_time[i][0])
			tsup      = str(list_time[i][1])
			list_cuts =[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_estimator["HEAT"] + ">2"]+ ["abs(EC1-EC2)<2"] + ["1E6*UT1+UT2>" + tinf + " && 1E6*UT1+UT2<" + tsup]
			fit_val  = get_heat_spec_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)
			heat_slope_file.write(str(tinf) + "," + str(fit_val) + "\n")
			del hist



# 			heat_slope_file.write(str(tinf) + "," + str(fit_val) + "\n")


def plot_fitted_results(bolo_name):

	array = np.loadtxt('../Analyse_' + bolo_name + "/Text_files/" + bolo_name + "_heat_slope_distribution.txt", delimiter=",")
	list_time  = array[:,0].astype("float")
	list_slope = array[:,1].astype("float")
	
	h = TH1F("Slope vs time", "Slope vs time", 10, min(list_time), max(list_time))
	h.SetStats(0)
	h.SetMinimum(-6)
	
	h.GetXaxis().SetTitle("Time")
	h.GetXaxis().CenterTitle(kTRUE)
	h.GetXaxis().SetTitleSize(0.06)
	h.GetXaxis().SetTitleOffset(0.8)
	
	h.GetYaxis().SetTitle("Slope of heat spectrum in time bin")
	h.GetYaxis().CenterTitle(kTRUE)
	h.GetYaxis().SetTitleSize(0.06)
	h.GetYaxis().SetTitleOffset(0.8)
	h.GetXaxis().SetTimeDisplay(1)

	h.Draw()
	gr = TGraph(len(array[:,0]), list_time, list_slope)
	gr.Draw("sameL*")
	raw_input()

# def launch_fit_manybins(bolo_name, data_dir, tree_name = "data"):

    
#     # Load start and end times and duration
#     path_name= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')  
#     file_name= bolo_name + "_polar_start_and_end_time.txt"  
#     file_start_and_end = script_utils.open_text_file(path_name, file_name , "r")

#     #Create lists to hold the contents of the .txt file + the list for each polar duration
#     list_tmin, list_tmax, list_string_polar, list_duration = [], [], [], []

#     #Read the file lines and fill the lists
#     list_lines_start_and_end = [elem.rstrip().split(",") for elem in file_start_and_end.readlines()]
    
#     for k in range(len(list_lines_start_and_end)):
#         list_string_polar.append(list_lines_start_and_end[k][0])
#         list_tmin.append(float(list_lines_start_and_end[k][1]))
#         list_tmax.append(float(list_lines_start_and_end[k][2]))
#         list_duration.append(float(list_lines_start_and_end[k][3]))
# 	file_start_and_end.close()


# 	#Load estimators
# 	estimator_path_name = script_utils.create_directory('../Analyse_' + bolo_name + "/Text_files/")  
# 	estimator_file_name =bolo_name + "_estimators.txt" 
# 	d_estimator         ={}
# 	assert(os.path.isfile(estimator_path_name + estimator_file_name) )
# 	with open(estimator_path_name + estimator_file_name, 'r') as estimator_file: 
# 		list_estimator_lines= [elem.rstrip().split(",") for elem in estimator_file.readlines()]
# 		for line in list_estimator_lines:
# 			d_estimator[line[0]] = line[1]

# 	#Load standard cuts
# 	TCut_path_name = script_utils.create_directory('../Cut_files/')  
# 	TCut_file_name ="TCuts.txt" 
# 	file_TCut      ="" 
# 	#Add an exception if the file does not exist
# 	try:
# 		file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
# 	except IOError:
# 		script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
# 		sys.exit()

# 	# Load the cut values. 
# 	list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
# 	standard_cuts        =""
# 	# Add a boolean flag to check if the bolo has its cuts in the file
# 	is_bolo_in_file      =False
# 	for line in list_file_TCut_lines:
# 		if bolo_name == line[0]:
# 			standard_cuts = line[1]
# 			is_bolo_in_file = True
# 	assert(is_bolo_in_file)

#     #Get the plots for given bolo
# 	fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X = "heat_spec", "heat_spec", "", "Rate",d_estimator["HEAT"], 30,  1,10
# 	tmin = list_tmin[0]
# 	tmax = list_tmax[0]
# 	bin_time = 200.
# 	bin_width= float(tmax-tmin)/bin_time
# 	list_time=[[tmin, tmin + bin_width]]
# 	for i in range(1,int(bin_time)+1):
# 		list_time.append([list_time[i-1][1], list_time[i-1][1]+bin_width])

#     with open('../Analyse_' + bolo_name + "/Text_files/" + bolo_name + "_heat_slope_distribution.txt", "w") as heat_slope_file:
# 		for i in range(len(list_time)):
# 			tinf= str(list_time[i][0])
# 			tsup= str(list_time[i][1])
# 			list_cuts =[standard_cuts] + ["EIA<1 && EIB<1 && EIC<1 && EID<1"] + [d_estimator["HEAT"] + ">1"]+ ["abs(EC1-EC2)<2"] + ["1E6*UT1+UT2>" + tinf + " && 1E6*UT1+UT2<" + tsup]
# 			fit_val = get_heat_spec_time(bolo_name, data_dir, tree_name, fig_title, h_name, channel_X_title, channel_Y_title, channel_X, bin_X, min_X, max_X, list_cuts)

bolo_list   =["FID827", "FID837", "FID838", "FID828"]
# bolo_list = ["FID846"]
data_dir    = "../Fond_ERA_merged/"
for bolo_name in bolo_list:
	launch_scaled_plot( bolo_name, data_dir)
	# launch_fit( bolo_name, data_dir)
	# plot_fitted_results(bolo_name)
	pass



