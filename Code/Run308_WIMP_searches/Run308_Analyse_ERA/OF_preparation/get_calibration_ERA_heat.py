from ROOT import *
import numpy as np
import sys,os
import new_conversion_ANA_SAMBA as conv
import script_utils as script_utils
import matplotlib.pylab as plt
import BDT_file_handler as BDT_fh
import PyROOTPlots as PyRPl

def get_calib_pop(bolo_name):

	"""Get calib data for given bolometer
	
	Detail:
		void

	Args:
		bolo_name = (str) bolo name 

	Returns:
		AssertionError
		IOError

	Raises:
		void
	"""
	
	file_tree   = TFile("../Fond_merged/"+bolo_name+"_fond.root")
	tree        = file_tree.Get("data")
	
	#Create background hist directory
	path_txt_files = script_utils.create_directory("./Text_files/" + bolo_name + "/")
	
	#Load standard cuts
	TCut_path_name = script_utils.create_directory('../Cut_files/')  
	TCut_file_name ="TCuts_ANA.txt" 
	file_TCut      ="" 
	#Add an exception if the file does not exist
	try:
		file_TCut = script_utils.open_text_file(TCut_path_name, TCut_file_name , "r")
	except IOError:
		script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
		sys.exit()

	# Load the cut values. 
	list_file_TCut_lines =[line.rstrip().split(",") for line in file_TCut.readlines()]
	standard_cuts =""
	# Add a boolean flag to check if the bolo has its cuts in the file
	is_bolo_in_file =False
	for line in list_file_TCut_lines:
		if bolo_name == line[0]:
			standard_cuts = line[1]
			is_bolo_in_file = True
	assert(is_bolo_in_file)
	
	#All events up to rather low energies
	l_calib_pop = TEventList("l_calib_pop")	
	tree.Draw(">>l_calib_pop",standard_cuts +  "&& 0.5*(EIA+EIB+EIC+EID) >4 && abs(EC1-EC2)<2 && 0.5*(EC1+EC2)>4 && 0.5*(EC1+EC2)<20")
	pop_len = l_calib_pop.GetN()
	pop_file_name = bolo_name + "_calib_pop.txt"
	pop_file = script_utils.open_text_file(path_txt_files, pop_file_name, 'w')
	for k in range(pop_len):
		counter = l_calib_pop.GetEntry(k)
		tree.GetEntry(counter)
		pop_file.write(str(int(tree.RUN)) + "," + str(int(tree.SN)) + "," + str(tree.EC1) + "," + str(tree.EC2) + "," + str(tree.EIA) + "," + str(tree.EIB) + "," + str(tree.EIC) + "," + str(tree.EID) + "\n")
	pop_file.close()

	del l_calib_pop




def get_calibration_coefficient(bolo_name):

	"""Calibrate bolometer
	
	Detail:
		Use calib data to calibrate

	Args:
		bolo_name = (str) bolo name 

		
	Returns:
		void

	Raises:
		void
	"""

	# pre_run_list = glob.glob("../Amp_files/" + bolo_name + "*Chal1*")
	# samba_run_list = [elem[30:38] for elem in pre_run_list]
	# ana_run_list = sorted([conv.samba2ana(elem) for elem in samba_run_list])
	# samba_run_list = [conv.ana2samba(elem) for elem in ana_run]
	
	data_types = {'names': ('RUN', 'SN', 'EC1', 'EC2', 'EIA', 'EIB', 'EIC', 'EID'), 'formats': ('i', 'i', 'f', 'f', 'f', 'f', 'f', 'f')}
	arr_jules = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name + "_calib_pop.txt", delimiter=',',  dtype=data_types)

	d_calib={}
	for index, elem in enumerate(arr_jules["RUN"]):
		d_calib[str(arr_jules["RUN"][index])+"_" +str(arr_jules["SN"][index])]=[arr_jules["EC1"][index], arr_jules["EC2"][index]]
	

	fChal1    = TFile("../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_Chal1.root", "read")
	fChal2    = TFile("../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_Chal2.root", "read")
	fbasic    = TFile("../Amp_files_merged/" + bolo_name + "/" + bolo_name + "_wiener_basic.root", "read")
	
	tChal1    = fChal1.Get("wienerntp_" + bolo_name + "_Chal1")
	tChal2    = fChal2.Get("wienerntp_" + bolo_name + "_Chal2")
	tbasic    = fbasic.Get("basicntp_" + bolo_name )

	tChal1.AddFriend(tChal2)
	tChal1.AddFriend(tbasic)

	list_EC1, list_EC2, list_Chal1, list_Chal2 = [], [], [], []

	nEntries = tChal1.GetEntries()

	for i in range(nEntries):

		sys.stdout.write('\r' + str(i)+' / '+str(-1+nEntries))
		sys.stdout.flush()
		
		tChal1.GetEntry(i)
		RunNum, SambaNum = str(conv.samba2ana(tbasic.Run[:-1])) , str(tbasic.SambaNum)
		try:
			EC1, EC2 = d_calib[RunNum+ "_" + SambaNum][0], d_calib[RunNum+ "_" + SambaNum][1]
			Chal1, Chal2 = tChal1.WienerAmpl, tChal2.WienerAmpl
			list_EC1.append(EC1)
			list_EC2.append(EC2)
			list_Chal1.append(Chal1)
			list_Chal2.append(Chal2)

		except KeyError:
			pass

	arr_EC1=np.array(list_EC1).astype(float)
	arr_EC2=np.array(list_EC2).astype(float)
	arr_Chal1=np.array(list_Chal1).astype(float)
	arr_Chal2=np.array(list_Chal2).astype(float)


	# #remove outliers far from the calib equivalent in ADU
	# med_Chal1 = np.median(arr_Chal1)
	# med_Chal2 = np.median(arr_Chal2)

	# lChal1 = np.where(np.abs(arr_Chal1 -med_Chal1) < np.sqrt(med_Chal1))[0]
	# lChal2 = np.where(np.abs(arr_Chal2 -med_Chal2) < np.sqrt(med_Chal2))[0]

	# arr_EC1 = arr_EC1[lChal1]
	# arr_EC2 = arr_EC2[lChal2]
	# arr_Chal1 = arr_Chal1[lChal1]
	# arr_Chal2 = arr_Chal2[lChal2]

	p1 = np.polyfit(arr_EC1, arr_Chal1, 1)
	p2 = np.polyfit(arr_EC2, arr_Chal2, 1)

	pol1 = np.poly1d(p1)
	pol2 = np.poly1d(p2)

	print p1
	print p2

	#Create figure path
	figure_path = script_utils.create_directory("./Figures/" + bolo_name)

	plt.plot(arr_EC1, arr_Chal1, 'o', arr_EC1, pol1(arr_EC1), '-')
	plt.show()
	plt.savefig(figure_path + bolo_name + "_Chal1_calibration.png")
	raw_input()

	plt.plot(arr_EC2, arr_Chal2, 'o', arr_EC2, pol2(arr_EC2), '-')
	plt.show()
	plt.savefig(figure_path + bolo_name + "_Chal2_calibration.png")
	raw_input()

	coeff1 = p1[0]
	coeff2 = p2[0]

	raw_input("Check fit...")

	#Write output bolo calib file
	with open("./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff.txt", "w") as fc:
		fc.write("bolo_name,calibcoeff1,calibcoeff2\n")
		fc.write(bolo_name+","+str(coeff1)+","+str(coeff2))


def refine_calibration_coefficient(bolo_name):

	"""Re Calibrate bolometer so that 10.37 keV peak IS as 10.37
	
	Detail:
		Use gamma to calibrate at 10 keV

	Args:
		bolo_name = (str) bolo name 

		
	Returns:
		void

	Raises:
		void
	"""

	# cut = "FWIA<2&&FWIB<2&&FWIC<2&&FWID<2&&FWC1_ERA<2&&FWC2_ERA<3&&(CHIA-RCIA)<0.3&&(CHIB-RCIB)<0.3&&(CHIC-RCIC)<0.3&&(CHID-RCIC)<0.3&&CHIC1_ERA<2&&CHIC2_ERA<2&&KTH>0&&KTH<1"
	cut = "FWIA<2&&FWIB<2&&FWIC<2&&FWID<2&&FWC1_ERA<2&&FWC2_ERA<3&&CHIA<4&&CHIB<4&&CHIC<4&&CHID<4&&CHIC1_ERA<2&&CHIC2_ERA<2&&KTH>0&&KTH<1"

	tree, ftree = PyRPl.open_ROOT_object(Ana_path + "Fond_ERA_merged/" + bolo_name + "_fond.root", "t_merged")
	hgammaEC1 = TH1F("hgammaEC1", "hgammaEC1", 70, 0, 15)
	hgammaEC2 = TH1F("hgammaEC2", "hgammaEC2", 70, 0, 15)
	hgammaEIA = TH1F("hgammaEIA", "hgammaEIA", 100, 0, 15)
	hgammaEIB = TH1F("hgammaEIB", "hgammaEIB", 100, 0, 15)
	hgammaEIC = TH1F("hgammaEIC", "hgammaEIC", 100, 0, 15)
	hgammaEID = TH1F("hgammaEID", "hgammaEID", 100, 0, 15)
	string_EI = "0.5*(EIA+EIB+EIC+EID)"
	string_ER = "((1+8./3)*0.5*(EC1_ERA+EC2_ERA) - 0.333*(1.5*EIA+4*EIB+1.5*EIC+4*EID))"
	string_Q = string_EI+"/"+string_ER
	tree.Project("hgammaEC1", "EC1_ERA", "abs(EIB-EID)<2 && abs(EC1-EIB)<2 && (EIB)/EC1>0.7 && EIB>5")	
	tree.Project("hgammaEC2", "EC2_ERA", "abs(EIB-EID)<2 && abs(EC2-EIB)<2 && (EIB)/EC1>0.7 && EIB>5")
	tree.Project("hgammaEIA", "EIA", cut + "&&" + string_Q + ">0.7 && EIB>3 && EIA>3 && EID<2 && EIC<2")	
	tree.Project("hgammaEIB", "EIB",  cut + "&&" + string_Q + ">0.7 && EIB>3 && EID>3 && EIA<2 && EIC<2")
	tree.Project("hgammaEIC", "EIC",  cut + "&&" + string_Q + ">0.7 && EIB<2 && EIA<2 && EIC>3 && EID>3")
	tree.Project("hgammaEID", "EID",  cut + "&&" + string_Q + ">0.7 && EIB>3 && EID>3 && EIA<2 && EIC<2")

	class Gamma:
		def __call__( self, x, par ):
			peak_10_4 = par[0]*TMath.Gaus(x[0], par[1], par[2])
			return peak_10_4

	# Call FidGamma function to get derivative
	fFidGamma = TF1( "FidGamma", Gamma(), 0, 14., 3 )
	fFidGamma.SetNpx(500)
	fFidGamma.SetParName(0,"A_10.4")
	fFidGamma.SetParName(1,"mu_10.4")
	fFidGamma.SetParName(2,"s_10.4")

	fFidGamma.SetParameter(0,1)
	fFidGamma.SetParameter(1,10.4)
	fFidGamma.SetParameter(2,0.4)

	list_hist_ion = [hgammaEIA, hgammaEIB, hgammaEIC, hgammaEID]
	list_corr = []
	cc = TCanvas("cc", "cc")
	cc.Divide(2,2)
	for i in range(4):
		cc.cd(i+1)
		list_hist_ion[i].Draw()
		list_hist_ion[i].Fit("FidGamma","","",9.5,12)
		list_corr.append(fFidGamma.GetParameter(1)/10.37)

	raw_input()
	sys.exit()

	#Get the old calibration coefficients
	calib_file = "./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff.txt"
	calibcoeff1, calibcoeff2=1,1
	assert(os.path.isfile(calib_file))
	with open(calib_file, "r") as fc:
		lines = fc.readlines()
		calibcoeff1, calibcoeff2 = float(lines[1].rstrip().split(",")[1]), float(lines[1].rstrip().split(",")[2]) 

	hgammaEC1.Fit("FidGamma","","",9.5, 12)	
	coeff1 = (fFidGamma.GetParameter(1)/10.37)*calibcoeff1
	A1_ratio = 10.37/fFidGamma.GetParameter(1)
	raw_input()

	hgammaEC2.Fit("FidGamma","","",9.5, 12)
	coeff2 = (fFidGamma.GetParameter(1)/10.37)*calibcoeff2
	A2_ratio = 10.37/fFidGamma.GetParameter(1)
	raw_input()
	# sys.exit()


	#Write output bolo calib file
	with open("./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff_refined.txt", "w") as fc:
		fc.write("bolo_name,calibcoeff1,calibcoeff2,A1/A1_old,A2/A2_old\n")
		fc.write(bolo_name+","+str(coeff1)+","+str(coeff2)+","+str(A1_ratio)+","+str(A2_ratio))

	# Write output bolo calib file
	with open("./Text_files/" + bolo_name + "/" + bolo_name + "_ion_calibcorrection.txt", "w") as fc:
		fc.write("bolo_name,corr_IA,corr_IB,corr_IC,corr_ID\n")
		fc.write(bolo_name+","+str(list_corr[0])+","+str(list_corr[1])+","+str(list_corr[2])+","+str(list_corr[3]))