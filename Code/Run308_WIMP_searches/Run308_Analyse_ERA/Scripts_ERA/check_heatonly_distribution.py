from ROOT import *
import PyROOTPlots as PyRPl
import script_utils as script_utils
import Analysis_utilities as Ana_ut

def check_heatonly_heat_distribution(bolo_name="FID837"):

	"""Check the simulation of the heatonly distribution
	
	Detail:
		Plot the 2D correlation plots and the 1D spectrum

	Args:
		bolo_name = (str) bolometer name
		
	Returns:
		void

	Raises:
		void
	"""

	ttrue,ftrue = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/FID837_fond.root", "data")
	tsimu,fsimu = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_min1_min4_5/Heatonly/ROOT_files/FID837_heatonly_tree.root", "t_new0")

	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	# 1D plots
	list_names1D      = ["EC1", "EC2", "Estimator"]
	list_tree_chann1D = ["EC1", "EC2", "0.5*(EC1+EC2)"]
	list_leg1D        = [TLegend(0.57, 0.57,0.9,0.9,"","brNDC") for i in range(len(list_names1D))]
	list_htrue1D    = [TH1F("htrue2D" + name, "htrue2D" + name, 100, 0, 5) for name in list_names1D]
	list_hsimu1D    = [TH1F("hsimu2D" + name, "hsimu2D" + name, 100, 0, 5) for name in list_names1D]
	for index, chan in enumerate(list_tree_chann1D):
		ttrue.Project(list_htrue1D[index].GetName(), chan, "0.5*(EIB+EID)<0 &&" + standard_cuts)
		tsimu.Project(list_hsimu1D[index].GetName(), chan, "0.5*(EIB+EID)<0")
		PyRPl.process_TH1(list_htrue1D[index], X_title = chan, color=kRed)
		PyRPl.process_TH1(list_hsimu1D[index], X_title = chan, color=kBlack)
		list_htrue1D[index].Scale(1./list_htrue1D[index].Integral())
		list_hsimu1D[index].Scale(1./list_hsimu1D[index].Integral())


	# 2D plots
	list_names1D      = ["EC1vsEC2"]
	list_tree_chann2D = ["EC1:EC2"]
	list_leg2D       = [TLegend(0.19, 0.67,0.44,0.9,"","brNDC") for i in range(len(list_names1D))]
	list_htrue2D    = [TH2F("htrue2D" + name, "htrue2D" + name, 1000, -2, 5, 1000, -2, 5) for name in list_names1D]
	list_hsimu2D    = [TH2F("hsimu2D" + name, "hsimu2D" + name, 1000, -2, 5, 1000, -2, 5) for name in list_names1D]
	for index, chan in enumerate(list_tree_chann2D):
		ttrue.Project(list_htrue2D[index].GetName(), chan, "0.5*(EIB+EID)<0 &&" + standard_cuts)
		tsimu.Project(list_hsimu2D[index].GetName(), chan, "0.5*(EIB+EID)<0")
		PyRPl.process_TH2(list_htrue2D[index], X_title = chan[:3],  Y_title = chan[4:], color=kRed)
		PyRPl.process_TH2(list_hsimu2D[index], X_title = chan[:3],  Y_title = chan[4:], color=kBlack)



	cc =TCanvas("cc", "cc")
	cc.Divide(2,2)
	for i in range(3):
		cc.cd(i+1)
		list_hsimu1D[i].SetMaximum(1.2*max(list_hsimu1D[i].GetMaximum(), list_htrue1D[i].GetMaximum()))
		list_hsimu1D[i].Draw()
		list_htrue1D[i].Draw("same")

		list_leg1D[i].AddEntry(list_hsimu1D[i].GetName(), "Simu" , "l")
		list_leg1D[i].AddEntry(list_htrue1D[i].GetName(), "True data" , "l")
		list_leg1D[i].Draw("same")
	cc.cd(4)
	list_hsimu2D[0].Draw()
	list_htrue2D[0].Draw("same")

	list_leg2D[0].AddEntry(list_hsimu2D[0].GetName(), "Simu" , "p")
	list_leg2D[0].AddEntry(list_htrue2D[0].GetName(), "True data" , "p")
	list_leg2D[0].Draw("same")
	
	raw_input()
	fig_dir = script_utils.create_directory("./Figures/")
	cc.Print(fig_dir + bolo_name + "_heatonly_heat_correlation.eps")


def check_heatonly_ion_distribution(bolo_name="FID837"):

	"""Check the simulation of the heatonly distribution
	
	Detail:
		Plot the 2D correlation plots and the 1D spectrum

	Args:
		bolo_name = (str) bolometer name
		
	Returns:
		void

	Raises:
		void
	"""

	ttrue,ftrue = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/FID837_ana_0.5_0_5_fond.root", "data")
	tsimu,fsimu = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_0.5_0_5/Heatonly/ROOT_files/FID837_heatonly_tree.root", "t_new0")


	#1D plot
	channel1D = "EIB"

	htrue1D = TH1F("htrue1D", "htrue1D", 300, -5, 5)
	hsimu1D = TH1F("hsimu1D", "hsimu1D", 300, -5, 5)

	ttrue.Project("htrue1D", channel1D, "0.5*(EIB+EID)<0.7")
	tsimu.Project("hsimu1D", channel1D, "0.5*(EIB+EID)<0.7")

	PyRPl.process_TH1(htrue1D, X_title = channel1D, color=kRed)
	PyRPl.process_TH1(hsimu1D, X_title = channel1D, color=kBlack)

	htrue1D.Scale(1./htrue1D.Integral())
	hsimu1D.Scale(1./hsimu1D.Integral())

	c1D =TCanvas("c1D", "c1D")
	htrue1D.Draw()
	hsimu1D.Draw("same")

	# 2D plots
	list_names      = ["AB", "AC", "AD", "BC", "BD", "CD"]
	list_tree_chann = ["EIA:EIB", "EIA:EIC", "EIA:EID", "EIB:EIC", "EIB:EID", "EIC:EID"]
	list_leg        = [TLegend(0.57, 0.57,0.9,0.9,"","brNDC") for i in range(6)]
	list_htrue2D    = [TH2F("htrue2D" + name, "htrue2D" + name, 1000, -2, 2, 1000, -2, 2) for name in list_names]
	list_hsimu2D    = [TH2F("hsimu2D" + name, "hsimu2D" + name, 1000, -2, 2, 1000, -2, 2) for name in list_names]
	for index, chan in enumerate(list_tree_chann):
		ttrue.Project(list_htrue2D[index].GetName(), chan, "0.5*(EIB+EID)<0.7")
		tsimu.Project(list_hsimu2D[index].GetName(), chan, "0.5*(EIB+EID)<0.7")
		PyRPl.process_TH2(list_htrue2D[index], X_title = chan[:3],  Y_title = chan[4:], color=kRed)
		PyRPl.process_TH2(list_hsimu2D[index], X_title = chan[:3],  Y_title = chan[4:], color=kBlack)


	l = TLine(-2,0,2,0)
	l2 = TLine(0,-2,0,2)
	l.SetLineWidth(3)
	l2.SetLineWidth(3)

	c2D =TCanvas("c2D", "c2D")
	c2D.Divide(3,2)
	for i in range(6):
		c2D.cd(i+1)
		list_hsimu2D[i].Draw()
		l.Draw("same")
		l2.Draw("same")
		list_htrue2D[i].Draw("same")

		# list_leg[i].AddEntry(list_hsimu2D[i].GetName(), "Simu" , "p")
		# list_leg[i].AddEntry(list_htrue2D[i].GetName(), "True data" , "p")
		# list_leg[i].Draw("same")

	raw_input()
	fig_dir = script_utils.create_directory("./Figures/")
	c2D.Print(fig_dir + bolo_name + "_heatonly_ionisation_correlation.eps")

def check_heatonly_ion_distribution_noioncut(bolo_name="FID837"):

	"""Check the simulation of the heatonly distribution
	
	Detail:
		Plot the 2D correlation plots and the 1D spectrum

	Args:
		bolo_name = (str) bolometer name
		
	Returns:
		void

	Raises:
		void
	"""

	ttrue,ftrue = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/FID837_fond.root", "data")
	tsimu,fsimu = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/ana_min1_min4_5/Heatonly/ROOT_files/FID837_heatonly_tree.root", "t_new0")

	#Load standard cuts
	standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")

	# 2D plots
	list_names      = ["AB", "AC", "AD", "BC", "BD", "CD"]
	list_tree_chann = ["EIA:EIB", "EIA:EIC", "EIA:EID", "EIB:EIC", "EIB:EID", "EIC:EID"]
	list_leg        = [TLegend(0.57, 0.57,0.9,0.9,"","brNDC") for i in range(6)]
	list_htrue2D    = [TH2F("htrue2D" + name, "htrue2D" + name, 300, -2, 2, 300, -2, 2) for name in list_names]
	list_hsimu2D    = [TH2F("hsimu2D" + name, "hsimu2D" + name, 300, -2, 2, 300, -2, 2) for name in list_names]
	for index, chan in enumerate(list_tree_chann):
		ttrue.Project(list_htrue2D[index].GetName(), chan, "0.5*(EC1+EC2)<0.5 &&" + standard_cuts)
		tsimu.Project(list_hsimu2D[index].GetName(), chan, "0.5*(EC1+EC2)<0.5")
		PyRPl.process_TH2(list_htrue2D[index], X_title = chan[:3],  Y_title = chan[4:], color=kRed)
		PyRPl.process_TH2(list_hsimu2D[index], X_title = chan[:3],  Y_title = chan[4:], color=kBlack)

	print ttrue.GetEntries("0.5*(EC1+EC2)<0.5 &&" + standard_cuts)
	print tsimu.GetEntries("0.5*(EC1+EC2)<0.5 ")

	l = TLine(-2,0,2,0)
	l2 = TLine(0,-2,0,2)
	l.SetLineWidth(3)
	l2.SetLineWidth(3)

	c2Dtrue =TCanvas("c2Dtrue", "c2Dtrue")
	c2Dtrue.Divide(3,2)
	for i in range(6):
		c2Dtrue.cd(i+1)
		list_htrue2D[i].Draw("colz")
		l.Draw("same")
		l2.Draw("same")

		list_leg[i].AddEntry(list_htrue2D[i].GetName(), "True data" , "p")
		# list_leg[i].Draw("same")

	c2Dsimu =TCanvas("c2Dsimu", "c2Dsimu")
	c2Dsimu.Divide(3,2)
	for i in range(6):
		c2Dsimu.cd(i+1)
		list_hsimu2D[i].Draw("")
		list_htrue2D[i].Draw("same")
		l.Draw("same")
		l2.Draw("same")

		list_leg[i].AddEntry(list_htrue2D[i].GetName(), "Simu data" , "p")
		# list_leg[i].Draw("same")

	raw_input()
	fig_dir = script_utils.create_directory("./Figures/")
	c2Dsimu.Print(fig_dir + bolo_name + "_heatonly_ionisation_correlation_noioncut.png")

# check_heatonly_ion_distribution()
check_heatonly_ion_distribution_noioncut()
# check_heatonly_heat_distribution()