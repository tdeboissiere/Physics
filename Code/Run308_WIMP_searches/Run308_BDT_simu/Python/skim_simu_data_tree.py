from ROOT import *
import BDT_file_handler as BDT_fh
import script_utils as script_utils
from array import array
import Analysis_utilities as Ana_ut

def skim_true_event_tree(bolo_name, analysis_type, nsimu):

	""" Skim simulated data tree so that surface events are rejected with a veto cut
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)
		nsimu = (int) number of simulated tree
	Returns:
		void

	Raises:
		AssertionError 
	"""

	#Load simulated data tree
	ori_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type + "/True_events/ROOT_files/"
	file_tree = TFile(ori_path + bolo_name + "_true_events_tree.root", "read")

	#Load FWHM
	d_std = BDT_fh.open_true_event_FWHM_file(bolo_name,)

 	file_skim = TFile(ori_path + bolo_name + "_true_events_nosurf_tree.root", "recreate")

	#Loop over trees
	for i in range(nsimu):
		tree = file_tree.Get("t_new" + str(i))

		#Event list
		l_event    = TEventList("l_event")
		tree.Draw(">>l_event","EIA<2*" + str(d_std["FWIA"]) + "&&EIC<2*" + str(d_std["FWIC"]) )
		pop_len = l_event.GetN()

		skim_tree = TTree("t_new" + str(i), "t_new" + str(i))

		EC1 =  array("f", [ 0 ] )
		EC2 =  array("f", [ 0 ] )
		EIA =  array("f", [ 0 ] )
		EIB =  array("f", [ 0 ] )
		EIC =  array("f", [ 0 ] )
		EID =  array("f", [ 0 ] )
		ENR =  array("f", [ 0 ] )

		skim_tree.Branch("EC1", EC1, "EC1/F")
		skim_tree.Branch("EC2", EC2, "EC2/F")
		skim_tree.Branch("EIA", EIA, "EIA/F")
		skim_tree.Branch("EIB", EIB, "EIB/F")
		skim_tree.Branch("EIC", EIC, "EIC/F")
		skim_tree.Branch("EID", EID, "EID/F")
		skim_tree.Branch("ENR", ENR, "ENR/F")

		for k in range(pop_len):
			counter = l_event.GetEntry(k)
			tree.GetEntry(counter)

			EC1[0] = tree.EC1
			EC2[0] = tree.EC2
			EIA[0] = tree.EIA
			EIB[0] = tree.EIB
			EIC[0] = tree.EIC
			EID[0] = tree.EID
			ENR[0] = 0

			skim_tree.Fill()

		skim_tree.Write()
		del l_event
		del skim_tree

	file_skim.Close()


bolo_name = "FID837"
#convention ana_u_v_w_x :  cut @ u keV Heat, v sigma heat only ion band width, w sigma veto
# analysis_type = "ana_1.5_min2_50"
analysis_type = "ana_1.5_0_5"
nsimu = 10
skim_true_event_tree(bolo_name, analysis_type, nsimu)