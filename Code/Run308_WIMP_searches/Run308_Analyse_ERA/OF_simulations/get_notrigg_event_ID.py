from ROOT import *
from array import array
import PyROOTPlots as PyRPl
import numpy as np
import sys
import conversion_ANA_SAMBA as conv
import os
import script_utils as script_utils
ROOT.gSystem.Load("/home/irfulx204/mnt/tmain/Desktop/New_ERA/ERA/lib/EraLib.so")
from ROOT import EdwPulse,EdwEvent,FitPulse,EdwTemplate,NoiseSpectrum

def get_notrigg_event_ID(bolo_name):
    
	"""Create run by run files with simulated low energy pulses
	
	Detail:
		Detailed description

	Args:
		bolo_name = (str) bolometer name
		
	Returns:
		void

	Raises:
		void
	"""
	
	bolo_name = "FID837"
	
	#Open tree with no-trigg events
	event_path = "../Event_files/" + bolo_name + "/"
	assert(event_path + bolo_name + "_all_notrigged_event_tree.root")
	t,f = PyRPl.open_ROOT_object(event_path + bolo_name + "_all_notrigged_event_tree.root", "out")
	nEntries = t.GetEntries()
	
	#Get the event of the event tree
	Evt=EdwEvent()
	t.SetBranchAddress("Event",Evt)
	
	print nEntries

	with open("./Text_files/" + bolo_name + "/" + bolo_name + "_notrigged_event_ID.txt", "w") as fID:		
		for i in range(nEntries):

			#Progress bar
			sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
			sys.stdout.flush()

			t.GetEntry(i)
			Run = Evt.Run()
			SambaNum = Evt.SambaNum()
			fID.write(str(Run) + "," + str(SambaNum) + "\n")


bolo_name = "FID837"
get_notrigg_event_ID(bolo_name)