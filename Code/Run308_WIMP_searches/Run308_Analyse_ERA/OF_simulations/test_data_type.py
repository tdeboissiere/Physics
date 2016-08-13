from ROOT import *
from array import array
import PyROOTPlots as PyRPl
import numpy as np
import sys
import conversion_ANA_SAMBA as conv

ROOT.gSystem.Load("/home/irfulx204/mnt/tmain/Desktop/New_ERA/ERA/lib/EraLib.so")
from ROOT import EdwPulse,EdwEvent,FitPulse,EdwTemplate,NoiseSpectrum



def test_data_type():
    
	"""Check the type of the pulse trace
	
	Detail:
		Detailed description

	Args:

		
	Returns:
		void

	Raises:
		void
	"""
	
	bolo_name = "FID837"
	
	#Open tree with no-trigg traces
	t,f = PyRPl.open_ROOT_object("../Output/" + bolo_name + "_all_notrigg_event_tree.root", "out")
	nEntries = t.GetEntries()
	
	#Get the event of the event tree
	Evt=EdwEvent()
	t.SetBranchAddress("Event",Evt)
	
	list_runs_ana = np.loadtxt("./Text_files/list_runs_ana_s2.txt").astype(int)
	list_runs_ana = sorted(list(list_runs_ana))

	#Loop over runs to fill the corresponding tree
	for run_ana in list_runs_ana:
		for i in range(nEntries):

			sys.stdout.write('\r' + str(i)+' / '+str(-1+nEntries))
			sys.stdout.flush()
			
			t.GetEntry(i)
			Run = conv.samba2ana(Evt.Run())

			if int(Run) == int(run_ana):

				tr = np.ones(1024).astype(int)
				pp =Evt.Pulse("chalA " + bolo_name)
				tA = array( "h", 1024*[ 1 ] )
				tA[0] = 100
				pp.SetTrace(1024, tA)

				print
				print pp.Trace()[0]
				print
				raw_input()




test_data_type()