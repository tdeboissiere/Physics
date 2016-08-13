import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from matplotlib import gridspec
import conversion_ANA_SAMBA as conv

ROOT.gSystem.Load("/home/irfulx204/mnt/tmain/Desktop/New_ERA/ERA/lib/EraLib.so")
from ROOT import EdwPulse,EdwEvent,FitPulse,EdwTemplate,NoiseSpectrum

def check_event_tree(run_name):

	"""Check the event pulse tree has been correctly built
	
	Detail:
		Detailed description

	Args:
		run_name = (str) run_name (ANA convention)
		
	Returns:
		void

	Raises:
		void
	"""

	bolo_name = "FID837"

	plt.ion()

	f = TFile("./ROOT_files/TraceDir/" + bolo_name + "_" + str(conv.ana2samba(int(run_name))) + "_simulated_event_tree.root", "read")
	t = f.Get("EdwTree")

	#Get the event of the event tree
	Evt=EdwEvent()
	t.SetBranchAddress("Event",Evt)

	print t.GetEntries()

	for i in range(t.GetEntries()):

		print i

		t.GetEntry(i)

		pulse_A = np.array(Evt.Pulse("chalA " + bolo_name).Trace())
		pulse_B = np.array(Evt.Pulse("chalB " + bolo_name).Trace())

		fig = plt.figure(figsize=(8, 6)) 
		gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
		ax0 = plt.subplot(gs[0])
		ax0.plot(pulse_A)
		ax1 = plt.subplot(gs[1])
		ax1.plot(pulse_B)
		plt.show()
		raw_input()
		plt.close()


run_name = "222502"
check_event_tree(run_name)
