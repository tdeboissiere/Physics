import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from matplotlib import gridspec


def check_pulse_tree(run_name):

	"""Check the simu pulse tree has been correctly built
	
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

	f = TFile("./ROOT_files/" + bolo_name + "_" + run_name + "_simulated_pulse_tree.root", "read")
	t = f.Get("tout")

	print t.GetEntries()

	for i in range(t.GetEntries()):

		print i

		t.GetEntry(i)

		pulse_A = np.zeros(1024)
		pulse_B = np.zeros(1024)

		for i in range(1024):
			pulse_A[i]=t.trace_A[i]
			pulse_B[i]=t.trace_B[i]

		fig = plt.figure(figsize=(8, 6)) 
		gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
		ax0 = plt.subplot(gs[0])
		ax0.plot(pulse_A)
		ax1 = plt.subplot(gs[1])
		ax1.plot(pulse_B)
		plt.show()
		raw_input()
		plt.close()


run_name = "222500"
check_pulse_tree(run_name)
