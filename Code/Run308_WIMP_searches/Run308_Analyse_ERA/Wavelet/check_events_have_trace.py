from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import conversion_ANA_SAMBA as conv
import sys

def check_events_have_trace(bolo_name):

	"""Check events have a trace in the trace tree
	
	Detail:
		Events come from the all events file (events for low mass analysis)

	Args:
		bolo_name = (type) bolo name
		
	Returns:
		void

	Raises:
		void
	"""

	data_types = {"names": ("RUN", "SN"), "formats": ("i", "i")}
	arr_events = np.loadtxt("./Populations/"+ bolo_name + "/" + bolo_name + "_all_events.txt", delimiter=",",  dtype=data_types)

	t, f = PyRPl.open_ROOT_object("./ROOT_files/" + bolo_name + "/" + bolo_name + "_all_trigged_trace_tree.root", "out")
	nEntries = t.GetEntries()

	d_events = {}

	# for i in arr_events.shape[0]:
	# 	d_events[str(arr_events["RUN"][i]) + "_" + str(arr_events["SN"][i])] = True
		
	for i in range(nEntries):

		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		t.GetEntry(i)
		Run, SN = int(conv.samba2ana(t.Run[:-1])), int(t.SambaNum)
		d_events[str(Run) + "_" + str(SN)] = True

	count = 0
	for i in range(arr_events.shape[0]):
		try :
			d_events[str(arr_events["RUN"][i]) + "_" + str(arr_events["SN"][i])]
		except KeyError:
			count+=1
	print
	print count



bolo_name = "FID837"
check_events_have_trace(bolo_name)