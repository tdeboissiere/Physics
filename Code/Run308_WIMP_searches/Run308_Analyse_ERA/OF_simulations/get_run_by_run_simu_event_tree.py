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

class Normalised_TemplateChalA_FID837:
    def __call__( self, x, par ):

        A1    = 0.691248887981/8.54005337642
        A2    =6.12435148311
        A3    = 53.8071234798
        t1    = 0.0932996328902
        t2    = 0.0085815845335
        t3    = 0.0314935265196
        t4    = 0.143222506399
        start = 512.742674442
      
        xnorm =x[0]-start
        pulse = par[0]*A1*(1-TMath.Exp(-xnorm*t1))*(TMath.Exp(-xnorm*t2)+A2*TMath.Exp(-xnorm*t3)+A3*TMath.Exp(-xnorm*t4))
        if x[0]<start :
            return 0
        else :
            return pulse

class Normalised_TemplateChalB_FID837:
    def __call__( self, x, par ):

        A1    = 0.990831812329/8.15133971563
        A2    = 4.15335723092
        A3    = 27.8726768628
        t1    = 0.105953931947
        t2    = 0.010041224142
        t3    = 0.033455593789
        t4    = 0.12666380985
        start = 512.675832237

        xnorm =x[0]-start
        pulse = par[0]*A1*(1-TMath.Exp(-xnorm*t1))*(TMath.Exp(-xnorm*t2)+A2*TMath.Exp(-xnorm*t3)+A3*TMath.Exp(-xnorm*t4))
        if x[0]<start :
            return 0
        else :
            return pulse

def create_event_simu_tree(d_channel_calib):
    
	"""Create run by run files with simulated low energy pulses
	
	Detail:
		Detailed description

	Args:
		d_channel_calib = (dict) calib coefficients
		
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
	
	data_types = {"names": ("RUN", "start", "end"), "formats": ("i", "i", "i")}
	file_path = "./Text_files/" + bolo_name + "/"
	file_start_end_name = file_path + bolo_name + "_list_runs_start_end_index.txt"
	assert(os.path.isfile(file_start_end_name))
	arr_start_end = np.loadtxt(file_start_end_name, delimiter=",",  dtype=data_types)
	d_start_end = {}
	for i in range(len(arr_start_end["RUN"])):
		d_start_end[int(arr_start_end["RUN"][i])] = [int(arr_start_end["start"][i]), int(arr_start_end["end"][i])]

	list_runs_ana =   [int(elem) for elem in arr_start_end["RUN"] ]

	#Get the template pulses
	fpulseA = TF1("fpulseA", Normalised_TemplateChalA_FID837(), 0, 1023, 1)
	fpulseA.SetParameter(0,1)
	standard_pulse_A = np.array([fpulseA.Eval(i) for i in range(1024)])
	
	fpulseB = TF1("fpulseB", Normalised_TemplateChalB_FID837(), 0, 1023, 1)
	fpulseB.SetParameter(0,1)
	standard_pulse_B = np.array([fpulseB.Eval(i) for i in range(1024)])

	event_path = script_utils.create_directory("../Event_files/" + bolo_name + "/" + "Event_simu/")

	#remove bad runs
	list_runs_ana.remove(231600)
	list_runs_ana.remove(232100)

	#Loop over runs to fill the corresponding tree
	for run_ana in list_runs_ana[18:]:
		print run_ana
		fout = TFile(event_path + bolo_name + "_" + str(conv.ana2samba(run_ana)) + "_simu_event_tree.root", "recreate")
		tout = TTree("EdwTree", "EdwTree")
		
		#Branches of the tree
		WienerAmplA = array("f", [0.])
		WienerAmplB = array("f", [0.])
		True_heat   = array("f", [0.])
	
		#Arrays for the trace modification
		trace_A = array( "h", 1024*[ 0 ] )
		trace_B = array( "h", 1024*[ 0 ] )	

		tout.Branch("WienerAmplA", WienerAmplA, "WienerAmplA/F")
		tout.Branch("WienerAmplB", WienerAmplB, "WienerAmplB/F")
		tout.Branch("True_heat", True_heat, "True_heat/F")
		tout.Branch("Event","EdwEvent",Evt)

		#increase start and end index for precaution
		start_index = d_start_end[run_ana][0]
		end_index = d_start_end[run_ana][1]

		for i in range(start_index, end_index):

			sys.stdout.write('\r' + str(i)+' / '+str(-1+nEntries))
			sys.stdout.flush()
			
			t.GetEntry(i)


			Run = conv.samba2ana(Evt.Run())

			if int(Run) == int(run_ana):

				True_heat[0]   = np.random.uniform(0,10,1)[0]
				WienerAmplA[0] = True_heat[0]*d_channel_calib["Chal1"]
				WienerAmplB[0] = True_heat[0]*d_channel_calib["Chal2"]

				base_A =Evt.Pulse("chalA " + bolo_name).Trace()
				base_B =Evt.Pulse("chalB " + bolo_name).Trace()

				base_A  = np.array([int(elem) for elem in base_A])
				base_B  = np.array([int(elem) for elem in base_B])

				pulse_A = standard_pulse_A*WienerAmplA[0]
				pulse_B = standard_pulse_B*WienerAmplB[0]

				pulse_A = [int(round(pulse_A[i])) for i in range(1024)]
				pulse_B = [int(round(pulse_B[i])) for i in range(1024)]

				pulse_A+=base_A 
				pulse_B+=base_B

				for i in range(1024):
					trace_A[i] = pulse_A[i]
					trace_B[i] = pulse_B[i]

				Evt.Pulse("chalA " + bolo_name).SetTrace(1024, trace_A)
				Evt.Pulse("chalB " + bolo_name).SetTrace(1024, trace_B)

				tout.Fill()

		tout.Write()
		fout.Close()
		del tout
		del fout
		print

