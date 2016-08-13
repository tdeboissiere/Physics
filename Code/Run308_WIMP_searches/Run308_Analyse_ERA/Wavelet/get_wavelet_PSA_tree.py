import numpy as np
import mlpy.wavelet as wave
from sklearn import preprocessing
from scipy import stats
from ROOT import *
from array import array
from mlpy import dtw_std
import PyROOTPlots as PyRPl
import conversion_ANA_SAMBA as conv
import sys
import glob
ROOT.gSystem.Load("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/lib/EraLib.so")
from ROOT import EdwPulse,EdwEvent,FitPulse,EdwTemplate,NoiseSpectrum

def get_wavelet_PSA_tree(file_name, bolo_name):

	"""Compute wavelet PSA for given event file name
	
	Detail:
		void

	Args:
		bolo_name = (str) bolo name 
		file_name = (str) event file name 
		
	Returns:
		void

	Raises:
		void
	"""

	#Get event file
	data_types = {"names": ("RUN", "SN"), "formats": ("i", "i")}
	pop_path = "./Populations/"+ bolo_name + "/"
	arr_events = np.loadtxt(pop_path + bolo_name + "_" + file_name + ".txt", delimiter=",",  dtype=data_types)

	d_events={}
	for i in range(arr_events.shape[0]):
		d_events[str(arr_events["RUN"][i]) + "_" + str(arr_events["SN"][i])] = True

	#Get reconstructed-wavelet-filtered pulse
	x_refA         = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_pulseA_ref_fit.txt")
	scales_refA    = wave.autoscales(N=x_refA.shape[0], dt=1, dj=0.05, wf='dog', p=2)
	X_refA         = wave.cwt(x=x_refA, dt=1, scales=scales_refA, wf='dog', p=2)
	X_ref_flatA    = X_refA.flatten("F")
	lA             = np.where(np.abs(X_ref_flatA) <4)[0]
	X_ref_flatA[lA] =0
	X_ref_flatA    = X_ref_flatA.reshape(X_refA.shape, order = "F")
	

	X_ref_flatA[100:200,:]=0
	Xi_refA = wave.icwt(X_ref_flatA, dt=1, scales=scales_refA, wf='dog', p=2)

	pearsA = array("f", [0.])
	dtw_valA = array("f", [0.])
	Run = array("i", [0])
	SambaNum = array("i", [0])

	trace_path = "../Trace_files/" + bolo_name + "/"
	t, f = PyRPl.open_ROOT_object(trace_path + bolo_name +"_ol_pa_trigged_trace_tree.root", "out")

	root_path = "./ROOT_files/" + bolo_name + "/"
	fout = TFile(root_path + bolo_name + "_" + "ol_pa_event_wavelet_PSA_tree.root", "recreate")
	tout = TTree("twavelet", "twavelet")

	tout.Branch("pearsA", pearsA, "pearsA/F")
	tout.Branch("dtw_valA", dtw_valA, "dtw_valA/F")

	tout.Branch("SambaNum", SambaNum, "SambaNum/I")
	tout.Branch("Run", Run, "Run/I")

	nEntries = t.GetEntries()

	for i in range(nEntries):

		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		t.GetEntry(i)

		SambaNum[0] = int(t.SambaNum)
		Run[0] = int(conv.samba2ana(t.Run[:-1]))

		try:
			# d_events[str(Run[0]) + "_" + str(SambaNum[0])]
			pulse_A = np.zeros(1024)
			for i in range(1024):
				pulse_A[i]=int(t.arr_A[i])

			#Wavelet action !
			scalesA = wave.autoscales(N=pulse_A.shape[0], dt=1, dj=0.05, wf='dog', p=2)
			XA = wave.cwt(x=pulse_A, dt=1, scales=scalesA, wf='dog', p=2)
			
			
			XflatA = XA.flatten("F")
			XflatA[lA]=0
			XflatA = XflatA.reshape(XA.shape, order = "F")
			XflatA[100:200,:]=0
			
			XiA = wave.icwt(XflatA, dt=1, scales=scalesA, wf='dog', p=2)
			
			pearsA[0] = stats.pearsonr(preprocessing.scale(XiA), preprocessing.scale(Xi_refA))[0]
			dtw_valA[0] =  dtw_std(preprocessing.scale(XiA), preprocessing.scale(Xi_refA), dist_only=True)
			
			tout.Fill()

		except KeyError:
			pass

	tout.Write()
	fout.Close()
	raw_input()

def get_wavelet_PSA_simu_tree(file_name, bolo_name):

	"""Compute wavelet PSA for given event file name
	
	Detail:
		Use the simulated data

	Args:
		bolo_name = (str) bolo name 
		file_name = (str) event file name 
		
	Returns:
		void

	Raises:
		void
	"""


	#Get reconstructed-wavelet-filtered pulse
	x_refA         = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name +"_pulseA_ref_fit.txt")
	scales_refA    = wave.autoscales(N=x_refA.shape[0], dt=1, dj=0.05, wf='dog', p=2)
	X_refA         = wave.cwt(x=x_refA, dt=1, scales=scales_refA, wf='dog', p=2)
	X_ref_flatA    = X_refA.flatten("F")
	lA             = np.where(np.abs(X_ref_flatA) <4)[0]
	X_ref_flatA[lA] =0
	X_ref_flatA    = X_ref_flatA.reshape(X_refA.shape, order = "F")
	

	X_ref_flatA[100:200,:]=0
	Xi_refA = wave.icwt(X_ref_flatA, dt=1, scales=scales_refA, wf='dog', p=2)

	pearsA = array("f", [0.])
	dtw_valA = array("f", [0.])
	Run = array("i", [0])
	SambaNum = array("i", [0])

	trace_path = "../Event_files/" + bolo_name + "/Event_simu/"
	list_run_event_file = sorted(glob.glob(trace_path + "/*.root*"))
	#keep a subselection to shorten processing
	short_run_file = [list_run_event_file[i] for i in [1,20,40,60,80]]
	list_short_run_name = [list_run_event_file[i][40:48] for i in [1,20,40,60,80]]

	np.savetxt("./Text_files/" + bolo_name + "/run_simu_used_for_PSA.txt", np.array(list_short_run_name), fmt="%s")

	c = TChain("EdwTree", "EdwTree")
	for run_file in short_run_file:
		c.AddFile(run_file)
	
	Evt=EdwEvent()
	c.SetBranchAddress("Event",Evt)
	
	root_path = "./ROOT_files/" + bolo_name + "/"
	fout = TFile(root_path + bolo_name + "_" + "simu_event_wavelet_PSA_tree.root", "recreate")
	tout = TTree("twavelet", "twavelet")

	tout.Branch("pearsA", pearsA, "pearsA/F")
	tout.Branch("dtw_valA", dtw_valA, "dtw_valA/F")

	tout.Branch("SambaNum", SambaNum, "SambaNum/I")
	tout.Branch("Run", Run, "Run/I")

	nEntries = c.GetEntries()

	for i in range(nEntries):

		#Progress bar
		sys.stdout.write("\r" + str(i)+" / "+str(-1+nEntries))
		sys.stdout.flush()

		c.GetEntry(i)

		SambaNum[0] = int(Evt.SambaNum())
		Run[0] = int(conv.samba2ana(Evt.Run()))
		
		PulseA=Evt.Pulse("chalA " + bolo_name)
		pulse_A=np.array(PulseA.Trace()).astype(int)
		
		#Wavelet action !
		scalesA = wave.autoscales(N=pulse_A.shape[0], dt=1, dj=0.05, wf='dog', p=2)
		XA = wave.cwt(x=pulse_A, dt=1, scales=scalesA, wf='dog', p=2)
			
			
		XflatA = XA.flatten("F")
		XflatA[lA]=0
		XflatA = XflatA.reshape(XA.shape, order = "F")
		XflatA[100:200,:]=0
			
		XiA = wave.icwt(XflatA, dt=1, scales=scalesA, wf='dog', p=2)
			
		pearsA[0] = stats.pearsonr(preprocessing.scale(XiA), preprocessing.scale(Xi_refA))[0]
		dtw_valA[0] =  dtw_std(preprocessing.scale(XiA), preprocessing.scale(Xi_refA), dist_only=True)
			
		tout.Fill()


	tout.Write()
	fout.Close()
	raw_input()

bolo_name = "FID837"
file_name = "all_events"
get_wavelet_PSA_tree(file_name, bolo_name)
# get_wavelet_PSA_simu_tree(file_name, bolo_name)