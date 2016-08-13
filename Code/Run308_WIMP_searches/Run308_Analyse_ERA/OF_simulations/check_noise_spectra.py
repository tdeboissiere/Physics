from ROOT import *
ROOT.gSystem.Load("/home/irfulx204/mnt/tmain/Desktop/New_ERA/ERA/lib/EraLib.so")
from ROOT import EdwPulse,EdwEvent,FitPulse,EdwTemplate,NoiseSpectrum

import sys, glob, os
import scipy.integrate as SI
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from matplotlib import gridspec
from array import array
import conversion_ANA_SAMBA as conv

import PyROOTPlots as PyRPl

def normalizenoise(SPEC):
    # Calcule puissance spectrale
    # fait un truc en ADU/sqrt(Hz)...
    spectre=SPEC.Spectrum()
    sampling_ms=SPEC.Sampling_ms()
    N=len(spectre)
    norm_noise=[spectre[0]] # frequence zero
    for i in range(1,N/2):
        norm_noise.append(spectre[i]+spectre[N-i])
    norm_noise.append(spectre[N/2]) # frequence Nyquist
    norm_noise=[sqrt(x) for x in norm_noise]
    delta_f_khz=1./(N*sampling_ms)
    frequencies_hz=[1000*delta_f_khz*i for i in range(N/2+1)] # N/2+1 frequencies

    return frequencies_hz,norm_noise

def compute_low_freq_integral(bolo_name):

    plt.ion()

    fana = np.loadtxt("../list_runs_ana_s2.txt", dtype = "str")
    fsamba = np.loadtxt("../list_runs_samba_s2.txt", dtype = "str")

    d_ana_to_samba = {}
    d_samba_to_ana = {}

    for i, elem in enumerate(fana) :
        d_ana_to_samba[fana[i]] = fsamba[i]
        d_samba_to_ana[fsamba[i]] = fana[i]

    #cd Desktop/New_ERA/ERA/
    list_pulsefiles=[]
    list_pulsefiles.extend(glob.glob("../Spectra_Files/*.root*"))
    list_pulsefiles=sorted(list_pulsefiles)

    Tree=TChain("SpectrumTree")
    for index, filou in enumerate(list_pulsefiles) :  
        sys.stdout.write('\r' + str(index)+' / '+str(-1+len(list_pulsefiles)))
        sys.stdout.flush()
        Tree.AddFile(filou)  

    print

    Spec=NoiseSpectrum()
    Tree.SetBranchAddress("Spectrum",Spec)
    nEntries=Tree.GetEntries()
    # nEntries =1000

    list_I=[]
    list_bad_periods = []

    for i in range(0,nEntries):

        # sys.stdout.write('\r' + str(i)+' / '+str(-1+nEntries))
        # sys.stdout.flush()

        Tree.GetEntry(i) 

        if (Spec.Channel() == "chalA FID837"):

            frequencies_hz,norm_noise = normalizenoise(Spec)
            I1 = SI.simps(norm_noise[2:11], frequencies_hz[2:11])
            list_I.append(I1)
            
            if (I1 > 1200) :
                list_bad_periods.append([Spec.StartValidity(), Spec.EndValidity()])

            # plt.plot(frequencies_hz, norm_noise)
            # plt.show()
            # raw_input()
        
        # plt.close()   
    np.savetxt("./Text_files/" + bolo_name + "_bad_noise_periods.txt", np.array(list_bad_periods), delimiter = ",")
    # np.savetxt("./Text_files/" + bolo_name + "_noise_integral.txt", np.array(list_I))
    # list_I = list(np.loadtxt("./Text_files/" + bolo_name + "_noise_integral.txt"))

    # h = TH1F("h", "h", 200, 0.5*min(list_I), 10000)
    # for elem in list_I:
    #     h.Fill(elem)

    # PyRPl.process_TH1(h, X_title = "Integral 2 to 10 Hz", Y_title= "# of periods")

    # cc = TCanvas("cc", "cc")
    # h.Draw()
    # print float(h.Integral(1,10))/float(h.Integral()), h.GetBinCenter(10)
    # print float(h.Integral(1,20))/float(h.Integral()), h.GetBinCenter(20)
    # print float(h.Integral(1,30))/float(h.Integral()), h.GetBinCenter(30)
    # print float(h.Integral(1,40))/float(h.Integral()), h.GetBinCenter(40)
    # print float(h.Integral(1,50))/float(h.Integral()), h.GetBinCenter(50)
    # print float(h.Integral(1,60))/float(h.Integral()), h.GetBinCenter(60)


    # raw_input()
    # cc.Print("./Figures/" + bolo_name + "_noise_integral.png")

def get_bad_periods(bolo_name):


    # plt.ion()

    # fana = np.loadtxt("../list_runs_ana_s2.txt", dtype = "str")
    # fsamba = np.loadtxt("../list_runs_samba_s2.txt", dtype = "str")

    # d_ana_to_samba = {}
    # d_samba_to_ana = {}

    # for i, elem in enumerate(fana) :
    #     d_ana_to_samba[fana[i]] = fsamba[i]
    #     d_samba_to_ana[fsamba[i]] = fana[i]

    # #cd Desktop/New_ERA/ERA/
    # list_pulsefiles=[]
    # list_pulsefiles.extend(glob.glob("../Spectra_Files/*.root*"))
    # list_pulsefiles=sorted(list_pulsefiles)

    # Tree=TChain("SpectrumTree")
    # for index, filou in enumerate(list_pulsefiles) :  
    #     sys.stdout.write('\r' + str(index)+' / '+str(-1+len(list_pulsefiles)))
    #     sys.stdout.flush()
    #     Tree.AddFile(filou)  

    # print

    # Spec=NoiseSpectrum()
    # Tree.SetBranchAddress("Spectrum",Spec)
    # nEntries=Tree.GetEntries()
    # # nEntries =1000

    # list_I=[]

    # for i in range(0,nEntries):

    #     # sys.stdout.write('\r' + str(i)+' / '+str(-1+nEntries))
    #     # sys.stdout.flush()

    #     Tree.GetEntry(i) 

    #     if (Spec.Channel() == "chalA FID837"):

    #         frequencies_hz,norm_noise = normalizenoise(Spec)
    #         temp_noise = norm_noise[2:11]
    #         temp_freq = frequencies_hz[2:11]
    #         # I1 = SI.simps(norm_noise[2:11], frequencies_hz[2:11])
    #         list_I.append(SI.simps(temp_noise, temp_freq))
            
    #         # plt.plot(frequencies_hz, norm_noise)
    #         # plt.show()
    #         # raw_input()
        
    #     # plt.close()   
    # # np.savetxt("./Text_files/" + bolo_name + "_noise_integral.txt", np.array(list_I))
    list_I = list(np.loadtxt("./Text_files/" + bolo_name + "_noise_integral.txt"))

    h = TH1F("h", "h", 200, 0.5*min(list_I), 10000)
    for elem in list_I:
        h.Fill(elem)

    PyRPl.process_TH1(h, X_title = "Integral 2 to 10 Hz", Y_title= "# of periods")

    cc = TCanvas("cc", "cc")
    h.Draw()
    print float(h.Integral(1,10))/float(h.Integral()), h.GetBinCenter(10)
    print float(h.Integral(1,20))/float(h.Integral()), h.GetBinCenter(20)
    print float(h.Integral(1,30))/float(h.Integral()), h.GetBinCenter(30)
    print float(h.Integral(1,40))/float(h.Integral()), h.GetBinCenter(40)
    print float(h.Integral(1,50))/float(h.Integral()), h.GetBinCenter(50)
    print float(h.Integral(1,60))/float(h.Integral()), h.GetBinCenter(60)


    raw_input()
    cc.Print("./Figures/" + bolo_name + "_noise_integral.png")

def plot_noise_spectra():

    plt.ion()

    #cd Desktop/New_ERA/ERA/
    list_pulsefiles=[]
    list_pulsefiles.extend(glob.glob("./ROOT_files/SpectraDir/*.root*"))
    list_pulsefiles=sorted(list_pulsefiles)

    Tree=TChain("SpectrumTree")
    for index, filou in enumerate(list_pulsefiles) :  
        sys.stdout.write('\r' + str(index)+' / '+str(-1+len(list_pulsefiles)))
        sys.stdout.flush()
        Tree.AddFile(filou)  

    print

    Spec=NoiseSpectrum()
    Tree.SetBranchAddress("Spectrum",Spec)
    nEntries=Tree.GetEntries()
    # nEntries =1000

    for i in range(0,nEntries):


        Tree.GetEntry(i) 

        if (Spec.Channel() == "chalA FID837"):

            frequencies_hz,norm_noise = normalizenoise(Spec)
            plt.plot(frequencies_hz, norm_noise)
            plt.show()
            # raw_input()
        
        # plt.close()          
    raw_input()



bolo_name = "FID837"

# compute_low_freq_integral(bolo_name)
plot_noise_spectra()