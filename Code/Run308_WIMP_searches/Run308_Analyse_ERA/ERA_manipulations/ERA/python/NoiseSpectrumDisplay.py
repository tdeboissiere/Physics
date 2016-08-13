#! /usr/bin/env python

# TODO!

from GlobalParams import ReadGlobalParams
import os, sys

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
eradir=gParams['eradir']+"/"
bolo=gParams["bolo"]
################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.figure import Figure
import random
from datetime import datetime
import RunParams
from ROOT import *
ROOT.gSystem.Load(eradir+'lib/EraLib.so')
from ROOT import NoiseSpectrum
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import Tkinter as Tk

class NoiseSpectrumDisplay:
    def __init__(self,tkroot,anadir,bolo):
        self.anadir=anadir
        self.bolo=bolo
        self.BoloDir=anadir+"/"+bolo+"/"
        self.spectradir=self.BoloDir+"Spectra/"
        self.noise=NoiseSpectrum()
        # Liste voies
        self.channel_list=RunParams.ReadBoloChannels(self.bolo,self.anadir,remove_none=1)
        self.channel_dict=RunParams.ReadBoloChannels(self.bolo,self.anadir,remove_none=1,dico=1)
        self.Channel=self.channel_list[0]
        self.TkVarVoie=Tk.StringVar()
        self.TkVarVoie.set(self.Channel)
        # Liste runs
        self.runlist=[run.Name for run in RunParams.ReadListeRuns(self.BoloDir)]
        self.Run=(self.runlist)[0]
        self.TkVarRun=Tk.StringVar()
        self.TkVarNbPeriods=Tk.StringVar()
        self.TkVarPeriodNum=Tk.StringVar()
        self.TkVarRun.set(self.Run)
        self.updaterun()

        # (un frame gauche et un frame droit qui en contient d'autres)
        # Frame gauche : le choix des voies, plots, runs
        ################################
        self.framegauche=Tk.Frame(tkroot)
        self.framegauche.pack(side=Tk.LEFT,padx=10,fill=Tk.BOTH)
        for voie in self.channel_list:
            radio_voie=Tk.Radiobutton(self.framegauche,text=voie.ljust(20),variable=self.TkVarVoie,value=voie,command=self.updatevoie)
            radio_voie.pack(side=Tk.TOP,anchor=Tk.W)
        self.button_plotallperiods=Tk.Button(master=self.framegauche,text="Plot all periods",command=self.plotallperiods,width=14)
        self.button_plotallperiods.pack(side=Tk.TOP)
        #        self.button_plotallchannels=Tk.Button(master=self.framegauche,text="Plot all channels",command=self.plotallchannels,width=14)
        #        self.button_plotallchannels.pack(side=Tk.TOP)
        self.button_plotevol=Tk.Button(master=self.framegauche,text="2D noise evolution",command=self.plotevolspectrum,width=14)
        self.button_plotevol.pack(side=Tk.TOP)
        # Liste runs
        self.button_run=Tk.Button(master=self.framegauche,text="Select run",command=self.choixrun)
        self.button_run.pack(side=Tk.BOTTOM,pady=10)
        self.frameruns=Tk.Frame(self.framegauche)
        self.frameruns.pack(side=Tk.BOTTOM)
        self.scrollbar_run = Tk.Scrollbar(self.frameruns,orient=Tk.VERTICAL)
        self.listbox_run=Tk.Listbox(self.frameruns,yscrollcommand=self.scrollbar_run.set,height=4,width=15)
        self.scrollbar_run.config(command=self.listbox_run.yview)
        self.listbox_run.pack(side=Tk.LEFT)
        self.scrollbar_run.pack(side=Tk.LEFT, fill=Tk.Y)
        for r in self.runlist : self.listbox_run.insert(Tk.END,r)

        # A droite de ce frame : Le plot
        ################################
        self.framedroit=Tk.Frame(tkroot)
        self.framedroit.pack(side=Tk.RIGHT,expand=Tk.YES,fill=Tk.BOTH)
        self.frameplot=Tk.Frame(self.framedroit)
        self.frameplot.pack(side=Tk.TOP,expand=Tk.YES,fill=Tk.BOTH)
        fig = Figure(figsize=(6,4), dpi=100)
        self.SubPlot=fig.add_subplot(111)
        self.Canvas=FigureCanvasTkAgg(fig, master=self.frameplot)
        self.Canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.plotsinglespectrum()
        toolbar=NavigationToolbar2TkAgg(self.Canvas,self.frameplot)
        toolbar.pack(side=Tk.TOP)
        #toolbar.configure(background="grey")
        toolbar.update()

        # Un frame en dessous
        ################################
        self.frameinferieur=Tk.Frame(self.framedroit)
        self.frameinferieur.pack(side=Tk.TOP,fill=Tk.BOTH,padx=5,pady=5)
        self.label_period=Tk.Label(self.frameinferieur,text="  Select period:")
        self.label_period.pack(side=Tk.LEFT)
        self.entry_period=Tk.Entry(self.frameinferieur,width=4,textvariable=self.TkVarPeriodNum)
        self.entry_period.pack(side=Tk.LEFT)
        self.button_period=Tk.Button(master=self.frameinferieur,text="Go",command=self.choixperiod)
        self.button_period.pack(side=Tk.LEFT)
        self.label_nbperiods=Tk.Label(self.frameinferieur,textvariable=self.TkVarNbPeriods)
        self.label_nbperiods.pack(side=Tk.LEFT)
        self.label_nbperiods2=Tk.Label(self.frameinferieur,text="period(s)")
        self.label_nbperiods2.pack(side=Tk.LEFT)

        self.button_quit = Tk.Button(master=self.frameinferieur, text='Quit', command=sys.exit)
        self.button_quit.pack(side=Tk.RIGHT)

    #######################################  
    # Fonctions utilitaires du clickodrome
    #######################################

    def updaterun(self):
        self.NbPeriods=0
        self.StartValidity_Periods=[]
        noisefilename=self.spectradir+"spectra_"+self.Run+"_"+self.bolo+".root"
        noisefile=TFile(noisefilename,"READ")
        noisetree=noisefile.Get("SpectrumTree")
        noise=NoiseSpectrum()
        noisetree.SetBranchAddress("Spectrum",noise)
        for i in range(noisetree.GetEntries()) :
            noisetree.GetEntry(i)
            if noise.Channel()==self.Channel :
                self.NbPeriods+=1
                self.StartValidity_Periods.append(noise.StartValidity())
        noisefile.Close()
        self.PeriodNum=0
        self.TkVarNbPeriods.set(str(self.NbPeriods))
        self.TkVarPeriodNum.set(str(self.PeriodNum))

    def choixrun(self):
        #self.Run=self.TkVarRun.get()
        selection=self.listbox_run.curselection()
        if len(selection)!=1 :
            print "Wrong run selection"
            return
        self.Run=self.runlist[int(selection[0])]
        self.updaterun()
        self.plotsinglespectrum()

    def choixperiod(self):
        thenum=int(self.TkVarPeriodNum.get())
        if thenum<0 or thenum>=self.NbPeriods :
            print "Wrong period number!"
        else :
            self.PeriodNum=thenum
            self.plotsinglespectrum()
        self.TkVarPeriodNum.set(str(self.PeriodNum))

    def updatevoie(self):
        self.Channel=self.TkVarVoie.get()
        self.plotsinglespectrum()
        
    def plotsinglespectrum(self):
        noisefilename=self.spectradir+"spectra_"+self.Run+"_"+self.bolo+".root"
        self.noise=RunParams.GetNoiseSpectrum(noisefilename,self.Channel,self.StartValidity_Periods[self.PeriodNum]+5)
        freqs,bruit=self.normalizenoise()
        xmin=min(freqs)
        xmax=max(freqs)
        self.SubPlot.clear()
        self.SubPlot.set_xlabel('Frequency [Hz]')
        self.SubPlot.set_ylabel('Noise [ADU/sqrt(Hz)]')
        self.SubPlot.set_yscale('log')
        self.SubPlot.set_xscale('log')
        self.SubPlot.plot(freqs,bruit,color="blue")
        self.SubPlot.set_xlim([xmin,xmax])
        self.Canvas.show()

    def normalizenoise(self):
        # Calcule puissance spectrale
        # fait un truc en ADU/sqrt(Hz)...
        spectre=self.noise.Spectrum()
        sampling_ms=self.noise.Sampling_ms()
        N=len(spectre)
        norm_noise=[spectre[0]] # frequence zero
        for i in range(1,N/2):
            norm_noise.append(spectre[i]+spectre[N-i])
        norm_noise.append(spectre[N/2]) # frequence Nyquist
        norm_noise=[sqrt(x) for x in norm_noise]
        delta_f_khz=1./(N*sampling_ms)
        frequencies_hz=[1000*delta_f_khz*i for i in range(N/2+1)] # N/2+1 frequencies
        return frequencies_hz,norm_noise

    def plotallperiods(self):
        noisefilename=self.spectradir+"spectra_"+self.Run+"_"+self.bolo+".root"
        self.SubPlot.clear()
        self.SubPlot.set_xlabel('Frequency [Hz]')
        self.SubPlot.set_ylabel('Noise [ADU/sqrt(Hz)]')
        self.SubPlot.set_yscale('log')
        self.SubPlot.set_xscale('log')
        for periodnum in range(self.NbPeriods):
            self.noise=RunParams.GetNoiseSpectrum(noisefilename,self.Channel,self.StartValidity_Periods[periodnum]+5)
            freqs,bruit=self.normalizenoise()
            xmin=min(freqs)
            xmax=max(freqs)
            self.SubPlot.plot(freqs,bruit,color=[random.random(),random.random(),random.random()]) # code couleur [r,g,b] random
        self.SubPlot.set_xlim([xmin,xmax])
        self.Canvas.show()

    def plotallchannels(self):
        # !! A implementer !!
        pass

    def plotevolspectrum(self,timebin_hour=0.1):
        # !!! Fonction a ameliorer largement encore... !!!
        # !!! Essaie a faire: utiliser plt.matshow(array,0) au lieu de imshow ??
        grid_spectra=[]
        current_time=0
        t0=datetime.fromtimestamp(0)
        timebin_sec=timebin_hour*3600
        for run in self.runlist:
            noisefilename=self.spectradir+"spectra_"+run+"_"+self.bolo+".root"
            noisefile=TFile(noisefilename,"READ")
            noisetree=noisefile.Get("SpectrumTree")
            noise=NoiseSpectrum()
            noisetree.SetBranchAddress("Spectrum",noise)
            for i in range(noisetree.GetEntries()) :
                noisetree.GetEntry(i)
                if noise.Channel()==self.Channel :
                    self.noise=noise
                    freqs,n=self.normalizenoise()
                    n=[log10(x) for x in n]
                    n.reverse()
                    nb_bins_delay=int((noise.StartValidity()-current_time)/timebin_sec)
                    if current_time==0 :
                        nb_bins_delay=0
                        t0=datetime.fromtimestamp(noise.StartValidity())
                    nb_bins=int((noise.EndValidity()-noise.StartValidity())/timebin_sec)
                    current_time=noise.EndValidity()
                    if nb_bins<0 or nb_bins_delay<0 : print "pbl tps..."
                    for p in range(nb_bins_delay) :
                        grid_spectra.append([0.0 for x in n])
                    for p in range(nb_bins) :
                        grid_spectra.append(n)
            noisefile.Close()
        self.SubPlot.clear()
        self.SubPlot.set_xlabel('Time since '+t0.isoformat(' ')+' [hours]')
        self.SubPlot.set_ylabel('Frequency [Hz]')
        self.SubPlot.imshow(np.transpose(grid_spectra),aspect='auto',interpolation='none',cmap=cm.jet,extent=(0,len(grid_spectra)*timebin_hour,freqs[0],max(freqs)))
        # self.SubPlot.pcolor(np.transpose(grid_spectra),cmap=cm.jet)
        #  self.SubPlot.set_yscale('log') # IL FAUT TROUVER QQUE CHOSE...
        self.Canvas.show()


        
################################
# Le "main" de la gui
################################

tkroot = Tk.Tk()
tkroot.wm_title("EDELWEISS noise spectrum")
disp=NoiseSpectrumDisplay(tkroot,anadir,bolo)
tkroot.mainloop() # ou Tk, semble pareil

