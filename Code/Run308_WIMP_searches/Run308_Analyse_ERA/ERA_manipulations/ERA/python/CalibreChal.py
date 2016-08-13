#! /usr/bin/env python

# TODO: le faire aussi en fct de la tension!!

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ROOT
import os,re,sys
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams,SambaUtils
from CalibreChalFct import *
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
eradir=gParams["eradir"]+"/"
voie=gParams["voie"]
if voie[0:4]!="Chal" : print "Voie incorrecte?",voie
saveplots=int(gParams["saveplots"]) if "saveplots" in gParams else 1
################################

channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=0)
bolodir=anadir+"/"+bolo
ampldir=bolodir+"/Amplitudes/"
figdir=bolodir+"/Figures/"
listruns=RunParams.ReadListeRuns(bolodir)
## NB: pour l'instant, tous les runs... (fond + gamma)
## il faudra pouvoir choisir (au minimum virer les neutrons)
gains=RunParams.ReadGains(bolodir,voie)

# d'abord on definit une grossiere coupure chi2 sur la voie consideree
# (idem CalibreIon et GainChalEvol... on pourrait faire une fct commune?)
t_chi2check=TChain("wienerntp_"+bolo+"_"+voie,"Edelweiss Wiener Data Tree")
for run in listruns :
    t_chi2check.AddFile(ampldir+"tmpwiener_"+run.Name+"_"+voie+"_"+bolo+".root")
fit = TF1("fit","[0]+[1]*x+[2]*x*x",10,1000)
hg = TH2F("hg","Chi2 distribution",1000,10,1000,1000,0,10)
fit.FixParameter(0,1)
fit.FixParameter(1,0)
t_chi2check.Draw("WienerChi2:WienerAmpl >> hg","WienerChi2<10")
hg.Fit("fit","L")	
chicutfct = TF1("cutfct","[0]+[1]*x+[2]*x*x",0,100000)
chicutfct.SetParameter(0,2)
chicutfct.SetParameter(1,0)
chicutfct.SetParameter(2,4*fit.GetParameter(2))
chicutfct.SetLineStyle(2)
chicutfct.Draw("same")
gPad.Update()
foo=raw_input("Press enter to continue")	
chicutcoefs=[chicutfct.GetParameter(0),chicutfct.GetParameter(2)]

efid=[]
achal_norm=[]
for run in listruns :
   fchal=TFile(ampldir+"tmpwiener_"+run.Name+"_"+voie+"_"+bolo+".root","READ")
   tchal=fchal.Get("wienerntp_"+bolo+"_"+voie)
   fion=TFile(ampldir+"eion_"+run.Name+"_"+bolo+".root","READ")
   tion=fion.Get("eionntp_"+bolo)
   fbasic=TFile(ampldir+"basicntp_"+run.Name+"_"+bolo+".root","READ")
   tbasic=fbasic.Get("basicntp_"+bolo)
   current_tmin=0
   current_tmax=0
   for i in range(tchal.GetEntries()) :
       tchal.GetEntry(i)
       tion.GetEntry(i)
       tbasic.GetEntry(i)
       if tion.EventFlag==2 and tchal.WienerAmpl>0 and tion.Chi2CutCol1==1 and tion.Chi2CutCol2==1 and tchal.WienerChi2<chicutcoefs[0]+chicutcoefs[1]*tchal.WienerAmpl*tchal.WienerAmpl:
           # "Lazy way" pour recuperer le gain... peut etre ameliore si trop lent en cpu...
           gain=RunParams.GetGain(gains,tbasic.DateSec)
	   if gain!=0 :
	       efid.append(tion.Eion)
               achal_norm.append(tchal.WienerAmpl/gain)

all_e_sur_a=[efid[i]/achal_norm[i] for i in range(len(efid))]
all_loga=[math.log10(a) for a in achal_norm]
# temporaire tant qu'on a pas des cuts (chi2 ...) plus efficaces:
e_sur_a=[all_e_sur_a[i] for i in range(len(efid)) if all_loga[i]>1 and all_e_sur_a[i]>0.6 and all_e_sur_a[i]<1.5]
loga=[all_loga[i] for i in range(len(efid)) if all_loga[i]>1 and all_e_sur_a[i]>0.6 and all_e_sur_a[i]<1.5]

plt.plot(loga,e_sur_a,'o')
plt.xlabel("Log10(Anorm)")
plt.ylabel("Efid/Anorm")
plt.title("Heat NL calib")
popt,pcov=curve_fit(ChalNLFunc,np.asarray(loga),np.asarray(e_sur_a),p0=np.array([1,-0.1,1.05,0.2]))
print "Fit results: {0:.3f} {1:.3f} {2:.3f} {3:.3f}".format(popt[0],popt[1],popt[2],popt[3])
x=np.arange(0,4,0.01)
z=ChalNLFunc(x,popt[0],popt[1],popt[2],popt[3])
plt.plot(x,z,color='r')
plt.savefig(figdir+"calib_"+voie+".png") # TODO: sauver plot meme si interactif=0
plt.show()

filename=bolodir+"/liste_calibs.txt"
if (os.path.exists(filename)==0) :
    outputfile=open(filename,"w")
    outputfile.write("# Voie type tinf tsup coefficients\n")
    outputfile.close()
outputfile=open(filename,"a")
outputfile.write('"'+channels[voie]+'" StandardChal 0 0 {0:.3f} {1:.3f} {2:.3f} {3:.3f}\n'.format(popt[0],popt[1],popt[2],popt[3]))
outputfile.close()
