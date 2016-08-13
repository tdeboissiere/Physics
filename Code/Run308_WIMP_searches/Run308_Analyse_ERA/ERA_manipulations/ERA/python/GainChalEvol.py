#!/usr/bin/env python

# TODO = interactivite (a definir) pour choisir les valeurs tabulees du gain...
# par exemple pr une "periode" donnee ca pourrait etre une fct lineaire du temps (2 params)

################################
# PARAMETRES A EDITER EN DUR
eion_min=50 # range en Eion (bosse du compton, sans les raies du Ba)
eion_max=280
################################

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams,SambaUtils
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
voie=gParams["voie"]
if voie[0:4]!="Chal" : print "Voie incorrecte?",voie
overwrite=int(gParams["overwrite"]) if "overwrite" in gParams else 1
interactif=int(gParams["interactif"]) if "interactif" in gParams else 1
saveplots=int(gParams["saveplots"]) if "saveplots" in gParams else 1
measure_gain_by_run=1 # parametre en dur a ce stade!
################################

#channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=0)
bolodir=anadir+"/"+bolo
ampldir=bolodir+"/Amplitudes/"
figdir=bolodir+"/Figures/"

listruns=RunParams.ReadListeRuns(bolodir)
listperiods=RunParams.ReadListePeriods(bolodir)

# Selection des periodes et runs (au besoin) et creation du fichier
filename=bolodir+"/gain_"+voie+".txt"
if overwrite==0 and os.path.exists(filename) :
    existing_gains=RunParams.ReadGains(bolodir,voie)    
    for gain in existing_gains:
        listperiods=[p for p in listperiods if (p.Tinf!=gain.Tinf or p.Tsup!=gain.Tsup)]
    runs_to_keep=[p.Run for p in listperiods]
    listruns=[r for r in listruns if r.Name in runs_to_keep]

if os.path.exists(filename)==0 or overwrite==1 :
    outfile=open(filename,"w")
    outfile.write("# Tinf Tsup Gain(ADUs/keV)\n")
    outfile.close()

# d'abord on definit une grossiere coupure chi2 sur la voie consideree
# (idem CalibreIon)
t_chi2check=TChain("wienerntp_"+bolo+"_"+voie,"Edelweiss Wiener Data Tree")
for run in listruns :
    t_chi2check.AddFile(ampldir+"tmpwiener_"+run.Name+"_"+voie+"_"+bolo+".root")
fit = TF1("fit","[0]+[1]*x+[2]*x*x",10,1000)
hg = TH2F("hg","Chi2 distribution",1000,10,1000,1000,0,10)
fit.FixParameter(0,1)
fit.FixParameter(1,0)
option_plot="goff" if interactif==0 else ""
t_chi2check.Draw("WienerChi2:WienerAmpl >> hg","WienerChi2<10",option_plot)
if interactif==1 :
    hg.Fit("fit","L")
else :
    hg.Fit("fit","LN")
chicutfct = TF1("cutfct","[0]+[1]*x+[2]*x*x",0,100000)
chicutfct.SetParameter(0,2)
chicutfct.SetParameter(1,0)
chicutfct.SetParameter(2,4*fit.GetParameter(2))
chicutfct.SetLineStyle(2)
if interactif==1 :
    chicutfct.Draw("same")
    gPad.Update()
    foo=raw_input("Press enter to continue")
chicutcoefs=[chicutfct.GetParameter(0),chicutfct.GetParameter(2)]

vect_gain=[]
vect_gain_nochicut_ion=[]
vect_gain_nochicut_chal=[]
print "Calcul gains pour runs",[run.Name for run in listruns]
for run in listruns :
    fchal=TFile(ampldir+"tmpwiener_"+run.Name+"_"+voie+"_"+bolo+".root","READ")
    tchal=fchal.Get("wienerntp_"+bolo+"_"+voie)
    fion=TFile(ampldir+"eion_"+run.Name+"_"+bolo+".root","READ")
    tion=fion.Get("eionntp_"+bolo)
    fbasic=TFile(ampldir+"basicntp_"+run.Name+"_"+bolo+".root","READ")
    tbasic=fbasic.Get("basicntp_"+bolo)
    for i in range(tchal.GetEntries()) :
        tchal.GetEntry(i)
        tion.GetEntry(i)
        tbasic.GetEntry(i)
        # on choisit uniquement les evts fiduciels bien reconstruits.. (chi2 cut)
        if tion.EventFlag==2 and tion.Eion>eion_min and tion.Eion<eion_max and tchal.WienerAmpl!=0 :
            if tion.Chi2CutCol1==1 and tion.Chi2CutCol2==1 :
                vect_gain_nochicut_chal.append((tbasic.DateSec,tchal.WienerAmpl/tion.Eion))
            if tchal.WienerChi2<chicutcoefs[0]+chicutcoefs[1]*tchal.WienerAmpl*tchal.WienerAmpl :
                vect_gain_nochicut_ion.append((tbasic.DateSec,tchal.WienerAmpl/tion.Eion))
            if tion.Chi2CutCol1==1 and tion.Chi2CutCol2==1 and tchal.WienerChi2<chicutcoefs[0]+chicutcoefs[1]*tchal.WienerAmpl*tchal.WienerAmpl :
                vect_gain.append((tbasic.DateSec,tchal.WienerAmpl/tion.Eion))

listmedians=[]
for period in listperiods :
    local_tinf=period.Tinf
    local_tsup=period.Tsup
    if measure_gain_by_run==1 : # dans ce cas on coupe par run
        local_periods=[pp for pp in listperiods if pp.Run==period.Run]
        local_tinf=min([pp.Tinf for pp in local_periods])
        local_tsup=min([pp.Tsup for pp in local_periods])
    ratios_select=[r for (t,r) in vect_gain if t>local_tinf and t<local_tsup]
    if len(ratios_select) != 0 :
        listmedians.append(np.median(ratios_select))
    else :
        listmedians.append(0) # convention

# Fichier gain : juste les medianes pour l'instant
outfile=open(filename,"a")
for i in range(len(listperiods)) :
    outfile.write(str(listperiods[i].Tinf)+" "+str(listperiods[i].Tsup)+" "+str(listmedians[i])+"\n")
outfile.close()

if interactif==1 :
    x=[t for (t,r) in vect_gain]
    y=[r for (t,r) in vect_gain]
    x_nochicut_chal=[t for (t,r) in vect_gain_nochicut_chal]
    y_nochicut_chal=[r for (t,r) in vect_gain_nochicut_chal]
    x_nochicut_ion=[t for (t,r) in vect_gain_nochicut_ion]
    y_nochicut_ion=[r for (t,r) in vect_gain_nochicut_ion]
    plt.ylim([-0.1,max(y)+1])
    p3,=plt.plot(x_nochicut_chal,y_nochicut_chal,'rs',ms=3)
    p2,=plt.plot(x_nochicut_ion,y_nochicut_ion,'gs',ms=3)
    p1,=plt.plot(x,y,'o',ms=5)
    plt.xlabel('Unix Time')
    plt.ylabel('Heat/Fiducial ion')
    plt.title('Gain chaleur')
    plt.legend((p1,p2,p3),('selected evts','no chi2 cut ion','no chi2 cut chal'),loc='lower right')
    for i in range(len(listmedians)) :
        plt.plot([listperiods[i].Tinf,listperiods[i].Tsup],[listmedians[i],listmedians[i]],color='r')
    plt.ylim([-0.1,max(y)+1])
    plt.savefig(figdir+"gainevol_"+voie+".png") # TODO: sauver plot meme si interactif=0
    plt.ylim([-0.1,max(y)+1])
    plt.show()
    plt.ylim([-0.1,max(y)+1])
