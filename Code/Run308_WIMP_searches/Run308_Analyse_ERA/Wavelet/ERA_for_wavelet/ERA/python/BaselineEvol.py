#! /usr/bin/env python

import sys,os
import numpy as np
from GlobalParams import ReadGlobalParams
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams,SambaUtils

################################
# PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
voie=gParams["voie"]
interactif=int(gParams["interactif"]) if "interactif" in gParams else 0
overwrite=int(gParams["overwrite"]) if "overwrite" in gParams else 1
saveplots=int(gParams["saveplots"]) if "saveplots" in gParams else 0
################################

bolodir=anadir+"/"+bolo
ampldir=bolodir+"/Amplitudes/"
figdir=bolodir+"/Figures/"

channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1)
coefs=RunParams.GetCalibCoef(bolodir,channels[voie])
# Pour la chaleur, on utilise le regime "basse energie" de la fct de calibration
if voie[0:4]=="Chal" : calibcoef=coefs[2]
else : calibcoef=coefs[0] # ionisation = 1 seul coef
liste_periods=RunParams.ReadListePeriods(bolodir)
liste_gains=[]
if voie[0:4]=="Chal" :
    liste_gains=RunParams.ReadGains(bolodir,voie=voie)

# Selection des periodes a suivre (au besoin) et creation du fichier
filename=bolodir+"/baseline_"+voie+".txt"
if overwrite==0 and os.path.exists(filename) :
    existing_baselines=RunParams.ReadBaselines(bolodir,voie)
    for base in existing_baselines:
        liste_periods=[p for p in liste_periods if (p.Tinf!=base.Tinf or p.Tsup!=base.Tsup)]
if os.path.exists(filename)==0 or overwrite==1 :
    outputfile=open(filename,"w")
    outputfile.write("# Tinf Tsup FWHM(keV) FWHM_from_spectrum(keV) \n")
    outputfile.close()
outputfile=open(filename,"a")

nbins=100 # ou 400??
min_histo=-5
max_histo=5
normal_law=TF1("f","gaus",0.5*min_histo,0.5*max_histo)
normal_law.FixParameter(1,0)

periodnum=0
gain=0
for theperiod in liste_periods :
    tinf=theperiod.Tinf
    tsup=theperiod.Tsup
    run=theperiod.Run
    periodnum+=1 # Ca sert pour le nom du plot si on le sauvegarde
    fbasic=TFile(ampldir+"basicntp_"+run+"_"+bolo+".root","READ")
    tbasic=fbasic.Get("basicntp_"+bolo)
    f=TFile(ampldir+"tmpwiener_"+run+"_"+voie+"_"+bolo+".root","READ")
    t=f.Get("wienerntp_"+bolo+"_"+voie)
    Baseline=TH1F("Baseline","Baseline",nbins,min_histo,max_histo)
    SigmaAmpl_List=[] # les sigmaampl du ntp => on prend la moyenne, simplement.
    for ievt in range(tbasic.GetEntries()) :
        tbasic.GetEntry(ievt)
        if tbasic.DateSec<tinf or tbasic.DateSec>tsup: continue
        t.GetEntry(ievt)
        if voie[0:4]=="Chal" : gain=RunParams.GetGain(liste_gains,tbasic.DateSec)
        # Remplissage de l'histo Baseline:
        if tbasic.Saturation==0 and t.WienerZeroAmpl!=0 and tbasic.IsBoloTrigger==0 and t.WienerChi2<2.0 and t.WienerChi2>0.5 and t.WienerZeroAmpl <10 and eval("tbasic.PileUp"+voie+"==0") :
            # NB: coupure assez mechante sur chi2 (0.5)
            if voie[0:4]=="Chal" :
                if gain!=0 :
		    zeroampl_kev=calibcoef*t.WienerZeroAmpl/gain
            	    Baseline.Fill(zeroampl_kev)
	    else :
		zeroampl_kev=t.WienerZeroAmpl/calibcoef
		Baseline.Fill(zeroampl_kev)
        # Remplissage de la liste des SigmaAmpl :
        if t.WienerSigmaAmpl!=0 :
            if voie[0:4]=="Chal" :
                if gain!=0 : SigmaAmpl_List.append(2.35482*calibcoef*t.WienerSigmaAmpl/gain)
            else :
                SigmaAmpl_List.append(2.35482*t.WienerSigmaAmpl/calibcoef)
    # Calcul du fwhm a partir de l'histo Baseline:
    if Baseline.GetEntries()==0 : 
	print run,periodnum," : No data to compute baseline.."
	fwhm=0
    else :
	if interactif==1 :
            Baseline.Fit("f","RQ")
            gPad.Update()
            gPad.SaveAs(figdir+"baseline_"+voie+"_"+str(periodnum)+".gif")
            # 1) sauver aussi si interactif = 0...
            # 2) ATTENTION: ca va pas marcher en mode "update" (numerotation des periodes)
            # mais a priori on voudra pas garder tous ces plots...
            foo=raw_input("Next..")
        else :
            Baseline.Fit("f","RQN")
        fwhm=2.35482*normal_law.GetParameter(2);
    # Le FWHM a partir de la liste des SigmaAmpl:
    fwhm_fromspectra=np.mean(SigmaAmpl_List) if len(SigmaAmpl_List)!=0 else 0
    outputfile.write(str(tinf)+" "+str(tsup)+" "+("{0:.3f}").format(fwhm)+" "+("{0:.3f}").format(fwhm_fromspectra)+"\n")
    f.Close()
    fbasic.Close()

outputfile.close()

