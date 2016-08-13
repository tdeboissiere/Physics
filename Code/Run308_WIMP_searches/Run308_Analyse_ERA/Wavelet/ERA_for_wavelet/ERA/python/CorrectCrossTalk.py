#!/usr/bin/env python

import sys,os
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams,SambaUtils
from array import array
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
overwrite=int(gParams["overwrite"]) if "overwrite" in gParams else 1
################################

bolodir=anadir+"/"+bolo
ampldir=bolodir+"/Amplitudes/"
listruns=RunParams.ReadListeRuns(bolodir)
listcrosstalks=RunParams.ReadListCrossTalks(bolodir)
listctruns=[ct.Run for ct in listcrosstalks]
channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=0,ion=1)
voies=[v for v in ["Col1","Col2","Vet1","Vet2","Gar1","Gar2"] if channels[v]!="NONE" ]

# Boucle par run et par voie
for run in listruns :
    print "Correcting CT for run",run.Name
    thectrun=RunParams.GetClosestRunFromList(run.Name,listctruns,bolodir) # recupere le run pertinent pour cts
    for voie in voies :
        filename=ampldir+"tmpwiener_"+run.Name+"_"+voie+"_"+bolo+".root"
        oldfile=TFile(filename,"READ")
        oldtree=oldfile.Get("wienerntp_"+bolo+"_"+voie)
        nentries=oldtree.GetEntries()
        # verifie si doit modifier fichier
        if not overwrite and oldtree.GetBranch("WienerAmplCTcorr") :
            oldfile.Close()
            continue
        # recupere coefs ct pertinents pour la voie
        thects=[ct for ct in listcrosstalks if ct.Run==thectrun and ct.Crosstalked==voie]
        # ouvre les fichiers necessaire en fct de ces coefs
        for cter in [ct.Crosstalker for ct in thects] :
            vars()["oldfile"+cter]=TFile(ampldir+"tmpwiener_"+run.Name+"_"+cter+"_"+bolo+".root","READ")
            vars()["oldtree"+cter]=vars()["oldfile"+cter].Get("wienerntp_"+bolo+"_"+cter)
        # fait nouveau fichier avec clone + nvlle branche
        tmpname=ampldir+"tmp_tmpwiener_"+run.Name+"_"+voie+"_"+bolo+".root"
        newfile=TFile(tmpname,"RECREATE")
        newtree=oldtree.CloneTree(0)
        ampl_corr=array('f',[0])
        newtree.Branch("WienerAmplCTcorr",ampl_corr,"WienerAmplCTcorr/F")
        for ievt in range(nentries) :
            oldtree.GetEntry(ievt)
            ampl_corr[0]=oldtree.WienerAmpl
            for ct in thects :
                vars()["oldtree"+ct.Crosstalker].GetEntry(ievt)
                ampl_corr[0] -= ct.Coef*vars()["oldtree"+ct.Crosstalker].WienerAmpl
            newtree.Fill()
        newfile.cd()
        newtree.Write('',TObject.kOverwrite)
        newfile.Close()
        # mv tmp file
        os.system("mv -f "+tmpname+" "+filename)
