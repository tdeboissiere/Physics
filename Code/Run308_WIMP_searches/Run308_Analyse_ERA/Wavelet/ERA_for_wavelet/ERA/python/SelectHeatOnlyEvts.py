#!/usr/bin/env python

import sys,os
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
################################

bolodir=anadir+bolo
ampldir=bolodir+"/Amplitudes/"
allruns=RunParams.ReadListeRuns(bolodir)
outputfile=open(bolodir+"/liste_evts.txt","w")
outputfile.write("# Run SambaNum \n")

tbasic=TChain("basicntp_"+bolo,"Edelweiss Basic Data Tree")
tion=TChain("eionntp_"+bolo,"Edelweiss Eion Tree")
tchal=TChain("eheatntp_"+bolo,"Edelweiss Eheat Tree")
# On prend TOUS les runs..
for run in allruns :
    tbasic.AddFile(ampldir+"basicntp_"+run.Name+"_"+bolo+".root")
    tion.AddFile(ampldir+"eion_"+run.Name+"_"+bolo+".root")
    tchal.AddFile(ampldir+"eheat_"+run.Name+"_"+bolo+".root")
nbevts=tbasic.GetEntries()
if nbevts!=tion.GetEntries() or nbevts!=tchal.GetEntries() : print "Mismatch nbevts!"
for i in range(nbevts) :
    tbasic.GetEntry(i)
    tion.GetEntry(i)
    tchal.GetEntry(i)
    # Critere (TRES SIMPLE) de selection des "chaleur seule" a ce stade:
    # Forcement tres tres imparfait...
    if tion.Eion==0 and tchal.EChalTot>10 :
        outputfile.write(str(tbasic.Run)+" "+str(tbasic.SambaNum)+"\n")
outputfile.close()
