#!/usr/bin/env python

from GlobalParams import ReadGlobalParams
import math
import numpy as np
import os,re,sys
import RunParams
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
overwrite=1
timestep=3600 # en seconds
if "overwrite" in gParams : overwrite=int(gParams["overwrite"])
if "timestep" in gParams : timestep=int(gParams["timestep"])
################################

bolodir=anadir+bolo
ampldir=bolodir+"/Amplitudes/"
runlist=RunParams.ReadListeRuns(bolodir)

# Selection des runs a lire (au besoin) et creation du fichier
filename=bolodir+"/liste_periods.txt"
if overwrite==0 and os.path.exists(filename) :
    existing_periods=RunParams.ReadListePeriods(bolodir)
    for p in existing_periods:
        runlist=[run for run in runlist if run.Name!=p.Run]

if os.path.exists(filename)==0 or overwrite==1 :
    outputfile=open(filename,"w")
    outputfile.write("# Tinf Tsup Run\n")
    outputfile.close()
outputfile=open(filename,"a")

# Boucle sur les runs
for therun in runlist :
    ntpfile=ampldir+"basicntp_"+therun.Name+"_"+bolo+".root"
    if os.path.exists(ntpfile)!=True :
        print "Error, no ntpfile for run ",therun.Name
    f=TFile(ntpfile,"READ")
    t=f.Get("basicntp_"+bolo)
    tinf,tsup=0,0 # Algo pour recuperer tinf,tsup meme si les evts sont pas ordonnes en temps...
    for i in range(t.GetEntries()) :
	t.GetEntry(i)
	thetime=t.DateSec
	if thetime < tinf or tinf==0 : tinf=thetime
	if thetime > tsup : tsup=thetime
    nbdivisions=float(tsup-tinf)/timestep
    if nbdivisions<1 : nbdivisions=1
    if nbdivisions-math.floor(nbdivisions) < 0.5 : nbdivisions=math.floor(nbdivisions)
    else : nbdivisions=math.ceil(nbdivisions)
    localperiodlength=float(tsup-tinf)/nbdivisions
    print therun.Name,"("+str(t.GetEntries())+" evts) :",nbdivisions," period[s] of length",localperiodlength,"secs"
    array_tinf=np.arange(nbdivisions)*localperiodlength+tinf
    array_tsup=(np.arange(nbdivisions)+1)*localperiodlength+tinf
    for p in range(0,int(nbdivisions)) :
        outputfile.write(str(long(array_tinf[p]))+" "+str(long(array_tsup[p]))+" "+therun.Name+"\n")

outputfile.close()
