#! /usr/bin/env python

import ROOT,os,sys
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
################################

tolerance_temps=2 # X secondes de tolerance dans l'ordre en temps des evts..
ntpdir=anadir+bolo+"/Amplitudes/"
toto=os.listdir(ntpdir)
listfiles=[x for x in toto if x[0:8]=="basicntp" and x[len(x)-5:]==".root"]
for file in listfiles :
    print "Check time ordering for",file,"tolerance:",tolerance_temps,"secs"
    fbasic=TFile(ntpdir+file,"READ")
    tbasic=fbasic.Get("basicntp_"+bolo)
    date_previous=0
    for i in range(tbasic.GetEntries()) :
        tbasic.GetEntry(i)
        if tbasic.DateSec < date_previous-tolerance_temps :
            print "Pbl time ordering:",file,"entry ",i
            print tbasic.DateSec,"previous:",date_previous
        date_previous=tbasic.DateSec
    fbasic.Close()
