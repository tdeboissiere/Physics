#!/usr/bin/env python

import sys,math,os
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams,SambaUtils
from CalibreChalFct import *
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

channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=0)
voies=[v for v in ["Chal1","Chal2"] if channels[v]!="NONE" ]
nbchals=len(voies)
bolodir=anadir+"/"+bolo
ampldir=bolodir+"/Amplitudes/"
listruns=RunParams.ReadListeRuns(bolodir)
calibcoefs=dict((voie,RunParams.GetCalibCoef(bolodir,channels[voie])) for voie in voies)
gains=dict((voie,RunParams.ReadGains(bolodir,voie=voie)) for voie in voies)
chi2coefs=dict( (voie,RunParams.GetChi2Coefs(bolodir,channels[voie])) for voie in voies)
if overwrite==0 :
    listruns=[run for run in listruns if not os.path.exists(ampldir+"eheat_"+run.Name+"_"+bolo+".root")]
 
for run in listruns :
    print "Writing Eheat tree for",run.Name
    fout=TFile(ampldir+"eheat_"+run.Name+"_"+bolo+".root","RECREATE")
    t=TTree("eheatntp_"+bolo,"Edelweiss Eheat Tree")
    fbasic=TFile(ampldir+"basicntp_"+run.Name+"_"+bolo+".root","READ")
    tbasic=fbasic.Get("basicntp_"+bolo)
    current_tinf=0
    current_tsup=0
    for v in voies :
        # pour declarer des variables a partir d'un string on utilise vars()...
        vars()["E"+v]=array('f',[0])
        vars()["Chi2"+v]=array('f',[0])
        vars()["Ldb"+v]=array('f',[0])
        vars()["Toffset"+v]=array('f',[0])
        vars()["Chi2Cut"+v]=array('i',[0])
        t.Branch("E"+v,vars()["E"+v],"E"+v+"/F")
        t.Branch("Chi2"+v,vars()["Chi2"+v],"Chi2"+v+"/F")
        t.Branch("Ldb"+v,vars()["Ldb"+v],"Ldb"+v+"/F")
        t.Branch("Toffset"+v,vars()["Toffset"+v],"Toffset"+v+"/F")
        t.Branch("Chi2Cut"+v,vars()["Chi2Cut"+v],"Chi2Cut"+v+"/I")
        vars()["file"+v]=TFile(ampldir+"tmpwiener_"+run.Name+"_"+v+"_"+bolo+".root","READ")
        vars()["tree"+v]=(vars()["file"+v]).Get("wienerntp_"+bolo+"_"+v)
        nbevts=vars()["tree"+v].GetEntries()
        # Select baselines of the corresponding run
        vars()["Baselines"+v]=RunParams.ReadBaselines(bolodir,v,run.Name)
        vars()["CurrentFWHM"+v]=0
    EChalTot=array('f',[0])
    LdbChalTot=array('f',[0])
    t.Branch("EChalTot",EChalTot,"EChalTot/F")
    t.Branch("LdbChalTot",LdbChalTot,"LdbChalTot/F")
    for i in range(nbevts) :
        tbasic.GetEntry(i)
        # Algo de mise a jour fwmh au besoin uniquement:
        if tbasic.DateSec<current_tinf or tbasic.DateSec>current_tsup :
            for v in voies :
                toto=[x for x in vars()["Baselines"+v] if x.Tinf<=tbasic.DateSec and x.Tsup>=tbasic.DateSec]
                if len(toto)!=1 : print "Not found a single fwhm for the evt:",i,toto
                current_tinf=toto[0].Tinf
                current_tsup=toto[0].Tsup
                vars()["CurrentFWHM"+v]=toto[0].FWHM
        for v in voies :
            vars()["tree"+v].GetEntry(i)
            #toto=[x.FWHM for x in vars()["Baselines"+v] if x.Tinf<=tbasic.DateSec and x.Tsup>=tbasic.DateSec]
            #(vars()["Ldb"+v])[0]=toto[0]
            (vars()["Ldb"+v])[0]=vars()["CurrentFWHM"+v]
            gain=RunParams.GetGain(gains[v],tbasic.DateSec)
            if gain!=0 and vars()["tree"+v].WienerAmpl!=0 :
		ampl_norm=vars()["tree"+v].WienerAmpl/gain
            	(vars()["E"+v])[0]=ampl_norm*ChalNLEnergy(log10(ampl_norm),calibcoefs[v])
            else :
		(vars()["E"+v])[0]=0
	    (vars()["Chi2"+v])[0]=vars()["tree"+v].WienerChi2
            (vars()["Toffset"+v])[0]=vars()["tree"+v].WienerTime
            chi2limit=chi2coefs[v][0]+chi2coefs[v][1]*vars()["tree"+v].WienerAmpl*vars()["tree"+v].WienerAmpl
            if (vars()["Chi2"+v])[0]>0.01 and (vars()["Chi2"+v])[0]<chi2limit :
                (vars()["Chi2Cut"+v])[0]=1
            else : (vars()["Chi2Cut"+v])[0]=0
        if nbchals==2 :
            invdenom=0
	    if LdbChal1[0]!=0 and LdbChal2[0]!=0 :
	        invdenom=1./(LdbChal1[0]*LdbChal1[0]+LdbChal2[0]*LdbChal2[0])
            LdbChalTot[0]=LdbChal1[0]*LdbChal2[0]*math.sqrt(invdenom)
            EChalTot[0]=invdenom*(LdbChal2[0]*LdbChal2[0]*EChal1[0]+LdbChal1[0]*LdbChal1[0]*EChal2[0])
        else :
            EChalTot[0]=EChal1[0]
            LdbChalTot[0]=LdbChal1[0]
        t.Fill()
        if floor(i/10000)*10000 == i : print i
    for v in voies :
        vars()["file"+v].Close()
    fout.cd()
    t.Write('',TObject.kOverwrite);
    fout.Close()

