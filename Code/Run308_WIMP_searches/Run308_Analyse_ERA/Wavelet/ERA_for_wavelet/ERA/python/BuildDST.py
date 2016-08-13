#!/usr/bin/env python

import sys,os
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
from ROOT import *
import RunParams
from array import array
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams["bolo"]
################################

epsilon=3 # pour la formule effet Luke
bolodir=anadir+bolo
ampldir=bolodir+"/Amplitudes/"
channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=1)
nbchals=len([v for v in ["Chal1","Chal2"] ])
nbions=len([v for v in ["Col1","Col2","Vet1","Vet2","Gar1","Gar2"] ])

runtypes=["Gamma","Bckgd","Neutron"] # pas neutrons pour l'instant
allruns=RunParams.ReadListeRuns(bolodir)
allpolars=RunParams.ReadListPolars(bolodir,check_order=1) # "vflag=i -> allpolars[i]"
i_polar=-1

for runtype in runtypes:
    runs=[x for x in allruns if x.Type==runtype]
    if len(runs)==0 : continue
    print "Building DST for",runtype,"runs"
    fout=TFile(ampldir+"dst_"+runtype+"_"+bolo+".root","RECREATE")
    tout=TTree("dst_"+bolo,"Edelweiss DST Tree")
    var_datesec=array('L',[0])
    var_tempe=array('f',[0])
    var_echal=array('f',[0])
    var_eion=array('f',[0])
    var_vfid=array('f',[0])
    var_evtflag=array('I',[0])
    var_erec=array('f',[0])
    var_q=array('f',[0])
    var_chi2flag=array('i',[0]) # NB TODO crosscheck les lettres maj/min...
    tout.Branch("DateSec",var_datesec,"DateSec/l")
    tout.Branch("Temperature",var_tempe,"Temperature/F")



    ## BRANCHES A RAJOUTER:
    # Run SambaNum IsBoloTrigger
    # E/Ldb/Chi2Ion x6 (?)
    # E/Ldb/Chi2Chal x2 (?)
    tout.Branch("EChal",var_echal,"EChal/F")
    tout.Branch("EIon",var_eion,"EIon/F")
    tout.Branch("Vfid",var_vfid,"Vfid/F")
    tout.Branch("EventFlag",var_evtflag,"EventFlag/I")
    tout.Branch("Chi2Flag",var_chi2flag,"Chi2Flag/I")
    tout.Branch("Erec",var_erec,"Erec/F")
    tout.Branch("Q",var_q,"Q/F")
    tbasic=TChain("basicntp_"+bolo,"Edelweiss Basic Data Tree")
    tion=TChain("eionntp_"+bolo,"Edelweiss Eion Tree")
    tchal=TChain("eheatntp_"+bolo,"Edelweiss Eheat Tree")
    for run in runs :
        tbasic.AddFile(ampldir+"basicntp_"+run.Name+"_"+bolo+".root")
        tion.AddFile(ampldir+"eion_"+run.Name+"_"+bolo+".root")
        tchal.AddFile(ampldir+"eheat_"+run.Name+"_"+bolo+".root")
    nbevts=tbasic.GetEntries()
    if nbevts!=tion.GetEntries() or nbevts!=tchal.GetEntries() : print "Mismatch nbevts!"
    for i in range(nbevts) :
        tbasic.GetEntry(i)
        tion.GetEntry(i)
        tchal.GetEntry(i)
        var_datesec[0]=tbasic.DateSec
        var_tempe[0]=tbasic.Temperature
        if nbchals==1 : var_echal[0]=tchal.EChal1
        else : var_echal[0]=tchal.EChalTot
        var_eion[0]=tion.Eion
        var_vfid[0]=(allpolars[tion.VFlag]).Vfid
        var_evtflag[0]=tion.EventFlag
        # Calcul effet Luke (incomplet pour l'instant...)
        Eluke=0
        if tion.EventFlag==2 :
            Eluke=allpolars[tion.VFlag].Vfid*tion.Eion/epsilon
        if tion.EventFlag>2 and nbions==2 :
            Eluke=abs(allpolars[tion.VFlag].Vcol1)*tion.ECol1*tion.IonFlags[0]+abs(allpolars[tion.VFlag].Vcol2)*tion.ECol2*tion.IonFlags[1]+abs(allpolars[tion.VFlag].Vvet1)*tion.EVet1*tion.IonFlags[2]+abs(allpolars[tion.VFlag].Vvet2)*tion.EVet2*tion.IonFlags[3]
            Eluke/=epsilon
        Erec=var_echal[0]*(1+allpolars[tion.VFlag].Vfid/epsilon)-Eluke
        Q=0
        if Erec!=0 : Q=tion.Eion/Erec
        var_erec[0]=Erec
        var_q[0]=Q
        # Flag chi2 ne marche que pour les fiduciels et seulement avec chal1 ... a ameliorer comme le reste
        var_chi2flag[0]=0
        if tion.Chi2CutCol1==1 and tion.Chi2CutCol2==1 and tchal.Chi2CutChal1==1 : var_chi2flag[0]=1
        # TODO: completer remplissage...
        tout.Fill()
        if floor(i/10000)*10000 == i : print i
    fout.cd()
    tout.Write('',TObject.kOverwrite);
    fout.Close()
