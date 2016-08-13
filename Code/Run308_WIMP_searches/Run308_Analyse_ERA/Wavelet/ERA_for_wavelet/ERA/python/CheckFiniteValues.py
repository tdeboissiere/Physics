#! /usr/bin/env python

## Todo = enrichir la liste des variables qu'on teste

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

# Fonction intermediaire: liste nans des variables donnees d'un arbre
def CheckFinite(tree,nom_vars) :
    nbnans=[0]*len(nom_vars)
    for i in range(tree.GetEntries()) :
        tree.GetEntry(i)
        for j in range(len(nom_vars)) :
            x=eval("tree."+nom_vars[j])
            if TMath.Finite(x)==0 : nbnans[j]+=1
    return nbnans

# Le "main"
ntpdir=anadir+bolo+"/Amplitudes/"
toto=os.listdir(ntpdir)
listfiles=[x for x in toto if x[0]!="." and x[len(x)-5:]==".root"] # virer les trucs genre .DS_store
for file in listfiles :
    # on definit la liste des variables a tester a la main
    # on teste pas les strings ni les tableaux
    treename=""
    listvars=[]
    if file[:8]=="basicntp" :
        treename="basicntp_"+bolo
        listvars=["SambaNum","TriggerBit1","TriggerBit2","IsBoloTrigger","DateSec","TimeStamp","GigaStamp","SimpleAmplChal1","SimpleAmplChal2","SimpleAmplVet1","SimpleAmplVet2","SimpleAmplCol1","SimpleAmplCol2","SimpleAmplColTot"]
    elif file[:4]=="eion" :
        treename="eionntp_"+bolo
        listvars=["ECol1","Chi2Col1","LdbCol1","ToffsetCol1","ECol2","Chi2Col2","LdbCol2","ToffsetCol2","EVet1","Chi2Vet1","LdbVet1","ToffsetVet1","EVet2","Chi2Vet2","LdbVet2","ToffsetVet2","Eion","Tion","EventFlag","VFlag"]
    elif file[:5]=="eheat" :
        treename="eheatntp_"+bolo
        listvars=["EChal1","Chi2Chal1","LdbChal1","ToffsetChal1","EChal2","Chi2Chal2","LdbChal2","ToffsetChal2","EChalTot","LdbChalTot"]
    elif file[:3]=="dst" :
        treename="dst_"+bolo
        listvars=["DateSec","EChal","EIon","Vfid","EventFlag","Erec","Q"]
    else :
        print "Fichier ntuple de type inconnu:",file
        continue
    f=TFile(ntpdir+file,"READ")
    t=f.Get(treename)
    vars_ok=[x for x in listvars if t.GetLeaf(x)]
    nbnans=CheckFinite(t,vars_ok)
    for i in range(len(vars_ok)) :
        if nbnans[i]!=0 : print nbnans[i],"non-finite values for",vars_ok[i],"in",file
    f.Close()

