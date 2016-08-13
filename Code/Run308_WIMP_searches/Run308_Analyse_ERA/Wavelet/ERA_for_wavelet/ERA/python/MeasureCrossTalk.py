#! /usr/bin/env python

import ROOT
import sys,os
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
eradir=gParams["eradir"]+"/"
saveplots=int(gParams["saveplots"]) if "saveplots" in gParams else 1
################################

bolodir=anadir+bolo
ampldir=bolodir+"/Amplitudes/"
figdir=bolodir+"/Figures/"

# feuille de style root
gStyle.SetLabelFont(72,"xyz")
gStyle.SetTitleFont(72,"xyz")

# On recupere voies et runs
channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=1)
allruns=RunParams.ReadListeRuns(bolodir)
gammaruns=[run.Name for run in allruns if run.Type=="Gamma"]

if "Col1" not in channels.keys() or "Col2" not in channels.keys() or "Vet1" not in channels.keys() or "Vet2" not in channels.keys() :
    print "Some channels not in the channel list!"

list_crosstalks=[("Vet1","Col1"),("Vet2","Col2")] # EN DUR POUR L'INSTANT
filename=bolodir+"/liste_crosstalks.txt"
if (os.path.exists(filename)==0) :
    outputfile=open(filename,"w")
    outputfile.write("# Crosstalked crosstalker run coefficient\n")
    outputfile.close()

for run in gammaruns :
    t=TChain("basicntp_"+bolo,"Edelweiss Basic Data Tree")
    t.AddFile(ampldir+"basicntp_"+run+"_"+bolo+".root")
    nbevts=t.GetEntries()
    tc1=TChain("wienerntp_"+bolo+"_Col1","Edelweiss Wiener Data Tree")
    tc2=TChain("wienerntp_"+bolo+"_Col2","Edelweiss Wiener Data Tree")
    tc1.AddFile(ampldir+"tmpwiener_"+run+"_Col1_"+bolo+".root")
    tc2.AddFile(ampldir+"tmpwiener_"+run+"_Col2_"+bolo+".root")
    tv1=TChain("wienerntp_"+bolo+"_Vet1","Edelweiss Wiener Data Tree")
    tv2=TChain("wienerntp_"+bolo+"_Vet2","Edelweiss Wiener Data Tree")
    tv1.AddFile(ampldir+"tmpwiener_"+run+"_Vet1_"+bolo+".root")
    tv2.AddFile(ampldir+"tmpwiener_"+run+"_Vet2_"+bolo+".root")
    if tc1.GetEntries()!=nbevts or tc2.GetEntries()!=nbevts or tv1.GetEntries()!=nbevts or tv2.GetEntries()!=nbevts : print "Wrong nb of evts"
    t.AddFriend(tc1)
    t.AddFriend(tc2)
    t.AddFriend(tv1)
    t.AddFriend(tv2)
    print "Run:",run,"    evts:",nbevts

    for x in list_crosstalks:
        crosstalked=x[0]
        crosstalker=x[1]
        cutampl_ctalked="(wienerntp_"+bolo+"_"+crosstalked+".WienerAmpl>-500 && wienerntp_"+bolo+"_"+crosstalked+".WienerAmpl<50)"
        cutampl_ctalker="(wienerntp_"+bolo+"_"+crosstalker+".WienerAmpl>-100 && wienerntp_"+bolo+"_"+crosstalker+".WienerAmpl<2000)"
        fit = TF1("fit","[0]+[1]*x",-100,3000)
        hg = TH2F("hg","Un exemple de fit",1000,-100,2000,2000,-500,100)
        cc=TCanvas("cc","cc")
        cc.Divide(2,1)
        cc.cd(1)
        t.Draw("wienerntp_"+bolo+"_"+crosstalked+".WienerAmpl:wienerntp_"+bolo+"_"+crosstalker+".WienerAmpl >> hg",cutampl_ctalked+" && "+cutampl_ctalker)
        fit.FixParameter(0,0)
        fit.SetLineColor(kBlue)
        hg.Fit(fit,"","",100,1500)
        pente=fit.GetParameter(1)
        offset=fit.GetParameter(0)
        print "Crosstalk:",pente*100,"%"
        line_ct=TLine(0,0,1500,1500*pente)
        line_ct.SetLineColor(kRed)
        line_ct.SetLineWidth(2)
        line_ct.Draw("same")
        gPad.Update()
        foo=raw_input("After moving the line, press enter..")
        pente=line_ct.GetY2()/line_ct.GetX2()
        print "Crosstalk:",pente*100,"%"
        cutfct = TF1("cutfct","[0]",0,2000)
        cutfct.SetParameter(0,offset)
        cutfct.SetLineStyle(2)
        cutfct.Draw("same")
        cc.cd(2)
        t.Draw("wienerntp_"+bolo+"_"+crosstalked+".WienerAmpl-("+ str(pente)+"*wienerntp_"+bolo+"_"+crosstalker+".WienerAmpl)-"+str(offset)+ ":wienerntp_"+bolo+"_"+crosstalker+".WienerAmpl",cutampl_ctalked+" && "+cutampl_ctalker)
        if saveplots==1 : cc.SaveAs(figdir+"crosstalk_"+crosstalked+"_"+run+".gif")
        foo=raw_input("Press enter..")
        outputfile=open(filename,"a")
        outputfile.write(crosstalked+" "+crosstalker+" "+run+" {0:.3f}\n".format(pente))
        outputfile.close()
