#! /usr/bin/env python

# TODO:
# ameliorer un peu la procedure par exemple
# - calibrer d'abord les trous fiduciels avant les electrons
# - faire cut fiduciel grossier avant de calibrer les fiduciels

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
print "Gamma runs for calib:",gammaruns

if "Col1" not in channels.keys() or "Col2" not in channels.keys() :
    print "Collectrodes not in the channel list!"

t=TChain("basicntp_"+bolo,"Edelweiss Basic Data Tree")
for run in gammaruns :
    t.AddFile(ampldir+"basicntp_"+run+"_"+bolo+".root")
nbevts=t.GetEntries()
tc1=TChain("wienerntp_"+bolo+"_Col1","Edelweiss Wiener Data Tree")
tc2=TChain("wienerntp_"+bolo+"_Col2","Edelweiss Wiener Data Tree")
for run in gammaruns :
    tc1.AddFile(ampldir+"tmpwiener_"+run+"_Col1_"+bolo+".root")
    tc2.AddFile(ampldir+"tmpwiener_"+run+"_Col2_"+bolo+".root")
t.AddFriend(tc1)
t.AddFriend(tc2)
if tc1.GetEntries()!=nbevts or tc2.GetEntries()!=nbevts : print "Wrong nb of evts in collectrodes!"
for voie_sup in ["Vet1","Vet2","Gar1","Gar2"] :
    if voie_sup in channels.keys() :
        tsup=TChain("wienerntp_"+bolo+"_"+voie_sup,"Edelweiss Wiener Data Tree")
        for run in gammaruns : tsup.AddFile(ampldir+"tmpwiener_"+run+"_"+voie_sup+"_"+bolo+".root")
        t.AddFriend(tsup)
        if tsup.GetEntries()!=nbevts : print "Wrong nb of evts for",voie_sup

calibcoefs=dict()

## Calibration primaire : collectrodes
for voie in ["Col1","Col2"] :
    print "**** Calibrating",voie,"****"
    
    fit = TF1("fit","[0]+[1]*x+[2]*x*x",10,1000)
    hg = TH2F("hg","Chi2 distribution",1000,10,1000,1000,0,10)
    fit.FixParameter(0,1)
    fit.FixParameter(1,0)
    t.Draw("wienerntp_"+bolo+"_"+voie+".WienerChi2:wienerntp_"+bolo+"_"+voie+".WienerAmpl >> hg","wienerntp_"+bolo+"_"+voie+".WienerChi2<10")
    hg.Fit("fit","L")	
    cutfct = TF1("cutfct","[0]+[1]*x+[2]*x*x",0,100000)
    cutfct.SetParameter(0,2)
    cutfct.SetParameter(1,0)
    cutfct.SetParameter(2,4*fit.GetParameter(2))
    cutfct.SetLineStyle(2)
    cutfct.Draw("same")
    gPad.Update()
    foo=raw_input("Press enter to continue")	
    cut_ChiFit= "wienerntp_"+bolo+"_"+voie+".WienerChi2 < wienerntp_"+bolo+"_"+voie+".WienerAmpl*wienerntp_"+bolo+"_"+voie+".WienerAmpl*"+str(cutfct.GetParameter(2))+"+"+str(cutfct.GetParameter(0))

    h=TH1F("spectre","spectre",100,30,3000)
    h.SetLineWidth(2)
    # On fait le spectre corrige du crosstalk (NB le cut chi2 est fait sur les ampls non corrigees!)
    t.Draw("wienerntp_"+bolo+"_"+voie+".WienerAmplCTcorr >> spectre",cut_ChiFit)
    #cutcal=TCut("wienerntp_"+bolo+"_"+voie+".WienerChi2 < 200") # A AMELIORER!
    #t.Draw("wienerntp_"+bolo+"_"+voie+".WienerAmplCTcorr",cutcal,"same")
    linecal=TLine(500,0,500,h.GetMaximum())
    linecal.SetLineWidth(2)
    linecal.SetLineColor(kRed)
    linecal.Draw("same")
    gPad.Update()
    foo=raw_input("After moving the line, press enter..")
    calibcoefs[voie]=linecal.GetX1()/356.
    if saveplots==1 : gPad.SaveAs(figdir+"calib_"+voie+".gif")
    gPad.Close()

print "**** Cross-check col1 vs col2 ****"
t.SetMarkerStyle(7)
cutrange=TCut("wienerntp_"+bolo+"_Col1.WienerAmplCTcorr/"+str(calibcoefs["Col1"])+">-200 && wienerntp_"+bolo+"_Col1.WienerAmplCTcorr/"+str(calibcoefs["Col1"])+"<3000 && wienerntp_"+bolo+"_Col2.WienerAmplCTcorr/"+str(calibcoefs["Col2"])+">-200 && wienerntp_"+bolo+"_Col2.WienerAmplCTcorr/"+str(calibcoefs["Col2"])+"<3000")
t.Draw("wienerntp_"+bolo+"_Col1.WienerAmplCTcorr/"+str(calibcoefs["Col1"])+":wienerntp_"+bolo+"_Col2.WienerAmplCTcorr/"+str(calibcoefs["Col2"]),cutrange)
linediag=TLine(0,0,3000,3000)
linediag.SetLineWidth(2)
linediag.SetLineColor(kRed)
linediag.Draw("same")
if saveplots==1 : gPad.SaveAs(figdir+"calib_collectrodes.gif")

foo=raw_input("Press enter to continue")

## Calibration veto a partir des collectrodes
calibsec=dict()
if "Vet1" in channels.keys() : calibsec["Vet1"]="Col1"
if "Vet2" in channels.keys() : calibsec["Vet2"]="Col2"
if "Gar1" in channels.keys() : calibsec["Gar1"]="Col1"
if "Gar2" in channels.keys() : calibsec["Gar2"]="Col2"
for voie,calibrator in calibsec.items() :
    print "**** Calibrating",voie,"****"
    cutrange=TCut("wienerntp_"+bolo+"_"+voie+".WienerAmplCTcorr > -1000 && wienerntp_"+bolo+"_"+voie+".WienerAmplCTcorr < 10000 && wienerntp_"+bolo+"_"+calibrator+".WienerAmplCTcorr/"+str(calibcoefs[calibrator])+">-100 && wienerntp_"+bolo+"_"+calibrator+".WienerAmplCTcorr/"+str(calibcoefs[calibrator])+"<1000")
    t.Draw("wienerntp_"+bolo+"_"+voie+".WienerAmplCTcorr:wienerntp_"+bolo+"_"+calibrator+".WienerAmplCTcorr/"+str(calibcoefs[calibrator]),cutrange)
    linecal=TLine(0,0,500,500)
    if voie[0:3]=="Gar" :
        linecal=TLine(0,300,356,0)
    linecal.SetLineColor(kRed)
    linecal.Draw("same")
    gPad.Update()
    foo=raw_input("After moving the line, press enter..")
    calibcoefs[voie]=(linecal.GetY2()-linecal.GetY1())/(linecal.GetX2()-linecal.GetX1())
    if voie[0:3]=="Gar" :
        calibcoefs[voie]=linecal.GetY1()
    if saveplots==1 : gPad.SaveAs(figdir+"calib_"+voie+".gif")
    gPad.Close()

# Toto: Check gar1 vs gar2

filename=bolodir+"/liste_calibs.txt"
if (os.path.exists(filename)==0) :
    outputfile=open(filename,"w")
    outputfile.write("# Voie type tinf tsup coefficients\n")
    outputfile.close()
outputfile=open(filename,"a")
for voie,coef in calibcoefs.items() :
    print voie,":",coef
    outputfile.write('"'+channels[voie]+'" StandardIon 0 0 {0:.3f} \n'.format(coef))    
outputfile.close()

