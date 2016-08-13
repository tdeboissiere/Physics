#! /usr/bin/env python

from GlobalParams import ReadGlobalParams
import sys,os

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
wait=1
saveplots=int(gParams["saveplots"]) if "saveplots" in gParams else 1
run=gParams['run'] if "run" in gParams else "All"
################################

bolodir=anadir+"/"+bolo
ampldir=bolodir+"/Amplitudes/"
figdir=bolodir+"/Figures/"

channels=RunParams.ReadBoloChannels(bolo,anadir,dico=1,remove_none=1)
list_runs=RunParams.ReadListeRuns(bolodir)
#if run!="All" : # ON VEUT PAS ENCORE METTRE CA EN SERVICE
#    list_runs=[r for r in list_runs if r.Name==run]
t=TChain("wienerntp_"+bolo+"_"+voie,"Edelweiss Wiener Data Tree")
for therun in list_runs :
    t.AddFile(ampldir+"tmpwiener_"+therun.Name+"_"+voie+"_"+bolo+".root")

	
## Calibration primaire : collectrodes
# premier fit
c1=TCanvas()
fit = TF1("fit","[0]+[1]*x+[2]*x*x",1,1000)
hg = TH2F("hg","Chi2 distribution",1000,1,1000,1000,0,10)
fit.FixParameter(0,1)
fit.FixParameter(1,0)  
t.Draw("WienerChi2:WienerAmpl >> hg","WienerChi2<10")
hg.Fit("fit","L")	
gPad.Update()
foo=raw_input("Press enter to continue")

# Fit gaussien sur le Chi2
c2=TCanvas()
GGfit = TH1F("GGfit","histo Chi2 distri",100,0.7,1.8)
#GGfit = TH1F("GGfit","histo Chi2 distri",1000,0.9,1.05)
t.Draw("WienerChi2 >> GGfit","WienerAmpl<100")
GGfit.Fit("gaus","","E1",0.7,1.3)
coupure=gaus.GetParameter(1)+ 2.56*gaus.GetParameter(2)
#foo=raw_input("Press enter to continue")
gPad.Update()
print "Chi2 cut at low energy:",coupure
coupure=float(raw_input("Enter value by hand\n"))
print "coupure:",coupure
if saveplots==1 : c2.SaveAs(figdir+"chi2distri_b_"+voie+".gif")

#99% d'acceptance -> 2.56sigma
Gcutfct = TF1("Gcutfct","[0]+[1]*x+[2]*x*x",0,1000)
Gcutfct.SetParameter(0,coupure)
Gcutfct.SetParameter(1,0)
Gcutfct.SetParameter(2,3.5*fit.GetParameter(2))
Gcutfct.SetLineStyle(2)
#t.Draw("WienerChi2:WenerAmpl >> GGfit","WienerAmpl>0")
c1.cd()
Gcutfct.Draw("same")
if saveplots==1 : c1.SaveAs(figdir+"chi2distri_a_"+voie+".gif")
foo=raw_input("Press enter to continue")

c3=TCanvas()
cut_ChiFit="WienerChi2<WienerAmpl*WienerAmpl*"+str(Gcutfct.GetParameter(2))+"+"+str(Gcutfct.GetParameter(0))
h=TH1F("spectre","spectre",200,30,2000)
h.SetLineWidth(2)
t.Draw("WienerAmpl >> spectre",cut_ChiFit)
gPad.Update()
foo=raw_input("Press enter to continue")

param0=Gcutfct.GetParameter(0)
param2=Gcutfct.GetParameter(2)
file=bolodir+"/liste_coefschi2.txt"
if (os.path.exists(file)==0) :
	outputfile=open(file,"w")
	outputfile.write("# Voie coefficients\n")
	outputfile.close()
outputfile=open(bolodir+"/liste_coefschi2.txt","a")
outputfile.write('"'+channels[voie]+'" {0:.3f} {1:.7f}\n'.format(param0,param2))
outputfile.close()

        
