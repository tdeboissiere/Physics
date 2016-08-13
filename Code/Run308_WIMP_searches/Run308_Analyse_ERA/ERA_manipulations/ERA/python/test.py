from numpy import *
import matplotlib.pyplot as plt
plt.ion()
import os
from ROOT import *
ROOT.gSystem.Load("/Users/armengau/Edelweiss/ERA/"+'lib/EraLib.so')
from ROOT import EdwPulse,EdwEvent,NoiseSpectrum,FitPulse,EdwTemplate
from RunParams import *

anadir="/Users/armengau/Edelweiss/Data_Edw/TestERANew/RunC/"
bolo="FID806"
voie="chalA FID806"

##########################
## Tests sur fit
##########################
# trace
pulsedir=anadir+bolo+"/Traces/"
list_pulsefiles=os.listdir(pulsedir)
pulsefilename=pulsedir+list_pulsefiles[0]
pulsefile=TFile(pulsefilename,"READ")
gTree=pulsefile.Get("EdwTree")
gEvt=EdwEvent()
gTree.SetBranchAddress("Event",gEvt)
gTree.GetEntry(9)
p=gEvt.Pulse(voie)
pp=FitPulse(p)
plt.plot(p.Trace())
#noise
noisedir=anadir+bolo+"/Spectra/"
noisefilename=noisedir+"spectra_mg19e000_FID806.root"
#noisefilename=noisedir+"spectra_la17e002_ID404.root"
noise=GetNoiseSpectrum(noisefilename,voie,gEvt.DateSec())
plt.yscale('log')
plt.plot(noise)
# combo = le fit wiener !..
template=EdwTemplate()
tstr=GetTemplate(anadir+bolo,voie)
template.Initialize(tstr[0],tstr[1],tstr[2],tstr[3],tstr[4],p.TraceLength(),p.Pretrigger(),pp.FitTimeMin(),pp.FitTimeMax(),noise)
result=pp.WienerLoop(template,noise)
sortie=template.ComputeTrace(result[1]/pp.ModulationLength(),1)
sortie=[ result[0]*sortie[k]*pp.Sign()*template.Sign() for k in range(p.TraceLength())]
plt.plot(arange(pp.TraceLength()),sortie,color="red")
plt.plot(arange(pp.TraceLength()),pp.ProcessedTrace(),color="blue")

#tmpltnoise=[x*x*result[0]*result[0] for x in template.GetNonIntegerOffsetFFT(result[1]/pp.ModulationLength())]
diff=[pp.Sign()*result[0]*(template.GetNonIntegerOffsetFFT(result[1]/pp.ModulationLength()))[i]-(pp.ProcessedTraceFFT())[i] for i in range(pp.TraceLength())]
diff2=[d*d for d in diff]
ratio=[diff2[i]/noise[i] for i in range(pp.TraceLength())]
plt.yscale('log')
plt.plot(arange(pp.TraceLength()),noise)
plt.plot(arange(pp.TraceLength()),diff2)
plt.plot(ratio) # ratio centre sur 1.

##########################
## Tests sur histos/scatterplots de calib et crosstalks..
## et fcts de calibrations/calcul FWHM/chi2 distris, etc...
##########################

ampldir=anadir+bolo+"/Amplitudes/"
run="la17e002"

fb=TFile(ampldir+"basicntp_"+run+"_"+bolo+".root")
t=fb.Get("basicntp_"+bolo)
fc1=TFile(ampldir+"tmpwiener_"+run+"_Col1_"+bolo+".root")
tc1=fc1.Get("wienerntp_"+bolo+"_Col1")
fc2=TFile(ampldir+"tmpwiener_"+run+"_Col2_"+bolo+".root")
tc2=fc2.Get("wienerntp_"+bolo+"_Col2")
fv1=TFile(ampldir+"tmpwiener_"+run+"_Vet1_"+bolo+".root")
tv1=fv1.Get("wienerntp_"+bolo+"_Vet1")
fv2=TFile(ampldir+"tmpwiener_"+run+"_Vet2_"+bolo+".root")
tv2=fv2.Get("wienerntp_"+bolo+"_Vet2")
t.AddFriend(tc1)
t.AddFriend(tv1)
t.AddFriend(tc2)
t.AddFriend(tv2)

print "**** Calibrating Col1 ****"
t.Draw("wienerntp_"+bolo+"_Col1.WienerAmpl")
lcal=TLine(800,0,800,500)

lcal.Draw("same")
foo=raw_input("After moving the line, press enter..")

print "The 356keV line is set at",lcal.GetX1(),"ADUs"

