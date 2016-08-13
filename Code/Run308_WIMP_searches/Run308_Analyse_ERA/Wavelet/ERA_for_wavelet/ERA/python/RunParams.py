
# Utilitaires python pour lecture des petit fichiers ascii
# Essai de structure a partir de l'outil (un peu trop "advanced" peut-etre..)
# namedtuple, suggere pour mimiquer une structure C

from GlobalParams import ReadGlobalParams
import sys,os
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
eradir=gParams['eradir']
import ROOT
if "ENVIRONMENT" in os.environ :
    if os.environ["ENVIRONMENT"].find("BATCH")!=-1 : ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load(eradir+'/lib/EraLib.so')
from ROOT import vector,TFile,TTree
from ROOT import NoiseSpectrum,EdwTemplate
from collections import namedtuple
RunStruct=namedtuple('RunStruct','Name Type ModChal1 ModChal2 Vcol1 Vcol2 Vvet1 Vvet2 Vgar1 Vgar2')
# quasi meme syntaxe que mon code C
BaselineStruct=namedtuple('BaselineStruct','Tinf Tsup FWHM FWHM_from_spectrum')
PeriodStruct=namedtuple('PeriodStruct','Tinf Tsup Run')
GainStruct=namedtuple('PeriodStruct','Tinf Tsup Gain')
BoloStruct=namedtuple('BoloStruct','Name Mac BB1 BB2 Chal1 Chal2 Col1 Col2 Vet1 Vet2 Gar1 Gar2')
PolarStruct=namedtuple('PolarStruct','Vflag Vfid Vcol1 Vcol2 Vvet1 Vvet2 Vgar1 Vgar2')
CrossTalkStruct=namedtuple('CrossTalkStruct','Crosstalked Crosstalker Run Coef')

def ReadListeBolos(anadir) :
    bololist=[]
    file=anadir+"/liste_bolos.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split('"')
        head=[x.strip() for x in toto[0].split(' ') if x!=''] 
        voies=[x.strip() for x in toto[1:] if x.strip()!='']
        if len(head)==3 and len(voies)==8 : # Old file format, only one BB
            thestr=BoloStruct(head[0],head[1],head[2],"NONE",voies[0],voies[1],voies[2],voies[3],voies[4],voies[5],voies[6],voies[7])
        elif len(head)==4 and len(voies)==8 :
            thestr=BoloStruct(head[0],head[1],head[2],head[3],voies[0],voies[1],voies[2],voies[3],voies[4],voies[5],voies[6],voies[7])
        else : print "Wrong nb of variables in ReadListeBolos"
        bololist.append(thestr)
    f.close()
    return bololist

###### Old version :
#def ReadListeBolos(anadir,with_mac_bb=0):
#    # option with_mac_bb => retourne tableau des [bolo,mac,bb]
#    bololist=[]
#    file=anadir+"/liste_bolos.txt"
#    f=open(file,"r")
#    lines=f.readlines()
#    lines=lines[1:]
#    for theline in lines :
#        toto=theline.split(' ')
#        if not with_mac_bb :
#            bololist.append(toto[0])
#        else :
#            bololist.append(toto[0:3])
#    f.close()
#    return bololist

def ReadBoloChannels(bolo,anadir,remove_none=0,dico=0,ion=0) :
    bololist=ReadListeBolos(anadir)
    bolonamelist=[bb.Name for bb in bololist]
    index=bolonamelist.index(bolo)
    bb=bololist[index]
    channel_list=[bb.Chal1,bb.Chal2,bb.Col1,bb.Col2,bb.Vet1,bb.Vet2,bb.Gar1,bb.Gar2]
    channel_dict={"Chal1":bb.Chal1,"Chal2":bb.Chal2,"Col1":bb.Col1,"Col2":bb.Col2, "Vet1":bb.Vet1,"Vet2":bb.Vet2,"Gar1":bb.Gar1,"Gar2":bb.Gar2}
    if ion==1 :
        channel_list=[bb.Col1,bb.Col2,bb.Vet1,bb.Vet2,bb.Gar1,bb.Gar2]
        channel_dict={"Col1":bb.Col1,"Col2":bb.Col2,"Vet1":bb.Vet1,"Vet2":bb.Vet2,"Gar1":bb.Gar1,"Gar2":bb.Gar2}
    if (remove_none==1) :
        channel_list=[x for x in channel_list if x!="NONE"]
        channel_dict=dict((it,val) for it,val in channel_dict.items() if val!="NONE")
    if dico : channel_list=channel_dict
    return channel_list

def ReadBaselines(bolodir,voie,run="NONE") :
    f=open(bolodir+"/baseline_"+voie+".txt","r")
    lines=f.readlines()
    lines=lines[1:]
    baselinelist=[]
    for theline in lines :
        toto=theline.split(' ')
        if len(toto)==3 :
            thestr=BaselineStruct(int(toto[0]),int(toto[1]),float(toto[2]),0.0)
        elif len(toto)==4 :
            thestr=BaselineStruct(int(toto[0]),int(toto[1]),float(toto[2]),float(toto[3]))
        else : print "Wrong nb of variables in ReadBaseline"
        baselinelist.append(thestr)
    f.close()
    if run!="NONE" :
        periodlist=ReadListePeriods(bolodir)
        if len(periodlist)!=len(baselinelist) : print "Mismatch period/baseline list."
        baselinelist=[b for i,b in enumerate(baselinelist) if periodlist[i].Run==run]
    return baselinelist

def ReadListePeriods(bolodir) :
    periodlist=[]
    file=bolodir+"/liste_periods.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split(' ')
        if len(toto)!=3 : print "Wrong nb of variables in ReadPeriods"
        thestr=PeriodStruct(int(toto[0]),int(toto[1]),toto[2].rstrip("\n"))
        periodlist.append(thestr)
    f.close()
    return periodlist

def ReadListCrossTalks(bolodir) :
    ctlist=[]
    file=bolodir+"/liste_crosstalks.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split(' ')
        if len(toto)!=4 : print "Wrong nb of variables in ReadCrossTalks"
        thestr=CrossTalkStruct(toto[0],toto[1],toto[2],float(toto[3]))
        ctlist.append(thestr)
    f.close()
    return ctlist

def ReadGains(bolodir,voie="Chal1") :
    gainlist=[]
    file=bolodir+"/gain_"+voie+".txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split(' ')
        if len(toto)!=3 : print "Wrong nb of variables in ReadGains"
        thestr=GainStruct(int(toto[0]),int(toto[1]),float(toto[2]))
        gainlist.append(thestr)
    f.close()
    return gainlist

def GetGain(vectgain,datesec) :
    gg=[g.Gain for g in vectgain if datesec>=g.Tinf and datesec<=g.Tsup]
    # il peut en trouver 2 si c'est entre 2 periodes.. on prend alors la 1ere
    if len(gg) not in [1,2] : print "Pbl finding gain for",datesec,"Nb gains candidates:",len(gg)
    return gg[0]

def ReadListeRuns(bolodir) :
    runlist=[]
    file=bolodir+"/liste_runs.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split(' ')
        if len(toto)!=10 : print "Wrong nb of variables for run",toto[0]
        thestr=RunStruct(toto[0],toto[1],float(toto[2]),float(toto[3]),float(toto[4]),float(toto[5]),float(toto[6]),float(toto[7]),float(toto[8]),float(toto[9]))
        runlist.append(thestr)
    f.close()
    return runlist

def ReadListPolars(bolodir,check_order=0) :
    polarlist=[]
    file=bolodir+"/liste_polars.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split(' ')
        if len(toto)!=8 : print "Wrong nb of variables for vflag",toto[0]
        thestr=PolarStruct(int(toto[0]),float(toto[1]),float(toto[2]),float(toto[3]),float(toto[4]),float(toto[5]),float(toto[6]),float(toto[7]))
        polarlist.append(thestr)
    f.close()
    if (check_order) :
        for i in range(len(polarlist)) :
            if polarlist[i].Vflag!=i :
                print "Vflags NOT in order in polarlist..."
    return polarlist

def GetPolarFlag(run,polarlist):
    flag=-1
    for polar in polarlist:
        if run.Vcol1==polar.Vcol1 and run.Vcol2==polar.Vcol2 and run.Vvet1==polar.Vvet1 and run.Vvet2==polar.Vvet2 and run.Vgar1==polar.Vgar1 and run.Vgar2==polar.Vgar2 : flag=polar.Vflag
    if flag==-1 : print "No polar flag found for run",run.Name
    return flag

def GetNoiseSpectrum(noisefilename,channel,datesec) :
    noisefile=TFile(noisefilename,"READ")
    noisetree=noisefile.Get("SpectrumTree")
    noise=NoiseSpectrum()
    noisetree.SetBranchAddress("Spectrum",noise)
    for i in range(noisetree.GetEntries()) :
        noisetree.GetEntry(i)
        if noise.Channel()==channel and noise.StartValidity()<=datesec and noise.EndValidity()>=datesec :
            break
    else :
        print "No noise spectrum found for {0} at t={1}".format(channel,datesec)
    return noise

def GetTemplate(bolodir,voie,type="Standard",date=0) :
    # Retourne la structure "EdwTemplate-compatible" pour l'instant => utilisation de vector stl.
    tmpltstr=[]
    file=bolodir+"/liste_templates.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split('"')
        thevoie=toto[1]
        toto=(toto[2]).split(' ')
        thetype=toto[1]
        thetinf=int(toto[2])
        thetsup=int(toto[3])
        thecoefs=[float(x) for x in toto[4:] ]
        vectcoefs=vector('float')(len(thecoefs),0)
        for i in range(len(thecoefs)) : vectcoefs[i]=thecoefs[i]
        if thevoie==voie and (type=="Standard" or type==thetype) and (date==0 or (date<=thetinf and date>=thetsup)) :
            tmpltstr=[thevoie,thetype,thetinf,thetsup,vectcoefs]
            break
    else :
        print "No template found.."
    return tmpltstr

def GetCalibCoef(bolodir,voie,type="Standard",date=0) :
    calibcoefs=[]
    file=bolodir+"/liste_calibs.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split('"')
        thevoie=toto[1]
        toto=(toto[2].rstrip("\n").strip()).split(' ')
        thetype=toto[0]
        thetinf=int(toto[1])
        thetsup=int(toto[2])
        if thevoie==voie and (type=="Standard" or type==thetype) and (date==0 or (date<=thetinf and date>=thetsup)) :
            calibcoefs=[float(x) for x in toto[3:] ]
            break
    else :
        print "No calib found.."
    return calibcoefs

def GetChi2Coefs(bolodir,voie) :
    chi2coefs=[]
    file=bolodir+"/liste_coefschi2.txt"
    f=open(file,"r")
    lines=f.readlines()
    lines=lines[1:]
    for theline in lines :
        toto=theline.split('"')
        thevoie=toto[1]
        toto=(toto[2].rstrip("\n").strip()).split(' ')
        if len(toto)!=2 : print "Pbl format chi2coef file? ",theline
        if thevoie==voie :
            chi2coefs=[float(x) for x in toto ]
            break
    else :
        print "No chi2 coef found.."
    return chi2coefs

def GetClosestRunFromList(run,listruns,bolodir) :
    # Pour le run "run" trouve le run dans "runlist"
    # qui soit le plus proche en temps, et qui le precede
    # (eg. la plus proche calib gamma faite avant un run de fond)
    listperiods=ReadListePeriods(bolodir)
    tinf_run=min([p.Tinf for p in listperiods if p.Run==run])
    closest_tinf=listperiods[0].Tinf
    for p in listperiods :
        if p.Run in listruns and p.Tinf<=tinf_run and p.Tinf>=closest_tinf :
            closest_tinf=p.Tinf
            closest_run=p.Run
    if "closest_run" not in vars() : print "GetClosestRunFromList:Not found a good run in list!"
    return closest_run

def ReadListeEvts(bolodir) :
    if not os.path.exists(bolodir+"/liste_evts.txt") :
        print "No liste_evt file!"
        return None
    f=open(bolodir+"/liste_evts.txt","r")
    lines=f.readlines()
    evtlist=[]
    for theline in lines :
        toto=theline.split(' ')
        toto=[x for x in toto if x!='']
        if len(toto)==2 and toto[0]!="#" :
            evtlist.append((toto[0],int(toto[1])))
    return evtlist

