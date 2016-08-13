#! /usr/bin/env python

# Processing jusqu'a wiener spectra = partie "automatique"
# avant de visualiser les evts

import sys,os
import RunParams
from GlobalParams import ReadGlobalParams

#####################################
## Lecture des parametres
#####################################

paramfile="params_processing.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
eradir=gParams['eradir']
anadir=gParams['anadir']
bolo=gParams["bolo"]
sambadir=gParams["sambadir"]
bolodir=anadir+"/"+bolo
logdir=anadir+"/logs"
scriptdir=eradir+"/scripts"

#####################################
## D'abord: construction listes runs et polars, en interactif
#####################################

os.system(eradir+"/python/BuildRunList.py "+paramfile)
os.system(eradir+"/python/BuildPolarList.py "+paramfile)

# Liste des runs qu'on va processer (tous les runs de la liste)
runlist=[run.Name for run in RunParams.ReadListeRuns(bolodir)]
print "Run list:",runlist

## Options des jobs
qsub_options_light="qsub -l ct=1:00:00,vmem=1G,sps=1 -j y"
qsub_options_lourd="qsub -l ct=4:00:00,vmem=4G,sps=1 -j y"


#####################################
## SambaToRoot
#####################################

# NB on pourrait faire un array job mais en fait c'est plus complique..
etape1="start_sambatoroot"
for run in runlist :
    nomjob=etape1+"_"+run
    commande=qsub_options_lourd+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_sambatoroot.sh "+eradir+" "+anadir+" "+sambadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n"+script+"\neof")


#####################################
## SimpleProcess
#####################################

etape2="start_simpleprocess"
for run in runlist :
    nomjob=etape2+"_"+run
    holdjob=etape1+"_"+run
    commande=qsub_options_lourd+" -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_simpleprocess.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## PeriodList
#####################################

etape3="start_periodlist"
holdjobs=",".join([etape2+"_"+r for r in runlist])
commande=qsub_options_light+" -hold_jid "+holdjobs+" -N "+etape3+" -o "+logdir+"/"+etape3+".log"
script=scriptdir+"/run_buildperiodlist.sh "+eradir+" "+anadir+" "+bolo
os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## NoiseSpectra
#####################################

etape4="start_noisespectra"
holdjob=etape3
for run in runlist :
    nomjob=etape4+"_"+run
    commande=qsub_options_lourd+" -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_buildnoisespectra.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")

print "Once jobs are done, you can visualize events and make templates..."
