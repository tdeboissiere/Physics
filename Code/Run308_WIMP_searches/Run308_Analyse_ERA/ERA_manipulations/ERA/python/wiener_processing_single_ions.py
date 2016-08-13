#! /usr/bin/env python

# Processing : uniquement fit wiener

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

# Liste des runs qu'on va processer
runlist=[run.Name for run in RunParams.ReadListeRuns(bolodir)]
## Options des jobs
qsub_options_light="qsub -l ct=1:00:00,vmem=1G,sps=1 -j y"
qsub_options_lourd="qsub -l ct=4:00:00,vmem=4G,sps=1 -j y"

#####################################
## WienerProcess
#####################################

# NB: liste des voies en dur pour l'instant (-t 3-6 soit les ionisations)
etape5="start_wienerprocess"
for run in runlist :
    nomjob=etape5+"_"+run
    commande=qsub_options_lourd+" -t 3-6 -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_wienerprocess.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")

print "Once jobs are done, you can calibrate channels, etc..."
