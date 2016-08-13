#! /usr/bin/env python

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
## D'abord: maj listes runs et polars, en interactif
#####################################

os.system(eradir+"/python/BuildRunList.py "+paramfile)
os.system(eradir+"/python/BuildPolarList.py "+paramfile)

# Liste des runs qu'on va processer (critere un peu adhoc)
runlist=[run.Name for run in RunParams.ReadListeRuns(bolodir)]
runs_to_update=[run for run in runlist if (not os.path.exists(bolodir+"/Amplitudes/basicntp_"+run+"_"+bolo+".root") or not os.path.exists(bolodir+"/Amplitudes/eion_"+run+"_"+bolo+".root") or not os.path.exists(bolodir+"/Amplitudes/eheat_"+run+"_"+bolo+".root"))]
print "Runs to update:",runs_to_update

## Options des jobs
qsub_options_light="qsub -l ct=1:00:00,vmem=1G,sps=1 -j y"
qsub_options_lourd="qsub -l ct=4:00:00,vmem=4G,sps=1 -j y"


#####################################
## SambaToRoot
#####################################

# NB on pourrait faire un array job mais en fait c'est plus complique..
etape1="update_sambatoroot"
for run in runs_to_update :
    nomjob=etape1+"_"+run
    commande=qsub_options_lourd+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_sambatoroot.sh "+eradir+" "+anadir+" "+sambadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n"+script+"\neof")


#####################################
## SimpleProcess
#####################################

etape2="update_simpleprocess"
for run in runs_to_update :
    nomjob=etape2+"_"+run
    holdjob=etape1+"_"+run
    commande=qsub_options_lourd+" -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_simpleprocess.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## PeriodList
#####################################

etape3="update_periodlist"
holdjobs=",".join([etape2+"_"+r for r in runs_to_update])
commande=qsub_options_light+" -hold_jid "+holdjobs+" -N "+etape3+" -o "+logdir+"/"+etape3+".log"
script=scriptdir+"/run_buildperiodlist.sh "+eradir+" "+anadir+" "+bolo
os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## NoiseSpectra
#####################################

etape4="update_noisespectra"
holdjob=etape3
for run in runs_to_update :
    nomjob=etape4+"_"+run
    commande=qsub_options_lourd+" -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_buildnoisespectra.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## WienerProcess
#####################################

# NB: liste des voies en dur pour l'instant (-t 1-6)
etape5="update_wienerprocess"
for run in runs_to_update :
    nomjob=etape5+"_"+run
    holdjob=etape4+"_"+run
    commande=qsub_options_lourd+" -t 1-4 -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_wienerprocess.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")

# Voies veto synchronisees sur collectrices
etape5bis="update_wienersyncprocess"
for run in runs_to_update :
    nomjob=etape5bis+"_"+run
    holdjob=etape5+"_"+run
    commande=qsub_options_lourd+" -t 1-2 -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_wienersyncprocess.sh "+eradir+" "+anadir+" "+bolo+" "+run
    os.system(commande+" <<eof\n "+script+" \neof")

#####################################
## CorrectCrossTalk
#####################################

etape5ter="update_correctcrosstalk"
holdjobs=",".join([etape5bis+"_"+r for r in runs_to_update])
commande=qsub_options_light+" -hold_jid "+holdjobs+" -N "+etape5ter+" -o "+logdir+"/"+etape5ter+".log"
script=scriptdir+"/run_correctcrosstalk.sh "+eradir+" "+anadir+" "+bolo
os.system(commande+" <<eof\n "+script+" \neof")

#####################################
## BaselineEvol (ionisations)
#####################################

etape6="update_baselineevol"
holdjob=etape5ter
listevoies_ion=["Col1","Col2","Vet1","Vet2"] # en dur pour l'instant
for voie in listevoies_ion :
    nomjob=etape6+"_"+voie
    commande=qsub_options_light+" -hold_jid "+holdjobs+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_baselineevol.sh "+eradir+" "+anadir+" "+bolo+" "+voie
    os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## BuildIonEnergy + ClassifyEvents
#####################################

etape7="update_ionenergy"
holdjobs=",".join([etape6+"_"+v for v in listevoies_ion])
commande=qsub_options_light+" -hold_jid "+holdjobs+" -N "+etape7+" -o "+logdir+"/"+etape7+".log"
script=scriptdir+"/run_buildionenergy.sh "+eradir+" "+anadir+" "+bolo
os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## BaselineEvol(chals) + GainChalEvol
#####################################

etape8="update_chalevol"
holdjob=etape7
listevoies_chal=["Chal1","Chal2"] # en dur pour l'instant
for voie in listevoies_chal :
    nomjob=etape8+"_"+voie
    commande=qsub_options_lourd+" -hold_jid "+holdjob+" -N "+nomjob+" -o "+logdir+"/"+nomjob+".log"
    script=scriptdir+"/run_chalevol.sh "+eradir+" "+anadir+" "+bolo+" "+voie
    os.system(commande+" <<eof\n "+script+" \neof")


#####################################
## BuildChalEnergy + DST
#####################################

etape9="update_chalenergyanddst"
holdjobs=",".join([etape8+"_"+v for v in listevoies_chal])
commande=qsub_options_light+" -hold_jid "+holdjobs+" -N "+etape9+" -o "+logdir+"/"+etape9+".log"
script=scriptdir+"/run_buildchalenergyanddst.sh "+eradir+" "+anadir+" "+bolo
os.system(commande+" <<eof\n "+script+" \neof")
