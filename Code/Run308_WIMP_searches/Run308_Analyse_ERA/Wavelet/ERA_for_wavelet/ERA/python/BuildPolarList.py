#! /usr/bin/env python

# Cree un fichier "liste_polars.txt" par bolo a partir de la liste des runs
import RunParams,sys,os
from GlobalParams import ReadGlobalParams

################################
# LECTURE DES PARAMETRES
paramfile="params_python.txt"
if len(sys.argv) == 2 : paramfile=sys.argv[1]
gParams=ReadGlobalParams(paramfile)
anadir=gParams['anadir']+"/"
bolo=gParams['bolo']
################################

bolodir=anadir+"/"+bolo+"/"
runs=RunParams.ReadListeRuns(bolodir)
liste_volts=[ [x.Vcol1,x.Vcol2,x.Vvet1,x.Vvet2,x.Vgar1,x.Vgar2] for x in runs]
liste_volts_uniqs=[]
for thevolt in liste_volts:
    if thevolt not in liste_volts_uniqs : liste_volts_uniqs.append(thevolt)

outputfile=open(bolodir+"/liste_polars.txt",'w')
outputfile.write("# v_flag v_fid v_col1 v_col2 v_vet1 v_vet2 v_gar1 v_gar2\n")
i=0
for volt in liste_volts_uniqs:
    vfid=abs(volt[0]-volt[1])
    ligne=[i,vfid,volt[0],volt[1],volt[2],volt[3],volt[4],volt[5]]
    ligne=[str(x) for x in ligne]
    outputfile.write(" ".join(ligne)+"\n")
    i+=1
outputfile.close()
