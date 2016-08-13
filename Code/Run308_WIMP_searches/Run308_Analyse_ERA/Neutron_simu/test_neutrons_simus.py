#! /usr/bin/env python

import ROOT
import sys,os,glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
plt.ion()


datadir="./"
list_files=glob.glob(datadir+"/*-Run308.root")

distris=[]
names=[]
for i,f in enumerate(list_files):
    ff=ROOT.TFile(f,"READ")
    names.append(f[len(datadir)+1:-5])
    h=ff.Get("htot")
    if i==0 : 
        nb=h.GetNbinsX()
        xvals=np.asarray([h.GetBinCenter(p+1) for p in range(nb)])
    else :
        if h.GetNbinsX()!=nb : print "Pbl",h.GetNbinsX(),f
        if h.GetBinCenter(1)!=xvals[0] : print "Pbl",h.GetBinCenter(1),f

    yvals=np.asarray([h.GetBinContent(p+1) for p in range(nb)])
    distris.append(yvals)


cNorm  = colors.Normalize(vmin=0, vmax=len(names))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('jet') )
# for i,thename in enumerate(names) :
#     plt.clf()
#     plt.plot(xvals,distris[i],color=scalarMap.to_rgba(i),label=thename)
#     plt.title(thename)
#     plt.yscale("log")
#     foo=raw_input("Press enter..")
    
# plt.legend(ncol=2)

vv=np.asarray(distris)
fulldistri=np.asarray([np.sum(vv[:,p]) for p in range(len(xvals))])
threshold=1.5

plt.plot(xvals,fulldistri,color='black')
plt.ylim([1E-8,1E2])
plt.yscale("log")
# plt.axvline(threshold,'--')

## Et la normalisation...
np.sum(fulldistri[15:]) 
# seuil 1.5keV == bin 14.5
# ==> 9.45 evts pour  8147 kg.j
np.sum(fulldistri[15:70])*200*24/(36.*365.)
# Entre 1.5 et 7 keVnr, pour 24 bolos a 200 jours: 1.7 neutrons
np.sum(fulldistri[15:30])*200*7/(36.*365.)
# Entre 1.5 et 3 keVnr, pour 7 bolos a 200 jours: 0.2 neutrons
raw_input()