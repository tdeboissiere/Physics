#! /usr/bin/env python

################################################################################
# Fonction utilitaires generiques
# (paliant des manques des librairies usuelles)
################################################################################

import random,string
import numpy as np
import matplotlib.pyplot as plt

################################################################################

def randomname(length=6) :
    # Astuce pour donner un 'nom' aleatoire aux objets ROOT
    # Permet de creer plusieurs fois un TH1 etc.. en boucle
    thename=(''.join(random.choice(string.ascii_uppercase+string.digits) for _ in range(length)))
    return thename

################################################################################

def randomvalues(xarr,yarr,nb) :
    # tirage d'une v.a. en interpolant un tableau
    xarr=np.asarray(xarr)
    yarr=np.asarray(yarr)
    nbarr=len(xarr)
    if not np.all(np.diff(xarr) > 0) :
        print "xarr not in increasing order!"
        return 0
    if not np.all(yarr>=0) :
        print "negative values in yarr!"
        return 0
    # Calcul de la cdf: on prefere utiliser notre methode (effet de bord)
    #cdf=np.cumsum(yarr)
    cdf=np.zeros(nbarr)
    for i in range(1,nbarr) :
        mean=(yarr[i-1]+yarr[i])/2.
        integ=mean*(xarr[i]-xarr[i-1])
        cdf[i]=cdf[i-1]+integ
    cdf/=np.max(cdf)
    aa=np.random.rand(nb)
    values=np.interp(aa,cdf,xarr)
    return values

################################################################################
    
def densityplot(x,y,binsx=400,binsy=200,scaleperdeg2=False,maxdensity=0,mindensity=0,pltclf=True) :
    #http://oceanpython.org/2013/02/25/2d-histogram/
    H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
    if scaleperdeg2==True :
        delta_ra=abs(xedges[1]-xedges[0])
        delta_dec=abs(yedges[1]-yedges[0])
        cosdelta=np.cos(yedges[:-1]*np.pi/180.) # yedges a (binsy+1) elements
        for k in range(binsy) : H[:,k] = H[:,k]/(delta_ra*delta_dec*cosdelta[k])
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    if maxdensity>0 :
        wcut=np.where( (H>maxdensity) )
        H[wcut]=maxdensity
    if mindensity>0 :
        wcut=np.where( (H<mindensity) )
        H[wcut]=0
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    if pltclf : plt.clf()
    plt.pcolormesh(xedges,yedges,Hmasked)
    plt.colorbar()
    if scaleperdeg2==True : plt.title(r'density / deg$^2$')

################################################################################

def contourplot(x,y,ncont=10,colors=None,pltclf=True,binsx=100,binsy=100) : # EN CONSTRUCTION
    H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    xcenters=(xedges[:-1]+xedges[1:])/2.
    ycenters=(yedges[:-1]+yedges[1:])/2.
    if pltclf : plt.clf()
    plt.contour(xcenters,ycenters,H,ncont,colors=colors,linewidths=2)

################################################################################

def scalarfieldplot(var,x,y,binsx=400,binsy=200,maxvalue=0,minvalue=0) :
    # 2D map of var as a fct of x,y
    H_nbobjs, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
    H_somme_var, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy),normed=False,weights=var)
    # H_somme_var contient la somme des (var) dans chaque bin. Reste a diviser par le nb d'obj pour avoir la moyenne
    ww=np.where( (H_nbobjs>0) )
    H=H_somme_var*0
    H[ww]=H_somme_var[ww]/H_nbobjs[ww] # Ensuite c'est comme densityplot :
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    if maxvalue!=0 :
        wcut=np.where( (Hmasked>maxvalue) )
        Hmasked[wcut]=maxvalue
        Hmasked = np.ma.masked_where(Hmasked==maxvalue,Hmasked) 
    if minvalue!=0 :
        wcut=np.where( (Hmasked<minvalue) )
        Hmasked[wcut]=minvalue
        Hmasked = np.ma.masked_where(Hmasked==minvalue,Hmasked) 
    plt.clf()
    plt.pcolormesh(xedges,yedges,Hmasked)
    plt.colorbar()

################################################################################

def efficiencyplot(d1,d2,bins=200,range=(17,24)) :
	# distri d'une variable d1=avant cut, d2=apres cut => mesure cut efficiency en fct de cette var.
	# TODO: ajout barre d'erreur
	plt.clf()
	plt.subplot(2,1,1)
	n1,b1,p1=plt.hist(d1,bins=bins,range=range,histtype='step')
	n2,b2,p2=plt.hist(d2,bins=bins,range=range,histtype='step')
	toto=plt.ylim(0,np.max(n1)*1.1)
	ratio=np.zeros(len(n1))
	w=np.where( (n1>0) )
	ratio[w]=n2[w]/(1.0*n1[w])
	left,right=b1[:-1],b1[1:]
	X = np.array([left,right]).T.flatten()
	Y = np.array([ratio,ratio]).T.flatten()
	plt.subplot(2,1,2)
	plt.plot(X,Y,color='r')
	plt.plot(np.asarray(range),[1,1],color='b',linestyle='--')
	toto=plt.ylim((0,1.1))

################################################################################

def profilehistogram(xarr,yarr,nbins=50,xmin=None,xmax=None,plot_meanerror=True,color='r') :
    xarr=np.asarray(xarr)
    yarr=np.asarray(yarr)
    if xmin is None : xmin=np.min(xarr)
    if xmax is None : xmax=np.max(xarr)
    binedges=xmin+((xmax-xmin)/float(nbins))*np.arange(nbins+1)
    bincenters = xmin + ((xmax-xmin)/float(nbins))*np.arange(nbins) + ((xmax-xmin)/(2.0*nbins))
    binsarr=np.digitize(xarr,binedges)
    ymean=np.zeros(nbins)-999.99
    y_standdev=np.zeros(nbins)
    y_meanerror=np.zeros(nbins)
    for i in np.arange(nbins) :
        ww=np.where(binsarr==i+1)
        if len(ww[0])>0 :
            ymean[i]=np.mean(yarr[ww])
            y_standdev[i]=np.std(yarr[ww])
            y_meanerror[i]=y_standdev[i]/np.sqrt(len(ww[0]))
    the_error=y_meanerror if plot_meanerror else y_standdev
    wc=np.where(ymean!=-999.99)
    plt.errorbar(bincenters[wc],ymean[wc],xerr=(xmax-xmin)/(2.0*nbins),yerr=the_error[wc],fmt='o',color=color)
    