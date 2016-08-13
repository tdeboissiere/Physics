
import math
import scipy.special

def TemplateIon(xvect,c0,c1,c2,offset,ampl) :
    # convention: les c0/c1/c2 = parametres du fichier tmplt.txt
    xx = [x-offset for x in xvect]
    fprepulse=[0 for x in xx if x<0]
    fpostpulse=[ ampl*((1+scipy.special.erf(x/c2))/2.0)*(math.exp(-x/c0)-c1) for x in xx if x>=0]
    return fprepulse+fpostpulse

def TemplateIonBB2(xvect,c0,offset,ampl) :
    # convention: c0 = parametre du fichier tmplt.txt
    xx = [x-offset for x in xvect]
    fprepulse=[0 for x in xx if x<0]
    fpostpulse=[ ampl*((1+scipy.special.erf(x/c0))/2.0) for x in xx if x>=0]
    return fprepulse+fpostpulse

def TemplateChal(xvect,c0,c1,c2,c3,c4,c5,offset,ampl) :
    # convention: les c0/../c5 = parametres du fichier tmplt.txt
    if c1>c3 or c3>c5 : return [0 for x in xvect] # on force l'ordre ds le fit..
    xx = [x-offset for x in xvect]
    fprepulse=[0 for x in xx if x<0]
    fpostpulse=[ ampl*(1-math.exp(-x/c0))*(math.exp(-x/c1)+c2*math.exp(-x/c3)+c4*math.exp(-x/c5)) for x in xx if x>=0 ]
    return fprepulse+fpostpulse;

def TemplateChalNTD(xvect,c0,c1,offset,ampl) :
    # convention: les c0,c1 = parametres du fichier tmplt.txt
    # On force le domaine des parametres..
    if c0<0.001 or c0>100 or c1<0.1 or c1>400 : return [0 for x in xvect]
    xx = [x-offset for x in xvect]
    fprepulse=[0 for x in xx if x<0]
    fpostpulse=[ ampl*(1-math.exp(-x/c0))*(math.exp(-x/c1)) for x in xx if x>=0 ]
    return fprepulse+fpostpulse;

class RootTemplateChal:
    def __call__(self,x,p) :
        # p=[c0,..,c5,offset,ampl]
        xx=x[0]-p[6]
        if p[1]>p[3] or p[3]>p[5] or xx<0 :
            y=0
        else :
            y=p[7]*(1-math.exp(-xx/p[0]))*(math.exp(-xx/p[1])+p[2]*math.exp(-xx/p[3])+p[4]*math.exp(-xx/p[5]))
        return y
