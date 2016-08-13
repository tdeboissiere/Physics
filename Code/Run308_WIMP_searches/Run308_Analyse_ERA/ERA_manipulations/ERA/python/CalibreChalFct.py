
def ChalNLFunc(xvect,c0,c1,c2,c3) :
    f=[c0 for x in xvect]
    if c1!=0 :
        x0=((c2/c0)-1)/c1
        for i,x in enumerate(xvect) :
            if x<x0 : f[i]=c2
            else : f[i]=c0*(1+x*c1)
            if c3!=0 and abs(x-x0)<c3 :
                z=x+c3
                f[i]=(c2*(x0-x+c3)+c0*(z-x0)+0.5*c0*c1*(z*z-x0*x0))/(2*c3)
    return f

def ChalNLEnergy(x,c) :
    # same as ChalNLFunc but on single value...
    f=c[0]
    if c[1]!=0 :
        x0=((c[2]/c[0])-1)/c[1]
        if x<x0 : f=c[2]                    
        else : f=c[0]*(1+x*c[1])
        if c[3]!=0 and abs(x-x0)<c[3] :
            z=x+c[3]
            f=(c[2]*(x0-x+c[3])+c[0]*(z-x0)+0.5*c[0]*c[1]*(z*z-x0*x0))/(2*c[3])
    return f

class RootChalNLFunc:
    def __call__( self, xx, par ):
        # par=[c0,c1,c2,c3]
        if len(par)!= 4 : print "RootChalNLFunc: pbl nb of params"
        x=float(xx[0])
        f=0
        if par[1]==0 : f=par[0]
        else :
            x0=((par[2]/par[0])-1)/par[1]
            if x<x0 : f=par[2]
            else : f=par[0]*(1+par[1]*x)
            if par[3]!=0 and abs(x-x0)<par[3] :
                z=x+par[3]
                f=(par[2]*(x0-x+par[3])+par[0]*(z-x0)+0.5*par[0]*par[1]*(z*z-x0*x0))/(2*par[3])
        return f
