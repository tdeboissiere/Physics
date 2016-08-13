from ROOT import *
import numpy as np


class NR_to_EE:
   def __call__( self, x, par ):
		ENR=x[0]
		Vf=par[0]
		Q=0.16*pow(ENR,0.18)

		return ENR*(1+Q*float(Vf)/3.)/(1+float(Vf)/3.)


class EE_to_NR:
	def __call__( self, x, par ):

		f1 = TF1("r", NR_to_EE(), 0,1000,1)
		f1.SetNpx(1000)
		f1.SetParameter(0,par[0])
		return f1.GetX(x[0])


## test formula NR_to_EE
# f = TF1( 'test', NR_to_EE(), 0,100,1 )
# f.SetParameter(0,8.)
# print (1+8./3*(0.16*pow(10,0.18)))*10/(1+8./3.), f.Eval(10)

 
#test formula EE_to_NR
f = TF1( 'test', EE_to_NR(), 0,100,1 )
f.SetParameter(0,8.)
print  f.Eval(1)

def get_FWHM_Erec(fwhm_Eh):

	arr_heat = np.random.normal(0,float(fwhm_Eh)/2.3548, 10000)
	arr_Er   =[]

	fEE_to_NR = TF1( 'test', EE_to_NR(), 0,100,1 )
	fEE_to_NR.SetParameter(0,8.)

	for i in range(len(arr_heat)):
		# print i
		if arr_heat[i]>0:
			arr_Er.append(fEE_to_NR.Eval(arr_heat[i]))
   
    
	h = TH1F("h","h", 100,0,10)
	for i in range(len(arr_Er)):
		h.Fill(arr_Er[i])

	h.Draw()
	h.Fit("gaus", "LL", "", 0, 10)
	raw_input()




def get_boundaries(Eh_inf, Eh_sup):

	fEE_to_NR = TF1( 'test', EE_to_NR(), 0,500,1 )
	fEE_to_NR.SetParameter(0,8.)

	print "bound low: ", fEE_to_NR(Eh_inf), "bound sup: ", fEE_to_NR(Eh_sup)

get_FWHM_Erec(0.327)
# get_boundaries(0,400)
pass