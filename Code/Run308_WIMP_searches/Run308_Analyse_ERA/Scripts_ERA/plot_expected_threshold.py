from ROOT import *
import matplotlib.pylab as plt
import numpy as np

#EDW III keVee
arr = np.linspace(0,5,100)
arr_erf_ee = [0.5*(1+TMath.Erf( (elem-0.65) / (TMath.Sqrt(2)*0.1648))) for elem in arr]


#EDW III keVNR
#0.65 keVee == 1.625 keVNR
#sigma_kevee @ 0.165 keVee  == sigma_kevNR @ 0.388 keVNR
arr = np.linspace(0,5,100)
arr_erf_NR = [0.5*(1+ TMath.Erf( (elem-1.625) / (TMath.Sqrt(2)*0.388))) for elem in arr]



plt.plot(arr, arr_erf_ee, 'r', label = "keVee")
plt.plot(arr, arr_erf_NR, 'b', label = "keVNR")
plt.legend(loc = "upper right")
plt.savefig("../../FID837_trigger_eff.png")
plt.show()
raw_input()
