from ROOT import *
import numpy as np
import matplotlib.pylab as plt
import PyROOTPlots as PyRPl

x = np.random.uniform(low=0.0, high=2000, size=10000)
y = [np.random.normal(0.25*elem,50) for elem in x]

x_noise = np.random.uniform(low=0.0, high=2000, size=100)
y_noise = np.random.uniform(low=0.0, high=2000, size=100)

x_gaussian_noise = np.random.uniform(low=0.0, high=2000, size=4000)
y_gaussian_noise = [np.random.normal(0.25*elem,400) for elem in x]

x_g, y_g = [], []

for elx, ely in zip(x_gaussian_noise, y_gaussian_noise):
	if ely < 0.25*elx :
		x_g.append(elx)
		y_g.append(ely)



h =TH2F("h", "h", 100, 0, 2000, 100, 0,2000)
h_noise =TH2F("h_noise", "h_noise", 100, 0, 2000, 100, 0,2000)
h_gaussian_noise =TH2F("h_gaussian_noise", "h_gaussian_noise", 100, 0, 2000, 100, 0,2000)
PyRPl.fill_TH2(h,x,y)
PyRPl.fill_TH2(h_noise,x_noise,y_noise)
PyRPl.fill_TH2(h_gaussian_noise,x_g,y_g)

PyRPl.process_TH2(h, X_title = "Ionisation D (ADU)", Y_title = "Ionisation C (ADU)", X_title_offset=1.14, Y_title_offset = 1.14)

h.Draw() 
h.Fit("pol1")
h_noise.Draw("same")
h_gaussian_noise.Draw("same")
raw_input()

plt.scatter(x,y)
plt.ylim([0,1000])
plt.xlim([0,2000])
plt.show() 
raw_input()