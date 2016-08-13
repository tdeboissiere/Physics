
from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import sys, os
#Quick code to check reconstruction bias

d_channel={"Chal1":6.73247, "Chal2":8.54526}

tchal1, file1 = PyRPl.open_ROOT_object("./ROOT_files/AmplDir_merged/FID837_wiener_Chal1.root", "wienerntp_FID837_Chal1")
tchal2, file2 = PyRPl.open_ROOT_object("./ROOT_files/AmplDir_merged/FID837_wiener_Chal2.root", "wienerntp_FID837_Chal2")
tbasic, filebasic = PyRPl.open_ROOT_object("./ROOT_files/AmplDir_merged/FID837_wiener_basic.root", "basicntp_FID837")

# tchal1.Draw("((WienerAmpl/6.73247-OriHeat)/OriHeat):OriHeat>>hist1(200,0,10,200,-2,2)", "", "colz")
# tchal2.Draw("WienerAmpl/6.73247:OriHeat>>hist2(200,0,10,200,0,10)", "", "colz")
# 
tchal1.Draw("WienerAmpl/6.73247:OriHeat>>hist1(2000,0,10,2000,0,10)", "", "colz")
tchal2.Draw("WienerAmpl/8.54526:OriHeat>>hist2(200,0,10,200,0,10)", "", "colz")

# PyRPl.process_TH2(hist1, Y_title = "(WienerCalib - Orig Heat)/ Orig Heat", X_title = "Original heat (keV)")
# PyRPl.process_TH2(hist2, Y_title = "Wiener Calib Chal2", X_title = "Original heat (keV)")

# partition_10keV = list(0.1*np.array(range(100)))

# d_entries = {}
# for elem in partition_10keV:
# 	d_entries[str(elem)]=[]

# def find_index(value):
# 	#return the element of the dict that is closest to value
# 	return str((value - value%0.1))

# nEntries = tchal1.GetEntries()

# for i in range(nEntries):

# 	sys.stdout.write('\r' + str(i)+' / '+str(-1+nEntries))
# 	sys.stdout.flush()

# 	tchal1.GetEntry(i)
# 	key = find_index(tchal1.OriHeat)
# 	if tchal1.WienerChi2<1.1:
# 		d_entries[key].append(tchal1.WienerAmpl/6.73247)

# list_average_ori = list(0.05+ 0.1*np.array(range(100)))
# list_average_wiener = []



# for elem in partition_10keV:
# 	list_average_wiener.append(np.average(d_entries[str(elem)]))


# list_average_wiener = np.array(list_average_wiener).astype(float)
# list_average_ori = np.array(list_average_ori).astype(float)

# list_average_relative = (list_average_ori- list_average_wiener)/list_average_ori


# ccb = TCanvas("ccb", "ccb")
# gr = TGraph(len(list_average_wiener), list_average_ori, list_average_wiener)
# PyRPl.process_TGraph(gr, Y_title = "OFAvg", X_title = "OriHeatAvg (keV)")
# gr.Draw("A*")

# ccb2 = TCanvas("ccb2", "ccb2")
# gr2 = TGraph(len(list_average_wiener), list_average_ori, list_average_relative)
# PyRPl.process_TGraph(gr2, Y_title = "(OFAvg - OriHeatAvg)/OriHeatAvg", X_title = "OriHeatAvg (keV)")
# gr2.Draw("A*")

cc = TCanvas("cc", "cc")
hist1.Draw("colz")

cc2 = TCanvas("cc2", "cc2")
hist2.Draw("colz")

raw_input()


# cc.Print("./Figures/FID837_reconstruction_bias_chal1_other_view.png")
# cc2.Print("./Figures/FID837_reconstruction_bias_chal2.png")

# ccb.Print("./Figures/FID837_avg_bias_chal1.png")
# ccb2.Print("./Figures/FID837_avg_bias_chal2_relative.png")