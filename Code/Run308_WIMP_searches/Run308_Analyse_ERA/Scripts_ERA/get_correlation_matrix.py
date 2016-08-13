import numpy as np
import matplotlib.pylab as plt

from ROOT import*
from root_numpy import root2array

def get_correlation_matrix(bolo_name, analysis_type, tree_name = "data"):

	""" Get ionisation correlation matrix
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion and veto)

	Returns:
		void

	Raises:
		AssertionError 
	"""

	filename = "../Fond_ERA_merged/" + bolo_name + "_for_ion_corr_fond.root"

	# Convert a TTree in a ROOT file into a NumPy structured array
	arr = root2array(filename, tree_name)

	arr_EIA = arr["EIA"]
	arr_EIB = arr["EIB"]
	arr_EIC = arr["EIC"]
	arr_EID = arr["EID"]


	arr_stacked = np.vstack((arr_EIA, arr_EIB))
	arr_stacked = np.vstack((arr_stacked, arr_EIC))
	arr_stacked = np.vstack((arr_stacked, arr_EID))

	# R = np.corrcoef(arr_stacked)
	# plt.pcolor(R)
	# plt.colorbar()
	# plt.yticks(np.arange(0,4),range(0,4))
	# plt.xticks(np.arange(0,4),range(0,4))
	# plt.show()9+
	# raw_input()

	print
	print np.cov(arr_stacked)
	print 

	print np.corrcoef(arr_stacked)
	print
	
	L = np.linalg.cholesky(np.cov(arr_stacked))
	# print L
	# print np.ndarray.flatten(L)
	# LT = np.linalg.cholesky(np.cov(arr_stacked)).T
	# print L 
	# print 
	# print LT

	# print np.dot(L, LT)

	gaus_sample = np.random.normal(size=(4,1000000))
	proper = np.dot(L, gaus_sample)

	# np.savetxt("../Analyse_" + bolo_name + "/Text_files/" + bolo_name + "_Lmatrix_coeff.txt", np.ndarray.flatten(L))

	print np.dot(np.linalg.cholesky(np.cov(arr_stacked))*np.linalg.cholesky(np.cov(arr_stacked)).T)

bolo_name = "FID837"
analysis_type = "ana_0.5.0_0_5"
get_correlation_matrix(bolo_name, analysis_type)
