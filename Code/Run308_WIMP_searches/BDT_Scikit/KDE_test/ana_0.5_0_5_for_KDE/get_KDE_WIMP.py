from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
import pickle
import Analysis_utilities as Ana_ut
import script_utils as script_utils




def get_KDE(bolo_name = "FID837"):

    """Fit a KDE to WIMP 6D
      
    Detail:

    Args:
    bolo_name
        
    Returns:
    void 
             
    Raises:
    void 
    """
    
    
    sampling   = np.loadtxt("./" + bolo_name + "_WIMP_mass_10.csv", delimiter = ",", skiprows=1)
    
    kde        = KernelDensity(bandwidth=0.1, kernel='gaussian', algorithm='ball_tree', rtol = 1E-8)
    sampling = sampling[:20000,:-1]
    kde.fit(sampling)
    print "fit complete"

    pickle_out = open("./" + bolo_name + "_KDE_heatonly.pkl", "wb")
    pickle.dump(kde, pickle_out)
    pickle_out.close()

    # kde_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Pickle_files/" )
    # kde_file = open(kde_dir + bolo_name + "_KDE_heatonly.pkl", "rb")
    # kde = pickle.load(kde_file)
    # kde_file.close()

    ################################################################
    # verification plots
    # #############################################################

    # sample = kde.sample(sampling.shape[0])
    # print sample.shape 
    # print sampling.shape

    # h = TH2F("h", "h", 2000, 0, 5, 2000, 0, 5)
    # PyRPl.fill_TH2(h, sample[:,0], sample[:,1])
    # h.Draw()
    # raw_input()

    # plt.scatter(sample[:,0], sample[:,3], color = "b", s=1)
    # plt.scatter(sampling[:,0], sampling[:,3], color = "r", s=1)
    # plt.show()
    # raw_input()

    # plt.hist(sampling[:,1], color = "r", bins= 200, histtype = "step", normed = True)
    # plt.hist(sample[:,1], color = "b", bins= 200, histtype = "step", normed = True)
    # plt.yscale("log")
    # plt.show()
    # raw_input()

    ################################################################
    # use grid search cross-validation to optimize the bandwidth
    # #############################################################
    # params = {'bandwidth': np.logspace(-5, -1, 20)}
    # params = {'bandwidth': np.array([0.5, 1, 5])}
    # grid = GridSearchCV(KernelDensity(kernel='gaussian', algorithm='ball_tree'), params, verbose=10, n_jobs=10)
    # grid.fit(sampling)

    # print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))
    # Find 0.5 with all points
    # raw_input()


get_KDE()