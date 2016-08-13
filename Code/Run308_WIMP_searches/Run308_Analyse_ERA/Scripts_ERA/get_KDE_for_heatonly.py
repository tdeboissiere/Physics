from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
import pickle
import Analysis_utilities as Ana_ut
import script_utils as script_utils



def save_data(bolo_name = "FID837"):

  """Get 3D EC1, EC2, JOUR heatonly distribution
  
  Detail:

  Args:
    bolo_name
    
  Returns:
    void

  Raises:
    void
  """

  tree, file_tree = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/"+bolo_name+"_fond.root", "data")
  
  #Load standard cuts
  standard_cuts = Ana_ut.open_cut_file(bolo_name, "TCuts.txt")  
  standard_cuts = standard_cuts + "&&KTH<1&&KTH>0&&0.5*(EIB+EID)<0" 
    
  l_heatonly = TEventList("l_heatonly")
  arr_heatonly = []
  tree.Draw(">>l_heatonly",standard_cuts)
  pop_len       = l_heatonly.GetN()
  for k in range(pop_len):
    counter  = l_heatonly.GetEntry(k)
    tree.GetEntry(counter)
    arr_heatonly.append([tree.EC1, tree.EC2, tree.JOUR]) #, 1E6*tree.UT1+tree.UT2])
  out_dir = script_utils.create_directory("../Analyse_" + bolo_name + "/Text_files/Pop_BDT/" )
  np.savetxt(out_dir  +  bolo_name + "_heatonly_2D_and_time.txt", np.array(arr_heatonly), delimiter = ",")

def get_KDE(bolo_name = "FID837"):

    """Fit a KDE to heat-only data EC1, EC2, time
      
    Detail:

    Args:
    bolo_name
        
    Returns:
    void 
             
    Raises:
    void 
    """
    
    
    in_dir    = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Text_files/Pop_BDT/" )
    sampling   = np.loadtxt(in_dir + bolo_name + "_heatonly_2D_and_time.txt", delimiter = ",")
    
    kde        = KernelDensity(bandwidth=0.06, kernel='gaussian', algorithm='ball_tree', rtol = 1E-8)
    sampling = sampling[:,:]
    kde.fit(sampling)
    print "fit complete"

    out_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Pickle_files/" )
    pickle_out = open(out_dir + bolo_name + "_KDE_heatonly.pkl", "wb")
    pickle.dump(kde, pickle_out)
    pickle_out.close()

    # kde_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Pickle_files/" )
    # kde_file = open(kde_dir + bolo_name + "_KDE_heatonly.pkl", "rb")
    # kde = pickle.load(kde_file)
    # kde_file.close()

    ################################################################
    # verification plots
    # #############################################################

    sample = kde.sample(sampling.shape[0])
    print sample.shape 
    print sampling.shape

    # h = TH2F("h", "h", 2000, 0, 5, 2000, 0, 5)
    # PyRPl.fill_TH2(h, sample[:,0], sample[:,1])
    # h.Draw()
    # raw_input()

    # plt.scatter(sample[:,0], sample[:,1], color = "b", s=1)
    # plt.scatter(sampling[:,0], sampling[:,1], color = "r", s=1)
    # plt.show()
    # raw_input()

    plt.hist(sampling[:,1], color = "r", bins= 200, histtype = "step", normed = True)
    plt.hist(sample[:,1], color = "b", bins= 200, histtype = "step", normed = True)
    plt.yscale("log")
    plt.show()
    raw_input()

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


def plot_simu_vs_data(bolo_name, analysis_type):

  """Plot simu versus data
  
  Detail:

  Args:
    bolo_name
    
  Returns:
    void

  Raises:
    void
  """

  tree, file_tree = PyRPl.open_ROOT_object("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Fond_ERA_merged/"+bolo_name+ "_" + analysis_type + "_fond.root", "data")
  heat_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_" + bolo_name + "/" + analysis_type + "/Heatonly/ROOT_files/"
  theat, file_heat = PyRPl.open_ROOT_object(heat_path + bolo_name + "_heatonly_tree.root", "t_new0")

  print theat.GetEntries()

  # cc = TCanvas("cc", "cc")
  # # cc.Divide(2,1)
  # # cc.cd(1)
  # theat.Draw("EC1:EC2>>hist(1000,0,5,1000,0,5)", "", "")
  # # cc.cd(2)
  # tree.Draw("EC1:EC2>>hist2(1000,0,5,1000,0,5)", "", "same")
  # hist2.SetMarkerColor(kRed)
  # hist2.Draw("same")

  bin = 200
  cc1 = TCanvas("cc1", "cc1")
  cc1.Divide(2,1)
  cc1.cd(1)
  tree.Draw("EC1>>hist1EC1(" + str(bin) + ",0,5)", "EC2>0.8 && EC2<1.1")
  theat.Draw("EC1>>hist2EC1(" + str(bin) + ",0,5)","EC2>0.8 && EC2<1.1", "same")
  hist1EC1.Scale(1./hist1EC1.Integral())
  hist2EC1.Scale(1./hist2EC1.Integral())
  hist2EC1.SetLineColor(kRed)
  hist1EC1.Draw()
  hist2EC1.Draw("same")
  cc1.cd(2)
  tree.Draw("EC2>>hist1EC2(" + str(bin) + ",0,5)")
  theat.Draw("EC2>>hist2EC2(" + str(bin) + ",0,5)","", "same")
  hist1EC2.Scale(1./hist1EC2.Integral())
  hist2EC2.Scale(1./hist2EC2.Integral())
  hist2EC2.SetLineColor(kRed)
  hist1EC2.Draw() 
  hist2EC2.Draw("same")

  raw_input()

# save_data()
get_KDE()
# plot_simu_vs_data("FID837", "ana_0.5_0_5")
