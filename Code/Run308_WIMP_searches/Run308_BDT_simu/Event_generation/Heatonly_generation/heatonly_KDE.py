from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
import pickle
import Analysis_utilities as Ana_ut
import script_utils as script_utils



def generate_heatonly_from_KDE(bolo_name, analysis_type, heatonly_num_event):

    """
    Detail:
        Generate EC1, EC2, 1E6*UT1+UT2 sample from Heatonly KDE

    Args:
        bolo_name     = (str) bolometer name
        analysis_type = (str) type of analysis (cuts on heat and ion and veto)
        heatonly_num_event = (int) events for heatonly

    Returns:
        void

    Raises:
        void     
    """

    kde_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Pickle_files/" )
    kde_file = open(kde_dir + bolo_name + "_KDE_heatonly.pkl", "rb")
    kde = pickle.load(kde_file)
    kde_file.close()

    sample = kde.sample(2*heatonly_num_event)

    sample0 = sample[:heatonly_num_event]
    sample1 = sample[heatonly_num_event:]

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
    out_dir = script_utils.create_directory(gen_path + "BDT_" + bolo_name + "/" + analysis_type + "/Heatonly/Text_files/")
    np.savetxt(out_dir + bolo_name + "_heatonly_2D_time_0.txt" , sample0, delimiter = " ", fmt = "%.3f")
    np.savetxt(out_dir + bolo_name + "_heatonly_2D_time_1.txt" , sample1, delimiter = " ", fmt = "%.3f")

def generate_heatonly_from_KDE_scaling(bolo_name, analysis_type, heatonly_num_event):

    """
    Detail:
        Generate EC1, EC2, 1E6*UT1+UT2 sample from Heatonly KDE for scaling computations

    Args:
        bolo_name     = (str) bolometer name
        analysis_type = (str) type of analysis (cuts on heat and ion and veto)
        heatonly_num_event = (int) events for heatonly

    Returns:
        void

    Raises:
        void     
    """

    kde_dir = script_utils.create_directory("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_" + bolo_name + "/Pickle_files/" )
    kde_file = open(kde_dir + bolo_name + "_KDE_heatonly.pkl", "rb")
    kde = pickle.load(kde_file)
    kde_file.close()

    sample = kde.sample(heatonly_num_event)

    gen_path = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
    out_dir = script_utils.create_directory(gen_path + "BDT_" + bolo_name + "/" + analysis_type + "/Heatonly/Text_files/")
    np.savetxt(out_dir + bolo_name + "_heatonly_2D_time_scaling.txt" , sample, delimiter = " ", fmt = "%.3f")

