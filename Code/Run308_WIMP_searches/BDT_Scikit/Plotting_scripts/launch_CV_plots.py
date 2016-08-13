import numpy as np
import scipy as sp
import sys
sys.path.append("../")
import xgboost as xgb
import sklearn.cross_validation as cv
import itertools
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Normalizer, Binarizer
import data_preparation as dp 
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pylab as plt
import BDT_file_handler as BDT_fh
import BDT_CV_xgboost_for_visual as xgCV
import plotting_CV_functions as plfunc

def launch_CV_plots():

    """
    Detail:
        Launch plots to visualise parameter tuning 

    Args:
        void

    Returns:
        void

    Raises:
        void     
    """

    plt.ion()

    bolo_name = "FID837"
    analysis_type = "ana_1.5_0_5"
    exposure = 65.
    # list_mass = ["3", "4", "5", "6", "7", "10", "25"]
    list_mass = ["5", "6", "7", "10", "25"]
    list_mass = ["5"]

    list_variables = ["EIA", "EIB", "EIC", "EID", "EC1", "EC2", "test"]
    #Loop over masses
    for WIMP_mass in list_mass:
        d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                         "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}


        ###############################
        # Influence of training size
        ###############################
        
        # # xgboost parameter setup
        # param = {}
        # folds = 5
        # list_splits = [0.01, 0.02, 0.05, 0.1, 0.5]
        # eta, subsample, depth, num_rounds = 0.1, 0.9, 4, 100
        # param['objective'], param['eval_metric'], param['silent'] = 'binary:logistic', 'auc', 1
        # param['bst:eta'], param['bst:subsample'], param['bst:max_depth'] = eta, subsample, depth

        # plfunc.plot_size_training_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_splits,True, list_variables)
        # plfunc.plot_size_training_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_splits,False, list_variables)


        ###############################
        # Influence of boosting rounds
        ###############################
        
        #xgboost parameter setup
        # param = {}
        # folds = 5
        # list_boosting_rounds = [5, 10, 20, 40, 50, 80, 100]
        # eta, subsample, depth = 0.1, 0.9, 4
        # param['objective'], param['eval_metric'], param['silent'] = 'binary:logistic', 'auc', 1
        # param['bst:eta'], param['bst:subsample'], param['bst:max_depth'] = eta, subsample, depth

        # plfunc.plot_boosting_round_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, folds, list_boosting_rounds,True, list_variables)
        # plfunc.plot_boosting_round_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, folds, list_boosting_rounds,False, list_variables)


        ###############################
        # Influence of tree depth
        ###############################
        
        # # xgboost parameter setup
        # param = {}
        # folds = 5
        # list_depth = range(1,10)
        # eta, subsample, num_rounds = 0.1, 0.9, 100
        # param['objective'], param['eval_metric'], param['silent'] = 'binary:logistic', 'auc', 1
        # param['bst:eta'], param['bst:subsample'] = eta, subsample

        # plfunc.plot_depth_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_depth,True, list_variables)


        ###############################
        # Influence of learning rate
        ###############################
        
        # # xgboost parameter setup
        param = {}
        folds = 5
        list_eta = [0.01,0.02,0.05,0.08,0.1,0.2,0.5]
        subsample, num_rounds = 0.9, 100
        param['objective'], param['eval_metric'], param['silent'] = 'binary:logistic', 'auc', 1
        param['bst:subsample'] = subsample

        plfunc.plot_eta_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_eta,True, list_variables)


        ####################################################
        # Influence of relative training sample sizes
        # Here, same number of training events for each class
        ####################################################
        
        # # xgboost parameter setup
        # param = {}
        # folds = 5
        # list_size = [100, 500, 1000, 5000, 10000, 30000]
        # eta, subsample, depth, num_rounds = 0.1, 0.9, 4, 100
        # param['objective'], param['eval_metric'], param['silent'] = 'binary:logistic', 'auc', 1
        # param['bst:eta'], param['bst:subsample'], param['bst:max_depth'] = eta, subsample, depth

        # plfunc.plot_equal_training_size_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_size,True, list_variables)
        # plfunc.plot_equal_training_size_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_size,False, list_variables)

        ###################################################
        # Influence of relative training sample sizes
        # Here, varying relative contribution for each class
        ###################################################
        #xgboost parameter setup
        # param = {}
        # folds = 5
        # #We will vary the number of training events for one class
        # list_size = [1000, 5000, 10000, 20000, 30000, 50000]
        # list_size = [100, 1000, 50000] #, 10000, 20000]
        # eta, subsample, depth, num_rounds = 0.1, 0.9, 4, 100
        # param['objective'], param['eval_metric'], param['silent'] = 'binary:logistic', 'auc', 1
        # param['bst:eta'], param['bst:subsample'], param['bst:max_depth'] = eta, subsample, depth
        # varying_class = "WIMP_mass_5"

        # plfunc.plot_variable_training_size_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, varying_class, list_size, True, list_variables)
        # plfunc.plot_variable_training_size_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, varying_class, list_size, False, list_variables)

launch_CV_plots()

#Conclusions
# weight code OK (when I manually set equal weights for signal and backgrounds, same results as the xgboost option)
# with weights: not much to no dependence on number of training samples in each category (as expected) as long as it is above a sufficient value