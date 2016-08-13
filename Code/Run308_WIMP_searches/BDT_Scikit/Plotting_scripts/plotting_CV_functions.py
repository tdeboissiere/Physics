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
import script_utils as script_utils


def plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_xaxis, list_auc_train, list_auc_test, d_fig):

    """
    Detail:
        Plot the AUC + std dev 
    
    Args:
        WIMP_mass (str)            = WIMP mass
        bolo_name (str)            = bolometer name 
        analysis_type (str)        = type of analysis (which box cut)
        fold (int)                 = number of CV folds
        list_xaxis (int)           = xaxis points for the plot
        list_auc_train/test (list) = list of points to plot (eg AUC as a function of nu)
        d_fig (dict)               = stores figure information (title, axis title etc)

    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """

    train_AUC = np.array(list_auc_train)
    test_AUC= np.array(list_auc_test)
    
    train_AUC_mean = np.mean(train_AUC, axis=1)
    train_AUC_std = np.std(train_AUC, axis=1)
    test_AUC_mean = np.mean(test_AUC, axis=1)
    test_AUC_std = np.std(test_AUC, axis=1)
    

    plt.plot(list_xaxis, train_AUC_mean, 'o-', color ="r", label ="Training score")
    plt.plot(list_xaxis, test_AUC_mean, 'o-', color="g", label =str(folds) + " -fold Cross-Validation score")
    plt.fill_between(list_xaxis, train_AUC_mean - train_AUC_std, train_AUC_mean + train_AUC_std, alpha =0.1, color="r")
    plt.fill_between(list_xaxis, test_AUC_mean - test_AUC_std, test_AUC_mean + test_AUC_std, alpha =0.1, color="g")

    plt.xlabel(d_fig["x_title"])
    plt.ylabel(d_fig["y_title"])
    plt.ylim([0.98,1])
    plt.legend(loc ="best")
    # plt.show()
    out_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/" + d_fig["weight"] + "/")
    plt.savefig(out_dir + bolo_name + "_" + d_fig["fig_title"] + "_WIMP_mass_" + WIMP_mass + ".png")
    plt.clf()
    # raw_input()


def plot_size_training_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_splits, bool_weighted, list_variables):

    """
    Detail:
        Plot influence of training size on the AUC
    
    Args:
        WIMP_mass (str)     = WIMP mass
        bolo_name (str)     = bolometer name 
        d_event_dir (dict)  = dict to get the proper directory of each event class 
        exposure            = exposure in days
        analysis_type (str) = type of analysis (which box cut)
        param (dict)        = list of xgboost parameters
        num_rounds (int)    = number of boosting rounds
        fold (int)          = number of CV folds
        list_splits (int)   = list of data fraction (1pct, 2pct, 10 pct etc...)
        bool_weighted (bool) = True if include weights in BDT
        list_variables (list) = list of BDT features
    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """
    
    list_auc_test  = []
    list_auc_train = []


    for split in list_splits:
        
        script_utils.print_utility("CV with " +str(split*100)+ " % of the data set")
        
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, "", d_event_dir.keys(), exposure, list_variables, datasplit=split)
        
        AUC_train, AUC_test = [], []
        if bool_weighted:
            AUC_train, AUC_test = xgCV.CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        else:
            AUC_train, AUC_test = xgCV.CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        list_auc_train.append(AUC_train)
        list_auc_test.append(AUC_test)
        
        script_utils.print_utility("Mean AUC_train " +str(np.mean(AUC_train)))
        script_utils.print_utility("Mean AUC_test " +str(np.mean(AUC_test)))

    d_fig = {}
    if bool_weighted:
        d_fig["weight"] = "With_weights" 
    else:
        d_fig["weight"] = "No_weights"
    d_fig["fig_title"] = "training_size_influence"
    d_fig["x_title"] = "Size of training sample in %"
    d_fig["y_title"] = "AUC Score"
    plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_splits, list_auc_train, list_auc_test, d_fig)


def plot_equal_training_size_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_size, bool_weighted, list_variables):

    """
    Detail:
        Plot influence of training size on the AUC (with same number of events for each class)
    
    Args:
        WIMP_mass (str)     = WIMP mass
        bolo_name (str)     = bolometer name 
        d_event_dir (dict)  = dict to get the proper directory of each event class 
        exposure            = exposure in days
        analysis_type (str) = type of analysis (which box cut)
        param (dict)        = list of xgboost parameters
        num_rounds (int)    = number of boosting rounds
        folds (int)         = number of CV folds
        list_size (int)     = list of training sample size (same size for each class)
        bool_weighted (bool) = True if include weights in BDT
        list_variables (list) = list of BDT features
    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """
    
    list_auc_test  = []
    list_auc_train = []

    d_num_events = {}
    list_evt = ["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly", "WIMP_mass_" + str(WIMP_mass)]

    for size in list_size:
        
        for evt_type in list_evt:
            d_num_events[evt_type] = size

        script_utils.print_utility("CV with size of " +str(size)+ " events for each class")
        
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, "", d_event_dir.keys(), exposure, list_variables, d_num_events = d_num_events)
        
        AUC_train, AUC_test = [], []
        if bool_weighted:
            AUC_train, AUC_test = xgCV.CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        else:
            AUC_train, AUC_test = xgCV.CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        list_auc_train.append(AUC_train)
        list_auc_test.append(AUC_test)
        
        script_utils.print_utility("Mean AUC_train " +str(np.mean(AUC_train)))
        script_utils.print_utility("Mean AUC_test " +str(np.mean(AUC_test)))

    d_fig = {}
    if bool_weighted:
        d_fig["weight"] = "With_weights" 
    else:
        d_fig["weight"] = "No_weights"
    d_fig["fig_title"] = "equal_training_size_influence"
    d_fig["x_title"] = "Number of training events for each class"
    d_fig["y_title"] = "AUC Score"
    plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_size, list_auc_train, list_auc_test, d_fig)


def plot_variable_training_size_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, varying_class, list_size, bool_weighted, list_variables):

    """
    Detail:
        Plot influence of training size on the AUC (varying the contribution of one class wrt the others)
    
    Args:
        WIMP_mass (str)     = WIMP mass
        bolo_name (str)     = bolometer name 
        d_event_dir (dict)  = dict to get the proper directory of each event class 
        exposure            = exposure in days
        analysis_type (str) = type of analysis (which box cut)
        param (dict)        = list of xgboost parameters
        num_rounds (int)    = number of boosting rounds
        folds (int)         = number of CV folds
        varying_class (str) = the class for which we vary the number of training samples
        list_size (int)     = list of training sample size (for the class which we vary)
        bool_weighted (bool) = True if include weights in BDT
        list_variables (list) = list of BDT features
    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """
    
    list_auc_test  = []
    list_auc_train = []

    d_num_events = {}
    list_evt = ["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly", "WIMP_mass_" + str(WIMP_mass)]
    list_evt.remove(varying_class)

    for size in list_size:
        
        for evt_type in list_evt:
            d_num_events[evt_type] = 10000 #standard value

        d_num_events[varying_class] = size

        script_utils.print_utility("CV with size of " +str(size)+ " for " + varying_class + " events")
        
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, "", d_event_dir.keys(), exposure, list_variables, d_num_events = d_num_events)
        
        AUC_train, AUC_test = [], []
        if bool_weighted:
            AUC_train, AUC_test = xgCV.CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        else:
            AUC_train, AUC_test = xgCV.CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        list_auc_train.append(AUC_train)
        list_auc_test.append(AUC_test)
        
        script_utils.print_utility("Mean AUC_train " +str(np.mean(AUC_train)))
        script_utils.print_utility("Mean AUC_test " +str(np.mean(AUC_test)))

    d_fig = {}
    if bool_weighted:
        d_fig["weight"] = "With_weights" 
    else:
        d_fig["weight"] = "No_weights"
    d_fig["fig_title"] = "varying_" + varying_class + "_training_size"
    d_fig["x_title"] = "Number of training events for " + varying_class + " events"
    d_fig["y_title"] = "AUC Score"
    plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_size, list_auc_train, list_auc_test, d_fig)



def plot_boosting_round_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, folds, list_boosting_rounds, bool_weighted, list_variables):

    """
    Detail:
        Plot influence of number of boosting rounds on the AUC
    
    Args:
        WIMP_mass (str)             = WIMP mass
        bolo_name (str)             = bolometer name 
        d_event_dir (dict)          = dict to get the proper directory of each event class 
        exposure                    = exposure in days
        analysis_type (str)         = type of analysis (which box cut)
        param (dict)                = list of xgboost parameters
        folds (int)                 = number of CV folds
        list_boosting_rounds (list) = list of number of boosting rounds
        bool_weighted (bool) = True if include weights in BDT
        list_variables (list) = list of BDT features
    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """
    
    list_auc_test  = []
    list_auc_train = []

    for num_rounds in list_boosting_rounds:
        
        script_utils.print_utility("CV with " +str(num_rounds)+ " boosting rounds")
        
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, "", d_event_dir.keys(), exposure, list_variables, datasplit=0.1)
        
        AUC_train, AUC_test = [], []
        if bool_weighted:
            AUC_train, AUC_test = xgCV.CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        else:
            AUC_train, AUC_test= xgCV.CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        list_auc_train.append(AUC_train)
        list_auc_test.append(AUC_test)

        script_utils.print_utility("Mean AUC_train " +str(np.mean(AUC_train)))
        script_utils.print_utility("Mean AUC_test " +str(np.mean(AUC_test)))

    d_fig = {}
    if bool_weighted:
        d_fig["weight"] = "With_weights" 
    else:
        d_fig["weight"] = "No_weights"
    d_fig["fig_title"] = "boosting_rounds_influence"
    d_fig["x_title"] = "Number of boosting rounds"
    d_fig["y_title"] = "AUC Score"
    plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_boosting_rounds, list_auc_train, list_auc_test, d_fig)


def plot_depth_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_depth, bool_weighted, list_variables):

    """
    Detail:
        Plot influence of number of boosting rounds on the AUC
    
    Args:
        WIMP_mass (str)     = WIMP mass
        bolo_name (str)     = bolometer name 
        d_event_dir (dict)  = dict to get the proper directory of each event class 
        exposure            = exposure in days
        analysis_type (str) = type of analysis (which box cut)
        param (dict)        = list of xgboost parameters
        num_rounds (int)    = number of boosting rounds
        folds (int)         = number of CV folds
        list_depth (list)   = list of tree depth
        bool_weighted (bool) = True if include weights in BDT
        list_variables (list) = list of BDT features
    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """
    
    list_auc_test  = []
    list_auc_train = []


    for depth in list_depth:
        
        script_utils.print_utility("CV with depth of " +str(depth))
        
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, "", d_event_dir.keys(), exposure, list_variables, datasplit=0.1)
        
        param["bst:max_depth"] = depth
        AUC_train, AUC_test = [], []
        if bool_weighted:
            AUC_train, AUC_test = xgCV.CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        else:
            AUC_train, AUC_test = xgCV.CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        list_auc_train.append(AUC_train)
        list_auc_test.append(AUC_test)
        
        script_utils.print_utility("Mean AUC_train " +str(np.mean(AUC_train)))
        script_utils.print_utility("Mean AUC_test " +str(np.mean(AUC_test)))

    d_fig = {}
    if bool_weighted:
        d_fig["weight"] = "With_weights" 
    else:
        d_fig["weight"] = "No_weights"
    d_fig["fig_title"] = "depth_influence"
    d_fig["x_title"] = "Tree Max depth"
    d_fig["y_title"] = "AUC Score"
    plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_depth, list_auc_train, list_auc_test, d_fig)


def plot_eta_influence(WIMP_mass, bolo_name, d_event_dir, exposure, analysis_type, param, num_rounds, folds, list_eta, bool_weighted, list_variables):

    """
    Detail:
        Plot influence of number of boosting rounds on the AUC
    
    Args:
        WIMP_mass (str)     = WIMP mass
        bolo_name (str)     = bolometer name 
        d_event_dir (dict)  = dict to get the proper directory of each event class 
        exposure            = exposure in days
        analysis_type (str) = type of analysis (which box cut)
        param (dict)        = list of xgboost parameters
        num_rounds (int)    = number of boosting rounds
        folds (int)         = number of CV folds
        list_eta (list)   = list of tree eta
        bool_weighted (bool) = True if include weights in BDT
        list_variables (list) = list of BDT features
    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """
    
    list_auc_test  = []
    list_auc_train = []


    for eta in list_eta:
        
        script_utils.print_utility("CV with eta of " +str(eta))
        
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, "", d_event_dir.keys(), exposure, list_variables, datasplit=0.1)
        
        param["bst:eta"] = eta
        AUC_train, AUC_test = [], []
        if bool_weighted:
            AUC_train, AUC_test = xgCV.CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)
        else:
            AUC_train, AUC_test = xgCV.CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)

        list_auc_train.append(AUC_train)
        list_auc_test.append(AUC_test)
        
        script_utils.print_utility("Mean AUC_train " +str(np.mean(AUC_train)))
        script_utils.print_utility("Mean AUC_test " +str(np.mean(AUC_test)))

    d_fig = {}
    if bool_weighted:
        d_fig["weight"] = "With_weights" 
    else:
        d_fig["weight"] = "No_weights"
    d_fig["fig_title"] = "eta_influence"
    d_fig["x_title"] = "Tree learning rate"
    d_fig["y_title"] = "AUC Score"
    plotting_utility(WIMP_mass, bolo_name, analysis_type, folds, list_eta, list_auc_train, list_auc_test, d_fig)