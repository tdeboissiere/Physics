import csv
import sys
sys.path.append('../')
import numpy as np
import scipy as sp
import xgboost as xgb
import sklearn.cross_validation as cv
import itertools
import pandas as pd
from pandas.tools.plotting import scatter_matrix
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Normalizer, Binarizer
from sklearn.preprocessing import Imputer
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, AdaBoostClassifier
import data_preparation as dp 
import data_visualisation as dv 
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pylab as plt
from sklearn.decomposition import PCA, KernelPCA, RandomizedPCA
from sklearn.lda import LDA
import matplotlib.gridspec as gridspec
import BDT_file_handler as BDT_fh
import script_utils as script_utils
import pickle
import prettify_matplotlib as pret

def CV_xgboost_unweighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_round, folds):

    '''
    Cross validation for XGBoost performance (unweighted)
    '''
    # Load training data
    X, weights, labels = dp.prepare_data_for_scikit(d_train)
    X_ref, weights_ref, labels_ref = dp.prepare_data_for_scikit(d_test)
    scaler = StandardScaler()
    scaler.fit(X_ref)

    # X =scaler.transform(X)


    # Cross validate
    kf = cv.KFold(labels.size, n_folds=folds)
    kf_index = 1
    # These are the cutoffs used for the XGBoost predictions
    list_AUC = []
    for train_indices, test_indices in kf:

        script_utils.print_utility("CV round " + str(kf_index))

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = labels[train_indices], labels[test_indices]
        w_train, w_test = weights[train_indices], weights[test_indices]

        # Rescale weights so that their sum is the same as for the entire training set
        w_train *= (sum(weights) / sum(w_train))
        w_test  *= (sum(weights) / sum(w_test))

        xgmat = xgb.DMatrix(X_train, label=y_train)

        # scale weight of positive examples
        plst = param.items()#+[('eval_metric', 'ams@0.15')]

        watchlist = [] #[(xgb.DMatrix(X_test, label=y_test, weight=w_test), 'eval'), (xgmat, "train")]
        bst = xgb.train(plst, xgmat, num_round, watchlist)

        # Construct matrix for test set
        xgmat_test = xgb.DMatrix(X_test)
        y_pred_test = bst.predict(xgmat_test)

        list_AUC.append(roc_auc_score(y_test, y_pred_test, sample_weight = w_test))
        kf_index+=1

        print bst.get_fscore()

    return np.mean(np.array(list_AUC)), np.std(np.array(list_AUC))

def CV_xgboost_weighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_round, folds):

    '''
    Cross validation for XGBoost performance (weighted)
    '''
    # Load training data
    X, weights, labels = dp.prepare_data_for_scikit(d_train)
    X_ref, weights_ref, labels_ref = dp.prepare_data_for_scikit(d_test)

    # Cross validate
    kf = cv.KFold(labels.size, n_folds=folds)
    kf_index = 1
    # These are the cutoffs used for the XGBoost predictions
    list_AUC_train = []
    list_AUC_test = []

    for train_indices, test_indices in kf:

        script_utils.print_utility("CV round " + str(kf_index))

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = labels[train_indices], labels[test_indices]
        w_train, w_test = weights[train_indices], weights[test_indices]

        # Rescale weights so that their sum is the same as for the entire training set
        w_train *= (sum(weights) / sum(w_train))
        w_test  *= (sum(weights) / sum(w_test))

        sum_wpos = sum(w_train[y_train == 1])
        sum_wneg = sum(w_train[y_train == 0])

        xgmat = xgb.DMatrix(X_train, label=y_train, weight=w_train)

        # scale weight of positive examples
        param['scale_pos_weight'] = sum_wneg / sum_wpos
        plst = param.items()#+[('eval_metric', 'ams@0.15')]

        watchlist = [] #[(xgb.DMatrix(X_test, label=y_test, weight=w_test), 'eval'), (xgmat, "train")]
        bst = xgb.train(plst, xgmat, num_round, watchlist)

        # Construct matrix for test set
        xgmat_test = xgb.DMatrix(X_test, weight = w_test)
        y_pred_test = bst.predict(xgmat_test)
        y_pred_train = bst.predict(xgmat)

        list_AUC_test.append(roc_auc_score(y_test, y_pred_test, sample_weight = w_test))
        list_AUC_train.append(roc_auc_score(y_train, y_pred_train, sample_weight = w_train))
        kf_index+=1

        print bst.get_fscore()

    return list_AUC_train, list_AUC_test


def launch_CV():

    """
    Detail:
        Main script to launch cross validation

    Args:

    Returns:
        void

    Raises:
        void     
    """

    bolo_name = "FID837"
    data_dir = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
    analysis_type = "ana_1.5_0_5"
    exposure = 65.
    bin_X, min_X, max_X = 100, -20, 20
    MVA_tag = ""
    list_mass = ["5"]
    d_list_variables = {}

    d_list_variables[0] = ["EC", "EFID"]
    d_list_variables[1] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID"] 
    d_list_variables[2] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "ER"] 
    d_list_variables[3] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test"] 
    d_list_variables[4] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test", "HR"] 
    d_list_variables[5] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test", "ER", "HR"] 

    #Loop over masses
    for WIMP_mass in list_mass:
        d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                         "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

        d_AUC = {}

        for i_var in range(6):

            #Load data
            d_train = dp.get_data_array(bolo_name, 0, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[i_var], datasplit=.4)
            d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, d_list_variables[i_var], datasplit=1)

            # setup parameters for xgboost
            param = {}
            param['objective'] = 'binary:logistic'
            param['eval_metric'] = 'auc'
            param['silent'] = 1
            param['bst:eta'] = 0.1
            param['bst:subsample'] = 0.9
            param['bst:max_depth'] = 3
            num_rounds = 100

            folds = 5

            d_AUC[str(i_var)+"train"], d_AUC[str(i_var)+"test"] = CV_xgboost_weighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, num_rounds, folds)

        pickle_dir = script_utils.create_directory("./Pickle_files/" + bolo_name + "/" + analysis_type + "/")
        with open(pickle_dir + bolo_name + "_AUC_CV_various_features_mass_" + str(WIMP_mass) + ".pickle", 'wb') as handle:
          pickle.dump(d_AUC, handle)

def launch_boxplot():

    """
    Detail:
        Plot the boxplot

    Args:

    Returns:
        void

    Raises:
        void     
    """

    # plt.ion()

    bolo_name = "FID837"
    analysis_type = "ana_1.5_0_5"
    list_mass = ["5"]
    d_list_variables = {}

    plt.ion()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    d_list_variables[0] = ["EC", "EFID"]
    d_list_variables[1] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID"] 
    d_list_variables[2] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "ER"] 
    d_list_variables[3] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test"] 
    d_list_variables[4] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test", "HR"] 
    d_list_variables[5] = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test", "ER", "HR"] 

    for WIMP_mass in list_mass:
        d_AUC = {}
        pickle_dir = script_utils.create_directory("./Pickle_files/" + bolo_name + "/" + analysis_type + "/")
        with open(pickle_dir + bolo_name + "_AUC_CV_various_features_mass_" + str(WIMP_mass) + ".pickle", 'rb') as handle:
            d_AUC = pickle.load(handle)    

        fig = plt.figure()
        ax = plt.axes()
        plt.hold(True)

        d_pos = {0:[1,2], 1:[4,5], 2:[7,8], 3:[10,11], 4:[13,14], 5:[16,17], 6:[19,20]}

        nchoice = 6
        iter_list = range(nchoice)
        iter_list = [0,1,3,4]

        for index, i in enumerate(iter_list):

            # Add some noise for better visualisation
            temp_l = d_AUC[str(i)+"train"]
            dev = np.std(temp_l)
            for k in range(len(temp_l)):
                temp_l[k]= temp_l[k] + 2*np.random.normal(0,dev)

            d_AUC[str(i)+"train"] = temp_l

            bp = plt.boxplot(d_AUC[str(i)+"train"],0, "", positions = [d_pos[index][0]], widths = 0.6, whis= "range", patch_artist=True)
            pret.prettify_boxplot(bp, '#1b9e77', "k", 1.5, "-")
            bp = plt.boxplot(d_AUC[str(i)+"test"],0, "", positions = [d_pos[index][1]], widths = 0.6, whis= "range", patch_artist=True)
            pret.prettify_boxplot(bp, 'r', "k", 1.5, "-")



        plt.xlim(0,3*len(iter_list))
        plt.ylabel("AUC score", fontsize = 16, labelpad = 20)
        # plt.ylim(0,1.2)
        ax.set_xticklabels(["Sel" + str(k) for k in range(len(iter_list))])
        ax.set_xticks([np.mean(d_pos[k]) for k in range(len(iter_list))])


        hB, = plt.plot([1,1],color ='#1b9e77', linestyle = "-", linewidth = 2)
        hR, = plt.plot([1,1],'r-', linewidth = 2)
        plt.legend((hB, hR),('Train', 'Test'), loc="best")
        hB.set_visible(False)
        hR.set_visible(False)

        plt.show()
        raw_input()

        fig_dir = script_utils.create_directory("./Figures/" + bolo_name + "/" + analysis_type + "/")
        plt.savefig(fig_dir + bolo_name + "_feature_selection.png")


# launch_CV()
launch_boxplot()

