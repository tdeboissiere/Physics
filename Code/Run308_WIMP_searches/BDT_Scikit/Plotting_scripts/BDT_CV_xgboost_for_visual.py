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
import script_utils as script_utils


def CV_xgboost_unweighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_round, folds):

    """
    Detail:
        Launch CV with xgboot, no weights

    Args:
        d_train (dict)                      = dict with training data
        WIMP_mass (str)                     = WIMP mass
        bolo_name (str)                     = bolometer name 
        analysis_type (str)                 = type of analysis (which box cut)
        param (dict)                        = list of xgboost parameters
        num_round (int)                     = list of boosting rounds
        fold (int)                          = number of CV folds

    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """


    # Load training data
    X, weights, labels = dp.prepare_data_for_scikit(d_train)

    # Initialise CV
    kf = cv.KFold(labels.size, n_folds=folds)
    kf_index = 1
    # Initialise AUC
    AUC_train, AUC_test, AUC_train2, AUC_test2 = [], [], [], []

    for train_indices, test_indices in kf:

        script_utils.print_utility("CV round " + str(kf_index))

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = labels[train_indices], labels[test_indices]

        # Process weights for AUC computation
        w_train, w_test = weights[train_indices], weights[test_indices]
        w_train *= (sum(weights) / sum(w_train))
        w_test  *= (sum(weights) / sum(w_test))

        # sum_wpos = sum(w_train[y_train == 1])
        # sum_wneg = sum(w_train[y_train == 0])

        xgmat = xgb.DMatrix(X_train, label=y_train) 
        xgmat_test = xgb.DMatrix(X_test)

        plst = param.items()

        watchlist = [] #[(xgb.DMatrix(X_test, label=y_test, weight=w_test), 'eval'), (xgmat, "train")]
        bst = xgb.train(plst, xgmat, num_round, watchlist)

        y_pred_test = bst.predict(xgmat_test)
        y_pred_train = bst.predict(xgmat)

        AUC_train.append(roc_auc_score(y_train, y_pred_train, sample_weight = w_train))
        AUC_test.append(roc_auc_score(y_test, y_pred_test, sample_weight = w_test))

        # AUC_train2.append(roc_auc_score(y_train, y_pred_train))
        # AUC_test2.append(roc_auc_score(y_test, y_pred_test))

        kf_index+=1

    return AUC_train, AUC_test# , AUC_train2, AUC_test2


def CV_xgboost_weighted(d_train, WIMP_mass, bolo_name, analysis_type, param, num_round, folds):

    """
    Detail:
        Launch CV with xgboot, with weights

    Args:
        d_train (dict)                      = dict with training data
        WIMP_mass (str)                     = WIMP mass
        bolo_name (str)                     = bolometer name 
        analysis_type (str)                 = type of analysis (which box cut)
        param (dict)                        = list of xgboost parameters
        num_round (int)                     = list of boosting rounds
        fold (int)                          = number of CV folds

    Returns:
        AUC_train, AUC_test = list of AUC score at each CV

    Raises:
        void     
    """


    # Load training data
    X, weights, labels = dp.prepare_data_for_scikit(d_train)

    # Initialise CV
    kf = cv.KFold(labels.size, n_folds=folds)
    kf_index = 1
    # Initialise AUC
    AUC_train, AUC_test = [], []

    for train_indices, test_indices in kf:

        script_utils.print_utility("CV round " + str(kf_index))

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = labels[train_indices], labels[test_indices]

        #Include weights in model
        w_train, w_test = weights[train_indices], weights[test_indices]
        w_train *= (sum(weights) / sum(w_train))
        w_test  *= (sum(weights) / sum(w_test))

        sum_wpos = sum(w_train[y_train == 1])
        sum_wneg = sum(w_train[y_train == 0])

        xgmat = xgb.DMatrix(X_train, label=y_train, weight=w_train)
        xgmat_test = xgb.DMatrix(X_test, weight = w_test)

        # scale weight of positive examples
        param['scale_pos_weight'] = sum_wneg / sum_wpos
        plst = param.items()

        watchlist = [] #[(xgb.DMatrix(X_test, label=y_test, weight=w_test), 'eval'), (xgmat, "train")]
        bst = xgb.train(plst, xgmat, num_round, watchlist)

        y_pred_test = bst.predict(xgmat_test)
        y_pred_train = bst.predict(xgmat)

        AUC_train.append(roc_auc_score(y_train, y_pred_train, sample_weight = w_train) )
        AUC_test.append(roc_auc_score(y_test, y_pred_test, sample_weight = w_test) )

        kf_index+=1

    return AUC_train, AUC_test
