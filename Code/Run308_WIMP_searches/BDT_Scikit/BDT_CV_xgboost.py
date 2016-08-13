import csv
import sys
sys.path.append('../../python/')
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

def AMS(s, b):
    '''
    Approximate median significance:
        s = true positive rate
        b = false positive rate
    '''
    assert s >= 0
    assert b >= 0
    bReg = 10.
    return np.sqrt(2.0 * ((s + b + bReg) * np.log(1 + s / (b + bReg)) - s))


def get_rates(prediction, solution, weights):
    '''
    Returns the true and false positive rates.
    This assumes that:
        label 's' corresponds to 1 (int)
        label 'b' corresponds to 0 (int)
    '''
    assert prediction.size == solution.size
    assert prediction.size == weights.size

    # Compute sum of weights for true and false positives
    truePos  = sum(weights[(solution == 1) * (prediction == 1)])
    falsePos = sum(weights[(solution == 0) * (prediction == 1)])

    return truePos, falsePos


def plot_corr_and_hist(d_train, WIMP_mass, exposure, bolo_name, analysis_type):

    # X, weights, labels = dp.prepare_data_for_scikit(d_train)
    # scaler = StandardScaler()
    # scaler.fit(X)


    #Get scaling dict for data visualisation
    d_scaling = BDT_fh.open_MVA_scaling_file(bolo_name, analysis_type, "")

    list_ev = ["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly", "WIMP_mass_" + WIMP_mass]
    list_bg = ["S1Pb", "S2Pb", "S1Beta", "S2Beta", "S1Gamma", "S2Gamma", "FidGamma", "heatonly"]
    df = pd.concat( [d_train[event_type] for event_type in d_train.keys() ], ignore_index = True )
    d_color = {"S1Pb":"peru", "S2Pb":"saddlebrown", "S1Beta":"olive", "S2Beta":"forestgreen", "S1Gamma":"royalblue", "S2Gamma":"dodgerblue", 
                "FidGamma":"deepskyblue", "heatonly":"red", "WIMP_mass_" + WIMP_mass:"dimgray", "neutron":"magenta"}
    d_index = {}
    d_weight = {}
    for event_type in d_train.keys():
        d_index[event_type] = df.ix[df["EventID"]==event_type]
        if "WIMP" not in event_type:
            d_weight[event_type] = float(d_scaling["prop_" + event_type])*float(d_scaling["exp_per_day"])*exposure*np.ones(len(d_index[event_type]))/len(d_index[event_type])
        else:
            d_weight[event_type] = 8000*np.ones(len(d_index[event_type]))/len(d_index[event_type])

    list_weights = [d_weight[event_type] for event_type in list_ev]
    list_color = [d_color[event_type] for event_type in list_ev]
    gs = gridspec.GridSpec(3,3)

    for index, val in enumerate(list(df.columns.values)[:-3]):
    # for index, val in enumerate(list(df.columns.values)[:4]):
        if val == "prob":
            plt.subplot(gs[index]).hist([d_index[event_type][val] for event_type in list_ev], weights=list_weights, color=list_color , stacked=True, histtype="stepfilled", bins=200000)
        else:
            plt.subplot(gs[index]).hist([d_index[event_type][val] for event_type in list_ev], weights=list_weights, color=list_color , stacked=True, histtype="stepfilled", bins=200, range = (0,5))
        plt.subplot(gs[index]).set_yscale("log")
        plt.subplot(gs[index]).set_xlabel(val)
    plt.tight_layout()
    plt.show()
    raw_input()
    plt.savefig("test.png")

    # gs2D = gridspec.GridSpec(3,4)
    # for ix, valx in enumerate(list(df.columns.values)[:-3]):
    #     fig_index = 0
    #     fig = plt.figure(figsize=(10, 6.25))

    #     for iy, valy in enumerate(list(df.columns.values)[:-3]):
    #         if ix!= iy:
    #             print fig_index
    #             for event_type in ["heatonly"]:
    #                 plt.subplot(gs2D[fig_index]).scatter(d_index[event_type][valx], d_index[event_type][valy], color="r", s=0.01 )
    #             plt.subplot(gs2D[fig_index]).scatter(d_index["WIMP_mass_" + WIMP_mass][valx], d_index["WIMP_mass_" + WIMP_mass][valy], color="k", s=0.1 )
    #             plt.subplot(gs2D[fig_index]).set_xlabel(valx)
    #             plt.subplot(gs2D[fig_index]).set_ylabel(valy)
    #             fig_index+=1

    #     fig.tight_layout()
    #     plt.savefig("./Figures/" + valx + ".png")
    #     plt.close(fig)


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
    scaler = StandardScaler()
    scaler.fit(X_ref)

    # X =scaler.transform(X)
    # 
    #Cutoff for AMS
    npoints = 40
    cutoffs = sp.linspace(0.05, 0.45, npoints)
    all_AMS = {}
    for curr in range(npoints):
        all_AMS[curr] = []

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

        list_AUC.append(roc_auc_score(y_test, y_pred_test, sample_weight = w_test))
        kf_index+=1

        print bst.get_fscore()

        # Compute AMS
        res  = [(i, y_pred_test[i]) for i in xrange(len(y_pred_test))]
        # print "res", res[:10]
        rorder = {}
        # print sorted(res[:10], key = lambda x:-x[1])
        for k, v in sorted(res, key = lambda x:-x[1]):
            rorder[k] = len(rorder) + 1

        # Explore changing threshold_ratio and compute AMS
        best_AMS = -1.
        best_thres = 0.0
        for curr, threshold_ratio in enumerate(cutoffs):
            y_pred = sp.zeros(len(y_pred_test))
            ntop = int(threshold_ratio * len(rorder))
            for k, v in res:
                if rorder[k] <= ntop:
                    y_pred[k] = 1

            truePos, falsePos = get_rates(y_pred, y_test, w_test)
            this_AMS = AMS(truePos, falsePos)
            all_AMS[curr].append(this_AMS)
            if this_AMS > best_AMS:
                best_AMS = this_AMS
                best_thres = threshold_ratio
        print "Best AMS = %f at %.2f"%(best_AMS,best_thres)

    print "------------------------------------------------------"
    max_cut, max_ams, max_std = 0,0,0 
    for curr, cut in enumerate(cutoffs):
        if sp.mean(all_AMS[curr]) > max_ams :
            max_cut = cut 
            max_ams = sp.mean(all_AMS[curr])
            max_std = sp.std(all_AMS[curr])
        # print "Thresh = %.2f: AMS = %.4f, std = %.4f" % \
        #     (cut, sp.mean(all_AMS[curr]), sp.std(all_AMS[curr]))
    print "Thresh = %.2f: AMS = %.4f, std = %.4f" % \
        (max_cut, max_ams, max_std)
    print "------------------------------------------------------"

        # fpr, tpr, thresholds = roc_curve(y_test, y_pred_test)
        # fpr_train, tpr_train, thresholds_train = roc_curve(y_train, y_pred_train)
        # plt.plot(fpr, tpr, "-", label = "Fold " + str(kf_index -1) + " test sample")
        # plt.plot(fpr_train, tpr_train, "--", label = "Fold " + str(kf_index -1) + " train sample")

    # plt.xlabel("False Postives", fontsize = 19)
    # plt.ylabel("True Postives", fontsize = 19)
    # plt.legend(loc = "lower right")
    # plt.show() 
    # raw_input()


    return np.mean(np.array(list_AUC)), np.std(np.array(list_AUC))


        # fpr, tpr, thresholds = roc_curve(y_test, y_out)
        # fpr_train, tpr_train, thresholds_train = roc_curve(y_train, y_out_train)
        # mean_tpr += interp(mean_fpr, fpr, tpr)
        # mean_tpr[0] = 0.0
        # roc_auc = auc(fpr, tpr)
        # plt.plot(fpr, tpr, "-", label = "Fold " + str(folds) + " Num_rounds " + str(num_round) + " max_depth " + str(param['bst:max_depth']))
        # plt.plot(fpr_train, tpr_train, "--", label = "Fold " + str(folds) + " Num_rounds " + str(num_round) + " max_depth " + str(param['bst:max_depth']))
        # plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))
        # plt.show() 
        # raw_input()

        # res  = [(i, y_out[i]) for i in xrange(len(y_out))]
        # rorder = {}
        # for k, v in sorted(res, key = lambda x:-x[1]):
        #     rorder[k] = len(rorder) + 1

        # print res
        # raw_input()





    # plt.xlim([-0.05, 1.05])
    # plt.ylim([-0.05, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title('Receiver operating characteristic example')
    # plt.legend(loc="lower right")
    # plt.show()
    # raw_input()

    # print "------------------------------------------------------"
    # for curr, cut in enumerate(cutoffs):
    #     print "Thresh = %.2f: AMS = %.4f, std = %.4f" % \
    #         (cut, sp.mean(all_AMS[curr]), sp.std(all_AMS[curr]))
    # print "------------------------------------------------------"



def launch_CV():

    """
    Detail:
        Main script to launch classification

    Args:
        bolo_name (str)        = bolometer name
        analysis_type (str)    = type of analysis (which cuts)
        bool_train (bool)      = boolean to pick test/training sample
        list_event_type (list) = list of event class

    Returns:
        void

    Raises:
        void     
    """

    # plt.ion()

    bolo_name = "FID837"
    data_dir = "/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/"
    analysis_type = "ana_0.5_0_5"
    exposure = 65.
    bin_X, min_X, max_X = 100, -20, 20
    MVA_tag = ""
    # list_mass = ["3", "4", "5", "6", "7", "10", "25"]
    list_mass = ["4"]
    # list_mass = ["10"]

    # #Prepare data
    # d_data_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
    #                    "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "realdata": "realdata"}
    # for WIMP_mass in list_mass:
    #   d_data_dir["WIMP_mass_" + WIMP_mass] = "WIMP"

    # for event_type in d_data_dir.keys():
    #   #O for training, 1 for test
    #   dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 0)
    #   dp.root_to_csv(bolo_name, data_dir, analysis_type, event_type, d_data_dir, 1)

    list_variables = ["EC", "ECdiff", "EFID", "EIA", "EIB", "EIC", "EID", "max_ion", "test", "ER"]
    list_variables = ["EC", "EFID", "EFIDdiff",  "max_vet", "max_ion", "test", "ER"]
    # list_variables = ["EC","ECdiff", "EFID", "EIA", "EIB", "EIC", "EID", "test", "ER"]
    list_variables = ["EIA", "EIB", "EIC", "EID", "EC1", "EC2", "test", "ER"] # 0.99855
    list_variables = ["EIA", "EIB", "EIC", "EID", "EC", "ECdiff", "test"] #, "EFIDprod"] #, "ER", "ECdiff"]
    list_variables = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID", "test"] #, "test", "HR"] #, "test"]
    # list_variables = ["EIA" , "EIC", "EC", "test"] #, "ER", "ECdiff"]
    list_variables = ["EC1","EC2", "EIA", "EIC", "EFID", "testEC", "prod"]
    list_variables = ["EC1","EC2", "EFID", "to1_3"]
    list_variables = ["EC1","EC2", "EIA", "EIC", "EFID", "testEC", "to1_3"]
    # list_variables = ["EC1", "EC2", "test"]
    # list_variables = ["EC1", "EC2", "EIA", "EIB", "EIC", "EID"] #, "test", "HR"] #, "test"]
    list_variables = ["EC1","EC2", "EIA", "EIC", "EFID"]
    list_variables = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test", "prod", "to1_3", "EFIDdiff", "testEC"]
    list_variables = ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "EFIDdiff", "ECdiff", "test"]
    list_variables = ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test"] #, "to1_3", "sum_ion"]

    # 3 GeV : 5, 100, 0.1, ["EC1", "EC2", "EIA", "EIC", "EFID"]
    # 4 GeV : 5, 100, 0.1, ["EC1", "EC2", "EIA", "EIC", "EFID"]
    # 5 GeV : 5, 100, 0.1,["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
    # 6 GeV : 5, 100, 0.1,["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
    # 7 GeV : 5, 100, 0.1, ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]
    # 10 GeV : 7, 100, 0.1, ["EC1","EC2", "EIA", "EIC", "EIB", "EID", "test"]  
    # 25 GeV : 7, 100, 0.1, ["EC1","EC2", "EIA", "EIB", "EIC", "EID", "test", "prod"]

    #Loop over masses
    for WIMP_mass in list_mass:
        d_event_dir = {"S1Pb":"Beta_and_Pb", "S2Pb":"Beta_and_Pb", "S1Beta":"Beta_and_Pb", "S2Beta":"Beta_and_Pb",
                         "S1Gamma":"Gamma", "S2Gamma":"Gamma", "FidGamma":"Gamma", "heatonly":"Heatonly", "WIMP_mass_" + WIMP_mass: "WIMP"}

        #Load data
        d_train = dp.get_data_array(bolo_name, 0, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit=0.3)
        d_test  = dp.get_data_array(bolo_name, 1, analysis_type, MVA_tag, d_event_dir.keys(), exposure, list_variables, datasplit=0.3)
        # d_eval  = dp.get_eval_array(bolo_name, analysis_type)

        #Plot var 
        # plot_corr_and_hist(d_test, WIMP_mass, exposure, bolo_name, analysis_type)

        # setup parameters for xgboost
        param = {}
        param['objective'] = 'binary:logistic'
        param['eval_metric'] = 'auc'
        param['silent'] = 1

        folds = 5
        AUC=-1

        all_etas = [0.1]
        all_subsamples = [0.9]
        all_depth = [5]
        nums_rounds = [100]
        e_s_m = list(itertools.product(all_etas,all_subsamples,all_depth,nums_rounds))
        for e,s,m,r in e_s_m:
            param['bst:eta'] = e
            param['bst:subsample'] = s
            param['bst:max_depth'] = m
            AUC_temp, AUC_std = CV_xgboost_weighted(d_train, d_test, d_event_dir, WIMP_mass, bolo_name, analysis_type, param, r, folds)
            if AUC_temp > AUC:
                AUC = AUC_temp
                print "AUC: ", AUC, "+/-", "%.2E" % AUC_std, "Eta:", param["bst:eta"], "Max_depth:", param["bst:max_depth"], "Boosting_rounds:", r, "Subsample:", param["bst:subsample"]


launch_CV()