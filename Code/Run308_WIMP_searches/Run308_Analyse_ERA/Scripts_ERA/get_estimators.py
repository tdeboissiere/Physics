#!/usr/bin/env python

from ROOT import TTree, TFile, TH2F, TH1D, TGraph, TH1F, TCanvas
import script_utils as script_utils
import sys
import numpy as np

def compute_coeff( FWHM1,  FWHM2) :
    """
    Compute the best combination of two channels
    The estimator MUST be w1 * Grandeur 1 + (1-w1) * Grandeur 2
    """ 
    s1=float(FWHM1)/2.3548
    s2=float(FWHM2)/2.3548

    w1=(s2*s2)/(s1*s1+s2*s2)
    return w1

def compute_resolution( FWHM1,  FWHM2) :
    """
    Compute the resolution (i.e FWHM not sigma) of the best combination of two channels
    """
    s1=float(FWHM1)/2.3548
    s2=float(FWHM2)/2.3548

    w1=(s2*s2)/(s1*s1+s2*s2)

    res=np.sqrt(pow(w1*s1,2)+pow((1-w1)*s2,2))
    return 2.3548*res


def get_estimators(bolo_name):

    """
    Creates a .txt file with the prelim cut list of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A .txt file '_prelim_TCuts.txt'
    it contains the prelim cuts fin TCut form for each bolometer and is updated when a new one is added
    """

    #Access polar and FWHM info. 
    FWHM_path_name =script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    FWHM_file_name = bolo_name + "_FWHM.txt" 
    file_FWHM = script_utils.open_text_file(FWHM_path_name, FWHM_file_name, "r")
    list_FWHM = file_FWHM.readlines()[1].rstrip().split(",")

    FWIA,FWIB,FWIC,FWID = list_FWHM[0], list_FWHM[1], list_FWHM[2], list_FWHM[3]
    FWC1_ERA,FWC2_ERA= list_FWHM[4], list_FWHM[5]
    VFID, VET = list_FWHM[6], list_FWHM[7]

    #Compute the estimators
    wfid=str(compute_coeff(FWIB,FWID))[:5]
    wS1=str(compute_coeff(FWIB,FWIA))[:5]
    wS2=str(compute_coeff(FWID,FWIC))[:5]
    wheat=str(compute_coeff(FWC1_ERA,FWC2_ERA))[:5]
    
    e_fid=wfid+"*EIB+(1-"+wfid+")*EID"
    e_S1=wS1+"*EIB+(1-"+wS1+")*EIA"
    e_S2=wS2+"*EID+(1-"+wS2+")*EIC"
    e_heat=wheat+"*EC1+(1-"+wheat+")*EC2"

    e_rec_fid="("+e_heat+")*(1+"+VFID+"/3.)-("+e_fid+")*0.3333*"+VFID
    e_Q_fid="("+e_fid+")/("+e_rec_fid+")"

    e_rec_S1="("+e_heat+")*(1+"+VFID+"/3.)-("+e_S1+")*0.3333*"+VET
    e_Q_S1="("+e_S1+")/("+e_rec_S1+")"

    e_rec_S2="("+e_heat+")*(1+"+VFID+"/3.)-("+e_S2+")*0.3333*"+VET
    e_Q_S2="("+e_S2+")/("+e_rec_S2+")"       


    # Create file directory if it does not exist and define the file name
    estimator_file_path= script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')
    estimator_file_name= bolo_name + "_estimators.txt"

    estimator_file = script_utils.open_text_file(estimator_file_path, estimator_file_name, "w")
    estimator_file.write("FID," + e_fid + '\n')
    estimator_file.write("S1," + e_S1 + '\n')
    estimator_file.write("S2," + e_S2 + '\n')
    estimator_file.write("HEAT," + e_heat + '\n')

    estimator_file.write("ER_FID," + e_rec_fid + '\n')
    estimator_file.write("Q_FID," + e_Q_fid + '\n')

    estimator_file.write("ER_S1," + e_rec_S1 + '\n')
    estimator_file.write("Q_S1," + e_Q_S1 + '\n')

    estimator_file.write("ER_S2," + e_rec_S2 + '\n')
    estimator_file.write("Q_S2," + e_Q_S2 + '\n')

    estimator_file.close()