from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl

def get_ps_histogram(h, bin_X, min_X, max_X):

    """Get the p_s histogram
    
    Detail:
        For each bin, compute the integral from this
        bin to the last.
        This gives the p-value 

    Args:
        h  = (TH1F) histogram of the test statistic
        
    Returns:
        void

    Raises:
        void
    """

    hcumu = TH1F("hcumu", "hcumu", bin_X, min_X, max_X)
    for i in range(1, bin_X+1):
    	hcumu.SetBinContent(i, h.Integral(i, bin_X)/float(h.Integral()))

    # cc = TCanvas("cc", "cc")
    # gPad.SetLogy()
    # hcumu.Draw() 
    # raw_input()
    return hcumu

def get_stat_test_histogram(file_name,  bin_X, min_X, max_X):

    """Get the test stat histogram
    
    Detail:

    Args:
        file_name  = (str) the file name from which to extract the test stat
        
    Returns:
        void

    Raises:
        void
    """

    arr = np.loadtxt(file_name)
    htest_stat = TH1F("htest_stat", "htest_stat", bin_X, min_X, max_X)
    PyRPl.fill_TH1(htest_stat, arr)

    # cc = TCanvas("cc", "cc")
    # gPad.SetLogy()
    # htest_stat.Draw() 
    # raw_input()
    return htest_stat


def get_stat_test_obs_median_and_mean(file_name):

    """
    
    Detail: Use simulations with 0 signal 
    Load values of test stat for sigma = XX and simu with 0 signal

    Args:
        file_name  = (str) the file name from which to extract the test stat
        
    Returns:
        void

    Raises:
        void
    """

    arr = np.loadtxt(file_name)
    print np.mean(arr), np.median(arr)

bin_X,min_X,max_X = 100, 0, 10
# Nsignal = "100"
# file_name1 = "./Text_files/test_stat_Nsignal" + Nsignal + ".txt"
# h1 = get_stat_test_histogram(file_name1, bin_X, min_X, max_X)
# hcumu1 = get_ps_histogram(h1, bin_X, min_X, max_X)

# Nsignal = "200"
# file_name2 = "./Text_files/test_stat_Nsignal" + Nsignal + ".txt"
# h2 = get_stat_test_histogram(file_name2, bin_X, min_X, max_X)
# hcumu2 = get_ps_histogram(h2, bin_X, min_X, max_X)

# hcumu1.SetLineColor(kRed)
# hcumu1.Draw()
# hcumu2.Draw("same")
# raw_input()

Nsignal = "125"
get_stat_test_obs_median_and_mean("./Text_files/test_stat_obs_Nsignal" + Nsignal + ".txt")