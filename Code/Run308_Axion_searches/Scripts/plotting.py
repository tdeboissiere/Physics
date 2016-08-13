import utils_ana_jules as ut
import script_utils as script_utils
from ROOT import *
import PyROOTPlots as PyRPl
import numpy as np
import os

def choose_elec_recoil_cut(bolo_name):
    """Make plots to choose elec recoil cuts
    
    Detail:
        void

    Args:
        bolo_name  = (str) bolo name
        
    Returns:
        void

    Raises:
        void
    """
    
    #Get data tree
    t, f = PyRPl.open_ROOT_object("../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root", "t_axion")
    #Create hist2 and hist file
    h1Dcut  = TH1F("h1Dcut" , "h1Dcut", 100, 0, 15)
    h1Dfull = TH1F("h1Dfull" , "h1Dfull", 100, 0, 15)
    h2Dcut  = TH2F("h2Dcut" , "h2Dcut", 100, 0, 15, 100, 0, 15)
    h2Dfull = TH2F("h2Dfull" , "h2Dfull", 100, 0, 15, 100, 0, 15)
    h1Dcut.SetLineColor(kRed)
    h2Dcut.SetMarkerColor(kRed)

    cut = "EIA<2*FWIA/2.3548 && EIC<2*FWIC/2.3548 && abs(EC-EFID)<1.3 && EC>0.9 && EFID>1"
    # cut = "EIA<2*FWIA/2.3548 && EIC<2*FWIC/2.3548 && abs(EC-EFID)<0.5"

    t.Project("h1Dcut", "Ebest", cut)
    t.Project("h1Dfull", "Ebest")
    t.Project("h2Dcut", "EFID:EC", cut)
    t.Project("h2Dfull", "EFID:EC")

    cc = TCanvas("cc", "cc")
    h2Dfull.Draw()
    h2Dcut.Draw("same")

    cc2 = TCanvas("cc2", "cc2")
    gPad.SetLogy()
    h1Dfull.Draw()
    h1Dcut.Draw("same")


    raw_input()

def check_elec_recoil_bckgmodel(list_bolo_name, mass="0"):
    """Check the background model works ok
    
    Detail: 
        void

    Args:
        list_bolo_name = (str) list of bolo names
        
    Returns:
        void

    Raises:
        void
    """

    #Individual plots
    for bolo_name in list_bolo_name :

        hdata, fdata = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
        fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_perkeV.root", bolo_name + "_FidGamma")
        fmodel_noeff = file_model.Get(bolo_name + "_FidGamma_noeff")
        fmodel_noeff.SetLineStyle(7)

        cc = TCanvas("cc", "cc")
        PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
        hdata.Draw()
        fmodel.Draw("same")
        fmodel_noeff.Draw("same")
        gPad.SetLogy()
        cc.Print("../Figures/Data_and_model/Log_scale/" + bolo_name + "_elec_rec_perkeV.png")
        del cc

        cc = TCanvas("cc", "cc")
        PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
        hdata.Draw()
        fmodel.Draw("same")
        fmodel_noeff.Draw("same")
        cc.Print("../Figures/Data_and_model/Linear_scale/" + bolo_name + "_elec_rec_perkeV.png")
        del cc

    #Stacked plt
    hdata, fdata = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/Stacked_spec_perkeV.root", "hstacked")
    fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/Stacked_model_perkeV.root", "model_stacked")

    cc = TCanvas("cc", "cc")
    PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
    hdata.Draw()
    fmodel.Draw("same")
    gPad.SetLogy()
    cc.Print("../Figures/Data_and_model/Stacked_elec_rec_perkeV_log_scale.png")
    del cc

    cc = TCanvas("cc", "cc")
    PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
    hdata.Draw()
    fmodel.Draw("same")
    cc.Print("../Figures/Data_and_model/Stacked_elec_rec_perkeV_linear_scale.png")
    del cc

    #Stacked plt with limit
    hdata, fdata = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/Stacked_spec_perkeV.root", "hstacked")
    fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/Stacked_model_perkeV.root", "model_stacked")
    fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/Stacked_flux_perkeV.root", "flux_stacked_mass_" + str(mass))

    #Histograms parameters
    bin_X, min_X, max_X = hdata.GetNbinsX(), hdata.GetBinCenter(1) - 0.5*hdata.GetBinWidth(0), hdata.GetBinCenter(hdata.GetNbinsX()) + 0.5 * hdata.GetBinWidth(0)

    hbckg = TH1F("hbckg", "hbckg", bin_X, min_X, max_X)
    for i in range(1, bin_X+1):
        hbckg.SetBinContent(i, fmodel.Eval(hbckg.GetBinCenter(i)))

    hsignal = TH1F("hsignal", "hsignal", bin_X, min_X, max_X)
    for i in range(1, bin_X+1):
        #Signal at limit
        hsignal.SetBinContent(i, np.power(1.3E-11,4)*fflux.Eval(hsignal.GetBinCenter(i)))

    cc = TCanvas("cc", "cc")
    PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
    PyRPl.process_TH1(hbckg, X_title = "Best estimator (keVee)", Y_title = "c/keV", color = kGreen-3)
    PyRPl.process_TH1(hsignal, X_title = "Best estimator (keVee)", Y_title = "c/keV", color = kRed)
    hsignal.SetMaximum(1.1*hdata.GetMaximum())
    hsignal.SetMinimum(0.1)
    hsignal.Draw()
    hdata.Draw("same")
    hbckg.Draw("same")
    gPad.SetLogy()
    cc.Print("../Figures/Data_and_model/Stacked_elec_rec_perkeV_log_scale_with_limit.png")
    del cc


def check_elec_recoil_bckgmodel_with_tritium(list_bolo_name, mass="0"):
    """Check the background model works ok (with tritium in bckg model)
    
    Detail: 
        void

    Args:
        list_bolo_name = (str) list of bolo names
        
    Returns:
        void

    Raises:
        void
    """

    #Individual plots
    for bolo_name in list_bolo_name :

        hdata, fdata = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/" + bolo_name + "_spec_perkeV.root", "h_" + bolo_name)
        fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/" + bolo_name + "_model_with_tritium_perkeV.root", bolo_name + "_FidGamma")
        ftritium = file_model.Get("tritium")
        ftritium.SetLineColor(kBlue)
        fmodel_noeff = file_model.Get(bolo_name + "_FidGamma_noeff")
        fmodel_noeff.SetLineStyle(7)

        cc = TCanvas("cc", "cc")
        PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
        hdata.SetMinimum(1)
        hdata.Draw()
        ftritium.Draw("same")
        fmodel.Draw("same")
        fmodel_noeff.Draw("same")
        gPad.SetLogy()
        cc.Print("../Figures/Data_and_model/Log_scale/" + bolo_name + "_elec_rec_with_tritium_perkeV.png")
        del cc

        cc = TCanvas("cc", "cc")
        PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
        hdata.Draw()
        ftritium.Draw("same")
        fmodel.Draw("same")
        fmodel_noeff.Draw("same")
        cc.Print("../Figures/Data_and_model/Linear_scale/" + bolo_name + "_elec_rec_with_tritium_perkeV.png")
        del cc

    #Stacked plt
    hdata, fdata = PyRPl.open_ROOT_object("../ROOT_files/Axion/Spec/Stacked_spec_perkeV.root", "hstacked")
    fmodel, file_model = PyRPl.open_ROOT_object("../ROOT_files/Axion/Model/Stacked_model_with_tritium_perkeV.root", "model_stacked")
    ftritium = file_model.Get("tritium_stacked")

    cc = TCanvas("cc", "cc")
    PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
    hdata.GetXaxis().SetRangeUser(0,19)
    hdata.Draw()
    fmodel.Draw("same")
    ftritium.SetLineColor(kBlue)
    ftritium.Draw("same")
    gPad.SetLogy()
    cc.Print("../Figures/Data_and_model/Stacked_elec_rec_with_tritium_perkeV_log_scale.png")
    del cc

    cc = TCanvas("cc", "cc")
    PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
    hdata.GetXaxis().SetRangeUser(0,19)
    hdata.SetMaximum(3*ftritium.GetMaximum())
    hdata.Draw()
    fmodel.Draw("same")
    ftritium.SetLineColor(kBlue)
    ftritium.Draw("same")
    cc.Print("../Figures/Data_and_model/Stacked_elec_rec_with_tritium_perkeV_linear_scale.png")
    del cc

    #Stacked plt with limit
    fflux, file_flux = PyRPl.open_ROOT_object("../ROOT_files/Axion/CBRD_convolved/Stacked_flux_perkeV.root", "flux_stacked_mass_" + str(mass))

    #Histograms parameters
    bin_X, min_X, max_X = hdata.GetNbinsX(), hdata.GetBinCenter(1) - 0.5*hdata.GetBinWidth(0), hdata.GetBinCenter(hdata.GetNbinsX()) + 0.5 * hdata.GetBinWidth(0)

    htritium = TH1F("htritium", "htritium", bin_X, min_X, max_X)
    hbckg = TH1F("hbckg", "hbckg", bin_X, min_X, max_X)
    hsignal = TH1F("hsignal", "hsignal", bin_X, min_X, max_X)
    for i in range(1, bin_X+1):
        #Tritium
        htritium.SetBinContent(i, ftritium.Eval(htritium.GetBinCenter(i)))
        #Bckg
        hbckg.SetBinContent(i, fmodel.Eval(hbckg.GetBinCenter(i)))
        #Signal at limit
        hsignal.SetBinContent(i, np.power(1.3E-11,4)*fflux.Eval(hsignal.GetBinCenter(i)))

    cc = TCanvas("cc", "cc")
    PyRPl.process_TH1(hdata, X_title = "Best estimator (keVee)", Y_title = "c/keV")
    PyRPl.process_TH1(htritium, X_title = "Best estimator (keVee)", Y_title = "c/keV", color = kBlue)
    PyRPl.process_TH1(hbckg, X_title = "Best estimator (keVee)", Y_title = "c/keV", color = kGreen-3)
    PyRPl.process_TH1(hsignal, X_title = "Best estimator (keVee)", Y_title = "c/keV", color = kRed)
    hsignal.GetXaxis().SetRangeUser(0,19)
    hsignal.SetMaximum(1.1*hdata.GetMaximum())
    hsignal.SetMinimum(0.1)
    hsignal.Draw()
    hdata.Draw("same")
    hbckg.Draw("same")
    htritium.Draw("same")
    gPad.SetLogy()
    cc.Print("../Figures/Data_and_model/Stacked_elec_rec_with_tritium_perkeV_log_scale_with_limit.png")
    del cc


def plot_effect_cut_on_signal(bolo_name, mass="0"):
    """Plot effect of cut on the signal
    
    Detail:
        void

    Args:
        bolo_name  = (str) bolo name
        mass = (str) axion mass

    Returns:
        void

    Raises:
        void
    """

    try :
        os.remove("../ROOT_files/Axion/Cut_eff_signal/" + bolo_name + "_cut_eff_signal.root")
    except OSError:
        pass
    script_utils.print_utility("Computing for bolo " + bolo_name)
    t, f = PyRPl.open_ROOT_object("./Event_generation/ROOT_files/" + bolo_name + "_signal_tree.root", "t_new" + str(mass))
    cut = ut.get_electronic_cut_line(bolo_name)

    hcut = TH1F("hcut" + bolo_name, "hcut" + bolo_name, 200, 0, 10)
    hnocut = TH1F("hnocut" + bolo_name, "hnocut" + bolo_name, 200, 0, 10)

    t.Project("hcut" + bolo_name, "Ebest", cut)
    t.Project("hnocut" + bolo_name, "Ebest")

    cc = TCanvas("cc", "cc")
    hcut.SetMaximum(1.1*hnocut.GetMaximum())
    PyRPl.process_TH1(hcut, X_title = "Best estimator (keVee)", Y_title = "A.U.", color = kRed)
    PyRPl.process_TH1(hnocut, X_title = "Best estimator (keVee)", Y_title = "A.U.")
    hcut.Draw()
    hnocut.Draw("same")
    raw_input()
    cc.Print("../Figures/Effect_of_cut_on_signal_linear_scale.png")
    del cc

    cc = TCanvas("cc", "cc")
    hcut.SetMaximum(1.1*hnocut.GetMaximum())
    PyRPl.process_TH1(hcut, X_title = "Best estimator (keVee)", Y_title = "A.U.", color = kRed)
    PyRPl.process_TH1(hnocut, X_title = "Best estimator (keVee)", Y_title = "A.U.")
    hcut.Draw()
    hnocut.Draw("same")
    gPad.SetLogy()
    raw_input()
    cc.Print("../Figures/Effect_of_cut_on_signal_log_scale.png")
    del cc


def plot_igex():
    """ Superimpose igex plot /w tritium model used for EDW
    """

    arr = np.loadtxt("spec_igex.csv", delimiter = ",")
    arr_E, arr_flux = arr[:,0].astype(float), arr[:,1].astype(float)
    gr = TGraph(arr_E.shape[0], arr_E, arr_flux)

    class Beta_tritium:
        def __call__( self, x, par ):
            m_e = 9.10938215E-31 #in kg
            c = 299792458 # in m / s
            Z = 2 # for tritium
            alpha = 1/137. # fine structure constant
            E_0 = 18.6 * 1E3 * 1.6E-19 # tritium endpoint in Joule
            T = x[0]*1E3*1.6E-19 #kinetic energy in Joule
            E = T + m_e * np.power(c,2)
            p = np.sqrt( np.power(T,2) + 2*T*m_e*np.power(c,2) )
            eta = 2 * np.pi * alpha * Z * E / (p*c+1E-10)
            F = 2 * np.pi * eta / (1 - np.exp(-2*np.pi * eta))
            if T <E_0 :
                return par[0]*(1/4.37353532216e-56)*F*p*E*np.power(E_0-T,2) #this constant to have unit integral when par[0] = 1
            else :
                return 0

    ftritium = TF1("tritium", Beta_tritium(), 0, 20, 1)
    ftritium.SetNpx(200)
    ftritium.SetParameter(0,300)
    ftritium.SetLineColor(kBlue)
    ftritium.SetLineStyle(7)

    h = TH1F("h", "h", 100, 0, 20)
    h.SetMaximum(40)
    PyRPl.process_TH1(h, X_title = "Energy (keV)", Y_title = "Counts/keV")
    gr.SetLineColor(kRed) 
    gr.SetName("gr")   
    cc = TCanvas("cc", "cc")
    h.Draw()
    gr.Draw("sameC")
    ftritium.Draw("same")
    l = TLegend(0.1,0.7,0.48,0.9)
    l.AddEntry("gr", "igex", "leg")
    l.AddEntry("tritium","my model" , "leg")
    l.Draw("same")
    raw_input()
    cc.Print("../Figures/igex_vs_edw.png")
    del cc

if __name__ == "__main__":

    bolo_name = "FID825"
    list_bolo_name = ["FID824", "FID825", "FID827", "FID837", "FID838", "FID839", "FID841", "FID842"]
    # list_bolo_name = ["FID825"]
    # get_skimmed_axion_tree(bolo_name)
    # choose_elec_recoil_cut(bolo_name)
    # check_elec_recoil_bckgmodel(list_bolo_name, mass="0")    
    check_elec_recoil_bckgmodel_with_tritium(list_bolo_name, mass="0")    
    # plot_effect_cut_on_signal(bolo_name, mass="0")
    # plot_igex()