#include <string>
#include <iostream>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <fstream>
#include <sstream>
#include <string>
#include "TMatrixD.h"
#include "TArrayD.h"

using namespace std;


double compute_coeff(double FWHM1, double FWHM2)
{

    double s1=FWHM1/2.3548;
    double s2=FWHM2/2.3548;

    double w1=(s2*s2)/(s1*s1+s2*s2);
    return w1;
}


double compute_resolution(double FWHM1, double FWHM2) 
{

    double s1=FWHM1/2.3548;
    double s2=FWHM2/2.3548;

    double w1=(s2*s2)/(s1*s1+s2*s2);

    double res=sqrt(pow(w1*s1,2)+pow((1-w1)*s2,2));
    return 2.3548*res;

}


int build_Gamma_tree(TString bolo_name)
{

    /*Build the training tree for Gamma events and the corresponding cut eff file
    
    Detail:
        Use 1D fit to data as the background model 
        Loop over events, select those which pass the cut and write to tree
        Also compute the efficiency of the event selection
        Write this efficiency to a .txt file

    Args:
        bolo_name     = (str) bolometer name
        analysis_type = (str) name of analysis (name indicates which ion cut, which resolution...)
        event_type    = (str) class of event (S1Gamma S1Pb etc...)
        d_cut         = (dict) dictionnary which states which cut to use on heat/ion/veto
        num_event       = (int) number of events to simulate


    Returns:
        void

    Raises:
        void
    */

    // C r e a t e     t r e e
    // ---------------------------------------
    TString outfilename;
    TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_axion/");

    // Sample the 1D Gamma bckg model
    TString Ana_path   = gen_path + TString("Scripts/Bolo_ana/");
    TFile * file_Gamma = new TFile(Ana_path+  bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Gamma_spectrum_deconv.root"), "read");
    TF1 * func_FidGamma = (TF1*)file_Gamma->Get("FidGamma_deconv");

    // Get the trigger eff function 
    TFile * feff = new TFile(gen_path + TString("ROOT_files/Axion/Trigger_eff/") + bolo_name + TString("_trigger_eff.root"), "read" );
    TF1 * ftrigger = (TF1*)feff->Get("trigger_eff");

    //Get the histogram of the baselines for event generation
    TFile * file_fwhm = new TFile(Ana_path + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_fwhm_time_hist.root"), "read"); 
    TH1F* hFWIA = (TH1F*)file_fwhm->Get("hfwia");
    TH1F* hFWIC = (TH1F*)file_fwhm->Get("hfwic");
    TH1F* hFWFID = (TH1F*)file_fwhm->Get("hfwfid");
    TH1F* hFWC = (TH1F*)file_fwhm->Get("hfwc");
    TH1F* hJOUR = (TH1F*)file_fwhm->Get("hist_jour");

    TString Gamma_out_path = gen_path + TString("Event_generation/ROOT_files/");
    TFile *ftree = new TFile(Gamma_out_path + bolo_name+"_FidGamma_tree.root", "recreate");

    //Get the correlation matrix
    TString corr_path = TString("./Text_files/");
    ifstream file_corr(corr_path + bolo_name + TString("_Lmatrix_coeff.txt"));
    TArrayD Lcoeff(16);
    TMatrixD Lmatrx(4,4);
    for (int k = 0; k<16; k++) {file_corr >> Lcoeff[k];}
    Lmatrx.SetMatrixArray(Lcoeff.GetArray());


    // F i l l     t r e e
    // ---------------------------------------
    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;

    //New tree to be filled
    TTree *t_new = new TTree("t_new", "t_new");

    float EIA, EIC, EIB, EID, FWIA, FWIC, EC, EFID, FWC, FWFID;
    t_new->Branch("EIA", &EIA, "EIA/F");
    t_new->Branch("EIC", &EIC, "EIC/F");

    t_new->Branch("FWIA", &FWIA, "FWIA/F");
    t_new->Branch("FWIC", &FWIC, "FWIC/F");

    t_new->Branch("EC", &EC, "EC/F");
    t_new->Branch("EFID", &EFID, "EFID/F");

    t_new->Branch("FWC", &FWC, "FWC/F");
    t_new->Branch("FWFID", &FWFID, "FWFID/F");


    for (int k=0; k<1000000; k++) 
    {

        float evt_heat = func_FidGamma->GetRandom(0,20);

        //Get a day
        float jour = hJOUR->GetRandom();

        FWC = hFWC->GetBinContent(hFWC->FindBin(jour));
        FWFID = hFWFID->GetBinContent(hFWFID->FindBin(jour));
        FWIA = hFWIA->GetBinContent(hFWIA->FindBin(jour));
        FWIC = hFWIC->GetBinContent(hFWIC->FindBin(jour));

        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = ion_matrx(0,0);
        EIC = ion_matrx(2,0);
        EIB = evt_heat + ion_matrx(1,0);
        EID = evt_heat + ion_matrx(3,0);
        EC = evt_heat + u.Gaus(0,FWC/2.3548);

        float EFID = 0.5*EIB+0.5*EID;
        float accept = ftrigger->Eval(EC);

        if (u.Binomial(1,0.9999*accept)==1) 
        {
            t_new->Fill();
        }
        
    }


    ftree->Close();


    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
    build_Gamma_tree(bolo_name);


    return 0;

}