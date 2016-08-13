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

int get_MVA_scaling(TString bolo_name, TString  analysis_type, float exposure)
{

    /*Build the training tree for Beta/Pb events and the corresponding cut eff file
    
    Detail:
        Use 1D fit to data as the background model 
        Loop over events, select those which pass the cut and write to tree
        Also compute the efficiency of the event selection
        Write this efficiency to a .txt file

    Args:
        bolo_name     = (str) bolometer name
        analysis_type = (str) name of analysis (name indicates which ion cut, which resolution...)

    Returns:
        void

    Raises:
        void
    */


    ///////////////////////////////////////////////////////////////
    ///
    /// No need to het heat rate here (no tree created)
    ///
    ////////////////////////////////////////////////////////////////


    // C r e a t e     t r e e
    // ---------------------------------------
    TString outfilename;
    TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");

    //Include NR computations directly
    TString util_path = TString("/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/");
    TFile * file_EE_to_NR = new TFile(util_path + TString("conv_EE_to_NR.root"), "read");
    TF1* fEE_to_NR = (TF1*)file_EE_to_NR->Get("conv");

    //Include wavelet efficiency
    TString swave_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/") + bolo_name + TString("/");
    TFile file_wave(swave_path + TString("PSA_cut_eff.root"), "read");
    TF1 * f_wave_eff = (TF1*)file_wave.Get("PSA_eff");

    int nFidGamma, nS1Gamma, nS2Gamma, nS1Beta, nS2Beta, nS1Pb, nS2Pb, nheatonly;
	nFidGamma=78568;
	nS1Gamma=15535;
	nS2Gamma=14729;
	nS1Beta=62245;
	nS2Beta=72891;
	nS2Beta=72891;
	nS1Pb=3790 ;
	nS2Pb=6931;
	nheatonly=32147707;
    cout <<  nFidGamma<< "," << nS1Gamma<< "," << nS2Gamma<< "," << nS1Beta<< "," << nS2Beta<< "," << nS1Pb<< "," << nS2Pb<< "," << nheatonly << endl;

    //Get the function that converts EI to EC for surface events
    TString Ana_path   = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/");
    TString convert_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_ion_heat_surface_relation.root");
    TFile * file_conv = new TFile(convert_path, "read");
    TF1* func_convS1Beta = (TF1*)file_conv->Get("S1Beta");
    TF1* func_convS2Beta = (TF1*)file_conv->Get("S2Beta");
    TF1* func_convS1Pb  = (TF1*)file_conv->Get("S1Pb");
    TF1* func_convS2Pb  = (TF1*)file_conv->Get("S2Pb");

    // Sample the  bckg model
    TFile * file_Gamma = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Gamma_spectrum_deconv.root"), "read");
    TF1 * func_FidGamma = (TF1*)file_Gamma->Get("FidGamma_deconv");
    TF1 * func_S1Gamma = (TF1*)file_Gamma->Get("S1Gamma_deconv");
    TF1 * func_S2Gamma = (TF1*)file_Gamma->Get("S2Gamma_deconv");

    TFile * file_Beta = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Beta_spectrum_extrapol.root"), "read");
    TF1 * func_S1Beta = (TF1*)file_Beta->Get("S1Beta_extra");
    TF1 * func_S2Beta = (TF1*)file_Beta->Get("S2Beta_extra");

    TFile * file_Pb = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Pb_spectrum_extrapol.root"), "read");
    TF1 * func_S1Pb = (TF1*)file_Pb->Get("S1Pb_extra");
    TF1 * func_S2Pb = (TF1*)file_Pb->Get("S2Pb_extra");

    TString heatonly_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_heatonly_2D.root");
    TFile * file_heatonly = new TFile(heatonly_path, "read");
    TH2F * heatonly2D = (TH2F*)file_heatonly->Get("heat2D");
    
    //Get the histogram of the baselines for event generation
    TFile * file_fwhm = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_fwhm_time_hist_for_BDT.root"), "read"); 
    TH1F* hFWIA = (TH1F*)file_fwhm->Get("hfwia");
    TH1F* hFWIB = (TH1F*)file_fwhm->Get("hfwib");
    TH1F* hFWIC = (TH1F*)file_fwhm->Get("hfwic");
    TH1F* hFWID = (TH1F*)file_fwhm->Get("hfwid");
    TH1F* hOWC1 = (TH1F*)file_fwhm->Get("howc1");
    TH1F* hOWC2 = (TH1F*)file_fwhm->Get("howc2");
    TH1F* hJOUR = (TH1F*)file_fwhm->Get("hist_jour");

    //Get the correlation matrix
    TString corr_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_") + bolo_name + TString("/Text_files/");
    ifstream file_corr(corr_path + bolo_name + TString("_Lmatrix_coeff.txt"));
    TArrayD Lcoeff(16);
    TMatrixD Lmatrx(4,4);
    for (int k = 0; k<16; k++) {file_corr >> Lcoeff[k];}
    Lmatrx.SetMatrixArray(Lcoeff.GetArray());

    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;
    float Luke_factor = (1+V_fid/3)/(1+V_surf/3);

    float heat_thresh = 0.65;

	float s_IA=0.38644470868;
	float s_IB=0.354170205538;
	float s_IC=0.304484457279;
	float s_ID=0.297265160523;

	float cut_vetA=1.9322235434;
	float cut_vetC=1.52242228639;

	float ECinf=0.5;
	float ECsup=15;
	float EIinf=0;
	float EIsup=15;


    int cFidGamma, cS1Gamma, cS2Gamma, cS1Beta, cS2Beta, cS1Pb, cS2Pb, cheatonly;
    cFidGamma=0;
    cS1Gamma=0;
    cS2Gamma=0;
    cS1Beta=0;
    cS2Beta=0;
    cS2Beta=0;
    cS1Pb=0;
    cS2Pb=0;
    cheatonly=0;


    float  EC1, EC2, EIA, EIB, EIC, EID, ENR;

    //////////////////////////////////////////////////
    // Heatonly events
    //////////////////////////////////////////////////
    for (int k=0; k<nheatonly; k++) 
    {
        double heat_EC1, heat_EC2;

        heatonly2D->GetRandom2(heat_EC1, heat_EC2);
        
        //Get JOUR
        float jour = hJOUR->GetRandom();
        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);

        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EC1 = heat_EC1 ;
        EC2 = heat_EC2 ; 
        EIA = ion_matrx(0,0);
        EIB = ion_matrx(1,0);
        EIC = ion_matrx(2,0);
        EID = ion_matrx(3,0);

		float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
		float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 &&  ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
                ENR= fEE_to_NR->Eval(EC);
                cheatonly++;
        }
    }

    //////////////////////////////////////////////////
    // Gamma events
    //////////////////////////////////////////////////


    for (int k=0; k<nFidGamma; k++) 
    {

        // float evt_heat = func_FidGamma->GetRandom(ECinf, ECsup);
        float evt_heat = func_FidGamma->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);

        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = ion_matrx(0,0);
        EIC = ion_matrx(2,0);
        EIB = evt_heat + ion_matrx(1,0);
        EID = evt_heat + ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0,s_EC1);
        EC2 = evt_heat + u.Gaus(0,s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
            if (u.Binomial(1,0.9999*accept)==1) 
            {
                ENR= fEE_to_NR->Eval(EC);
                cFidGamma++;
            }
        }
    }
    
    for (int k=0; k<nS1Gamma; k++) 
    {

        // float evt_heat = func_S1Gamma->GetRandom(ECinf, ECsup);
        float evt_heat = func_S1Gamma->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);


        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = Luke_factor*evt_heat + ion_matrx(0,0);
        EIC = ion_matrx(2,0);
        EIB = Luke_factor*evt_heat + ion_matrx(1,0);
        EID = ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0,s_EC1);
        EC2 = evt_heat + u.Gaus(0,s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
            if (u.Binomial(1,0.9999*accept)==1) 
            {
                ENR= fEE_to_NR->Eval(EC);
                cS1Gamma++;
            }
        }
    }
    
    for (int k=0; k<nS2Gamma; k++) 
    {

        float evt_heat = func_S2Gamma->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);


        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = ion_matrx(0,0);
        EIC = Luke_factor*evt_heat + ion_matrx(2,0);
        EIB = ion_matrx(1,0);
        EID = Luke_factor*evt_heat + ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0,s_EC1);
        EC2 = evt_heat + u.Gaus(0,s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
            if (u.Binomial(1,0.9999*accept)==1) 
            {
                ENR= fEE_to_NR->Eval(EC);
                cS2Gamma++;
            }
        }
    }
    


    //////////////////////////////////////////
    // Beta
    //////////////////////////////////////////

    for (int k=0; k<nS1Beta; k++) 
    {

        float evt_heat = func_S1Beta->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);


        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = func_convS1Beta->Eval(evt_heat) + ion_matrx(0,0);
        EIC = ion_matrx(2,0);
        EIB = func_convS1Beta->Eval(evt_heat) + ion_matrx(1,0);
        EID = ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0,s_EC1);
        EC2 = evt_heat + u.Gaus(0,s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
                if (u.Binomial(1,0.9999*accept)==1) 
                {
                    ENR= fEE_to_NR->Eval(EC);
                    cS1Beta++;
                }
        }
    }

    for (int k=0; k<nS2Beta; k++) 
    {

        float evt_heat = func_S2Beta->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);


        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = ion_matrx(0,0);
        EIC = func_convS2Beta->Eval(evt_heat) + ion_matrx(2,0);
        EIB = ion_matrx(1,0);
        EID = func_convS2Beta->Eval(evt_heat) + ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0,s_EC1);
        EC2 = evt_heat + u.Gaus(0,s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
            if (u.Binomial(1,0.9999*accept)==1) 
            {
                ENR= fEE_to_NR->Eval(EC);
                cS2Beta++;
            }
        }
    }

    //////////////////////////////////////////
    // Pb
    //////////////////////////////////////////

    for (int k=0; k<nS1Pb; k++) 
    {

        float evt_heat = func_S1Pb->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);


        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = func_convS1Pb->Eval(evt_heat) + ion_matrx(0,0);
        EIC = ion_matrx(2,0);
        EIB = func_convS1Pb->Eval(evt_heat) + ion_matrx(1,0);
        EID = ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0, s_EC1);
        EC2 = evt_heat + u.Gaus(0, s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
                if (u.Binomial(1,0.9999*accept)==1) 
                {
                    ENR= fEE_to_NR->Eval(EC);
                    cS1Pb++;
                }
        }
    }

    for (int k=0; k<nS2Pb; k++) 
    {

        float evt_heat = func_S2Pb->GetRandom(0,20);

        //Get JOUR
        float jour = hJOUR->GetRandom();

        float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
        float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
        float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
        float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);


        TArrayD reso_coeff(4);
        for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
        TMatrixD reso_matrx(4,1);
        reso_matrx.SetMatrixArray(reso_coeff.GetArray());
        TMatrixD ion_matrx = Lmatrx*reso_matrx;

        EIA = ion_matrx(0,0);
        EIC = func_convS2Pb->Eval(evt_heat) + ion_matrx(2,0);
        EIB = ion_matrx(1,0);
        EID = func_convS2Pb->Eval(evt_heat) + ion_matrx(3,0);
        EC1 = evt_heat + u.Gaus(0, s_EC1);
        EC2 = evt_heat + u.Gaus(0, s_EC2);

			float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
			float EI =0.413*EIB+0.587*EID;
        float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
        {
                if (u.Binomial(1,0.9999*accept)==1) 
                {
                    ENR= fEE_to_NR->Eval(EC);
                    cS2Pb++;
                }
        }
    }

    double sum_counts = cFidGamma+ cS1Gamma+ cS2Gamma+ cS1Beta+ cS2Beta+ cS1Pb+ cS2Pb+ cheatonly;
    cout << "total_counts  " << sum_counts << endl;
    cout << cheatonly << "," << cFidGamma << "," << cS1Gamma << "," << cS2Gamma << "," << cS1Beta << "," << cS2Beta << "," << cS1Pb << "," << cS2Pb << endl;
    cout << cheatonly/sum_counts <<"," << cFidGamma/sum_counts <<"," << cS1Gamma/sum_counts <<"," << cS2Gamma/sum_counts <<"," << cS1Beta/sum_counts <<"," << cS2Beta/sum_counts <<"," << cS1Pb/sum_counts <<"," << cS2Pb/sum_counts << endl; 

    ofstream scale_file(Ana_path + TString("Analyse_") + bolo_name + TString("/Text_files/") + bolo_name + TString("_MVA_scaling_") + analysis_type + TString(".txt"));
    scale_file << "prop_heatonly,prop_FidGamma,prop_S1Gamma,prop_S2Gamma,prop_S1Beta,prop_S2Beta,prop_S1Pb,prop_S2Pb,exp_per_day" << endl; 
    scale_file << cheatonly/sum_counts <<"," << cFidGamma/sum_counts <<"," << cS1Gamma/sum_counts <<"," << cS2Gamma/sum_counts <<"," << cS1Beta/sum_counts <<"," << cS2Beta/sum_counts <<"," << cS1Pb/sum_counts <<"," << cS2Pb/sum_counts << "," << sum_counts/exposure <<endl; 
    scale_file.close();
    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_0.5_0_5");
	float exposure =10000;

    get_MVA_scaling(bolo_name,  analysis_type, exposure);

    return 0;

}