#include <string>
#include <iostream>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
#include "TPaletteAxis.h"
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


int build_heatonly_tree(TString bolo_name, TString  analysis_type, int num_event)
{

    /* 
        Build the training tree for heatonly events and the corresponding cut eff file

       Detail:
        Use a fitted 1D x3 exponential spectrum to generate the data
        Loop over events, select those which pass the cut and write to tree
        Also compute the efficiency of the event selection
        Write this efficiency to a .txt file

    Args:
        bolo_name     = (str) bolometer name
        analysis_type = (str) name of analysis (name indicates which ion cut, which resolution...)
        num_event     = (int) number of events to simulate
    Returns:
        void

    Raises:
        void
    */

    // C r e a t e     t r e e
    // ---------------------------------------
    TString outfilename;
    TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");
    TString txt_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_") + bolo_name + TString("/Text_files/") + bolo_name + TString("_bckg_cuteff_") + analysis_type + TString(".txt") ;
    TString heat_txt_path = gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + TString("/Heatonly/Text_files/");

    //Include NR computations directly
    TString util_path = TString("/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/");
    TFile * file_EE_to_NR = new TFile(util_path + TString("conv_EE_to_NR.root"), "read");
    TF1* fEE_to_NR = (TF1*)file_EE_to_NR->Get("conv");

    //Include wavelet efficiency
    TString swave_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/") + bolo_name + TString("/");
    TFile file_wave(swave_path + TString("PSA_cut_eff.root"), "read");
    TF1 * f_wave_eff = (TF1*)file_wave.Get("PSA_eff");

    //Sample the 2D heat bckg model
    TString Ana_path   = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/");
    TString heatonly_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_heatonly_2D.root");
    // TFile * file_heatonly = new TFile(heatonly_path, "read");
    // TH3F * heatonly2D = (TH2F*)file_heatonly->Get("heat2D");

    //Get the histogram of the heatonly rate over time
    TFile * file_heatrate = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_heatonly_rate_over_time_") + analysis_type + TString(".root"), "read"); 
    TH1F* heatrate = (TH1F*)file_heatrate->Get("hheat");


    //Get the histogram of the baselines for event generation
    TFile * file_fwhm = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_fwhm_time_hist_for_BDT.root"), "read"); 
    TH1F* hOWC1 = (TH1F*)file_fwhm->Get("howc1");
    TH1F* hOWC2 = (TH1F*)file_fwhm->Get("howc2");

    TString heatonly_out_path = gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + TString("/Heatonly/ROOT_files/");
    TFile *ftree = new TFile(heatonly_out_path + bolo_name+ TString("_heatonly_tree.root"), "recreate");

    //Get the correlation matrix
    TString corr_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_") + bolo_name + TString("/Text_files/");
    ifstream file_corr(corr_path + bolo_name + TString("_Lmatrix_coeff.txt"));
    TArrayD Lcoeff(16);
    TMatrixD Lmatrx(4,4);
    for (int k = 0; k<16; k++) {file_corr >> Lcoeff[k];}
    // for (int k = 0; k<16; k++) {cout << Lcoeff[k] << endl;}
    Lmatrx.SetMatrixArray(Lcoeff.GetArray());
    // for (int k = 0; k<4; k++)
    // {
    //     for (int i=0; i<4; i++)
    //     {
    //         cout << Lmatrx(k,i) << "," ;
    //     }
    //     cout << endl;
    // } 


    // F i l l     t r e e
    // ---------------------------------------
    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;

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

    int cheatonly = 0;

    //2 simulations: one for training, one for application
    for (int i =0 ; i<2 ; i++)
    {

        int cheatonly, cheatonlybis;
        cheatonly=0;
        cheatonlybis=0;

        //New tree to be filled
        TString s_index = Form("%d", i);
        TTree *t_new = new TTree("t_new" + s_index,"t_new" + s_index);

        float  EC1, EC2, EIA, EIB, EIC, EID, HR;
        t_new->Branch("EC1", &EC1, "EC1/F");
        t_new->Branch("EC2", &EC2, "EC2/F");
        t_new->Branch("EIA", &EIA, "EIA/F");
        t_new->Branch("EIB", &EIB, "EIB/F");
        t_new->Branch("EIC", &EIC, "EIC/F");
        t_new->Branch("EID", &EID, "EID/F");
        t_new->Branch("HR", &HR, "HR/F");


        // Read .txt file with heat only simulated by KDE
        ifstream fheatonly(heat_txt_path + bolo_name + TString("_heatonly_2D_time_") + s_index +TString(".txt"));
        double heat_EC1, heat_EC2, jour;

        //////////////////////////////////////////////////
        // Heatonly events
        //////////////////////////////////////////////////
        for (int k=0; k<num_event; k++) 
        {

            // fheatonly >> heat_EC1 >> heat_EC2;
            fheatonly >> heat_EC1 >> heat_EC2 >> jour;

            // jour = heatrate->GetRandom();
            HR = heatrate->GetBinContent(heatrate->FindBin(jour));

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
            float EI = 0.5*EIB+0.5*EID;
            float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cheatonly++;
            // if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup) cheatonly++;
            if (abs(EC1-EC2)<1 &&  ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            // if ( 1==1) 
            {
                        cheatonlybis++;
                        t_new->Fill();
            }
        }

        // cout << "cheatonlybis" << "," << cheatonlybis << endl;
        string tmpstr;
        stringstream convert;
        convert << cheatonly/float(num_event);
        TString s_cheatonly = convert.str();
        string pattern = "heatonly";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cheatonly);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) {feffout << vec_str[u] << endl;} //; cout << vec_str[u] << endl;}
        feffout.close();

        t_new->Write();

        fheatonly.close();
    }

    ftree->Close();


    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_0.5_0_5");
	int num_event =800000;

    build_heatonly_tree(bolo_name,  analysis_type, num_event);


    return 0;

}