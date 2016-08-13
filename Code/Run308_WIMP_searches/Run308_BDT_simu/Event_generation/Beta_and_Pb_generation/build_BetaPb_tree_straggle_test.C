#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
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
#include <string>

using namespace std;

int build_BetaPb_tree(TString bolo_name, TString  analysis_type, TString event_type,  int num_event)
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
    TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");
    TString txt_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_") + bolo_name + TString("/Text_files/") + bolo_name + TString("_bckg_cuteff_") + analysis_type + TString(".txt") ;

    //Include NR computations directly
    TString util_path = TString("/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/");
    TFile * file_EE_to_NR = new TFile(util_path + TString("conv_EE_to_NR.root"), "read");
    TF1* fEE_to_NR = (TF1*)file_EE_to_NR->Get("conv");

    //Include wavelet efficiency
    TString swave_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/") + bolo_name + TString("/");
    TFile file_wave(swave_path + TString("PSA_cut_eff.root"), "read");
    TF1 * f_wave_eff = (TF1*)file_wave.Get("PSA_eff");

    // Sample the 1D eta/Pb bckg model
    TString Ana_path   = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/");
    TFile * file_Beta = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Beta_spectrum_extrapol.root"), "read");
    TF1 * func_S1Beta = (TF1*)file_Beta->Get("S1Beta_extra");
    TF1 * func_S2Beta = (TF1*)file_Beta->Get("S2Beta_extra");
    TFile * file_Pb = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Pb_spectrum_extrapol.root"), "read");
    TF1 * func_S1Pb = (TF1*)file_Pb->Get("S1Pb_extra");
    TF1 * func_S2Pb = (TF1*)file_Pb->Get("S2Pb_extra");


    //Get the histogram of the baselines for event generation
    TFile * file_fwhm = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_fwhm_time_hist_for_BDT.root"), "read"); 
    TH1F* hFWIA = (TH1F*)file_fwhm->Get("hfwia");
    TH1F* hFWIB = (TH1F*)file_fwhm->Get("hfwib");
    TH1F* hFWIC = (TH1F*)file_fwhm->Get("hfwic");
    TH1F* hFWID = (TH1F*)file_fwhm->Get("hfwid");
    TH1F* hOWC1 = (TH1F*)file_fwhm->Get("howc1");
    TH1F* hOWC2 = (TH1F*)file_fwhm->Get("howc2");
    TH1F* hJOUR = (TH1F*)file_fwhm->Get("hist_jour");

    //Get the function that converts EI to EC for surface events
    TString convert_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_ion_heat_surface_relation.root");
    TFile * file_conv = new TFile(convert_path, "read");
    TF1* func_convS1Beta = (TF1*)file_conv->Get("S1Beta");
    TF1* func_convS2Beta = (TF1*)file_conv->Get("S2Beta");
    TF1* func_convS1Pb  = (TF1*)file_conv->Get("S1Pb");
    TF1* func_convS2Pb  = (TF1*)file_conv->Get("S2Pb");

    //Create TFile + TTree to store the training data

    TString Surf_out_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  "_" + event_type + TString("_tree_strag_test.root");
    TFile *ftree = new TFile(Surf_out_path , "recreate");

    // F i l l     t r e e
    // ---------------------------------------
    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;

    float heat_thresh = 0.65;

	float cut_vetA=1.9322235434;
	float cut_vetC=1.52242228639;

	float ECinf=0.5;
	float ECsup=15;
	float EIinf=0;
	float EIsup=15;

    //2 simulations: one for training, one for application
    for (int i =0 ; i<1 ; i++)
    {


        //New tree to be filled
        TTree *t_new = new TTree("t_new" ,"t_new" );

        float  EC1, EC2, EIA, EIB, EIC, EID, ENR;
        t_new->Branch("EC1", &EC1, "EC1/F");
        t_new->Branch("EC2", &EC2, "EC2/F");
        t_new->Branch("EIA", &EIA, "EIA/F");
        t_new->Branch("EIB", &EIB, "EIB/F");
        t_new->Branch("EIC", &EIC, "EIC/F");
        t_new->Branch("EID", &EID, "EID/F");
        t_new->Branch("ENR", &ENR, "ENR/F");

        //////////////////////////////////////////
        // Beta
        //////////////////////////////////////////
        if (event_type == TString("S1Beta"))
        {
            for (int k=0; k<num_event; k++) 
            {

                // float evt_heat = func_S1Beta->GetRandom(ECinf, ECsup);
                float evt_heat = func_S1Beta->GetRandom(0,20);
                float strag=evt_heat*0.0765+0.1;

                float jour = hJOUR->GetRandom();

                float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
                float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
                float s_IA  = hFWIA->GetBinContent(hFWIA->FindBin(jour))/2.3548;                
                float s_IB  = hFWIB->GetBinContent(hFWIB->FindBin(jour))/2.3548;
                float s_IC  = hFWIC->GetBinContent(hFWIC->FindBin(jour))/2.3548;
                float s_ID  = hFWID->GetBinContent(hFWID->FindBin(jour))/2.3548;

                float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EIB = compute_coeff(s_IB*2.3548, s_ID*2.3548);
                
                EIA = func_convS1Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) + pow(strag,2)));
                EIC = u.Gaus(0,s_IC);
                EIB = func_convS1Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) + pow(strag,2)));
                EID = u.Gaus(0,s_ID);
                EC1 = evt_heat + u.Gaus(0,s_EC1);
                EC2 = evt_heat + u.Gaus(0,s_EC2);

                float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
                float EI =coeff_EIB*EIB+(1-coeff_EIB)*EID;
                float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

                if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
                {
                        if (u.Binomial(1,0.9999*accept)==1) 
                        {
                            ENR= fEE_to_NR->Eval(EC);
                            t_new->Fill();
                        }
                }
            }



        }

        if (event_type == TString("S2Beta"))
        {        
            for (int k=0; k<num_event; k++) 
            {

                // float evt_heat = func_S2Beta->GetRandom(ECinf, ECsup);
                float evt_heat = func_S2Beta->GetRandom(0,20);
                float strag=evt_heat*0.0765+0.1;

                float jour = hJOUR->GetRandom();

                float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
                float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
                float s_IA  = hFWIA->GetBinContent(hFWIA->FindBin(jour))/2.3548;                
                float s_IB  = hFWIB->GetBinContent(hFWIB->FindBin(jour))/2.3548;
                float s_IC  = hFWIC->GetBinContent(hFWIC->FindBin(jour))/2.3548;
                float s_ID  = hFWID->GetBinContent(hFWID->FindBin(jour))/2.3548;

                float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EIB = compute_coeff(s_IB*2.3548, s_ID*2.3548);

                EIA = u.Gaus(0,s_IA);
                EIC = func_convS2Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) + pow(strag,2)));
                EIB = u.Gaus(0,s_IB);
                EID = func_convS2Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) + pow(strag,2)));
                EC1 = evt_heat + u.Gaus(0,s_EC1);
                EC2 = evt_heat + u.Gaus(0,s_EC2);

                float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
                float EI =coeff_EIB*EIB+(1-coeff_EIB)*EID;
                float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

                if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
                {
                        if (u.Binomial(1,0.9999*accept)==1) 
                        {
                            ENR= fEE_to_NR->Eval(EC);
                            t_new->Fill();
                        }
                }
            }

        }


        //////////////////////////////////////////
        // Pb
        //////////////////////////////////////////


        if (event_type == TString("S1Pb"))
        {
            for (int k=0; k<num_event; k++) 
            {

                // float evt_heat = func_S1Pb->GetRandom(ECinf, ECsup);
                float evt_heat = func_S1Pb->GetRandom(0,20);
                float strag=evt_heat*0.05;

                float jour = hJOUR->GetRandom();

                float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
                float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
                float s_IA  = hFWIA->GetBinContent(hFWIA->FindBin(jour))/2.3548;                
                float s_IB  = hFWIB->GetBinContent(hFWIB->FindBin(jour))/2.3548;
                float s_IC  = hFWIC->GetBinContent(hFWIC->FindBin(jour))/2.3548;
                float s_ID  = hFWID->GetBinContent(hFWID->FindBin(jour))/2.3548;

                float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EIB = compute_coeff(s_IB*2.3548, s_ID*2.3548);

                EIA = func_convS1Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) + pow(strag,2)));
                EIC = u.Gaus(0,s_IC);
                EIB = func_convS1Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) + pow(strag,2)));
                EID = u.Gaus(0,s_ID);
                EC1 = evt_heat + u.Gaus(0, s_EC1);
                EC2 = evt_heat + u.Gaus(0, s_EC2);

                float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
                float EI =coeff_EIB*EIB+(1-coeff_EIB)*EID;
                float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

                if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
                {
                        if (u.Binomial(1,0.9999*accept)==1) 
                        {
                            ENR= fEE_to_NR->Eval(EC);
                            t_new->Fill();
                        }
                }
            }
        }

        if (event_type == TString("S2Pb"))
        {
            for (int k=0; k<num_event; k++) 
            {

                // float evt_heat = func_S2Pb->GetRandom(ECinf, ECsup);
                float evt_heat = func_S2Pb->GetRandom(0,20);
                float strag=evt_heat*0.05;

                float jour = hJOUR->GetRandom();

                float s_EC1 = hOWC1->GetBinContent(hOWC1->FindBin(jour))/2.3548;
                float s_EC2 = hOWC2->GetBinContent(hOWC2->FindBin(jour))/2.3548;
                float s_IA  = hFWIA->GetBinContent(hFWIA->FindBin(jour))/2.3548;                
                float s_IB  = hFWIB->GetBinContent(hFWIB->FindBin(jour))/2.3548;
                float s_IC  = hFWIC->GetBinContent(hFWIC->FindBin(jour))/2.3548;
                float s_ID  = hFWID->GetBinContent(hFWID->FindBin(jour))/2.3548;

                float s_heat = compute_resolution(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EC1 = compute_coeff(s_EC1*2.3548, s_EC2*2.3548);
                float coeff_EIB = compute_coeff(s_IB*2.3548, s_ID*2.3548);

                EIA = u.Gaus(0,s_IA);
                EIC = func_convS2Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) + pow(strag,2)));
                EIB = u.Gaus(0,s_IB);
                EID = func_convS2Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) + pow(strag,2)));
                EC1 = evt_heat + u.Gaus(0, s_EC1);
                EC2 = evt_heat + u.Gaus(0, s_EC2);

                float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
                float EI =coeff_EIB*EIB+(1-coeff_EIB)*EID;
                float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

                if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
                {
                        if (u.Binomial(1,0.9999*accept)==1) 
                        {
                            ENR= fEE_to_NR->Eval(EC);
                            t_new->Fill();
                        }
                }
            }
        }


    t_new->Write();
    }

    ftree->Close();
    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_0.5_0_5");
    int num_event[4];

	num_event[0] =320;
	num_event[1] =320;
	num_event[2] =20;
	num_event[3] =130;

    build_BetaPb_tree(bolo_name,  analysis_type, TString("S1Beta"),  num_event[0]);
    build_BetaPb_tree(bolo_name,  analysis_type, TString("S2Beta"),  num_event[1]);
    build_BetaPb_tree(bolo_name,  analysis_type, TString("S1Pb"),  num_event[2]);
    build_BetaPb_tree(bolo_name,  analysis_type, TString("S2Pb"),  num_event[3]);


    return 0;

}

        
