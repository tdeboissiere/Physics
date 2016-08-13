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

using namespace std;

int build_true_event_tree(TString bolo_name, TString  analysis_type, int nsimu)
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
	nFidGamma=563;
	nS1Gamma=164;
	nS2Gamma=123;
	nS1Beta=364;
	nS2Beta=438;
	nS2Beta=438;
	nS1Pb=19 ;
	nS2Pb=35;
	nheatonly=57440;
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

    TFile * file_Beta = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Beta_from_bolo_spectrum_extrapol.root"), "read");
    TF1 * func_S1Beta = (TF1*)file_Beta->Get("S1Beta_extra");
    TF1 * func_S2Beta = (TF1*)file_Beta->Get("S2Beta_extra");

    TFile * file_Pb = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Pb_from_bolo_spectrum_extrapol.root"), "read");
    TF1 * func_S1Pb = (TF1*)file_Pb->Get("S1Pb_extra");
    TF1 * func_S2Pb = (TF1*)file_Pb->Get("S2Pb_extra");

    TString heatonly_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_heatonly_2D.root");
    TFile * file_heatonly = new TFile(heatonly_path, "read");
    TH2F * heatonly2D = (TH2F*)file_heatonly->Get("heat2D");

    //Create TFile + TTree to store the training data
    TString out_path = gen_path + TString("/BDT_") + bolo_name + TString("/") + analysis_type + TString("/True_events/ROOT_files/");
    TFile *ftree = new TFile(out_path + bolo_name+"_true_events_tree.root", "recreate");

    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;
    float Luke_factor = (1+V_fid/3)/(1+V_surf/3);

    float heat_thresh = 0.65;

	float s_heat=0.18727705113;
	float s_EC1_ERA=0.237812128419;
	float s_EC2_ERA=0.276031934772;
	float s_IA=0.38899269577;
	float s_IB=0.361814166808;
	float s_IC=0.307881773399;
	float s_ID=0.304484457279;

	float cut_vetA=19.4496347885;
	float cut_vetC=15.39408867;

	float ECinf=1.5;
	float ECsup=15;
	float EIinf=-2;
	float EIsup=15;

    //Loop over number of times we want to simulate the experiment
    for (int i =0 ; i<nsimu ; i++)
    {

        cout << "Starting simulation " << i << endl;
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

        //New tree to be filled
        TString s_index = Form("%d", i);
        TTree *t_new = new TTree("t_new" + s_index,"t_new" + s_index);

        float  EC1, EC2, EIA, EIB, EIC, EID, ENR;
        t_new->Branch("EC1", &EC1, "EC1/F");
        t_new->Branch("EC2", &EC2, "EC2/F");
        t_new->Branch("EIA", &EIA, "EIA/F");
        t_new->Branch("EIB", &EIB, "EIB/F");
        t_new->Branch("EIC", &EIC, "EIC/F");
        t_new->Branch("EID", &EID, "EID/F");
        t_new->Branch("ENR", &ENR, "ENR/F");

        //////////////////////////////////////////////////
        // Heatonly events
        //////////////////////////////////////////////////
        for (int k=0; k<nheatonly; k++) 
        {
            double heat_EC1, heat_EC2;

            heatonly2D->GetRandom2(heat_EC1, heat_EC2);

            EIA = u.Gaus(0,s_IA);
            EIC = u.Gaus(0,s_IC);
            EIB = u.Gaus(0,s_IB);
            EID = u.Gaus(0,s_ID);
            EC1 = heat_EC1 ;
            EC2 = heat_EC2 ; 

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 &&  ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cheatonly++;
                    }
            }
        }

        //////////////////////////////////////////////////
        // Gamma events
        //////////////////////////////////////////////////


        for (int k=0; k<nFidGamma; k++) 
        {

            // float evt_heat = func_FidGamma->GetRandom(ECinf, ECsup);
            float evt_heat = func_FidGamma->GetRandom(0,20);

            EIA = u.Gaus(0,s_IA);
            EIC = u.Gaus(0,s_IC);
            EIB = evt_heat + u.Gaus(0,s_IB);
            EID = evt_heat + u.Gaus(0,s_ID);
            EC1 = evt_heat + u.Gaus(0,s_EC1_ERA);
            EC2 = evt_heat + u.Gaus(0,s_EC2_ERA);

            // EIA = u.Gaus(0,s_IA);
            // EIC = u.Gaus(0,s_IC);
            // EIB = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) -pow(s_EC1_ERA,2)));
            // EID = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) -pow(s_EC1_ERA,2)));
            // EC1 = evt_heat ;
            // EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) -pow(s_EC1_ERA,2)));

            // EIA = u.Gaus(0,TMath::Sqrt(pow(s_IA,2) ));
            // EIC = u.Gaus(0,TMath::Sqrt(pow(s_IC,2) ));
            // EIB = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) ));
            // EID = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) ));
            // EC1 = evt_heat ;
            // EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) ));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cFidGamma++;
                    }
            }
        }
    
        for (int k=0; k<nS1Gamma; k++) 
        {

            // float evt_heat = func_S1Gamma->GetRandom(ECinf, ECsup);
            float evt_heat = func_S1Gamma->GetRandom(0,20);

            EIA = Luke_factor*evt_heat + u.Gaus(0,s_IA);
            EIC = u.Gaus(0,s_IC);
            EIB = Luke_factor*evt_heat + u.Gaus(0,s_IB);
            EID = u.Gaus(0, s_ID);
            EC1 = evt_heat + u.Gaus(0,s_EC1_ERA);
            EC2 = evt_heat + u.Gaus(0,s_EC2_ERA);

            // EIA = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) -pow(s_EC1_ERA,2)));
            // EIC = u.Gaus(0,s_IC);
            // EIB = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) -pow(s_EC1_ERA,2)));
            // EID = u.Gaus(0,s_ID);
            // EC1 = evt_heat ;
            // EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) -pow(s_EC1_ERA,2)));

            // EIA = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) ));
            // EIC = u.Gaus(0,s_IC);
            // EIB = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) ));
            // EID = u.Gaus(0,s_ID);
            // EC1 = evt_heat ;
            // EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) ));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cS1Gamma++;
                    }
            }
        }
    
        for (int k=0; k<nS2Gamma; k++) 
        {

            // float evt_heat = func_S2Gamma->GetRandom(ECinf, ECsup);
            float evt_heat = func_S2Gamma->GetRandom(0,20);

            EIA = u.Gaus(0,s_IA);
            EIC = Luke_factor*evt_heat + u.Gaus(0,s_IC);
            EIB = u.Gaus(0,s_IB);
            EID = Luke_factor*evt_heat + u.Gaus(0,s_ID);
            EC1 = evt_heat + u.Gaus(0,s_EC1_ERA);
            EC2 = evt_heat + u.Gaus(0,s_EC2_ERA);

            // EIA = u.Gaus(0,s_IA);
            // EIC = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) -pow(s_EC1_ERA,2)));
            // EIB = u.Gaus(0,s_IB);
            // EID = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) -pow(s_EC1_ERA,2)));
            // EC1 = evt_heat ;
            // EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) -pow(s_EC1_ERA,2)));

            // EIA = u.Gaus(0,s_IA);
            // EIC = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) ));
            // EIB = u.Gaus(0,s_IB);
            // EID = Luke_factor*evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) ));
            // EC1 = evt_heat ;
            // EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) ));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cS2Gamma++;
                    }
            }
        }
    


        //////////////////////////////////////////
        // Beta
        //////////////////////////////////////////

        for (int k=0; k<nS1Beta; k++) 
        {

            // float evt_heat = func_S1Beta->GetRandom(ECinf, ECsup);
            float evt_heat = func_S1Beta->GetRandom(0,20);
            float strag=evt_heat*0.1;

            EIA = func_convS1Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIC = u.Gaus(0,s_IC);
            EIB = func_convS1Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EID = u.Gaus(0,s_ID);
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cS1Beta++;
                    }
            }
        }

        for (int k=0; k<nS2Beta; k++) 
        {

            // float evt_heat = func_S2Beta->GetRandom(ECinf, ECsup);
            float evt_heat = func_S2Beta->GetRandom(0,20);
            float strag=evt_heat*0.1;

            EIA = u.Gaus(0,s_IA);
            EIC = func_convS2Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIB = u.Gaus(0,s_IB);
            EID = func_convS2Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cS2Beta++;
                    }
            }
        }

        //////////////////////////////////////////
        // Pb
        //////////////////////////////////////////

        for (int k=0; k<nS1Pb; k++) 
        {

            // float evt_heat = func_S1Pb->GetRandom(ECinf, ECsup);
            float evt_heat = func_S1Pb->GetRandom(0,20);
            float strag=evt_heat*0.05;

            EIA = func_convS1Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIC = u.Gaus(0,s_IC);
            EIB = func_convS1Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EID = u.Gaus(0,s_ID);
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cS1Pb++;
                    }
            }
        }

        for (int k=0; k<nS2Pb; k++) 
        {

            // float evt_heat = func_S2Pb->GetRandom(ECinf, ECsup);
            float evt_heat = func_S2Pb->GetRandom(0,20);
            float strag=evt_heat*0.05;

            EIA = u.Gaus(0,s_IA);
            EIC = func_convS2Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIB = u.Gaus(0,s_IB);
            EID = func_convS2Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.574*EC1+0.426*EC2;
			float EI =0.414*EIB+0.586*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                        cS2Pb++;
                    }
            }
        }

        double sum_counts = cFidGamma+ cS1Gamma+ cS2Gamma+ cS1Beta+ cS2Beta+ cS1Pb+ cS2Pb+ cheatonly;
        cout << "total_counts  " << sum_counts << endl;
        cout << cheatonly << "," << cFidGamma << "," << cS1Gamma << "," << cS2Gamma << "," << cS1Beta << "," << cS2Beta << "," << cS1Pb << "," << cS2Pb << endl;
        cout << cheatonly/sum_counts <<"," << cFidGamma/sum_counts <<"," << cS1Gamma/sum_counts <<"," << cS2Gamma/sum_counts <<"," << cS1Beta/sum_counts <<"," << cS2Beta/sum_counts <<"," << cS1Pb/sum_counts <<"," << cS2Pb/sum_counts << endl; 

        t_new->Write();
    }
    ftree->Close();
    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_1.5_min2_50");
	int nsimu=10;

    build_true_event_tree(bolo_name,  analysis_type, nsimu);

    return 0;

}