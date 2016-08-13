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
    TFile * file_Beta = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Beta_from_bolo_spectrum_extrapol.root"), "read");
    TF1 * func_S1Beta = (TF1*)file_Beta->Get("S1Beta_extra");
    TF1 * func_S2Beta = (TF1*)file_Beta->Get("S2Beta_extra");
    TFile * file_Pb = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Pb_from_bolo_spectrum_extrapol.root"), "read");
    TF1 * func_S1Pb = (TF1*)file_Pb->Get("S1Pb_extra");
    TF1 * func_S2Pb = (TF1*)file_Pb->Get("S2Pb_extra");

    //Get the function that converts EI to EC for surface events
    TString convert_path = Ana_path + TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_ion_heat_surface_relation.root");
    TFile * file_conv = new TFile(convert_path, "read");
    TF1* func_convS1Beta = (TF1*)file_conv->Get("S1Beta");
    TF1* func_convS2Beta = (TF1*)file_conv->Get("S2Beta");
    TF1* func_convS1Pb  = (TF1*)file_conv->Get("S1Pb");
    TF1* func_convS2Pb  = (TF1*)file_conv->Get("S2Pb");

    //Create TFile + TTree to store the training data
    TString Surf_out_path = gen_path + TString("/BDT_") + bolo_name + TString("/") + analysis_type + TString("/Beta_and_Pb/ROOT_files/");
    TFile *ftree = new TFile(Surf_out_path + bolo_name+"_"+ event_type + "_tree.root", "recreate");
    TTree *t_new = new TTree("t_new","t_new");

    float  EC1, EC2, EIA, EIB, EIC, EID, ENR;
    t_new->Branch("EC1", &EC1, "EC1/F");
    t_new->Branch("EC2", &EC2, "EC2/F");
    t_new->Branch("EIA", &EIA, "EIA/F");
    t_new->Branch("EIB", &EIB, "EIB/F");
    t_new->Branch("EIC", &EIC, "EIC/F");
    t_new->Branch("EID", &EID, "EID/F");
    t_new->Branch("ENR", &ENR, "ENR/F");

    // F i l l     t r e e
    // ---------------------------------------
    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;

    float heat_thresh = 0.65;

	float s_heat=0.164769831833;
	float s_EC1_ERA=0.193647018855;
	float s_EC2_ERA=0.354170205538;
	float s_IA=0.403431289281;
	float s_IB=0.367759470019;
	float s_IC=0.329539663666;
	float s_ID=0.314251741125;

	float cut_vetA=1.21029386784;
	float cut_vetC=0.988618990997;

	float ECinf=1.5;
	float ECsup=15;
	float EIinf=0.0;
	float EIsup=15;

    //Fill tree differently (it depends on the event type)

    if (event_type == TString("S1Beta"))
    {

        int cS1Beta = 0;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_S1Beta->GetRandom(ECinf, ECsup);
            float evt_heat = func_S1Beta->GetRandom(0,20);
            float strag=evt_heat*0.1;

            EIA = func_convS1Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIC = u.Gaus(0,TMath::Sqrt(pow(s_IC,2) - pow(s_EC1_ERA, 2)));
            EIB = func_convS1Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EID = u.Gaus(0,TMath::Sqrt(pow(s_ID,2) - pow(s_EC1_ERA, 2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cS1Beta++;
            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                    }
            }
        }

        string tmpstr;
        stringstream convert;
        convert << cS1Beta/float(num_event);
        TString s_cS1Beta = convert.str();
        string pattern = "S1Beta";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cS1Beta);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) {feffout << vec_str[u] << endl;} //; cout << vec_str[u] << endl;}
        feffout.close();
    }

    if (event_type == TString("S2Beta"))
    {

        int cS2Beta;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_S2Beta->GetRandom(ECinf, ECsup);
            float evt_heat = func_S2Beta->GetRandom(0,20);
            float strag=evt_heat*0.1;

            EIA = u.Gaus(0,TMath::Sqrt(pow(s_IA,2) - pow(s_EC1_ERA, 2)));
            EIC = func_convS2Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIB = u.Gaus(0,TMath::Sqrt(pow(s_IB,2) - pow(s_EC1_ERA, 2)));
            EID = func_convS2Beta->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cS2Beta++;
            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                    }
            }
        }

        string tmpstr;
        stringstream convert;
        convert << cS2Beta/float(num_event);
        TString s_cS2Beta = convert.str();
        string pattern = "S2Beta";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cS2Beta);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) feffout << vec_str[u] << endl;
        feffout.close();

    }

    if (event_type == TString("S1Pb"))
    {

        int cS1Pb;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_S1Pb->GetRandom(ECinf, ECsup);
            float evt_heat = func_S1Pb->GetRandom(0,20);
            float strag=evt_heat*0.05;

            EIA = func_convS1Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIC = u.Gaus(0,TMath::Sqrt(pow(s_IC,2) - pow(s_EC1_ERA, 2)));
            EIB = func_convS1Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EID = u.Gaus(0,TMath::Sqrt(pow(s_ID,2) - pow(s_EC1_ERA, 2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cS1Pb++;
            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                    }
            }
        }

        string tmpstr;
        stringstream convert;
        convert << cS1Pb/float(num_event);
        TString s_cS1Pb = convert.str();
        string pattern = "S1Pb";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cS1Pb);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) {feffout << vec_str[u] << endl;}
        feffout.close();

    }

    if (event_type == TString("S2Pb"))
    {

        int cS2Pb;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_S2Pb->GetRandom(ECinf, ECsup);
            float evt_heat = func_S2Pb->GetRandom(0,20);
            float strag=evt_heat*0.05;

            EIA = u.Gaus(0,TMath::Sqrt(pow(s_IA,2) - pow(s_EC1_ERA, 2)));
            EIC = func_convS2Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EIB = u.Gaus(0,TMath::Sqrt(pow(s_IB,2) - pow(s_EC1_ERA, 2)));
            EID = func_convS2Pb->Eval(evt_heat) + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) + pow(strag,2) - pow(s_EC1_ERA, 2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) - pow(s_EC1_ERA, 2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cS2Pb++;
            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= fEE_to_NR->Eval(EC);
                        t_new->Fill();
                    }
            }
        }

        string tmpstr;
        stringstream convert;
        convert << cS2Pb/float(num_event);
        TString s_cS2Pb = convert.str();
        string pattern = "S2Pb";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cS2Pb);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) feffout << vec_str[u] << endl;
        feffout.close();

    }

    t_new->Write();
    ftree->Close();


    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_1.5_0_3");
    int num_event[4];

	num_event[0] =200000;
	num_event[1] =200000;
	num_event[2] =40000;
	num_event[3] =30000;

    build_BetaPb_tree(bolo_name,  analysis_type, TString("S1Beta"),  num_event[0]);
    build_BetaPb_tree(bolo_name,  analysis_type, TString("S2Beta"),  num_event[1]);
    build_BetaPb_tree(bolo_name,  analysis_type, TString("S1Pb"),  num_event[2]);
    build_BetaPb_tree(bolo_name,  analysis_type, TString("S2Pb"),  num_event[3]);


    return 0;

}

        
