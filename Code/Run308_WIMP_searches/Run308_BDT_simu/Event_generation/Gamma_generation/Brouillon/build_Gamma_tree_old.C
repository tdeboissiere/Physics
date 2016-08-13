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


using namespace std;

int build_Gamma_tree(TString bolo_name, TString  analysis_type, TString event_type,  int num_event)
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

    // Sample the 1D Gamma bckg model
    TString Ana_path   = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/");
    TFile * file_Gamma = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_Gamma_spectrum_extrapol.root"), "read");
    TF1 * func_FidGamma = (TF1*)file_Gamma->Get("FidGamma_extra");
    TF1 * func_S1Gamma = (TF1*)file_Gamma->Get("S1Gamma_extra");
    TF1 * func_S2Gamma = (TF1*)file_Gamma->Get("S2Gamma_extra");


    //Luke factor
    double Luke_factor = (1+8./3)/(1+5.5/3);

    TString Gamma_out_path = gen_path + TString("/BDT_") + bolo_name + TString("/") + analysis_type + TString("/Gamma/ROOT_files/");
    TFile *ftree = new TFile(Gamma_out_path + bolo_name+"_"+ event_type + "_tree.root", "recreate");
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

    if (event_type == TString("FidGamma"))
    {

        int cFidGamma=0;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_FidGamma->GetRandom(ECinf, ECsup);
            float evt_heat = func_FidGamma->GetRandom(0,20);

            EIA = u.Gaus(0,TMath::Sqrt(pow(s_IA,2) -pow(s_EC1_ERA,2)));
            EIC = u.Gaus(0,TMath::Sqrt(pow(s_IC,2) -pow(s_EC1_ERA,2)));
            EIB = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) -pow(s_EC1_ERA,2)));
            EID = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) -pow(s_EC1_ERA,2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) -pow(s_EC1_ERA,2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );
            
            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cFidGamma++;
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
        convert << cFidGamma/float(num_event);
        TString s_cFidGamma = convert.str();
        string pattern = "FidGamma";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cFidGamma);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) {feffout << vec_str[u] << endl;} //; cout << vec_str[u] << endl;}
        feffout.close();

    }

    if (event_type == TString("S1Gamma"))
    {

        int cS1Gamma = 0;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_S1Gamma->GetRandom(ECinf, ECsup);
            float evt_heat = func_S1Gamma->GetRandom(0,20);

            EIA = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IA,2) -pow(s_EC1_ERA,2)));
            EIC = u.Gaus(0,TMath::Sqrt(pow(s_IC,2) -pow(s_EC1_ERA,2)));
            EIB = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IB,2) -pow(s_EC1_ERA,2)));
            EID = u.Gaus(0,TMath::Sqrt(pow(s_ID,2) -pow(s_EC1_ERA,2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) -pow(s_EC1_ERA,2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cS1Gamma++;
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
        convert << cS1Gamma/float(num_event);
        TString s_cS1Gamma = convert.str();
        string pattern = "S1Gamma";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cS1Gamma);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) {feffout << vec_str[u] << endl;} //; cout << vec_str[u] << endl;}
        feffout.close();

    }

    if (event_type == TString("S2Gamma"))
    {

        int cS2Gamma=0;

        for (int i=0; i<num_event; i++) 
        {

            // float evt_heat = func_S2Gamma->GetRandom(ECinf, ECsup);
            float evt_heat = func_S2Gamma->GetRandom(0,20);

            EIA = u.Gaus(0,TMath::Sqrt(pow(s_IA,2) -pow(s_EC1_ERA,2)));
            EIC = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_IC,2) -pow(s_EC1_ERA,2)));
            EIB = u.Gaus(0,TMath::Sqrt(pow(s_IB,2) -pow(s_EC1_ERA,2)));
            EID = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_ID,2) -pow(s_EC1_ERA,2)));
            EC1 = evt_heat ;
            EC2 = evt_heat + u.Gaus(0,TMath::Sqrt(pow(s_EC2_ERA,2) -pow(s_EC1_ERA,2)));

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
            float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cS2Gamma++;
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
        convert << cS2Gamma/float(num_event);
        TString s_cS2Gamma = convert.str();
        string pattern = "S2Gamma";
        vector <TString> vec_str;
        ifstream feffin(txt_path);
        while (getline(feffin, tmpstr))
        {
            size_t found = tmpstr.find(pattern);
            if (found!=string::npos) vec_str.push_back(TString(pattern) + TString(",") + s_cS2Gamma);
            else vec_str.push_back(TString(tmpstr));
        }
        feffin.close();
        ofstream feffout(txt_path);
        for (int u = 0; u<vec_str.size(); u++) {feffout << vec_str[u] << endl;} //; cout << vec_str[u] << endl;}
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
    int num_event[3];

	num_event[0] =10000;
	num_event[1] =10000000;
	num_event[2] =250000000;

    build_Gamma_tree(bolo_name,  analysis_type, TString("FidGamma"),  num_event[0]);
    build_Gamma_tree(bolo_name,  analysis_type, TString("S1Gamma"),  num_event[1]);
    build_Gamma_tree(bolo_name,  analysis_type, TString("S2Gamma"),  num_event[2]);


    return 0;

}