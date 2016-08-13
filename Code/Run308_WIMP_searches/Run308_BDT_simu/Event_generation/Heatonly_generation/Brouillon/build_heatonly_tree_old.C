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
#include "TPaletteAxis.h"
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

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
    TFile * file_heatonly = new TFile(heatonly_path, "read");
    TH2F * heatonly2D = (TH2F*)file_heatonly->Get("heat2D");

    TString heatonly_out_path = gen_path + TString("/BDT_") + bolo_name + TString("/") + analysis_type + TString("/Heatonly/ROOT_files/");
    TFile *ftree = new TFile(heatonly_out_path + bolo_name+ TString("_heatonly_tree.root"), "recreate");
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

    int cheatonly = 0;

    for (int i=0; i<num_event; i++) 
    {
        double heat_EC1, heat_EC2;

        heatonly2D->GetRandom2(heat_EC1, heat_EC2);

        EIA = u.Gaus(0,s_IA);
        EIC = u.Gaus(0,s_IC);
        EIB = u.Gaus(0,s_IB);
        EID = u.Gaus(0,s_ID);
        EC1 = heat_EC1 ;
        EC2 = heat_EC2 ; 

			float EC =0.769*EC1+0.231*EC2;
			float EI =0.422*EIB+0.578*EID;
        float accept = (f_wave_eff->Eval(EC))*0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

        if (abs(EC1-EC2)<1 && ECinf < EC && EC < ECsup) cheatonly++;
        if (abs(EC1-EC2)<1 &&  ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
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
    ftree->Close();

    return 0;

}

int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_1.5_0_3");
	int num_event =200000;

    build_heatonly_tree(bolo_name,  analysis_type, num_event);


    return 0;

}