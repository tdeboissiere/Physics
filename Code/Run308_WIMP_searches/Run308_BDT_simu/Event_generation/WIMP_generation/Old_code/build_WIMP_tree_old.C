
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
#include "WIMP_recoil_generation.C"
#include <string>

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


int fill_root_tree (int Wmass, int num_events, TString bolo_name, TString analysis_type)
{

    ///--------------------------------------------------------------------------------------------------///
    ///
    /// Build a tree 
    /// Arguments:
    /// Wmass       => int      => Wimp mass
    /// num_events  => int      => Number of events in the .txt file
    /// V_fid       => float    => Fiduvial bias
    /// s_**        => float    => resolution (FWHM/2.3548) for the corresponding channel
    ///
    /// Ouput:
    ///     a .root file containing EIA/B/C/D and EC1/EC2 simulated from a list of Er
    ///--------------------------------------------------------------------------------------------------///

    // R e a d   f i l e    a n d   c r e a t e     v e c t o r s
    // --------------------------------------------------------------------
    TString  s_Wmass= Form("%d", Wmass); 
    TString filename;
    
    TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");


    // C r e a t e     t r e e
    // ---------------------------------------
    TString outfilename;

    //Include NR computations directly
    TString util_path = TString("/home/irfulx204/mnt/tmain/Desktop/Miscellaneous/Python/Useful_scripts/Utilities/");
    TFile * file_EE_to_NR = new TFile(util_path + TString("conv_EE_to_NR.root"), "read");
    TF1* fEE_to_NR = (TF1*)file_EE_to_NR->Get("conv");

    // //Include wavelet efficiency
    // TString swave_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/") + bolo_name + TString("/");
    // TFile file_wave(swave_path + TString("PSA_cut_eff.root"), "read");
    // f_wave_eff = (TF1*)file_wave.Get("PSA_eff");

    //Get the histogram of the baselines for event generation
    TString Ana_path   = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/");
    TFile * file_fwhm = new TFile(Ana_path+ TString("Analyse_") + bolo_name + TString("/ROOT_files/") + bolo_name +  TString("_fwhm_time_hist_for_BDT.root"), "read"); 
    TH1F* hFWIA = (TH1F*)file_fwhm->Get("hfwia");
    TH1F* hFWIB = (TH1F*)file_fwhm->Get("hfwib");
    TH1F* hFWIC = (TH1F*)file_fwhm->Get("hfwic");
    TH1F* hFWID = (TH1F*)file_fwhm->Get("hfwid");
    TH1F* hOWC1 = (TH1F*)file_fwhm->Get("howc1");
    TH1F* hOWC2 = (TH1F*)file_fwhm->Get("howc2");
    TH1F* hJOUR = (TH1F*)file_fwhm->Get("hist_jour");

    outfilename=gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + TString("/WIMP/ROOT_files/") + bolo_name + TString("_WIMP_mass_")+s_Wmass+TString("_tree.root");
    TFile *ftree = new TFile(outfilename, "recreate");

    // F i l l     t r e e
    // ---------------------------------------
    //Prepare random generator
    gRandom->SetSeed(0);
    TRandom3 u(0);

    //Voltages
    float V_fid=8;
    float V_surf=5.5;

	float cut_vetA=1.9322235434;
	float cut_vetC=1.52242228639;

	float ECinf=0.5;
	float ECsup=15;
	float EIinf=0.0;
	float EIsup=15;

    // C r e a t e   W I M P    f u n c t i o n
    // ---------------------------------------
    TF1* frecoil = new TF1("frecoil", computerate_sfg, 0, 30, 1);
    frecoil->SetParameter(0,Wmass);
    float ERsimu, Qsimu;

    //2 simulations: one for training, one for application
    for (int i =0 ; i<2 ; i++)
    {
        TString s_index = Form("%d", i);

        // //Load simulated events
        // filename=gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + TString("/WIMP/ROOT_files/recoils_mass_")+s_Wmass+TString("_GeV_") +s_index + TString(".txt");

        // ifstream file(filename);    
        // float tempflt;

        // vector<float> recoil_vec;
        // vector<float> Q_vec;

        // for (int k = 0; k < num_events; k += 1)
        // {
        //     file >> tempflt;
        //     recoil_vec.push_back(tempflt); 
        //     Q_vec.push_back(0.16*pow(tempflt,0.18));   
        // }   


        //New tree to be filled
        TTree *t_new = new TTree("t_new" + s_index,"t_new" + s_index);

        float  EC1, EC2, EIA, EIB, EIC, EID, ENR;
        t_new->Branch("EC1", &EC1, "EC1/F");
        t_new->Branch("EC2", &EC2, "EC2/F");
        t_new->Branch("EIA", &EIA, "EIA/F");
        t_new->Branch("EIB", &EIB, "EIB/F");
        t_new->Branch("EIC", &EIC, "EIC/F");
        t_new->Branch("EID", &EID, "EID/F");
        t_new->Branch("ENR", &ENR, "ENR/F");

        //Loop over the events
        for (int k=0; k<num_events; k++) 
        {

            ERsimu = frecoil->GetRandom(0,30);
            Qsimu = 0.16*pow(ERsimu,0.18);

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
            EIC = u.Gaus(0,s_IC);
            // EC1 = (1+Q_vec[k]*V_fid/3.)/(1+V_fid/3.)*recoil_vec[k] + u.Gaus(0,s_EC1);
            // EC2 = (1+Q_vec[k]*V_fid/3.)/(1+V_fid/3.)*recoil_vec[k] + u.Gaus(0,s_EC2);     
            // EIB = recoil_vec[k]*Q_vec[k] + u.Gaus(0,s_IB);
            // EID = recoil_vec[k]*Q_vec[k] + u.Gaus(0,s_ID);

            EC1 = (1+Qsimu*V_fid/3.)/(1+V_fid/3.)*ERsimu + u.Gaus(0,s_EC1);
            EC2 = (1+Qsimu*V_fid/3.)/(1+V_fid/3.)*ERsimu + u.Gaus(0,s_EC2);     
            EIB = ERsimu*Qsimu + u.Gaus(0,s_IB);
            EID = ERsimu*Qsimu + u.Gaus(0,s_ID);

            float EC =coeff_EC1*EC1+(1-coeff_EC1)*EC2;
            float EI =coeff_EIB*EIB+(1-coeff_EIB)*EID;
            double accept = 0.5*(1+TMath::Erf((EC-0.65)/(TMath::Sqrt(2)*s_heat)));

            if (abs(EC1-EC2)< 1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        // ENR= recoil_vec[k];
                        ENR= ERsimu;
                        t_new->Fill();
                    }
            }

            
        }

        t_new->Write();
    } 
    ftree->Close();


    return 0;
}


int intermediate(vector<int>& Wmass_vec, vector<int>& num_events, TString bolo_name, TString analysis_type) 
{

    ///--------------------------------------------------------------------------------------------------///
    ///
    /// Intermediate function 
    /// Arguments:
    /// Wmass_vec   => vector   => Wimp mass vector
    /// num_events  => int      => Number of events in the .txt file
    /// V_fid       => float    => Fiduvial bias
    /// s_**        => float    => resolution (FWHM/2.3548) for the corresponding channel
    ///
    /// Ouput:
    ///     a .root file containing EIA/B/C/D and EC1/EC2 simulated from a list of Er
    ///--------------------------------------------------------------------------------------------------///

    for(std::vector<int>::iterator it_mass = Wmass_vec.begin(); it_mass != Wmass_vec.end(); it_mass++) 
    {
        cout << "**********************************************************************" << endl;
        cout << "Now creating tree for WIMPS with " << *it_mass << " GeV mass" << endl;
        cout << "**********************************************************************" << endl;
        fill_root_tree (*it_mass, num_events[it_mass-Wmass_vec.begin()], bolo_name, analysis_type);
    }

    return 0;
}


int main()
{


	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_0.5_0_5");

    ///--------------------------------------------------------------------------------------------------///
    ///
    /// Build a tree for each mass stored in Wmass_vec
    /// Arguments:
    ///     num_events  => int      => the number of events to be simulated
    /// Ouput:
    ///     a .root file containing EIA/B/C/D and EC1/EC2 simulated from a list of Er
    ///--------------------------------------------------------------------------------------------------///

    //------------------------------///
    // Other parameters
    //------------------------------///

    vector<int> Wmass_vec;
    vector<int> num_events;

	Wmass_vec.push_back(5);num_events.push_back(200000);
	Wmass_vec.push_back(6);num_events.push_back(100000);
	Wmass_vec.push_back(7);num_events.push_back(50000);
	Wmass_vec.push_back(10);num_events.push_back(30000);
	Wmass_vec.push_back(25);num_events.push_back(20000);

    intermediate(Wmass_vec, num_events, bolo_name, analysis_type);


    return 0;

}