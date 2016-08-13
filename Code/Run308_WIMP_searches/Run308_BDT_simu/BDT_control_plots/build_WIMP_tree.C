
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

    outfilename= TString("./ROOT_files/") + bolo_name + TString("/") + analysis_type + TString("/") + bolo_name + TString("_WIMP_mass_")+s_Wmass+TString("_tree_atlimit.root");
    TFile *ftree = new TFile(outfilename, "recreate");

    //Get the correlation matrix
    TString corr_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Analyse_") + bolo_name + TString("/Text_files/");
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

    float heat_thresh = 0.65;

	float s_heat=0.138865296416;
	float s_EC1=0.157125870562;
	float s_EC2=0.180907083404;
	float s_IA=0.38644470868;
	float s_IB=0.354170205538;
	float s_IC=0.304484457279;
	float s_ID=0.297265160523;

	float cut_vetA=1.9322235434;
	float cut_vetC=1.52242228639;

	float ECinf=0.5;
	float ECsup=15;
	float EIinf=-4;
	float EIsup=15;

    // C r e a t e   W I M P    f u n c t i o n
    // ---------------------------------------
    TF1* frecoil = new TF1("frecoil", computerate_sfg, 0, 30, 1);
    frecoil->SetNpx(1000);
    frecoil->SetParameter(0,Wmass);
    float ERsimu, Qsimu;

    //2 simulations: one for training, one for application
    for (int i =0 ; i<2 ; i++)
    {
        TString s_index = Form("%d", i);


        //New tree to be filled
        TTree *t_new = new TTree("t_new" + s_index,"t_new" + s_index);
        TTree *t_newnocut = new TTree("t_newnocut" + s_index,"t_newnocut" + s_index);

        float  EC1, EC2, EIA, EIB, EIC, EID, ENR;
        t_new->Branch("EC1", &EC1, "EC1/F");
        t_new->Branch("EC2", &EC2, "EC2/F");
        t_new->Branch("EIA", &EIA, "EIA/F");
        t_new->Branch("EIB", &EIB, "EIB/F");
        t_new->Branch("EIC", &EIC, "EIC/F");
        t_new->Branch("EID", &EID, "EID/F");
        t_new->Branch("ENR", &ENR, "ENR/F");

        t_newnocut->Branch("EC1", &EC1, "EC1/F");
        t_newnocut->Branch("EC2", &EC2, "EC2/F");
        t_newnocut->Branch("EIA", &EIA, "EIA/F");
        t_newnocut->Branch("EIB", &EIB, "EIB/F");
        t_newnocut->Branch("EIC", &EIC, "EIC/F");
        t_newnocut->Branch("EID", &EID, "EID/F");
        t_newnocut->Branch("ENR", &ENR, "ENR/F");

        //Loop over the events
        for (int k=0; k<num_events; k++) 
        {

            ERsimu = frecoil->GetRandom(0,30);
            Qsimu = 0.16*pow(ERsimu,0.18);

            TArrayD reso_coeff(4);
            for (int k = 0; k<4; k++) {reso_coeff[k] = u.Gaus(0,1);}
            TMatrixD reso_matrx(4,1);
            reso_matrx.SetMatrixArray(reso_coeff.GetArray());
            TMatrixD ion_matrx = Lmatrx*reso_matrx;

            EIA = ion_matrx(0,0);
            EIC = ion_matrx(2,0);
            EC1 = (1+Qsimu*V_fid/3.)/(1+V_fid/3.)*ERsimu + u.Gaus(0,s_EC1);
            EC2 = (1+Qsimu*V_fid/3.)/(1+V_fid/3.)*ERsimu + u.Gaus(0,s_EC2);     
            EIB = ERsimu*Qsimu + ion_matrx(1,0);
            EID = ERsimu*Qsimu + ion_matrx(3,0);

			float EC =0.570*EC1+0.43*EC2;
			float EI =0.413*EIB+0.587*EID;
            float accept = 0.5*(1 + TMath::Erf( (EC-heat_thresh) / (TMath::Sqrt(2)*s_heat) ) );

            ENR=ERsimu;
            if (abs(EC1-EC2)< 1 && ECinf < EC && EC < ECsup && EIinf<EI  && EI<EIsup && EIA < cut_vetA && EIC < cut_vetC) 
            {
                    if (u.Binomial(1,0.9999*accept)==1) 
                    {
                        ENR= ERsimu;
                        t_new->Fill();
                    }
            }

            
        }

        t_new->Write();    } 
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
	TString analysis_type = TString("ana_0.5_min4_5");

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



    // Signal at limit for ana_0.5_min4_5
    Wmass_vec.push_back(3);num_events.push_back(int(1128*10000000./11169));
    Wmass_vec.push_back(4);num_events.push_back(int(1360*1000000/15328));
    Wmass_vec.push_back(5);num_events.push_back(int(286*200000/12356));
    Wmass_vec.push_back(6);num_events.push_back(int(122*100000/13749));
    Wmass_vec.push_back(7);num_events.push_back(int(48*50000/11166));
    Wmass_vec.push_back(10);num_events.push_back(int(9*30000/13522));
    Wmass_vec.push_back(25);num_events.push_back(int(7*20000/16355));

    // For ana_0.5_min4_5
    // Wmass_vec.push_back(3);num_events.push_back(10000000);
    // Wmass_vec.push_back(4);num_events.push_back(1000000);
    // Wmass_vec.push_back(5);num_events.push_back(200000);
    // Wmass_vec.push_back(6);num_events.push_back(100000);
    // Wmass_vec.push_back(7);num_events.push_back(50000);
    // Wmass_vec.push_back(10);num_events.push_back(30000);
    // Wmass_vec.push_back(25);num_events.push_back(20000);


    // We obtain the following number of events  below ana_0.5_min4_5
    // mass3 = 11169
    // mass4 = 15328
    // mass5 = 12356
    // mass6 = 13749
    // mass7 = 11166
    // mass10 = 13522
    // mass25 = 16355



    // Signal at limit for ana_0.5_0_5
	// Wmass_vec.push_back(3);num_events.push_back(int(684*10000000./8463));
	// Wmass_vec.push_back(4);num_events.push_back(int(1018*1000000/13418));
	// Wmass_vec.push_back(5);num_events.push_back(int(236*200000/11293));
	// Wmass_vec.push_back(6);num_events.push_back(int(106*100000/12872));
	// Wmass_vec.push_back(7);num_events.push_back(int(43*50000/10574));
	// Wmass_vec.push_back(10);num_events.push_back(int(9*30000/13112));
	// Wmass_vec.push_back(25);num_events.push_back(int(7*20000/16235));

    // For the num_events below ana_0.5_0_5
    // Wmass_vec.push_back(3);num_events.push_back(10000000);
    // Wmass_vec.push_back(4);num_events.push_back(1000000);
    // Wmass_vec.push_back(5);num_events.push_back(200000);
    // Wmass_vec.push_back(6);num_events.push_back(100000);
    // Wmass_vec.push_back(7);num_events.push_back(50000);
    // Wmass_vec.push_back(10);num_events.push_back(30000);
    // Wmass_vec.push_back(25);num_events.push_back(20000);

    // We obtain the following number of events  below ana_0.5_0_5
    // mass3 =8463
    // mass4 =13418
    // mass5 =11293
    // mass6 =12872
    // mass7 =10574
    // mass10 =13112
    // mass25 =16235


    intermediate(Wmass_vec, num_events, bolo_name, analysis_type);


    return 0;

}