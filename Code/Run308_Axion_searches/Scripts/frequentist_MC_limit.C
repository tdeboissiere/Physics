
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
#include "TF2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <fstream>
#include <sstream>
#include <string>
#include "TMatrixD.h"
#include "TArrayD.h"
using namespace std;

int Nobs_true;
vector<double> pred_points;
vector<double> s_array;
vector<double> b_array;

double approx_factorial(int n)
{

    if (n<5) 
    {
        return TMath::Log(TMath::Factorial(n)); 
    } 
        
    else 
    {
        double fact = 0.5*TMath::Log(TMath::Pi()) + n*TMath::Log(n/TMath::Exp(1)) + (1./6.)*TMath::Log(8*pow(n,3)+4*pow(n,2)+n+1./30.);
        return fact;
    }
        
}

double likelihood(double * x, double * par)
{


}




int main() 
{

    //Initialize Random generation
    gRandom->SetSeed(0);
    TRandom3 u(0);

    TString axion_mass = TStrin("0");

    TFile * file_data = new TFile("../ROOT_files/Axion/Spec/Stacked_spec_perkeV.root", "read");
    TFile * file_model = new TFile("../ROOT_files/Axion/Model/Stacked_model_perkeV.root", "read");
    TFile * file_flux = new TFile("../ROOT_files/Axion/CBRD_convolved/Stacked_flux_perkeV.root", "read");

    hdata = (TH1F*)file_data->Get("hstacked");
    fmodel = (TF1*)file_model->Get("model_stacked");
    fflux = (TF1*)file_flux->Get(TString("flux_stacked_mass_") + axion_mass);

    //Number of simulations
    int simu_number = 10;

    // Read bckg and signal PDF
    TFile f("./ROOT_files/FID837_BDT_bckg_and_sig_4GeV.root", "read");
    TF1* f_bckg_norm = (TF1*)f.Get("bckg_norm");
    TF1* f_signal_norm = (TF1*)f.Get("signal_norm");

    //Output best fit value TH1F 
    TH1F * hbestfit = new TH1F("hbestfit", "hbestfit", 100, Nsignal-5*sqrt(Nsignal), Nsignal + 5*sqrt(Nsignal)); 
    TH1F * hstat = new TH1F("hstat", "hstat", 200, 0,5); 

    vector<double > test_stat_vec;

    for (int nsimu =0; nsimu<simu_number; nsimu++)
    {

        TString s_simu = Form("%d", nsimu);
        TH1F * hsimu = new TH1F(TString("hsimu") + s_simu, TString("hsimu") + s_simu, 300, 0, 15)


        //Fill s_array and b_array 
        for (int k = 0; k<pred_points.size(); k++)
        {
            s_array.push_back(f_signal_norm->Eval(pred_points[k]));
            b_array.push_back(f_bckg_norm->Eval(pred_points[k]));
        }



    }

    ofstream file("./Text_files/test_stat_obs_Nsignal125.txt");
    for (int k = 0; k<test_stat_vec.size(); k++)
    {
        file << test_stat_vec[k] << endl;
    }

    return 0;

}

